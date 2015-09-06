#' Predict data using randomForest object, tracing variable usage along the way
#' This is essentially randomForest::predict with some modifications for
#' monitoring variable frequency, abundance, and deficiency
#' @export
trace.forest <-
function (object, newdata, response=NULL, type = "response", norm.votes = TRUE, proximity = FALSE, nodes=FALSE, cutoff, ...)
{
    predict.all<-TRUE
    if (!inherits(object, "randomForest"))
    stop("object not of class randomForest")
    if (is.null(object$forest)) stop("No forest component in the object")
    out.type <- 1
    if (is.na(out.type))
    stop("type must be one of 'response', 'prob', 'vote'")
    if (out.type != 1 && object$type == "regression")
    stop("'prob' or 'vote' not meaningful for regression")
    if (out.type == 2)
    norm.votes <- TRUE
    if (missing(newdata)) {
        p <- if (! is.null(object$na.action)) {
            napredict(object$na.action, object$predicted)
        } else {
            object$predicted
        }
        if (object$type == "regression") return(p)
        if (proximity & is.null(object$proximity))
        warning("cannot return proximity without new data if random forest object does not already have proximity")
        if (out.type == 1) {
            if (proximity) {
                return(list(pred = p,
                proximity = object$proximity))
            } else return(p)
        }
        v <- object$votes
        if (!is.null(object$na.action)) v <- napredict(object$na.action, v)
        if (norm.votes) {
            t1 <- t(apply(v, 1, function(x) { x/sum(x) }))
            class(t1) <- c(class(t1), "votes")
            if (proximity) return(list(pred = t1, proximity = object$proximity))
            else return(t1)
        } else {
            if (proximity) return(list(pred = v, proximity = object$proximity))
            else return(v)
        }
    }
    if (missing(cutoff)) {
        cutoff <- object$forest$cutoff
    } else {
        if (sum(cutoff) > 1 || sum(cutoff) < 0 || !all(cutoff > 0) ||
        length(cutoff) != length(object$classes)) {
            stop("Incorrect cutoff specified.")
        }
        if (!is.null(names(cutoff))) {
            if (!all(names(cutoff) %in% object$classes)) {
                stop("Wrong name(s) for cutoff")
            }
            cutoff <- cutoff[object$classes]
        }
    }

    if (inherits(object, "randomForest.formula")) {
        newdata <- as.data.frame(newdata)
        rn <- row.names(newdata)
        Terms <- delete.response(object$terms)
        x <- model.frame(Terms, newdata, na.action = na.omit)
        keep <- match(row.names(x), rn)
    } else {
        if (is.null(dim(newdata)))
        dim(newdata) <- c(1, length(newdata))
        x <- newdata
        if (nrow(x) == 0)
        stop("newdata has 0 rows")
        if (any(is.na(x)))
        stop("missing values in newdata")
        keep <- 1:nrow(x)
        rn <- rownames(x)
        if (is.null(rn)) rn <- keep
    }
    vname <- if (is.null(dim(object$importance))) {
        names(object$importance)
    } else {
        rownames(object$importance)
    }
    if (is.null(colnames(x))) {
        if (ncol(x) != length(vname)) {
            stop("number of variables in newdata does not match that in the training data")
        }
    } else {
        if (any(! vname %in% colnames(x)))
        stop("variables in the training data missing in newdata")
        x <- x[, vname, drop=FALSE]
    }
    if (is.data.frame(x)) {
        isFactor <- function(x) is.factor(x) & ! is.ordered(x)
        xfactor <- which(sapply(x, isFactor))
        if (length(xfactor) > 0 && "xlevels" %in% names(object$forest)) {
            for (i in xfactor) {
                if (any(! levels(x[[i]]) %in% object$forest$xlevels[[i]]))
                stop("New factor levels not present in the training data")
                x[[i]] <-
                factor(x[[i]],
                levels=levels(x[[i]])[match(levels(x[[i]]), object$forest$xlevels[[i]])])
            }
        }
        cat.new <- sapply(x, function(x) if (is.factor(x) && !is.ordered(x))
        length(levels(x)) else 1)
        if (!all(object$forest$ncat == cat.new))
        stop("Type of predictors in new data do not match that of the training data.")
    }
    mdim <- ncol(x)
    ntest <- nrow(x)
    ntree <- object$forest$ntree
    maxcat <- max(object$forest$ncat)
    nclass <- object$forest$nclass
    nrnodes <- object$forest$nrnodes
    ## get rid of warning:
    op <- options(warn=-1)
    on.exit(options(op))
    x <- t(data.matrix(x))

    if (predict.all) {
        treepred <- matrix(integer(ntest * ntree), ncol=ntree)
    } else {
        treepred <- numeric(ntest)
    }
    proxmatrix <- if (proximity) matrix(0, ntest, ntest) else numeric(1)
    nodexts <- if (nodes) integer(ntest * ntree) else integer(ntest)

    # frequency of each feature during correct classification of each class
    abundant = integer(mdim*nclass)
    deficient = integer(mdim*nclass)
    frequency = integer(mdim*nclass)

    countts <- matrix(0, ntest, nclass)
    # use my own C prediction routine to get frequency information
    t1 <- .C("_predict",
        mdim = as.integer(mdim),
        ntest = as.integer(ntest),
        nclass = as.integer(object$forest$nclass),
        maxcat = as.integer(maxcat),
        nrnodes = as.integer(nrnodes),
        jbt = as.integer(ntree),
        xts = as.double(x),
        response = as.integer(response),
        xbestsplit = as.double(object$forest$xbestsplit),
        pid = object$forest$pid,
        cutoff = as.double(cutoff),
        countts = as.double(countts),
        treemap = as.integer(aperm(object$forest$treemap,
        c(2, 1, 3))),
        nodestatus = as.integer(object$forest$nodestatus),
        cat = as.integer(object$forest$ncat),
        nodepred = as.integer(object$forest$nodepred),
        treepred = as.integer(treepred),
        jet = as.integer(numeric(ntest)),
        bestvar = as.integer(object$forest$bestvar),
        nodexts = as.integer(nodexts),
        ndbigtree = as.integer(object$forest$ndbigtree),
        predict.all = as.integer(predict.all),
        prox = as.integer(proximity),
        proxmatrix = as.double(proxmatrix),
        nodes = as.integer(nodes),
        frequency = as.integer(frequency),
        abundant = as.integer(abundant),
        deficient = as.integer(deficient),
        DUP=FALSE,
        PACKAGE = "leaves")

    out.class <- factor(rep(NA, length(rn)),
    levels=1:length(object$classes),
    labels=object$classes)
    out.class[keep] <- object$classes[t1$jet]
    names(out.class)[keep] <- rn[keep]
    res <- out.class
    f <- matrix(t1$frequency, nrow=mdim)
    a <- matrix(t1$abundant, nrow=mdim)
    d <- matrix(t1$deficient, nrow=mdim)
    dimnames(f) <- list(NULL, object$classes)
    dimnames(a) <- list(NULL, object$classes)
    dimnames(d) <- list(NULL, object$classes)

    if (predict.all) {
        treepred <- matrix(object$classes[t1$treepred],
        nrow=length(keep), dimnames=list(rn[keep], NULL))
        # this is where the magic happens
        res <- list(aggregate=res, individual=treepred,
                    frequency=f,
                    abundant=a,
                    deficient=d)
    }
    if (proximity)
    res <- list(predicted = res, proximity = structure(t1$proxmatrix,
    dim = c(ntest, ntest),
    dimnames = list(rn[keep], rn[keep])))
    if (nodes) attr(res, "nodes") <- matrix(t1$nodexts, ntest, ntree,
    dimnames=list(rn[keep], 1:ntree))

    #for(i in seq(0, 479)){
    #    cat(sprintf("[%d]\t%d\t%d\t%d\n", i, t1$frequency[i], t1$abundant[i], t1$deficient[i]));
    #}
    res
}

#' @useDynLib leaves _traceTree
# C code to trace through tree, counting frequency of variables for positive cases
traceTreeC <- function(tree, x, bestvar, bestsplit, noderep, endnode, prediction,
frequency, abundancy, deficiency){
    .C("_traceTree",
    tree= as.integer(t(tree[,c("left daughter", "right daughter", "status")])),
    x= as.double(x),
    bestvar = as.integer(tree[,3]),
    bestsplit = as.double(tree[,4]),
    nodepred = as.integer(rep(0,nrow(tree))),
    endnode = as.integer(0),    # currently unused
    prediction = as.integer(0), # currently unused
    frequency = as.integer(frequency),
    abundance = as.integer(abundancy),
    deficiency = as.integer(deficiency),
    DUP=FALSE, PACKAGE="leaves")
}

# R code to trace through tree, counting frequency of variables for positive cases
traceTreeR <- function(tree, x, noderep, endnode, prediction,
frequency, abundancy, deficiency){
    # count frequency during tree traversal
    j <- 1
    repeat{
        #cat("#k =", j, tree[j,1], tree[j,2], tree[j,5], tree[j,3], "\n")
        bestvar <- tree[j,3]
        frequency[bestvar] <- frequency[bestvar] + 1
        if(x[bestvar]<=tree[j,4]){
            #cat("\t", as.double(x[bestvar]), "<=", tree[j,4], "->",tree[j,1], "\n" )
            deficiency[bestvar] <- deficiency[bestvar] + 1
            j <- tree[j,1]
        } else {
            #cat("\t", as.double(x[bestvar]), ">", tree[j,4], "->", tree[j,2], "\n")
            abundancy[bestvar] <- abundancy[bestvar] + 1
            j <- tree[j,2]
        }

    #     break if terminal node
        if( tree[j,5] == -1 ) break
    }
    result <- list(frequency=frequency, abundance=abundancy, deficiency=deficiency)
    return(result)
}

generateForest <- function(df, ntree=500, mtry=50,
                keep.forest=TRUE, importance=TRUE, keep.inbag=TRUE){
    predictors <- df[,!(names(df) %in% c("class"))]
    response <- df$class
    randomForest::randomForest(x = predictors,
                 y = response,
                ntree=ntree,
                mtry=mtry,
                keep.forest=keep.forest,
                importance=importance,
                keep.inbag=keep.inbag)
}