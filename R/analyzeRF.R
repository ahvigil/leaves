# trace through an individual RF classifier
analyzeRF <- function(model, df, seed=5, verbose=FALSE,
                        rf_data_path=file.path(getwd(), "data", "rf"),
                        feature_data_path=file.path(getwd(), "data", "feature"),
                        output_path=file.path(getwd(), "output", model_name))
{
    mtry <- model$mtry
    ntree      <- model$ntree
    seed_mult <- model$seed_mult
    nForests <- model$forests
    model_name <- model$name

    # ggplot2 for generating graphics
    library(ggplot2)

    rf_params <- sprintf("%s_n%04d_m%02d_%04d", model_name, ntree, mtry, seed*seed_mult)
    model_file <- file.path(rf_data_path, sprintf("%s.model", rf_params))

    if(file.access(model_file)==0){
        load(model_file)
        if(verbose) cat("Loaded RF model from", model_file, "\n")
    } else {
        if(verbose) cat("RF model file", model_file, "does not exist!\n")
        set.seed( seed * seed_mult )

        # Train RF model
        if(verbose) cat("Constructing model...\n")
        model <- generateForest(df, ntree=ntree, mtry=mtry)

        # Save model file for future use.
        if(verbose) cat("Saving model...\n")
        save( model, file = model_file )
    }

    save_file <- file.path(output_path, "binary", sprintf("%s.save", rf_params))
    if(file.access(save_file)==-1){
        if(verbose) cat("No saved trace at", save_file, "\n")
        prediction <- trace.forest(model, df, response=df$class)

        if(verbose) cat(sprintf("prediction is %d by %d\n", nrow(prediction$individual), ncol(prediction$individual)))

        frequency <- prediction$frequency[,"1"]
        deficient <- prediction$deficient[,"1"]
        abundant <- prediction$abundant[,"1"]

        if(verbose) cat("\nRF trace complete\n")

        # get variable importance measures
        # positive, negative, mean decrease accuracy, mean decrease gini
        importance <- randomForest::importance(model)
        # variable forest-wide counts
        used <- randomForest::varUsed(model)

        # output matrix
        O <- cbind(1:480, frequency, importance, used, deficient, abundant)
        # assign column names
        dimnames(O)[[2]][1] <- "index"
        dimnames(O)[[2]][2] <- "frequency"
        # importance measures:
        dimnames(O)[[2]][3] <- "pVI"
        dimnames(O)[[2]][4] <- "nVI"
        dimnames(O)[[2]][5] <- "mdaVI"
        dimnames(O)[[2]][6] <- "mdgVI"

        dimnames(O)[[2]][7] <- "used"
        dimnames(O)[[2]][8] <- "deficient"
        dimnames(O)[[2]][9] <- "abundant"

        save(O,file = save_file)
    } else{
        load(save_file)
        if(verbose) cat("Loaded previous trace from", save_file, "\n")
    }

    # write output files
    labels<-names(df)

    sink(sprintf("%s/text/%s.trace.txt", output_path, rf_params))
    cat(model_name, "\n", sep = "")
    print(model)
    cat("\n", sep = "")
    cat("index\tfrequency\tabundant\tdeficient\tpVI\tnVI\tmdaVI\tmdgVI\tprevalence\tproperty\n", sep = "")
    for(i in 1:480){
        cat(O[i,"index"], "\t",
            O[i,"frequency"], "\t",
            O[i,"abundant"], "\t",
            O[i,"deficient"],"\t",
            O[i,"pVI"],"\t",
            O[i,"nVI"], "\t",
            O[i,"mdaVI"], "\t",
            O[i,"mdgVI"], "\t",
            O[i,"used"], "\t",
            labels[i], "\n", sep = "")
    }
    sink()

    result.frame <- as.data.frame(O)
    x_label <- rownames(result.frame)[order(-result.frame$frequency)]

    jpeg(filename=sprintf("%s/graphics/%s.jpeg", output_path, rf_params), width=600, height=600)
    p<-qplot(order(order(-frequency)), frequency, data=result.frame, geom="point", xlab="rank", ylab="frequency")
    p<-p+ geom_area(aes(y=deficient), colour = "red", fill="red") +
    geom_ribbon(aes(ymin=deficient, ymax=frequency), colour = "green", fill="green")
    suppressWarnings(print(p))
    invisible(dev.off())
    jpeg(filename=sprintf("%s/graphics/%s.top.jpeg", output_path, rf_params), width=600, height=600)
    p<-qplot(data=result.frame, geom="point", xlab="property", ylab="frequency") + aes(x=x_label, y=frequency[order(-frequency)]) + xlim(x_label[1:10]) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    p<-p+geom_area(aes(x=x_label,y=deficient[order(-frequency)]), colour = "red") +
    geom_ribbon(aes(x=x_label,ymin=deficient[order(-frequency)], ymax=frequency[order(-frequency)]), colour = "green") + xlab("property") + ylab("frequency")
    suppressWarnings(print(p))
    invisible(dev.off())
    return(O)
}

# analyze multiple RF classifiers for a dataset
#' @export
analyzeModel <- function(model, nForests, verbose=FALSE,
                        rf_data_path=file.path(getwd(), "data", "rf"),
                        feature_data_path=file.path(getwd(), "data", "feature"),
                        output_path=file.path(getwd(), "output", model_name)) {
    # ggplot2 for generating graphics
    library(ggplot2)

    mtry <- model$mtry
    ntree      <- model$ntree
    seed_mult <- model$seed_mult
    nForests <- model$forests
    model_name <- model$name

    # sensible defaults for data and output paths
    results <- array(0,c(480, 9, nForests))

    data_file <- file.path(feature_data_path, sprintf("%s.arff.gz", model_name))
    rf_params <- sprintf("%s_n%04d_m%02d", model_name, ntree, mtry)

    # create output paths if they don't already exist
    dir.create(file.path(output_path, "binary"), showWarnings = FALSE, recursive=TRUE)
    dir.create(file.path(output_path, "graphics"), showWarnings = FALSE, recursive=TRUE)
    dir.create(file.path(output_path, "text"), showWarnings = FALSE, recursive=TRUE)

    if(file.access(file.path(output_path, "binary", paste(rf_params, "save", sep=".")))==-1){
        if(verbose) cat("Reading data from", data_file, "\n")
        df <- foreign::read.arff(data_file)
        if(verbose) cat("done reading data\n")
        dimnames(results)[[1]] <- dimnames(df)[[2]][1:480]

        for(seed in 0:(nForests-1)){
            if(verbose) cat("Generate forest", seed, "\n")
            results[,,seed+1] <- analyzeRF(model, df, seed, verbose=verbose,
                                        feature_data_path=feature_data_path,
                                        rf_data_path=rf_data_path,
                                        output_path=output_path)
            colnames(results)<-c('index','frequency','pVI','nVI','mdaVI','mdgVI','prevalence', 'deficient', 'abundant')
        }

        save(results, file=file.path(output_path, "binary", paste(rf_params, "save", sep=".")))
    } else{
        load(sprintf("%s/binary/%s.save", output_path, rf_params))
    }

    # ranks is m x n matrix of frequency rank of variable m in forest n
    ranks <- array(0, c(480, 10))
    for(i in 1:nForests){
        ranks[,i] <- order(order(-results[,2,i]))
    }

    # mean, standard deviation, and standard error of frequency rank across forests
    rank_mean <- rowMeans(ranks)
    rank_sd <- matrixStats::rowSds(ranks)
    rank_sem <- rank_sd/sqrt(nForests)

    # mean, standard deviation, and standard error of frequency across forests
    frequency_mean <- rowMeans(results[,2,1:nForests])
    frequency_sd <- matrixStats::rowSds(results[,2,1:nForests])
    frequency_sem <- frequency_sd/sqrt(nForests)

    # mean/sd of abundance and deficiency measures
    deficiency_mean <- rowMeans(results[,8, 1:nForests])
    deficiency_sd <- matrixStats::rowSds(results[,8, 1:nForests])
    abundance_mean <- rowMeans(results[,9, 1:nForests])
    abundance_sd <- matrixStats::rowSds(results[,9, 1:nForests])

    # variables ranked by summing frequency across all forests
    frequency_sum <- rowSums(results[,2,1:nForests])
    frequency_rank <- order(order(-frequency_sum))

    # mean and standard deviation of variable importance measures across forests
    pVI_mean <- rowMeans(results[,3,1:nForests])
    pVI_rank <- order(order(-pVI_mean))
    nVI_mean <- rowMeans(results[,4,1:nForests])
    nVI_rank <- order(order(-nVI_mean))
    mdaVI_mean <- rowMeans(results[,5,1:nForests])
    mdaVI_rank <- order(order(-mdaVI_mean))
    mdgVI_mean <- rowMeans(results[,6,1:nForests])
    mdgVI_rank <- order(order(-mdgVI_mean))

    pVI_sd <- matrixStats::rowSds(results[,3,1:nForests])
    nVI_sd <- matrixStats::rowSds(results[,4,1:nForests])
    mdaVI_sd <- matrixStats::rowSds(results[,5,1:nForests])
    mdgVI_sd <- matrixStats::rowSds(results[,6,1:nForests])

    # put calculated statistics into a data frame for plotting by ggplot2
    result.frame <- data.frame(rank_mean, rank_sd, rank_sem, frequency_rank,
                               frequency_mean, frequency_sd, frequency_sem,
                               deficiency_mean, deficiency_sd, abundance_mean, abundance_sd,
                               pVI_mean, nVI_mean, mdaVI_mean, mdgVI_mean, pVI_rank,
                               nVI_rank, mdaVI_rank, mdgVI_rank, pVI_sd, nVI_sd, mdaVI_sd,
                               mdgVI_sd, row.names = rownames(results)[1:480])

    sink(sprintf("%s/text/%s_summary.txt", output_path, model_name))
    cat("index\trank mean\trank sd\tmdaVI mean\tmdaVI sd\tmdgVI mean\tmdgVI sd\n", sep = "")
    for(i in 1:480){
        cat(i, rank_mean[i], rank_sd[i], mdaVI_mean[i], mdaVI_sd[i], mdgVI_mean[i], mdgVI_sd[i], "\n", sep = "\t")
    }
    sink()

    # wilcoxon test of rank distributions
    breaks <- {}
    f_rank <- order(-frequency_mean)
    w <- vector(length=479)
    for(i in 1:479){
        suppressWarnings(r<-wilcox.test(ranks[f_rank[i],], y=ranks[f_rank[i+1],], paired=TRUE))
        w[i] <- r$p.value
        if(w[i]<=.05){
            breaks[length(breaks)+1] <- (rank_mean[f_rank[i]]+rank_mean[f_rank[i+1]])/2
        }
    }

    ## create graphics
    jpeg(filename=sprintf("%s/graphics/%s_rank_wilcoxon_by_index.jpeg", output_path, model_name), width=600, height=600)
    plot(rank_mean[1:479], w, xlab="index", ylab="p-value", main=sprintf("%s\nWilcoxon signed rank test", model_name))
    abline(h=.05, lty="dotted")
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_rank_wilcoxon_by_rank.jpeg", output_path, model_name), width=600, height=600)
    plot(1:50, w[1:50], xlab="rank", ylab="p-value", main=sprintf("%s\nWilcoxon signed rank test", model_name))
    abline(h=.05, lty="dotted")
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_rank_wilcoxon_summary.jpeg", output_path, model_name), width=600, height=600)
    plot(rank_mean, frequency_mean, xlab="rank", ylab="frequency", main=sprintf("%s\nWilcoxon signed rank test", model_name))
    for(b in breaks){
        abline(v=b, lty="dotted")
    }
    invisible(dev.off())

    #wilcoxon test of frequency distribution
    breaks <- {}
    f_rank <- order(-frequency_mean)
    ranked_f <-results[f_rank,2, 1:nForests]
    w <- vector(length=479)
    for(i in 1:479){
        suppressWarnings(r<-wilcox.test(ranked_f[i,], y=ranked_f[i+1,], paired=TRUE, alternative="greater"))
        w[i] <- r$p.value
        if(w[i]<=.05){
            breaks[length(breaks)+1] <- (frequency_mean[f_rank[i]]+frequency_mean[f_rank[i+1]])/2
        }
    }

    jpeg(filename=sprintf("%s/graphics/%s_wilcoxon_by_index.jpeg", output_path, model_name), width=600, height=600)
    plot(f_rank[1:479], w, xlab="index", ylab="p-value", main=sprintf("%s\nWilcoxon signed rank test", model_name))
    abline(h=.05, lty="dotted")
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_wilcoxon_by_rank.jpeg", output_path, model_name), width=600, height=600)
    plot(1:479, w, xlab="rank", ylab="p-value", col=ifelse(w<=.05, "red", "black"), main=sprintf("%s\nWilcoxon signed rank test", model_name))
    abline(h=.05, lty="dotted")
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_wilcoxon_summary.jpeg", output_path, model_name), width=600, height=600)
    plot(1:480, frequency_mean[f_rank], xlab="rank", ylab="frequency", main=sprintf("%s\nWilcoxon signed rank test (p<.05)", model_name))
    for(b in breaks){
        abline(h=b, lty="dotted")
    }
    invisible(dev.off())

    ##

    jpeg(filename=sprintf("%s/graphics/%s_summary_sd.jpeg", output_path, model_name), width=600, height=600)
    plot(rank_mean, rank_sd, xlab="mean rank", ylab="standard deviation", main=sprintf("%s\nvariable rank vs. standard deviation", model_name))
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_summary.jpeg", output_path, model_name), width=600, height=600)
    plot(rank_mean, frequency_mean, xlab="mean rank", ylab="mean frequency", main=sprintf("%s\nvariable rank vs. frequency", model_name))
    invisible(dev.off())

    jpeg(filename=sprintf("%s/graphics/%s_summary_participation.jpeg", output_path, model_name), width=600, height=600)
    plot(frequency_rank, frequency_mean, xlab="mean rank", ylab="mean frequency", main=sprintf("%s\nvariable rank vs. frequency", model_name))
    invisible(dev.off())

    #variable importance comparison
    rank <- order(-result.frame$frequency_mean)
    jpeg(filename=sprintf("%s/graphics/%s_summary_vi_comparison.jpeg", output_path, model_name), width=600, height=600)
    # Basic Scatterplot Matrix
    pairs(~frequency_rank+pVI_rank+nVI_rank+mdgVI_rank+mdaVI_rank,data=result.frame,
          main="Comparison of variable importance metrics")
    invisible(dev.off())

    # Plot showing top 10 properties for averaged forests by mean frequency
    rank <- order(-result.frame$frequency_mean)
    x_label <- rownames(results)[rank]
    jpeg(filename=sprintf("%s/graphics/%s_summary_top_by_frequency.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(data=result.frame, geom="point") +
    aes(x=x_label[1:10], y=frequency_mean[rank[1:10]]) +
    geom_errorbar(aes(ymin=frequency_mean[rank[1:10]]-frequency_sd[rank[1:10]], ymax=frequency_mean[rank[1:10]]+frequency_sd[rank[1:10]]), width=.1) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face="bold"), axis.title.x = element_blank()) +
    ylab("frequency") + xlim(x_label[1:10]) + ggtitle(sprintf("%s\nTop 10 Properties by Average Frequency", model_name))
    suppressWarnings(print(p))
    invisible(dev.off())

    # Plot showing top 10 properties for averaged forests by mean frequency rank
    rank <- order(result.frame$rank_mean)
    x_label <- rownames(results)[rank]
    jpeg(filename=sprintf("%s/graphics/%s_summary_top_by_rank.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(data=result.frame, geom="point") +
    aes(x=x_label[1:10], y=rank_mean[rank[1:10]]) +
    geom_errorbar(aes(ymin=rank_mean[rank[1:10]]-rank_sd[rank[1:10]], ymax=rank_mean[rank[1:10]]+rank_sd[rank[1:10]]), width=.1) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face="bold"), axis.title.x = element_blank()) +
    ylab("mean rank") + xlim(x_label[1:10]) + ggtitle(sprintf("%s\nTop 10 Properties by Average Rank", model_name))
    suppressWarnings(print(p))
    invisible(dev.off())

    # Plot showing top 10 properties for averaged forests by mda VI
    rank <- order(-result.frame$mdaVI_mean)
    x_label <- rownames(results)[rank]
    jpeg(filename=sprintf("%s/graphics/%s_summary_top_by_mdaVI.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(data=result.frame, geom="point") +
    aes(x=x_label[1:10], y=mdaVI_mean[rank[1:10]]) +
    geom_errorbar(aes(ymin=mdaVI_mean[rank[1:10]]-mdaVI_sd[rank[1:10]], ymax=mdaVI_mean[rank[1:10]]+mdaVI_sd[rank[1:10]]), width=.1) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face="bold"), axis.title.x = element_blank()) +
    ylab("MDA VI") + xlim(x_label[1:10]) + ggtitle(sprintf("%s\nTop 10 Properties by mdaVI", model_name))
    suppressWarnings(print(p))
    invisible(dev.off())

    # Plot showing top 10 properties for averaged forests by mdgVI
    rank <- order(-result.frame$mdgVI_mean)
    x_label <- rownames(results)[rank]
    jpeg(filename=sprintf("%s/graphics/%s_summary_top_by_mdgVI.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(data=result.frame, geom="point") +
    aes(x=x_label[1:10], y=mdgVI_mean[rank[1:10]]) +
    geom_errorbar(aes(ymin=mdgVI_mean[rank[1:10]]-mdgVI_sd[rank[1:10]], ymax=mdgVI_mean[rank[1:10]]+mdgVI_sd[rank[1:10]]), width=.1) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, face="bold"), axis.title.x = element_blank()) +
    ylab("MDG VI") + xlim(x_label[1:10]) + ggtitle(sprintf("%s\nTop 10 Properties by mdgVI", model_name))
    suppressWarnings(print(p))
    invisible(dev.off())

    # plots for abundance and deficiency
    jpeg(filename=sprintf("%s/graphics/%s.purity.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(order(order(-frequency_mean)), frequency_mean, data=result.frame, geom="point", xlab="rank", ylab="frequency")
    p<-p+ geom_area(aes(y=deficiency_mean), colour = "red", fill="red") +
    geom_ribbon(aes(ymin=deficiency_mean, ymax=frequency_mean), colour = "green", fill="green")
    suppressWarnings(print(p))
    invisible(dev.off())
    jpeg(filename=sprintf("%s/graphics/%s.top.purity.jpeg", output_path, model_name), width=600, height=600)
    p<-qplot(data=result.frame, geom="point", xlab="property", ylab="frequency") + aes(x=x_label, y=frequency_mean[order(-frequency_mean)]) + xlim(x_label[1:10]) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
    p<-p+geom_area(aes(x=x_label,y=deficiency_mean[order(-frequency_mean)]), colour = "red") +
    geom_ribbon(aes(x=x_label,ymin=deficiency_mean[order(-frequency_mean)], ymax=frequency_mean[order(-frequency_mean)]), colour = "green") + xlab("property") + ylab("frequency")
    suppressWarnings(print(p))
    invisible(dev.off())

    return(results)
}
##################################################
