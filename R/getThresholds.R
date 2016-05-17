#' Get variable threshold values for a forest
#' @export
getThresholds <- function(rf) {
  thresholds <- list()

  features <- rownames(rf$importance)

  i <- 1
  for(feature in features){
    thresholds[[feature]] <- rf$forest$xbestsplit[which(rf$forest$bestvar==i)]
    i<-i+1
  }

  thresholds
}

#' Calculate mean and standard deviation of thresholds in an RF model
#' @export
getThresholdStats <- function(rf) {
  thresholds <- list()

  features <- rownames(rf$importance)

  i <- 1
  for(feature in features){
    thresholds[[feature]] <- rf$forest$xbestsplit[which(rf$forest$bestvar==i)]
    i<-i+1
  }

  averages <- NULL
  sds <- NULL
  for(feature in features){
    averages <- rbind(averages, mean(thresholds[[feature]]))
    sds <- rbind(sds, sd(thresholds[[feature]]))
  }

  out <- cbind(averages, sds)

  colnames(out) <- c("mean", "SD")
  rownames(out) <- rownames(rf$importance)

  out
}

#' Retrieve the list of thresholds used to make a prediction on a set of observations
#' @export
getPredictionThresholds <- function(rf, newData, y_name="class") {
  classes <- levels(newData[,y_name])
  response <- as.numeric(newData[,y_name])

  n <- nrow(newData)
  mdim <- nrow(rf$importance)

  forest <- rf$forest

  prediction <- integer(n)
  pred <- integer(rf$ntree)

  abundant_thresholds <- vector("list", mdim)
  deficient_thresholds <- vector("list", mdim)

  # predict ith case
  for(i in 1:n){
    reference <- response[i]

    # nth tree prediction
    for(n in 1:rf$ntree){
      local_abundant_thresholds <- vector("list", mdim)
      local_deficient_thresholds <- vector("list", mdim)

      # start at node 1
      k <- 1

      # while node status is not terminal
      while(forest$nodestatus[k,n] != -1){
        # split on variable 'm'
        m <- forest$bestvar[k,n]

        threshold <- forest$xbestsplit[k,n]

        #deficient
        if(newData[i,m] <= threshold){
          local_deficient_thresholds[[m]] <- c(local_deficient_thresholds[[m]], threshold)
          k <- forest$treemap[k,1,n]
        #abundant
        } else{
          local_abundant_thresholds[[m]] <- c(local_abundant_thresholds[[m]], threshold)
          k <- forest$treemap[k,2,n]
        }
      }

      # record prediction
      pred[n] <- forest$nodepred[k,n]
      if(pred[n] == 0) stop("Prediction ended on an invalid node (prediction==0)")

      if(pred[n] == reference){
        for(j in 1:mdim){
          if(!is.null(local_deficient_thresholds[[j]])){
            deficient_thresholds[[j]] <- c(deficient_thresholds[[j]], local_deficient_thresholds[[j]])
          }
          if(!is.null(local_abundant_thresholds[[j]])){
            abundant_thresholds[[j]] <- c(abundant_thresholds[[j]], local_abundant_thresholds[[j]])
          }
        }
      }
    }
    prediction[i] <- as.integer(names(sort(-table(pred)))[1])

  }
  out <- list(abundant=abundant_thresholds,
              deficient=deficient_thresholds)

  out
}
