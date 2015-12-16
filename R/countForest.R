#' Add frequency, abundance, and deficiency information to a forest object's importance data
#' @export
countForest <- function(x, data, y_name="class") {
  classes <- levels(data[,y_name])
  results <- leaves::trace.forest(x, data)
  
  for(cl in classes){
    x$importance <- cbind(x$importance, results$frequency[,cl])
    colnames(x$importance)[ncol(x$importance)] <- sprintf("%s_frequency", cl)
    x$importance <- cbind(x$importance, results$abundant[,cl])
    colnames(x$importance)[ncol(x$importance)] <- sprintf("%s_abundant", cl)
    x$importance <- cbind(x$importance, results$deficient[,cl])
    colnames(x$importance)[ncol(x$importance)] <- sprintf("%s_deficient", cl)
  }
  
  x
}