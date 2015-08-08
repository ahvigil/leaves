# extract importance from multiple saved forests and return as 3d array

aggregateImportance <- function(names){

	importance <- array(0, dim=c(480, 4, length(names)))

	for(i in 1:length(names)){
		name <- sub(".model", "", names[i])
		load(sprintf("%s.model", name))
		importance[,,i] <- model$importance
		if(i==1) dimnames(importance) <- c(dimnames(model$importance), NULL)
	}

	importance
}

# select top n features using Mean Decrease in Gini
topGini <- function(importance, n){
	i<-rowMeans(importance, dims=2)
	order(-i[,"MeanDecreaseGini"])[1:n]
}

# select top n features using Mean Decrease in Accuracy
topAccuracy <- function(importance, n){
    i<-rowMeans(importance, dims=2)
    order(-i[,"MeanDecreaseAccuracy"])[1:n]
}
