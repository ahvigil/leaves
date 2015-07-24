context("traceRF tests")

## BEGIN TESTS
test_that("RF traces are reproducible", {
	for(n in c(0, 1)){
		arff.data <- foreign::read.arff("ASP_PROTEASE.4.ASP.OD1_test.arff.gz")
        load(file.path(file.path(getwd()), sprintf("baseline_0%d.model", n)))
		prediction <- analyzeRF::trace.forest(model, arff.data, response=arff.data$class, predict.all=TRUE)
		mdim <- ncol(arff.data)-1 # don't include classification value
        ntest <- nrow(arff.data)
        
        trees <- list(randomForest::getTree(model, 1), randomForest::getTree(model,2))
        for(k in 3:model$ntree){
            trees[[k]] <- randomForest::getTree(model, k)
        }
        
        # positive cases
        L<- which(arff.data$class==1)

        # matrix of tree predictions for positive cases
        P_tree <- prediction$individual[L,]
		
		frequency <- prediction$frequency[,"1"]
        c_frequency <- integer(480)
        r_frequency <- integer(480)
        deficient <- prediction$deficient[,"1"]
        c_deficient <- integer(480)
        r_deficient <- integer(480)
        abundant <- prediction$abundant[,"1"]
        c_abundant <- integer(480)
        r_abundant <- integer(480)
        
        #for each positive case
        for(i in L){
            feature <- arff.data[i,]
            expect_equal(as.numeric(arff.data$class[i]), 2)

            # trees with true positive predictions
            TP_tree <- which(prediction$individual[i,]==1)

            for(k in TP_tree){
                r_result <- analyzeRF::traceTreeR(trees[[k]], feature,
                                      rep(0,nrow(trees[[k]])), 0, 0, r_frequency,
                                      r_abundant, r_deficient)
                r_frequency <- r_result[["frequency"]]
                r_deficient <- r_result[["deficiency"]]
                r_abundant <- r_result[["abundance"]]

                # Same thing but in C
                c_result <- analyzeRF::traceTreeC(trees[[k]], feature, trees[[k]][,3], trees[[k]][,4],
                                      rep(0,nrow(trees[[k]])), 0, 0, c_frequency,
                                      c_abundant, c_deficient)

                c_frequency <- c_result[["frequency"]]
                c_deficient <- c_result[["deficiency"]]
                c_abundant <- c_result[["abundance"]]
            }
        }
        expect_equal(frequency, c_frequency)
        expect_equal(frequency, r_frequency)
        expect_equal(abundant, c_abundant)
        expect_equal(abundant, r_abundant)
        expect_equal(deficient, c_deficient)
        expect_equal(deficient, r_deficient)	
	}
})