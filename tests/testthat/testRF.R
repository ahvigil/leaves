context("RF tests")

## BEGIN TESTS
test_that("RF produces results identical to randomForest module", {
    set.seed(0)
    model <- leaves::rf(Species~., iris, importance=T)
    
    for(n in seq(5,model$ntree, 5)){
        set.seed(0)
        reference <- randomForest::randomForest(Species~., iris, importance=T, ntree=n)
        
        expect_equal(model$impmat[,,n], reference$importance)
    }
    
    expect_equal(model$importance, reference$importance)
    
    expect_equal(model$forest, reference$forest)
})