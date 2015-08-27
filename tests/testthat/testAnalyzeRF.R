context("leaves tests")

## BEGIN TESTS
test_that("Tests run from scratch produce correct results", {
    test.model <- {}
    test.model$name <- "ASP_PROTEASE.4.ASP.OD1_test"
    test.model$ntree      <- 501
    test.model$mtry 	  <- 50
    test.model$seed_mult <- 43
    test.model$forests <- 2

    output_path <- file.path(getwd(), "output")
    rf_data_path <- file.path(getwd())
    feature_data_path <- file.path(getwd())

    invisible(leaves::analyzeModel(test.model, nForests, rf_data_path=file.path(getwd()),
                            feature_data_path=file.path(getwd()),
                            output_path = output_path) )

    for(n in c(0,1)){
        # random forest model
        load(file.path(file.path(getwd()), sprintf("baseline_0%d.model", n)))
        baseline <- model
        rm(model)
        load(file.path(file.path(getwd()), sprintf("ASP_PROTEASE.4.ASP.OD1_test_n%04d_m%02d_%04d.model", test.model$ntree, test.model$mtry, n*test.model$seed_mult)))
        sprintf("ASP_PROTEASE.4.ASP.OD1_test_n%04d_m%02d_%04d.model", test.model$ntree, test.model$mtry, n*test.model$seed_mult)

        # test results should be equivalent for these properties
        properties <- c("type","predicted","err.rate","confusion","votes",
                        "oob.times","classes","importance","ntree","mtry","forest","inbag")
        for(property in properties){
            expect_equal(model[property], baseline[property])
        }

        # saved results
        load(file.path(file.path(getwd()), sprintf("baseline_0%d.save", n)))
        baseline <- O
        rm(O)
        load(file.path(file.path(getwd()), "output", "binary", sprintf("%s_n%04d_m%02d_%04d.save", test.model$name, test.model$ntree, test.model$mtry, n*test.model$seed_mult)))
        expect_equal(O, baseline)
    }

    # summary
    load(file.path(file.path(getwd()), "baseline_02.save"))
    baseline <- results
    rm(results)
    load(file.path(file.path(getwd()), "output", "binary", "ASP_PROTEASE.4.ASP.OD1_test_n0501_m50.save"))
    expect_equal(results, baseline)

    unlink(c(file.path(output_path),
             file.path(getwd(), "ASP_PROTEASE.4.ASP.OD1_test_n0501_m50_0000.model"),
             file.path(getwd(), "ASP_PROTEASE.4.ASP.OD1_test_n0501_m50_0043.model")),
           recursive=TRUE)
})