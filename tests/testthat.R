library(testthat)
library(miso)

if (requireNamespace("xml2")) {
    test_check(
        "miso",
        reporter = MultiReporter$new(
            reporters = list(
                JunitReporter$new(file = "test-results.xml"),
                CheckReporter$new()))
    )
} else {
    test_check("miso")
}
