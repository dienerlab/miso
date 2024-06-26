context("EM counting")

flog.threshold(WARN)

test_that("configuration is good", {
    config <- config_count()
    expect_type(config, "list")
    expect_s3_class(config, "config")
    expect_true(length(config) > 3)
    config <- config_count(method = "naive")
    expect_equal(config$method, "naive")
})

alns <- system.file("extdata/aln", package = "miso") %>%
      find_alignments()
ref <- system.file("extdata/genomes/zymo_mock.fna.gz",
                   package = "miso")

test_that("EM counting works", {
    cn <- count_references(alns, reference = ref)

    expect_named(cn, c("alignments", "counts", "steps"))
    expect_s3_class(cn$counts, "data.table")
    expect_true(cn$counts[, all(counts > 0)])
    expect_equal(cn$counts[, uniqueN(reference, na.rm = TRUE)], 10)
})

test_that("weighted EM counting works", {
    cn <- count_references(alns, reference = ref, weights = TRUE)

    expect_named(cn, c("alignments", "counts", "steps"))
    expect_s3_class(cn$counts, "data.table")
    expect_true(cn$counts[, all(counts > 0)])
    expect_equal(cn$counts[, uniqueN(reference, na.rm = TRUE)], 10)
})

test_that("naive counting works", {
    cn <- count_references(alns, method = "naive")

    expect_named(cn, c("alignments", "counts", "steps"))
    expect_s3_class(cn$counts, "data.table")
    expect_true(cn$counts[, all(counts > 0)])
    expect_equal(cn$counts[, uniqueN(reference, na.rm = TRUE)], 10)
})

flog.threshold(INFO)
