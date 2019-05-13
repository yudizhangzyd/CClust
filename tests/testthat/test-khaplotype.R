context("test-khaplotype")

# Name a fastq file
file <- "../data/sim.fastq"

# test error messages
test_that("check-khaplotype", {
  expect_error(khaplotype(K = 0, datafile = file))
  expect_error(khaplotype(n_init = 0, datafile = file))
  expect_error(khaplotype(K = 3, datafile = file, algorithm = "KMODES"))
  expect_error(khaplotype(K = 3, datafile = "sim.fasta"))
})
