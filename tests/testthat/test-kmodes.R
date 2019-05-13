context("test-kmdoes")

# Name a file
file <- "../data/zoo.int.data"

# test error messages
test_that("check-khaplotype", {
  expect_error(kmodes(K = 0, datafile = file))
  expect_error(kmodes(n_init = 0, datafile = file))
  expect_error(kmodes(K = 3, datafile = file, algorithm = "KMODES"))
  expect_error(kmodes(K = 3, datafile = "sim.fasta"))
})
