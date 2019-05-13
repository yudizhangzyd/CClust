context("test-read_fastq")

test_that("check-read_fastq", {
  expect_error(read_fastq(datafile = "sim.fasta"))
  expect_error(read_fastq(datafile = "../sim.fastq"))
})
