context("test-read_fastq(a)")

test_that("check-read_fastq", {
  expect_error(read_fastq(datafile = "sim.fasta"))
  expect_error(read_fastq(datafile = "../sim.fastq"))
})

test_that("check-read_fasta", {
  expect_error(read_fasta(datafile = "sim.fastq"))
  expect_error(read_fasta(datafile = "../sim.fasta"))
})
