context("test-ARI")

# Read the true assignments
true_assignments <- as.numeric(read.table("../data/assignment.txt", header = F, sep = ""))

# Clustering
res_khap <- khaplotype (K = 5, datafile = "../data/sim.fastq", n_init = 3)

test_that("check-ARI", {
  expect_error(ARI(est = rep(1:10), truth = true_assignments))
  expect_error(ARI(est = res_khap, truth = "A"))
})
