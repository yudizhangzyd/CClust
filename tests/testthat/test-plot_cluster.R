context("test-plot_cluster")

# Read the data and do the clustering
data <- read.table("../../data/zoo.int.data")
res_kmodes <- kmodes(K = 5, datafile = "../data/zoo.int.data",
                     algorithm = "KMODES_HARTIGAN_WONG", init_method = "KMODES_INIT_AV07_GREEDY",
                     n_init = 10)

test_that("check-plot_cluster", {
  expect_error(plot_cluster(data, res_kmodes$best_cluster_id))
})
