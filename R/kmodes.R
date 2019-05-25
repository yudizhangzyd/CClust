#' Clustering categorical datasets.
#'
#' @description
#' Implement three unsupervised clustering algorithms on categorical datasets.
#'
#' @usage kmodes(K = 1, datafile = NULL, n_init = 1, algorithm = "KMODES_HUANG",
#' init_method = "KMODES_INIT_RANDOM_SEEDS", seed = 1, shuffle = FALSE)
#'
#' @param K Number of clusters. Default is 1.
#' @param datafile Path to a data file.
#' @param n_init Number of initializations.
#' @param algorithm Algorithm to implement clustering. Default is "KMODES_HUANG". See details for the options available.
#' @param init_method Initialization methods. Default is "KMODES_INIT_RANDOM_SEEDS". See details for the options available.
#' @param seed Random number seed. Default is 1.
#' @param shuffle Incidate if shuffle the input order. Default is FALSE.
#'
#' @details
#' Algorithms avaiable:
#'
#' \itemize{
#'    \item \code{"KMODES_HUANG"}: MacQueen's algorithm
#'    \item \code{"KMODES_HARTIGAN_WONG"}: Hartigan and Wong algorithm
#'    \item \code{"KMODES_LLOYD"}: Lloyd's algorithm
#'  }
#'
#' Initialization methods avaiable:
#' \itemize{
#'    \item \code{"KMODES_INIT_RANDOM_SEEDS"}: Random sampling.
#'    \item \code{"KMODES_INIT_H97_RANDOM"}: Huang1997, randomized version.
#'    \item \code{"KMODES_INIT_HD17"}: Huang1997 interpretted by Python author de Vos.
#'    \item \code{"KMODES_INIT_CLB09_RANDOM"}: Cao2009, randomized version.
#'    \item \code{"KMODES_INIT_AV07"}: K-means++ adapted.
#'    \item \code{"KMODES_INIT_AV07_GREEDY"}: K-means++ greedy adapted.
#'  }
#'
#' @return Returns a list of clustering results.
#' @details
#' Value:
#' \itemize{
#'    \item \code{"best_cluster_size"}: Number of observations in each cluster of the best initialization.
#'    \item \code{"best_criterion"}: Optimized criterion in each cluster of the best initialization.
#'    \item \code{"best_cluster_id"}: Cluster assignment of each observation of the best initialization.
#'    \item \code{"best_modes"}: Estimated modes for each cluster of the best initialization.
#'    \item \code{"best_seed_index"}: Seed index of the best initialization.
#'    \item \code{"total_best_criterion"}: Total optimized criterion of the best initialization.
#'    \item \code{"clsuter_size"}: Number of clusters.
#'    \item \code{"data_dim"}: Dimension of input data.
#' }
#'
#' @useDynLib CClust r_kmodes
#' @importFrom checkmate expect_file_exists expect_string expect_choice
#' @importFrom Rdpack reprompt
#' @export kmodes
#'
#' @references {
#' \itemize{
#'     \item \insertRef{Lloyd1982}{CClust}
#'     \item \insertRef{MacQueen1967}{CClust}
#'     \item \insertRef{Huang1998}{CClust}
#'     \item \insertRef{Huang1997}{CClust}
#'     \item \insertRef{Hartigan1975}{CClust}
#'     }
#' }
#'
#' @examples
#' # Clustering with three initializations with default algorithm ("KMODES_HUANG")
#' datFile <- system.file("extdata", "zoo.int.data", package = "CClust")
#' res_kmodes <- kmodes(K = 5, datafile = datFile, n_init = 3, shuffle = TRUE)
#'
#' # Clustering with Harigan and Wong and K-means++ greedy adapted initialization method.
#' res_kmodes <- kmodes(K = 5, datafile = datFile,
#' algorithm = "KMODES_HARTIGAN_WONG", init_method = "KMODES_INIT_AV07_GREEDY")
#'

kmodes <- function(K = 1,
                   datafile = NULL,
                   n_init = 1,
                   algorithm = "KMODES_HUANG",
                   init_method = "KMODES_INIT_RANDOM_SEEDS",
                   seed = 1,
                   shuffle = FALSE)
{

  #Error checks
  if (K == 0)
    stop("The number of cluster K should be more than 0!")

  if (n_init == 0)
    stop ("The number of initialization should be more than 0!")

  checkmate::expect_file_exists(datafile, access = "r")
  checkmate::expect_string(algorithm, info = "The input algorithm must be a string.")
  checkmate::expect_string(init_method, info = "The input initialization method must be a string.")
  checkmate::expect_choice(algorithm, choices = c("KMODES_HUANG", "KMODES_LLOYD", "KMODES_HARTIGAN_WONG"),
                           info = "The algorithm specificed is not avaiable in kmodes clustering.")
  checkmate::expect_choice(init_method, choices = c("KMODES_INIT_RANDOM_SEEDS", "KMODES_INIT_H97_RANDOM", "KMODES_INIT_HD17",
                                                    "KMODES_INIT_CLB09_RANDOM", "KMODES_INIT_AV07", "KMODES_INIT_AV07_GREEDY"),
                           info = "The initialization method specificed is not avaiable in kmodes clustering.")

  if (!is.loaded("r_kmodes", PACKAGE = "CClust")) {
    dyn.load("../src/CClust.so")
  }

  K <- as.integer(K)
  n_init <- as.integer(n_init)
  seed <- as.integer(seed)
  shuffle <- as.integer(shuffle)

  if (algorithm == "KMODES_HUANG") {
    algorithm = as.integer(0)
  } else if (algorithm == "KMODES_LLOYD") {
    algorithm = as.integer(1)
  } else if (algorithm == "KMODES_HARTIGAN_WONG") {
    algorithm = as.integer(2)
  }

  if(init_method == "KMODES_INIT_RANDOM_SEEDS") {
    init_method = as.integer(0)
  } else if (init_method == "KMODES_INIT_H97_RANDOM") {
    init_method = as.integer(4)
  } else if (init_method == "KMODES_INIT_HD17") {
    init_method = as.integer(5)
  } else if (init_method == "KMODES_INIT_CLB09_RANDOM") {
    init_method = as.integer(6)
  } else if (init_method == "KMODES_INIT_AV07") {
    init_method = as.integer(7)
  } else if (init_method == "KMODES_INIT_AV07_GREEDY") {
    init_method = as.integer(8)
  } else if (init_method == "KMODES_INIT_H97") {
    init_method = as.integer(10)
  } else if (init_method == "KMODES_INIT_CLB09") {
    init_method = as.integer(11)
  }

  res <- .Call("r_kmodes", K , datafile, n_init, algorithm,
               init_method, seed, shuffle)
  names(res) <- c("best_cluster_size", "best_seed_idx",
                  "best_criterion", "best_cluster_id", "best_modes")
  res[[5]] <- matrix(res[[5]], nrow = K, byrow = TRUE)
  res$total_best_criterion <- sum(res$best_criterion)
  res$cluster_size <- K
  res$data_dim <- c(length(res[[4]]), dim(res[[5]])[2])

  return(res)
} # kmodes

