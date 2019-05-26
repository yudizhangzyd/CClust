#' Clustering amplicon datasets with quality scores.
#'
#' @description
#' Implement three unsupervised clustering algorithms on amplicon datasets with quality scores.
#'
#' @usage khaplotype(K = 1, datafile = NULL, n_init = 1, algorithm = "FASTQ_HW_EFFICIENT",
#' seed = 0, shuffle = FALSE)
#'
#' @param K Number of clusters. Default is 1.
#' @param datafile Path to a data file. Has to be a fastq file if want to conduct clustering on amplicon data.
#' @param n_init Number of initializations.
#' @param algorithm Algorithm to implement clustering. Default is "FASTQ_LLOYDS_EFFICIENT". See details for the options available.
#' @param seed Random number seed. Default is 1.
#' @param shuffle Incidate if shuffle the input order. Default is FALSE.
#'
#' @details
#' Algorithms avaiable:
#'
#' \itemize{
#'    \item \code{"FASTQ_LLOYDS_EFFICIENT"}: Efficient Lloyds algorithm
#'    \item \code{"FASTQ_HW_EFFICIENT"}: Efficient Hartigan and Wong algorithm
#'    \item \code{"FASTQ_MACQUEEN"}: MacQueen's algorithm
#'    \item \code{"FASTQ_LLOYDS"}: Lloyds algorithm
#'    \item \code{"FASTQ_HW"}: Hartigan and Wong algorithm
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
#'    \item \code{"total_best_criterion"}: Total optimized criterion of the best initialization.
#'    \item \code{"clsuter_size"}: Number of clusters.
#'    \item \code{"data_dim"}: Dimension of input data.
#'    \item \code{"data"}: Reads of the input data.
#' }
#'
#' @useDynLib CClust r_khaplotype
#' @importFrom checkmate expect_file_exists expect_string expect_choice
#' @importFrom Rdpack reprompt
#' @importFrom utils tail
#' @export khaplotype
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
#' # Clustering an amplicon dataset and run three initializations with default
#' # algorithm ("FASTQ_HW_EFFICIENT")
#' datFile <- system.file("extdata", "sim.fastq", package = "CClust")
#' res_khap <- khaplotype(K = 5, datafile = datFile, n_init = 3)
#'
#' # Clustering an amplicon dataset and run three initializations with
#' # MacQueen's algorithm (shuffle the data)
#' res_khap <- khaplotype(K = 5, datafile = datFile, n_init = 3,
#' algorithm = "FASTQ_MACQUEEN", shuffle = TRUE)
#'
#' # Clustering an amplicon dataset provide a different seed
#' res_khap <- khaplotype(K = 5, datafile = datFile, seed = 1)
#'

khaplotype <- function(K = 1,
                       datafile = NULL,
                       n_init = 1,
                       algorithm = "FASTQ_HW_EFFICIENT",
                       seed = 0,
                       shuffle = FALSE)
{

  #Error checks
  checkmate::expect_file_exists(datafile, access = "r")

  if (K == 0)
    stop ("The number of cluster K should be more than 0!")

  if (n_init == 0)
    stop ("The number of initialization should be more than 0!")

  if (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fastq")
    stop ("The input datafile has to be fastq file!")

  checkmate::expect_string(algorithm, info = "The input algorithm must be a string.")
  checkmate::expect_choice(algorithm, choices = c("FASTQ_LLOYDS_EFFICIENT", "FASTQ_LLOYDS", "FASTQ_HW_EFFICIENT",
                                                  "FASTQ_MACQUEEN", "FASTQ_HW"),
                           info = "The algorithm specificed is not avaiable in khaplotype clustering.")

  if (!is.loaded("r_khaplotype", PACKAGE = "CClust")) {
    dyn.load("../src/CClust.so")
  }

  K <- as.integer(K)
  n_init <- as.integer(n_init)
  seed <- as.integer(seed)
  shuffle <- as.integer(shuffle)
  run_with_quals <- as.integer(1)

  ## Varieties of algorithm for fastq data.
  if(algorithm == "FASTQ_LLOYDS") {
    algorithm = as.integer(0)
  } else if (algorithm == "FASTQ_LLOYDS_EFFICIENT") {
    algorithm = as.integer(1)
  } else if (algorithm == "FASTQ_HW_EFFICIENT") {
    algorithm = as.integer(2)
  } else if (algorithm == "FASTQ_MACQUEEN") {
    algorithm = as.integer(3)
  } else if (algorithm == "FASTQ_HW") {
    algorithm = as.integer(4)
  }

  res <- .Call("r_khaplotype", K, datafile, n_init, algorithm,
                seed, shuffle, run_with_quals)
  res[[4]] <- matrix(res[[4]], nrow = K, byrow = TRUE)
  res[[5]] <- matrix(res[[5]], nrow = length(res[[3]]), byrow = TRUE)

  names(res) <- c("best_cluster_size", "best_criterion",
                  "best_cluster_id", "best_modes", "data")

  res$total_best_criterion <- sum(res$best_criterion)
  res$cluster_size <- K
  res$data_dim <- c(length(res[[3]]), dim(res[[4]])[2])

  return(res)
} # khaplotype

