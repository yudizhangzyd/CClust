#' Plotting clusters after using \code{kmodes} or \code{khaplotype} function.
#'
#' @description
#' Visulize clusters after using the cluster algorithms implemented in the functions.
#' This cluster plot make use of the function \code{dapc} from package \code{adegenet}
#' and function \code{scatter} from package \code{ade4}. The plot is based on observation
#' assignment and discriminant analysis of principal components, therefore, when running
#' the code, it will asks you to choose the number PCs to retain and choose the number
#' discriminant functions to retain.
#'
#' @param res Results returned from \code{kmodes} or \code{khaplotype}.
#' @param isGene Indicate if the clustering data is gene sequences, default if FALSE.
#'
#' @importFrom adegenet dapc
#' @importFrom ade4 scatter
#' @importFrom checkmate expect_class
#' @importFrom stringr str_split
#' @export plot_cluster
#'
#' @return  A plot reflecting clustering results.
#'
#' @examples
#' # use function \code{kmodes}
#' \dontrun{
#' data <- system.file("extdata", "zoo.int.data", package = "CClust")
#' res_kmodes <- kmodes(K = 5, datafile = data, algorithm = "KMODES_HARTIGAN_WONG",
#' init_method = "KMODES_INIT_AV07_GREEDY", n_init = 10)
#' plot_cluster(res_kmodes)
#' }
#'
#' # use function \code{khaplotype}
#' \dontrun{
#' data <- system.file("extdata", "sim_small.fastq", package = "CClust")
#' res_khap <- khaplotype(K = 5, datafile = data, n_init = 10)
#' plot_cluster(res_khap, isGene = TRUE)
#' }

plot_cluster <- function(res, isGene = FALSE) {
  checkmate::expect_class(res, "list")

  if(isGene == FALSE) {
    dat <- res$data
  }
  else {
    if (is.null(res$best_seed_idx) == FALSE)
      stop ("The clustering data has to be gene sequences.")
    reads <- res$data
    dat <- matrix(0, nrow = res$data_dim[1], ncol = 4 * res$data_dim[2])
    for (i in 1:res$data_dim[1]) {
      reads[i, ][reads[i, ] == "A"] = c('1000')
      reads[i, ][reads[i, ] == "C"] = c('0100')
      reads[i, ][reads[i, ] == "G"] = c('0010')
      reads[i, ][reads[i, ] == "T"] = c('0001')
      dat[i, ] <- as.numeric(unlist(stringr::str_split(reads[i, ], "")))
    }
  }

  dapc <- adegenet::dapc(dat, res$best_cluster_id)
  ade4::scatter(dapc, posi.da = "bottomright", bg = "white",
            cstar = 0, scree.pca = TRUE, posi.pca = "bottomleft")
}


