#' Plotting clusters after using \code{kmodes} or \code{khaplotype} function.
#'
#' @description
#' Visulize clusters after using the cluster algorithms implemented in the functions.
#' This cluster plot make use of the function \code{dapc} from package \code{adegenet}
#' and function \code{scatter.dapc} from package \code{adegenet}. The plot is based on
#' observation assignment and discriminant analysis of principal components, therefore,
#' when running the code, it will asks you to choose the number of PCs and the number
#' of discriminant functions to retain.
#'
#' @param res Results returned from \code{kmodes} or \code{khaplotype}.
#' @param xax Integer specifying which principal components should be shown in x axes,
#' default is 1.
#' @param yax Integer specifying which principal components should be shown in y axes,
#' default is 2. Notice both xax and yax should be at least smaller than the number of
#' discriminant functions choosed.
#' @param isGene Indicate if the clustering data is gene sequences, default if FALSE.
#'
#' @importFrom adegenet dapc scatter.dapc
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
#'
#' #Plot the clusters by using other principle components.
#' plot_cluster(res_kmodes, xax = 2, yax = 3)
#' }
#'
#' # use function \code{khaplotype}
#' \dontrun{
#' data <- system.file("extdata", "sim_small.fastq", package = "CClust")
#' res_khap <- khaplotype(K = 5, datafile = data, n_init = 10)
#' plot_cluster(res_khap, isGene = TRUE)
#' }

plot_cluster <- function(res, xax = 1, yax = 2, isGene = FALSE) {
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
  adegenet::scatter.dapc(x = dapc, xax = xax, yax = yax, scree.pca = TRUE)
}


