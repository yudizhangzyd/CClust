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
#' @param dat Input data path.
#' @param res Results returned from \code{kmodes} or \code{khaplotype}.
#' @param isFastq Indicate if the input data is in fastq format, if TRUE,
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
#' data <- system.file("inst/extdata/zoo.int.data", package = "CClust")
#' res_kmodes <- kmodes(K = 5, datafile = data, algorithm = "KMODES_HARTIGAN_WONG",
#' init_method = "KMODES_INIT_AV07_GREEDY", n_init = 10)
#' plot_cluster(data, res_kmodes)
#'
#' # use function \code{khaplotype}
#' data <- system.file("inst/extdata/sim.fastq", package = "CClust")
#' res_khap <- khaplotype(K = 5, datafile = data, n_init = 10)
#' plot_cluster(data, res_khap, isFastq = TRUE)

plot_cluster <- function(dat, res, isFastq = FALSE) {
  checkmate::expect_class(res, "list")
  checkmate::expect_class(dat, "character")

  if(isFastq == FALSE) {
    if (tail(unlist(strsplit(dat, "[.]")), 1) != "fastq")
      stop ("The input datafile has to be fastq file when isFastq = TRUE.")
    dat <- read.table(dat)
  }
  else {
    if (tail(unlist(strsplit(dat, "[.]")), 1) != "fastq")
      stop ("The input datafile has to be fastq file.")
    data <- read_fastq(dat)
    reads <- data$reads
    dat <- matrix(0, nrow = data$dim[1], ncol = 4 * data$dim[2])
    for (i in 1:data$dim[1]) {
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


