#' Compute the adjusted rand index
#'
#' @description Compute the adjusted rand index given the estimated assignments and the true assiganments.
#'
#' @param est A list of results returned from \code{khaplotype} or \code{kmodes}.
#' @param truth A numeric or character vector of true assignemts.
#'
#' @importFrom mclust adjustedRandIndex
#' @importFrom checkmate expect_class checkClass
#' @importFrom Rdpack reprompt
#' @export ARI
#'
#' @references {
#' \itemize{
#'     \item \insertRef{Hubert1985}{CClust}
#'     }
#' }
#'
#' @return A numberic value between 0 and 1, which indicates agreement between two partitions.
#'
#' @examples
#' # Estimate cluster assignments by function `khaplotype`.
#' datFile <- system.file("inst/extdata/sim.fastq", package = "CClust")
#' res_khap <- khaplotype (K = 5, datafile = datFile, n_init = 3)
#' ass <- system.file("inst/extdata/assignment.txt", package = "CClust")
#' true_assignments <- as.numeric(read.table(ass, header = F, sep = ""))
#' # Compute the Adjusted Rand Index
#' ARI(res_khap, true_assignments)

ARI <- function(est, truth)
{
  checkmate::expect_class(est, "list")
  checkmate::expect_class(est$best_cluster_id, "integer")
  checkmate::assert(
    checkmate::checkClass(truth, "vector"),
    checkmate::checkClass(truth, "numeric"),
    checkmate::checkClass(truth, "integer")
  )
  mclust::adjustedRandIndex(est$best_cluster_id, truth)
}
