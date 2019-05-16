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
#' @param dat Input data.
#' @param res Results returned from \code{kmodes} or \code{khaplotype}.
#'
#' @importFrom adegenet dapc
#' @importFrom ade4 scatter
#' @importFrom checkmate expect_class
#' @export plot_cluster
#'
#' @return  A plot reflecting clustering results.
#'
#' @examples
#' # use function \code{kmodes}
#' data <- read.table("../data/zoo.int.data")
#' res_kmodes <- kmodes(K = 5, datafile = "../data/zoo.int.data",
#' algorithm = "KMODES_HARTIGAN_WONG", init_method = "KMODES_INIT_AV07_GREEDY",
#' n_init = 10)
#' plot_cluster(data, res_kmodes)
#'

plot_cluster <- function(dat, res) {
    checkmate::expect_class(res, "list")

    dapc <- adegenet::dapc(dat, res$best_cluster_id)
    ade4::scatter(dapc, posi.da = "bottomright", bg = "white",
            cstar = 0, scree.pca = TRUE, posi.pca = "bottomleft")
}
