#' Simulate clustered categorical datasets.
#'
#' @description
#' Simulate clustered categorical datasets by using continous time Marcov chain.
#'
#' @usage simulator(simK = 5, n_coordinates = 10, n_observations = 100,
#' n_categories = 4, sim_between_t = 2, sim_within_t = FALSE,
#' use_dirichlet = 0, sim_pi = c(0.1, 0.1, 0.2, 0.3, 0.3))
#'
#' @param simK Number of clusters.
#' @param n_coordinates Number of coordinate of the simulated dataset.
#' @param n_observations Number of observation of the simulated dataset.
#' @param sim_between_t Between cluster variation.
#' @param sim_within_t Within cluster variation.
#' @param use_dirichlet Indicate if cimulate datasets with dirichlet prior. Defalut is FALSE.
#' @param sim_pi Mixing proportions, a vector with the same length of specified number of clusters
#' and the sum of the values in this vector has to be 1.
#'
#' @return Returns a list of simulation dataset results.
#' @details
#' Value:
#' \itemize{
#'    \item \code{"CTMC_probabilities"}: Number of observations in each cluster of the best initialization.
#'    \item \code{"modes"}: Simulated modes.
#'    \item \code{"cluster_assignments"}: Simulated cluster assignments.
#'    \item \code{"cluster_sizes"}: Simulated cluster sizes.
#'    \item \code{"data"}: Simulated data.
#' }
#'
#' @useDynLib CClust r_simulate_data
#' @importFrom checkmate expect_class
#' @export simulator
#'
#' @examples
#' \dontrun{
#' #Simulate data with dim 100 * 10, 4 different categories and there are 5 true clusters.
#' data <- simulator(simK = 5, n_coordinates = 10, n_observations = 100, n_categories = 4,
#' sim_between_t = 2, sim_within_t = 1, use_dirichlet = TRUE, sim_pi = c(0.1, 0.1, 0.2, 0.3, 0.3))
#' }

simulator <- function(simK,
                      n_coordinates,
                      n_observations,
                      n_categories,
                      sim_between_t,
                      sim_within_t,
                      use_dirichlet = FALSE,
                      sim_pi)
{
  if (is.vector(sim_pi) == FALSE)
    stop ("The input sim_pi has to be a vector.")
  if (sum(sim_pi) != 1)
    stop ("The sum of sim_pi has to be 1.")
  checkmate::expect_class(sim_pi, "numeric")

  if (!is.loaded("r_simulate_data", PACKAGE = "CClust")) {
    dyn.load(dir(system.file("libs", package = "CClust")))
  }

  simK <- as.integer(simK)
  n_coordinates <- as.integer(n_coordinates)
  n_observations <- as.integer(n_observations)
  n_categories <- as.integer(n_categories)
  use_dirichlet <- as.integer(use_dirichlet)
  sim_between_t <- as.double(sim_between_t)
  sim_within_t <- as.double(sim_within_t)

  res <- .Call(getNativeSymbolInfo("r_simulate_data"), simK,
               n_coordinates, n_observations,
               n_categories, sim_between_t, sim_within_t,
               use_dirichlet, sim_pi)

  names(res) <- c("CTMC_probabilities", "modes", "cluster_assignments",
                  "cluster_sizes", "data")

  res$modes <- matrix(res$modes, ncol = n_coordinates, byrow = TRUE)
  res$data <- matrix(res$data, ncol = n_coordinates, byrow = TRUE)
  return(res)
}


