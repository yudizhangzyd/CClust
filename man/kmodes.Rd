% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmodes.R
\name{kmodes}
\alias{kmodes}
\title{Clustering categorical datasets.}
\usage{
kmodes(K = 1, datafile = NULL, n_init = 1, algorithm = "KMODES_HUANG",
init_method = "KMODES_INIT_RANDOM_SEEDS", seed = 1, shuffle = FALSE)
}
\arguments{
\item{K}{Number of clusters. Default is 1.}

\item{datafile}{Path to a data file.}

\item{n_init}{Number of initializations.}

\item{algorithm}{Algorithm to implement clustering. Default is "KMODES_HUANG". See details for the options available.}

\item{init_method}{Initialization methods. Default is "KMODES_INIT_RANDOM_SEEDS". See details for the options available.}

\item{seed}{Random number seed. Default is 1.}

\item{shuffle}{Incidate if shuffle the input order. Default is FALSE.}
}
\value{
Returns a list of clustering results.
}
\description{
Implement three unsupervised clustering algorithms on categorical datasets.
}
\details{
Algorithms avaiable:

\itemize{
   \item \code{"KMODES_HUANG"}: MacQueen's algorithm
   \item \code{"KMODES_HARTIGAN_WONG"}: Hartigan and Wong algorithm
   \item \code{"KMODES_LLOYD"}: Lloyd's algorithm
 }

Initialization methods avaiable:
\itemize{
   \item \code{"KMODES_INIT_RANDOM_SEEDS"}: Random sampling.
   \item \code{"KMODES_INIT_H97_RANDOM"}: Huang1997, randomized version.
   \item \code{"KMODES_INIT_HD17"}: Huang1997 interpretted by Python author de Vos.
   \item \code{"KMODES_INIT_CLB09_RANDOM"}: Cao2009, randomized version.
   \item \code{"KMODES_INIT_AV07"}: K-means++ adapted.
   \item \code{"KMODES_INIT_AV07_GREEDY"}: K-means++ greedy adapted.
 }

Value:
\itemize{
   \item \code{"best_cluster_size"}: Number of observations in each cluster of the best initialization.
   \item \code{"best_criterion"}: Optimized criterion in each cluster of the best initialization.
   \item \code{"best_cluster_id"}: Cluster assignment of each observation of the best initialization.
   \item \code{"best_modes"}: Estimated modes for each cluster of the best initialization.
   \item \code{"best_seed_index"}: Seed index of the best initialization.
   \item \code{"total_best_criterion"}: Total optimized criterion of the best initialization.
   \item \code{"clsuter_size"}: Number of clusters.
   \item \code{"data_dim"}: Dimension of input data.
   \item \code{"data"}: The input data.
}
}
\examples{
# Clustering with three initializations with default algorithm ("KMODES_HUANG")
datFile <- system.file("extdata", "zoo.int.data", package = "CClust")
res_kmodes <- kmodes(K = 5, datafile = datFile, n_init = 3, shuffle = TRUE)

# Clustering with Harigan and Wong and K-means++ greedy adapted initialization method.
res_kmodes <- kmodes(K = 5, datafile = datFile,
algorithm = "KMODES_HARTIGAN_WONG", init_method = "KMODES_INIT_AV07_GREEDY")

}
\references{
{
\itemize{
    \item \insertRef{Lloyd1982}{CClust}
    \item \insertRef{MacQueen1967}{CClust}
    \item \insertRef{Huang1998}{CClust}
    \item \insertRef{Huang1997}{CClust}
    \item \insertRef{Hartigan1975}{CClust}
    }
}
}
