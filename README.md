`CClust` is an R package developed to help with clustering categorical datasets. The package can be installed from GitHub using devtools and then loaded in the usual way.

Installation
------------

`CClust` can be installed from the GitHub repository using the devtools package.

``` r
# Installs CClust
devtools::install_github("yudizhangzyd/CClust")
```

Once installed, the package can be loaded as usual.

``` r
# Loads the library
library(CClust)
```

Functions
---------

The package contains the following functions.

#### read\_fastq 

`read_fastq` obtains amplicon datasets from a fastq file, extract and return a list of data information: reads, quality scores and dimension of data.

#### kmodes

`kmodes` is for clustering categorical datasets without quality information of the data, different from function, it has six different random initialization methods and three k-means algorithms (Lloyd's; MacQueen's; Hartigan and Wong's algorithm) were adapted to do clustering.

#### khaplotype

`khaplotype` is for clustering the amplicon datasets with quality scores, only random initialization is avaiable, also three k-means algorithms (Lloyd's; MacQueen's; Hartigan and Wong's algorithm) were adapted to do clustering.

#### ARI

`ARI` is for computing the adjusted rand index given the estimated assignments and the true assiganments.

#### plot\_cluster

`plot_cluster` 

Visulize clusters after using the cluster algorithms implemented in the functions. This cluster plot make use of the function \code{dapc} from package \code{adegenet} and function \code{scatter} from package \code{ade4}. The plot is based on observation assignment and discriminant analysis of principal components, therefore, when running
the code, it will asks you to choose the number PCs to retain and choose the number discriminant functions to retain.

Learn More
----------

To learn more about `CClust`, read through the [vignette](https://yudizhangzyd.github.io/CClust/articles/CClust-vignette.html) for `CClust` which contains:

-   details on the three algorithms
-   explanations of how to use the functions
-   interpretations of the clustering asscessment

Additionally, more information is available at the package [website](https://yudizhangzyd.github.io/CClust/index.html).



