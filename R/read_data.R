#' Read a fastq file
#' @description
#' Read a fastq file and encode Phred quality score from 0 to 93 using ASCII 33 to 126.
#' (Reads need to be the same length)
#'
#' @param datafile Path to a fastq file.
#'
#' @useDynLib CClust r_read_fastq
#'
#' @importFrom checkmate expect_file_exists
#' @importFrom utils tail
#' @export read_fastq
#'
#' @return Return a list contains reads, quality socres and dim of the data.
#'
#' @examples
#' # Read a fastq file
#' datFile <- system.file("extdata", "sim.fastq", package = "CClust")
#' dat <- read_fastq(datafile = datFile)
#' # Read a fasta file
#' datFile <- system.file("extdata", "sim.fasta", package = "CClust")
#' dat <- read_fasta(datafile = datFile)

read_fastq <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fastq")
    stop("The input datafile has to be fastq file!")

  if (!is.loaded("r_read_fastq", PACKAGE = "CClust")) {
    dyn.load("../src/CClust.so")
  }

  res <- .Call("r_read_fastq", datafile)

  names(res) <- c("reads", "quality", "dim")
  res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)
  res$quality <- matrix(res$quality, ncol = res$dim[2], byrow = TRUE)
  res$reads[res$reads == 0] <- "A"
  res$reads[res$reads == 1] <- "C"
  res$reads[res$reads == 2] <- "T"
  res$reads[res$reads == 3] <- "G"

  return(res)
} #read_fastq

#' Read a fasta file
#' @description
#' Read a fasta file (Reads need to be the same length).
#'
#' @param datafile Path to a fasta file.
#'
#' @useDynLib CClust r_read_fasta
#'
#' @importFrom checkmate expect_file_exists
#' @importFrom utils tail
#' @export read_fasta
#'
#' @return Return a list contains reads and dim of the data.
#'
#' @examples
#' # Read a fasta file
#' datFile <- system.file("extdata", "sim.fasta", package = "CClust")
#' dat <- read_fasta(datafile = datFile)

read_fasta <- function(datafile = NULL) {
  checkmate::expect_file_exists(datafile, access = "r")
  if (utils::tail(unlist(strsplit(datafile, "[.]")), 1) != "fasta")
    stop("The input datafile has to be fasta file!")

  if (!is.loaded("r_read_fasta", PACKAGE = "CClust")) {
    dyn.load("../src/CClust.so")
  }

  res <- .Call("r_read_fasta", datafile)
  names(res) <- c("reads", "dim")
  res$reads <- matrix(res$reads, ncol = res$dim[2], byrow = TRUE)

  res$reads[res$reads == 1] <- "A"
  res$reads[res$reads == 2] <- "C"
  res$reads[res$reads == 8] <- "T"
  res$reads[res$reads == 4] <- "G"
  res$reads[res$reads == 5] <- "R"
  res$reads[res$reads == 10] <- "Y"
  res$reads[res$reads == 6] <- "S"
  res$reads[res$reads == 9] <- "W"
  res$reads[res$reads == 12] <- "K"
  res$reads[res$reads == 3] <- "M"
  res$reads[res$reads == 15] <- "N"
  res$reads[res$reads == 0] <- "X"
  res$reads[res$reads == 7] <- "V"
  res$reads[res$reads == 14] <- "B"
  res$reads[res$reads == 13] <- "D"
  res$reads[res$reads == 11] <- "H"

  return(res)
}

