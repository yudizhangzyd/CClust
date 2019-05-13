#' Read a fastq file
#' @description
#' Read a fastq file and encode Phred quality score from 0 to 93 using ASCII 33 to 126.
#'
#' @param datafile Path to a fastq file.
#' @export read_fastq
#'
#' @return Return a list contains reads, quality socres and dim of the data.
#'
#' @examples
#' # Read a fastq file
#' dat <- read_fastq(datafile = "../data/sim.fastq")

read_fastq <- function(datafile = NULL)
{
  checkmate::expect_file_exists(datafile, access = "r")
  if (tail(unlist(strsplit(datafile, "[.]")), 1) != "fastq")
    stop("The input datafile has to be fastq file when setting run_with_quals = TRUE!")

  if (!is.loaded("r_read_fastq", PACKAGE = "khaplotype")) {
    dyn.load("../src/khaplotype.so")
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

