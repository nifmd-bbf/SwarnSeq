#' @title A test scRNA-seq Real Dataset for SwarnSeq
#' @description
#' A test dataset containing a single-cell RNA-seq (scRNA-seq) read counts matrix, ERCC spike-in counts, and spike-in concentration.
#'
#' @format A list of 3 element R objects.
#' @name Test_Dataset
#'
#' @usage data(TestData)
#'
#' \describe{
#' \item{CountData}{This is a gene-by-cell raw count matrix.}
#' \item{SpikeCounts}{This is the observed count matrix for spike-in transcripts, where each row is a spike-in and each column is a cell.}
#' \item{SpikeConc}{This is a vector of theoretical counts for each spike-in transcript in a single cell.}
#' }
#'
#'@source Petropoulos, S., Edsgärd, D., Reinius, B., Deng, Q., Panula, S. P., Codeluppi, S., Plaza Reyes, A., Linnarsson, S., Sandberg, R., & Lanner, F. (2016). Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos. Cell, 165(4), 1012–1026. https://doi.org/10.1016/j.cell.2016.03.023
#'
#'
