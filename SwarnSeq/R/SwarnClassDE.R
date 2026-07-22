#' @title ScRNA-seq Influential Gene Classification Function
#' @description
#' This function is used to clasify the influential genes of single-cell RNA-seq (scRNA-seq) data.
#' @param results This is obtained from 'SwarnUnadjLRT' or 'SwarnAdjLRT' function.
#' @param alpha Level of Significance
#' @returns
#' It returns the list with one more element, SwarnClassDE, that represents the SwarnClass.
#' @export
#' @examples
#' library(SingleCellExperiment)
#' # Load the test data.
#' data(SwarnSeqToyData); data(SpikeInData)
#' data <- assays(SwarnSeqToyData)[[1]][1:20, c(1:50, 350:399)]
#' groups <- SwarnSeqToyData$groups[c(1:50, 350:399)]
#' clusters <- SwarnSeqToyData$clusters[c(1:50, 350:399)]
#'
#' X <- data.frame(clusters = clusters, groups = groups)
#' testData <- SingleCellExperiment(assays = list(counts = data), colData = X)
#' # SpikeInData <- SingleCellExperiment(assays=list(SpikeCounts),rowData=SpikeConc)
#' # Make the spike-in single cell experiment object like this.
#'
#' res <- swarnAdjLrt(sce=testData,norm.method="log1p",RNAspike.use=TRUE,spike_in_sce=SpikeInData)
#' SwarnClass <- swarnClassDe(results = res, alpha = 0.01)
swarnClassDe <- function(results, alpha) {
    var <- length(results)
    if (!is.list(results) & !is.matrix(results[[var]])) {
        warning("Invalid input of the wrong data type of the results.")
        return(invisible(NULL))
    }
    if (ncol(results[[var]]) != 12) {
        warning("Invalid input of the wrong number of columns of the results. Must be the same object from 'SwarnSeqLRT' or 'SwarnUnadjLRT.'")
        return(invisible(NULL))
    }
    if (!is.numeric(alpha)) {
        warning("Invalid input of the wrong data type for m (number of tags).")
        return(invisible(NULL))
    }
    if (alpha <= 0 | alpha > nrow(results[[var]])) {
        warning("Invalid input of the wrong value of m.")
        return(invisible(NULL))
    }
    class <- ifelse(results[[var]][,9] < alpha & results[[var]][,11] < alpha, "DE&DZI", ifelse(results[[var]][,9] < alpha, "DE", ifelse(results[[var]][,11] < alpha, "DZI", "NonDE")))
    results$SwarnClassDE <- class
    return(results)
}
