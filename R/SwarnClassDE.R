#' @title ScRNA-seq Influential Gene Classification Function
#' @description
#' This function is used to clasify the influential genes of single-cell RNA-seq (scRNA-seq) data.
#' @param results This is obtained from 'SwarnUnadjLRT' or 'SwarnAdjLRT' function.
#' @param alpha Level of Significance
#' @returns
#' It returns the list with one more element, SwarnClassDE, that represents the SwarnClass.
#' @export
#' @examples
#' # Paste the example into your R console and run it.
#' # Load the test data.
#' data("TestData")
#' CountData <- as.matrix(TestData$CountData)
#' group <- c(rep(1,200), rep(2,200))
#' CellCluster <- c(rep(1,60), rep(2,90), rep(3, 80), rep(4, 95), rep(5, 75))
#' group <- as.factor(group)
#' CellCluster <- as.factor(CellCluster)
#' CellAuxil <- as.factor(c(rep("A",55), rep("B",90), rep("C", 150), rep("E", 105)))
#' results <- SwarnSeq::SwarnAdjLRT(CountData = CountData, norm.method = "DEseq.norm", group = group, CellCluster = CellCluster, RNAspike.use = F)
#' SwarnClass <- SwarnSeq::SwarnClassDE(results = results, alpha = 0.01)
SwarnClassDE <- function(results, alpha) {
    show.custom.warning <- function(x) warning(x)
    var <- length(results)
    if (!is.list(results) & !is.matrix(results[[var]])) {
        show.custom.warning("Invalid input of the wrong data type of the results.")
        return(invisible(NULL))
    }
    if (ncol(results[[var]]) != 12) {
        show.custom.warning("Invalid input of the wrong number of columns of the results. Must be the same object from 'SwarnSeqLRT' or 'SwarnUnadjLRT.'")
        return(invisible(NULL))
    }
    if (!is.numeric(alpha)) {
        show.custom.warning("Invalid input of the wrong data type for m (number of tags).")
        return(invisible(NULL))
    }
    if (alpha <= 0 | alpha > nrow(results[[var]])) {
        show.custom.warning("Invalid input of the wrong value of m.")
        return(invisible(NULL))
    }
    class <- ifelse(results[[var]][,9] < alpha & results[[var]][,11] < alpha, "DE&DZI", ifelse(results[[var]][,9] < alpha, "DE", ifelse(results[[var]][,11] < alpha, "DZI", "NonDE")))
    results$SwarnClassDE <- class
    return(results)
}
