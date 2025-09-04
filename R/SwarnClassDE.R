#' @title ScRNA-seq Influential Gene Classification Function
#' @description
#' This function is used to clasify the influential genes of single-cell RNA-seq (scRNA-seq) data.
#' @param results This is obtained from 'SwarnUnadjLRT' or 'SwarnAdjLRT' function.
#' @param alpha Level of Significance
#' @returns
#' It returns the list with one more element, SwarnClassDE, that represents the SwarnClass.
#' @export
#' @examples
#' # Do not run.
#' library(SwarnSeq)
#' # Load the test data.
#' data(TestData)
#' CountData <- as.matrix(TestData$CountData[1:100,1:50])
#' X <- cbind(CountData, TestData$CountData[1:100,300:349])
#' group <- c(rep(1,50), rep(2,50))
#' CellCluster <- c(rep(1,30), rep(2,20), rep(3, 10), rep(4, 15), rep(5, 25))
#' group <- as.factor(group)
#' clust <- as.factor(CellCluster)
#' # CellAuxil <- as.factor(c(rep("A",14), rep("B",26), rep("C", 40), rep("E", 20)))  # Optional
#' res <- SwarnSeq::SwarnAdjLRT(CountData=X,norm.method="DEseq.norm",group=group,CellCluster=clust)
#' SwarnClass <- SwarnSeq::SwarnClassDE(results = res, alpha = 0.01)
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
