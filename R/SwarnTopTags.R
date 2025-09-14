#' @title Most Statistically Significant Gene Filtering Function
#' @description
#' This function selects the most statistically significant genes in your scRNA-seq data.
#' @param results This is obtained from 'SwarnUnadjLRT' or 'SwarnAdjLRT' function.
#' @param m This is an integer that represents the number of statistically significant genes to be selected from the downstream analysis.
#' @returns
#' A list of the top genes along with their statistics.
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
#' top_genes <- SwarnSeq::SwarnTopTags(results = res, m = 50)
SwarnTopTags <- function(results, m){
    show.custom.warning <- function(x) warning(x);var <- length(results)
    if (!is.list(results) & !is.matrix(results[[var]])) {
        show.custom.warning("Invalid input of the wrong data type of the results.")
        return(invisible(NULL))
    }
    if (ncol(results[[var]]) != 12) {
        show.custom.warning("Invalid input of the wrong number of columns of the results. Must be the same object from 'SwarnUnadjLRT' or 'SwarnAdjLRT.'")
        return(invisible(NULL))
    }
    if (!is.numeric(m)) {
        show.custom.warning("Invalid input of the wrong data type for m (number of tags).")
        return(invisible(NULL))
    }
    if (m <= 0 | m > nrow(results[[var]])) {
        show.custom.warning("Invalid input of the wrong value of m.")
        return(invisible(NULL))
    }
    # Top Tags
    p.DE <- results[[var]][,9];p.DZI <- results[[var]][,11]
    id <- sort(p.DE, decreasing = FALSE, index.return = TRUE)$ix
    id <- id[seq_len(length.out = m)]
    top.DE <- results[[var]][id,] ##
    id2 <- sort(p.DZI, decreasing = FALSE, index.return = TRUE)$ix
    id2 <- id2[seq_len(length.out = m)]
    top.DZI <- results[[var]][id2,] ##
    out <- list(Top.DEG = top.DE, Top.DZI = top.DZI)
    rm(results, id, id2, p.DE, p.DZI)
    return(out)
}
