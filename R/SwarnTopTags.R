#' @title Most Statistically Significant Gene Filtering Function
#' @description
#' This function selects the most statistically significant genes in your scRNA-seq data.
#' @param results This is obtained from 'SwarnUnadjLRT' or 'SwarnAdjLRT' function.
#' @param m This is an integer that represents the number of statistically significant genes to be selected from the downstream analysis.
#' @returns
#' A list of the top genes along with their statistics.
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
#' # SpikeInData <- SingleCellExperiment(assays = list(SpikeCounts), rowData = SpikeConc)
#' # Make the spike-in single cell experiment object like this.
#'
#' res <- swarnAdjLrt(sce=testData,norm.method="DEseq.norm", RNAspike.use = FALSE)
#' top_genes <- swarnTopTags(results = res, m = 50)
swarnTopTags <- function(results, m){
    var <- length(results)
    if (!is.list(results) & !is.matrix(results[[var]])) {
        warning("Invalid input of the wrong data type of the results.")
        return(invisible(NULL))}
    if (ncol(results[[var]]) != 12) {
        warning("Invalid input of the wrong number of columns of the results. Must be the same object from 'SwarnUnadjLRT' or 'SwarnAdjLRT.'")
        return(invisible(NULL))}
    if (!is.numeric(m)) {
        warning("Invalid input of the wrong data type for m (number of tags).")
        return(invisible(NULL))}
    if (m <= 0 | m > nrow(results[[var]])) {
        warning("Invalid input of the wrong value of m.")
        return(invisible(NULL))}
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
