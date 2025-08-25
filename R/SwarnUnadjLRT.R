#' @title Differential gene expression analysis for single-cell data without Capture Efficiency Correction
#' @description
#' This function performs differential gene expression (DEG) analysis on a single-cell RNA-sequencing (scRNA-seq) raw read counts matrix.
#' The input is a non-negative integer matrix, where rows correspond to genes and columns to individual cells.
#' The primary objective is to identify genes with statistically significant changes in expression levels between teo predefined groups of cells.
#' This analysis is a downstream step following the initial read mapping and gene counting, which generates the input matrix. By contrasting the expression profiles of the two groups, the function identifies biomarkers or genes that characterize the biological differences between them.
#' @param CountData This is a gene-by-cell raw count matrix.
#' @param norm.method There are three distinct statistical normalization methods for scRNA-seq count data:DESeq2, a maximum likelihood approach described by Ye et al. (2027); TMM, a trimmed mean method introduced by Robinson et al. (2010), both designed to correct for library size and compositional biases; and a third method, which is the log1p standard normalization method for scRNA-seq data.
#' @param group This is a vector of factors with 2 levels that serves as a grouping variable that assigns each column (cell) of the count matrix to one of two predefined experimental conditions.
#' @param CellCluster This is a vector of factors with the optimal number of cluster levels, with each entry assigning a cell from the count matrix's columns to its computationally inferred cluster identity.
#' @param CellAuxil This is a vector of factors of cell-level metadata that provides additional information about each cell, such as its batch, donor, or quality control metrics, with each entry aligning to a column in the count matrix.
#' @param maxit This is the maximum number of iterations for the Expected-Maximization (EM) algorithm.
#' @param eps This is the convergence criteria for the Expected-Maximization (EM) algorithm.
#' @param muoffset This is the offset parameter for mean (mu), with a default value of NULL.
#' @param phioffset This is the offset parameter for the zero inflation (phi) parameter, with a default value of NULL.
#' @param weights This is the observation-wise weight for the cells, with a default value of NULL.
#' @param method A string that specifies the method for computing capture efficiencies. ("ML" or "")
#' @param p.value.adj.method A string that specifies the method for multiple hypothesis correction, which can be any of the following values: "holm," "hochberg," "hommel," "bonferroni," "BH," "BY," "fdr," or "none").
#' @returns
#' This function analyzes single-cell RNA-sequencing unadjusted normalized data to find differentially expressed genes and returns a list of the results.
#' @importFrom stats model.matrix
#' @importFrom stats glm.fit
#' @importFrom MASS glm.nb
#' @importFrom edgeR calcNormFactors
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
#' results <- SwarnSeq::SwarnUnadjLRT(CountData = CountData, norm.method = "DEseq.norm", group = group, CellCluster = CellCluster, CellAuxil = CellAuxil)
SwarnUnadjLRT <- function(CountData, norm.method = c("DEseq.norm", "TMM", "log1p"), group, CellCluster, CellAuxil = NULL, maxit = 100, eps = 1E-4, muoffset = NULL, phioffset = NULL, weights = NULL, method = NA, p.value.adj.method = "BH"){
    show.custom.warning <- function(x) warning(x);show.custom.metod <- function(x) message(x)
    if (!is.matrix(CountData)) {
        show.custom.warning("Wrong input data type of count data...")
        return(invisible(NULL))
    }
    if (sum(is.na(CountData)) > 0) {
        show.custom.warning("NAs are detected in the input count data...")
        return(invisible(NULL))
    }
    if (sum(CountData < 0) > 0) {
        show.custom.warning("Negative values are detected in the input count data...")
    }
    if (all(CountData == 0)) {
        show.custom.warning("All elements of the input count data are zeros...")
        return(invisible(NULL))
    }
    if (length(unique(group)) != 2) {
        show.custom.warning("Factor levels of group is not two...")
        return(invisible(NULL))
    }
    if (table(group)[1] < 2 | table(group)[2] < 2) {
        show.custom.warning("Too few samples (< 2) in a group...")
        return(invisible(NULL))
    }
    if (ncol(CountData) != length(group) | ncol(CountData) != length(CellCluster)) {
        show.custom.warning("The length of 'group' & 'CellCluster' must be equal to the number of columns of the count data...")
        return(invisible(NULL))
    }
    if (!is.null(CellAuxil) & ncol(CountData) != length(CellAuxil)) {
        show.custom.warning("The length of the cell-level covariate factor should be the same as the number of cells in your count data...")
        return(invisible(NULL))
    }
    if (!is.numeric(c(maxit, eps))) {
        show.custom.warning("The data type of maxit and eps is not numeric...")
        return(invisible(NULL))
    }
    if (length(maxit) != 1) {
        show.custom.warning("The length of 'maxit' is not one...")
        return(invisible(NULL))
    }
    if (length(eps) != 1) {
        show.custom.warning("The length of the convergence criterion is not one...")
        return(invisible(NULL))
    }
    CountData <- CountData[rowSums(CountData) > 0,]
    if (!is.null(CellAuxil)) {
        group <- group[colSums(CountData) > 0]
        CellCluster <- CellCluster[colSums(CountData) > 0]
        CellAuxil <- CellAuxil[colSums(CountData) > 0]
    }else {
        group <- group[colSums(CountData) > 0]
        CellCluster <- CellCluster[colSums(CountData) > 0]
    }
    CountData <- CountData[,colSums(CountData) > 0]
    if (is.null(dim(CountData))) {
        show.custom.warning("There may be an error in the input data dimensions. This could be due to having fewer than two genes with at least one read across all or some of the cells.")
        return(invisible(NULL))
    }
    if (norm.method == "DEseq.norm") {
        GM <- function(x) exp(mean(log(x[x > 0])))
        geomMean <- apply(CountData, 1, GM)
        fx <- function(x) x/geomMean
        samp <- apply(CountData, 2, fx)
        fxx <- function(xx) median(xx[xx != 0])
        size <- apply(samp, 2, fxx)
        fxxx <- function(xxx) xxx/size
        CountData <- t(apply(CountData, 1, fxxx))
        CountData <- ceiling(CountData);rm(geomMean, samp, size, GM, fx, fxx, fxxx)
    }else if (norm.method == "TMM") {
        m <- size_factor <- edgeR::calcNormFactors(CountData)
        fxxx <- function(xxx) xxx/size_factor
        CountData <- t(apply(CountData, 1, fxxx))
        CountData <- ceiling(CountData);rm(m, fxxx)
    }else {CountData <- log1p(CountData);CountData <- ceiling(CountData)}
    # Covariates...
    group <- as.factor(group);CellCluster <- as.factor(CellCluster)
    # Catch total mean and foldchange for each gene.
    totalMeans_1 <- rowMeans(CountData[row.names(CountData), group == levels(group)[1]])
    totalMeans_2 <- rowMeans(CountData[row.names(CountData), group == levels(group)[2]])
    All_Means <- cbind(totalMeans_1, totalMeans_2)
    colnames(All_Means) <- c("All_Control_Means", "All_Other_Means")
    Fold_Change <- totalMeans_2 / totalMeans_1
    Log2_Fold_Change <- log2(Fold_Change)
    results <- ZINBEM(counts_data = CountData, group, CellCluster, CellAuxil, weights, muoffset, phioffset, maxit, eps)
    # Adjusting the p-values...
    show.custom.metod(c("Correction for multiple hypothesis testing."))
    adj.method <- p.adjust.methods[match(p.value.adj.method, p.adjust.methods)]
    if (is.na(adj.method)) {adj.method <- "BH"}
    p.DE <- as.vector(results$DE.stat[,2])
    adj.pval.DE <- p.adjust(p.DE, method = adj.method, n = length(p.DE))
    FDR.DE <- p.adjust(p.DE, method = "fdr", n = length(p.DE)) ##
    p.DZI <- as.vector(results$DE.stat[,4])
    adj.pval.DZI <- p.adjust(p.DZI, method = adj.method, n = length(p.DZI))
    FDR.DZI <- p.adjust(p.DZI, method = "fdr", n = length(p.DZI)) ##
    p.adj <- cbind(adj.pval.DE, FDR.DE, adj.pval.DZI, FDR.DZI)
    colnames(p.adj) <- c("DE.Adj.Pvalue", "DE.FDR", "DZI.Adj.Pval", "DZI.FDR")
    rm(p.DE, p.DZI, adj.pval.DE, adj.pval.DZI, FDR.DE, FDR.DZI)
    Final_Results <- cbind(All_Means, Fold_Change, Log2_Fold_Change, results$DE.stat, p.adj)
    results$SwarnUnadjLRT <- Final_Results
    return(results)
}
