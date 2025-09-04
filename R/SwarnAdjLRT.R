#' @title Differential gene expression analysis for single-cell data
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
#' @param RNAspike.use This is a logical parameter. If this parameter is set to TRUE, then it requires you to provide spike counts (spike) and spike concentration (spike.conc) information.
#' @param spikes This is the observed count matrix for spike-in transcripts, where each row is a spike-in and each column is a cell. It's only required if you set RNAspike.use to TRUE.
#' @param spike.conc This is a vector of theoretical counts for each spike-in transcript in a single cell, and it's only required if you set RNAspike.use to TRUE.
#' @param CE.range This is a two-element vector that sets the lower and upper limits for the estimated range of capture efficiencies.
#' @returns
#' This function analyzes single-cell RNA-sequencing adjusted normalized data to find differentially expressed genes and returns a list of the results.
#' @importFrom stats model.matrix
#' @importFrom stats glm.fit
#' @importFrom MASS glm.nb
#' @importFrom edgeR calcNormFactors
#' @importFrom stats median
#' @importFrom stats lm
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
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
SwarnAdjLRT <- function(CountData, norm.method = c("DEseq.norm", "TMM", "log1p"), group, CellCluster, CellAuxil = NULL, maxit = 100, eps = 1E-4, muoffset = NULL, phioffset = NULL, weights = NULL, method = "ML", p.value.adj.method = "BH", RNAspike.use = FALSE, spikes, spike.conc, CE.range = c(.01,.5)){
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
    # Covariates...
    group <- as.factor(group);CellCluster <- as.factor(CellCluster)
    if (RNAspike.use) {
        if (method == "ML") {
            capeff.spike <- apply(spikes, 2, sum) / sum(spike.conc)
            CE <- capeff.spike;names(CE) <- colnames(spikes)
            rm(spikes, capeff.spike)
        }else {
            CE <- vector(mode = "numeric", length = ncol(spikes))
            for (i in seq_len(length.out = ncol(spikes))) {
                spik <- as.numeric(spikes[,i])
                spike.conc <- as.vector(spike.conc)
                mod.spike <- stats::lm(spik ~ spike.conc)
                CE[i] <- mod.spike$coefficients[2];names(CE) <- colnames(spikes)
                rm(spik, mod.spike)
            }
            CE <- ifelse(CE < 0, 0, CE)
            if (any(CE >= 1)) {
                show.custom.metod(c("CE can not be more than 1, please carefully check the inputs!"))
            }
        }
    }else {
        if (is.null(CE.range)) {
            CE.range <- c(0.01, 0.2)
        }else {
            if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) {
                warning("CE.range is invalid!")
                return(invisible(NULL))
            }
        }
        l.sz <- log10(colSums(CountData))
        l.max <- max(l.sz);l.min <- min(l.sz)
        ls.wt <- (l.sz - l.min) / (l.max - l.min)
        rand.CE <- CE.range[1] + (CE.range[2] - CE.range[1]) * ls.wt
        CE <- rand.CE;names(CE) <- colnames(CountData)
    }
    CountData <- sweep(CountData, 2, CE, '/')
    if (norm.method == "DEseq.norm") {
        GM <- function(x) exp(mean(log(x[x > 0])))
        geomMean <- apply(CountData, 1, GM)
        fx <- function(x) x/geomMean
        samp <- apply(CountData, 2, fx)
        fxx <- function(xx) median(xx[xx != 0])
        size <- apply(samp, 2, fxx)
        fxxx <- function(xxx) xxx/size
        CountData <- t(apply(CountData, 1, fxxx));CountData <- ceiling(CountData)
        rm(geomMean, samp, size, GM, fx, fxx, fxxx)
    }else if (norm.method == "TMM") {
        m <- size_factor <- edgeR::calcNormFactors(CountData)
        fxxx <- function(xxx) xxx/size_factor
        CountData <- t(apply(CountData, 1, fxxx))
        CountData <- ceiling(CountData);rm(m, fxxx)
    }else {CountData <- log1p(CountData);CountData <- ceiling(CountData)}
    # Catch total mean and foldchange for each gene.
    totalAdjMeans_1 <- rowMeans(CountData[row.names(CountData), group == levels(group)[1]])
    totalAdjMeans_2 <- rowMeans(CountData[row.names(CountData), group == levels(group)[2]])
    All_Adj_Means <- cbind(totalAdjMeans_1, totalAdjMeans_2)
    colnames(All_Adj_Means) <- c("All_Control_Adj_Means", "All_Other_Adj_Means")
    Adj_Fold_Change <- totalAdjMeans_2 / totalAdjMeans_1
    Adj_Log2_Fold_Change <- log2(Adj_Fold_Change)
    results <- ZINBEM(counts_data = CountData, group, CellCluster, CellAuxil, weights, muoffset, phioffset, maxit, eps)
    # Adjusting the p-values...
    show.custom.metod(c("\nCorrection for multiple hypothesis testing."))
    adj.method <- stats::p.adjust.methods[match(p.value.adj.method, stats::p.adjust.methods)]
    if (is.na(adj.method)) {adj.method <- "BH"}
    p.DE <- as.vector(results$DE.stat[,2])
    adj.pval.DE <- stats::p.adjust(p.DE, method = adj.method, n = length(p.DE))
    FDR.DE <- stats::p.adjust(p.DE, method = "fdr", n = length(p.DE))
    p.DZI <- as.vector(results$DE.stat[,4])
    adj.pval.DZI <- stats::p.adjust(p.DZI, method = adj.method, n = length(p.DZI))
    FDR.DZI <- stats::p.adjust(p.DZI, method = "fdr", n = length(p.DZI))
    p.adj <- cbind(adj.pval.DE, FDR.DE, adj.pval.DZI, FDR.DZI)
    colnames(p.adj) <- c("DE.Adj.Pvalue", "DE.FDR", "DZI.Adj.Pval", "DZI.FDR")
    rm(p.DE, p.DZI, adj.pval.DE, adj.pval.DZI, FDR.DE, FDR.DZI)
    Final_Results <- cbind(All_Adj_Means = All_Adj_Means, Adj_Fold_Change = Adj_Fold_Change, Adj_Log2_Fold_Change = Adj_Log2_Fold_Change, results$DE.stat, p.adj)
    results$SwarnAdjLRT <- Final_Results
    return(results)
}
