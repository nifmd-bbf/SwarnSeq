#' @title Two-Factor Normalized Counts Matrix Splitting Function
#' @description
#' This function adjusts the scRNA-seq data by capture efficiency for each cell, then normalizes the adjusted matrix, and finally splits the matrix into two distinct matrices based on the number of factors in the group vector.
#' @param CountData This is a gene-by-cell raw count matrix.
#' @param norm.method There are three distinct statistical normalization methods for scRNA-seq count data:DESeq2, a maximum likelihood approach described by Ye et al. (2027); TMM, a trimmed mean method introduced by Robinson et al. (2010), both designed to correct for library size and compositional biases; and a third method, which is the log1p standard normalization method for scRNA-seq data.
#' @param group This is a vector of factors with 2 levels that serves as a grouping variable that assigns each column (cell) of the count matrix to one of two predefined experimental conditions.
#' @param method A string that specifies the method for computing capture efficiencies. ("ML" or "")
#' @param RNAspike.use This is a logical parameter. If this parameter is set to TRUE, then it requires you to provide spike counts (spike) and spike concentration (spike.conc) information.
#' @param spikes This is the observed count matrix for spike-in transcripts, where each row is a spike-in and each column is a cell. It's only required if you set RNAspike.use to TRUE.
#' @param spike.conc This is a vector of theoretical counts for each spike-in transcript in a single cell, and it's only required if you set RNAspike.use to TRUE.
#' @param CE.range This is a two-element vector that sets the lower and upper limits for the estimated range of capture efficiencies.
#' @returns
#' This function returns a list of two distinct capture efficiency adjusted normalized counts matrices, based on the number of factors in the group vector for the scRNA-seq data.
#' @importFrom stats median
#' @importFrom stats lm
#' @export
#' @examples
#' # Load the test data.
#' data("TestData")
#' X <- as.matrix(TestData$CountData); spike_counts <- TestData$SpikeCounts
#' spike_conc <- TestData$SpikeConc; N.M <- "DEseq.norm"
#' group <- c(rep(1,200), rep(2,200))
#' group <- as.factor(group)
#' ExtData <- SwarnSeq::ExtAdjNormData(CountData=X,norm.method=N.M,group=group,CE.range=c(0.01,0.5))
ExtAdjNormData <- function(CountData, norm.method = c("DEseq.norm", "TMM", "log1p"), group, method = "ML", RNAspike.use = FALSE, spikes, spike.conc, CE.range = c(.01,.5)) {
    show.custom.metod <- function(x) message(x);show.custom.warning <- function(x) warning(x)
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
    if (ncol(CountData) != length(group)) {
        show.custom.warning("The length of 'group' & 'CellCluster' must be equal to the number of columns of the count data...")
        return(invisible(NULL))
    }
    CountData <- CountData[rowSums(CountData) > 0,]
    group <- group[colSums(CountData) > 0]
    CountData <- CountData[,colSums(CountData) > 0]
    if (is.null(dim(CountData))) {
        show.custom.warning("There may be an error in the input data dimensions. This could be due to having fewer than two genes with at least one read across all or some of the cells.")
        return(invisible(NULL))
    }
    if (RNAspike.use) {
        if (method == "ML") {
            capeff.spike <- apply(spikes, 2, sum) / sum(spike.conc)
            CE <- capeff.spike;names(CE) <- colnames(spikes);rm(spikes, capeff.spike)
        }else {
            CE <- vector(mode = "numeric", length = ncol(spikes))
            for (i in seq_len(length.out = ncol(spikes))) {
                spik <- as.numeric(spikes[,i])
                spike.conc <- as.vector(spike.conc)
                mod.spike <- stats::lm(spik ~ spike.conc)
                CE[i] <- mod.spike$coefficients[2]
                names(CE) <- colnames(spikes);rm(spik, mod.spike)
            }
            CE <- ifelse(CE < 0, 0, CE)
            if (any(CE >= 1)) {
                show.custom.metod("CE can not be more than 1, please carefully check the inputs!")
            }
        }
    }else {
        if (is.null(CE.range)) {
            CE.range <- c(0.01, 0.2)
        }else {
            if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) {
                show.custom.warning("CE.range is invalid!")
                return(invisible(NULL))
            }
        }
        l.sz <- log10(colSums(CountData));l.max <- max(l.sz);l.min <- min(l.sz)
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
        fxx <- function(xx) stats::median(xx[xx != 0])
        size <- apply(samp, 2, fxx)
        fxxx <- function(xxx) xxx/size
        CountData <- t(apply(CountData, 1, fxxx));CountData <- ceiling(CountData)
        rm(geomMean, samp, size, GM, fx, fxx, fxxx)
    }else if (norm.method == "TMM") {
        m <- size_factor <- edgeR::calcNormFactors(CountData)
        fxxx <- function(xxx) xxx/size_factor
        CountData <- t(apply(CountData, 1, fxxx));CountData <- ceiling(CountData)
        rm(m, fxxx)
    }else {CountData <- log1p(CountData);CountData <- ceiling(CountData)}
    control_data <- CountData[,group == levels(group)[1]]
    infected_data <- CountData[,group == levels(group)[2]]
    return(list(control = control_data, infected = infected_data))
}
