#' @title Cell Capture Efficiency Estimating Function
#' @description
#' This function estimates the capture efficiencies of cells from single-cell RNA-seq studies.
#' This function takes ERCC spike-in transcript and molecular concentration data, if available.
#' If spike-ins are not available, it uses count expression data.
#'
#' @param sce This is a gene-by-cell raw count matrix single cell experiment object with colData contains 'clusters', 'groups' and 'auxil' as a data frame object.
#' @param CE.range This is a two-element vector that sets the lower and upper limits for the estimated range of capture efficiencies.
#' @param RNAspike.use This is a logical parameter. If this parameter is set to TRUE, then it requires you to provide spike counts (spike) and spike concentration (spike.conc) information.
#' @param spike_in_sce This is a single cell experiment object which contains the spike-in counts data as the main assay and spike concentrations as its row data.
#' # spike-in counts: This is the observed count matrix for spike-in transcripts, where each row is a spike-in and each column is a cell. It's only required if you set RNAspike.use to TRUE.
#' # spike concentrations: This is a vector of theoretical counts for each spike-in transcript in a single cell, and it's only required if you set RNAspike.use to TRUE.
#' @param method A string that specifies the method for computing capture efficiencies. ("ML" or "")
#' @returns
#' This function returns the capture efficiencies of all the cells as a vector.
#' @importFrom stats lm
#' @importFrom SummarizedExperiment assays
#' @importFrom SingleCellExperiment rowData
#' @export
#' @examples
#' library(SingleCellExperiment)
#' # Load the test data.
#' data(SwarnSeqToyData); data(SpikeInData)
#' data <- assays(SwarnSeqToyData)[[1]][1:20, c(1:50, 350:399)]
#' groups <- SwarnSeqToyData$groups[c(1:50, 350:399)]
#' clusters <- SwarnSeqToyData$clusters[c(1:50, 350:399)]
#' X <- data.frame(clusters = clusters, groups = groups)
#' testData <- SingleCellExperiment(assays = list(counts = data), colData = X)
#'
#' capeff <- capEff(sce = testData,CE.range = c(0.01, 0.05),RNAspike.use = FALSE)
capEff <- function(sce, CE.range, RNAspike.use = FALSE, spike_in_sce, method = "ML")
    {
    CountData <- assays(sce)[[1]]
    if (RNAspike.use == TRUE) {
        spikes <- assays(spike_in_sce)[[1]]
        spike.conc <- SingleCellExperiment::rowData(spike_in_sce)[[1]]
    }
    if (RNAspike.use) {
        if (method == "ML") {
            capeff.spike <- apply(spikes, 2, sum) / sum(spike.conc)
            CE <- capeff.spike;names(CE) <- colnames(spikes);rm(spikes, capeff.spike)
        }else {
            CE <- vector(mode = "numeric", length = ncol(spikes))
            for (i in seq_len(length.out = ncol(spikes))) {
                spik <- as.numeric(spikes[,i]);spike.conc <- as.vector(spike.conc)
                mod.spike <- stats::lm(spik ~ spike.conc);CE[i] <- mod.spike$coefficients[2]
                names(CE) <- colnames(spikes);rm(spik, mod.spike)
            }
            CE <- ifelse(CE < 0, 0, CE)
            if (any(CE >= 1)) {
                message("CE can not be more than 1, please carefully check the inputs!")
            }
        }
    }else {
        if (is.null(CE.range)) {
            CE.range <- c(0.01, 0.2)
        }else {
            if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) {
                warning("CE.range is invalid!");return(invisible(NULL))
            }
        }
        l.sz <- log10(colSums(CountData));l.max <- max(l.sz);l.min <- min(l.sz)
        ls.wt <- (l.sz - l.min) / (l.max - l.min)
        rand.CE <- CE.range[1] + (CE.range[2] - CE.range[1]) * ls.wt
        CE <- rand.CE;names(CE) <- colnames(CountData)
    }
    return(CE)
}
