#' @title Cell Capture Efficiency Estimating Function
#' @description
#' This function estimates the capture efficiencies of cells from single-cell RNA-seq studies.
#' This function takes ERCC spike-in transcript and molecular concentration data, if available.
#' If spike-ins are not available, it uses count expression data.
#' @param CountData This is a gene-by-cell raw count matrix.
#' @param CE.range This is a two-element vector that sets the lower and upper limits for the estimated range of capture efficiencies.
#' @param RNAspike.use This is a logical parameter. If this parameter is set to TRUE, then it requires you to provide spike counts (spike) and spike concentration (spike.conc) information.
#' @param spikes This is the observed count matrix for spike-in transcripts, where each row is a spike-in and each column is a cell. It's only required if you set RNAspike.use to TRUE.
#' @param spike.conc This is a vector of theoretical counts for each spike-in transcript in a single cell, and it's only required if you set RNAspike.use to TRUE.
#' @param method A string that specifies the method for computing capture efficiencies. ("ML" or "")
#' @returns
#' This function returns the capture efficiencies of all the cells as a vector.
#' @importFrom stats lm
#' @export
#' @examples
#' # Load the data
#' data("TestData")
#' X <- as.matrix(TestData$CountData); spike.counts <- TestData$SpikeCounts
#' spike_conc <- TestData$SpikeConc
#'
#' capeff <- SwarnSeq::CapEff(CountData = X,CE.range = c(0.01, 0.05),RNAspike.use = FALSE,method = "")
CapEff <- function(CountData, CE.range, RNAspike.use, spikes, spike.conc, method){
    show.custom.metod <- function(x) message(x);show.custom.warning <- function(x) warning(x)
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
                show.custom.metod("CE can not be more than 1, please carefully check the inputs!")
            }
        }
    }else {
        if (is.null(CE.range)) {
            CE.range <- c(0.01, 0.2)
        }else {
            if (CE.range[1] < 0 | CE.range[1] > CE.range[2] | CE.range[2] > 1) {
                show.custom.warning("CE.range is invalid!");return(invisible(NULL))
            }
        }
        l.sz <- log10(colSums(CountData));l.max <- max(l.sz);l.min <- min(l.sz)
        ls.wt <- (l.sz - l.min) / (l.max - l.min)
        rand.CE <- CE.range[1] + (CE.range[2] - CE.range[1]) * ls.wt
        CE <- rand.CE;names(CE) <- colnames(CountData)
    }
    return(CE)
}
