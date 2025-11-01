library(SwarnSeq)
library(SingleCellExperiment)
library(SummarizedExperiment)
# Load the test data.
data(SwarnSeqToyData); data(SpikeInData)
data <- assays(SwarnSeqToyData)[[1]][1:20, c(1:50, 350:399)]
groups <- SwarnSeqToyData$groups[c(1:50, 350:399)]
clusters <- SwarnSeqToyData$clusters[c(1:50, 350:399)]
X <- data.frame(clusters = clusters, groups = groups)
testData <- SingleCellExperiment(assays = list(counts = data), colData = X)

test_that("This function returns a list", {
    result <- swarnAdjLrt(sce=testData,norm.method="log1p",RNAspike.use=FALSE)
    expect_type(result, "list")
})
