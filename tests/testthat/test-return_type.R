library(SwarnSeq)
library(SingleCellExperiment)
# Load the SwarnSeq test datasets.
data(SwarnSeqToyData); data(SpikeInData)
data <- assays(SwarnSeqToyData)[[1]][1:20, c(1:50, 350:399)]
groups <- SwarnSeqToyData$groups[c(1:50, 350:399)]
clusters <- SwarnSeqToyData$clusters[c(1:50, 350:399)]
X <- data.frame(clusters = clusters, groups = groups)
testData <- SingleCellExperiment(assays = list(counts = data), colData = X)

test_that("This function returns a list", {
    result <- swarnUnadjLrt(sce = testData, norm.method = "DEseq.norm")
    expect_type(result, "list")
})
test_that("This function returns a list", {
    result <- swarnAdjLrt(sce=testData,norm.method="log1p",RNAspike.use=FALSE)
    expect_type(result, "list")
})
test_that("Output vector length matches input cells", {
    cap_eff <- capEff(sce = testData, CE.range = c(0.05, 0.60), RNAspike.use = FALSE, spike_in_sce = SpikeInData)
    expect_length(cap_eff, dim(assays(testData)[[1]])[2])
})
test_that("all data frames in the list have the correct row count", {
    extractedAdjNormData <- extAdjNormData(sce=testData,norm.method = "log1p",CE.range = c(0.01,0.5))
    expect_equal(nrow(extractedAdjNormData$control), nrow(extractedAdjNormData$infected))
})
test_that("Output contains only allowed categories", {
    allowed_values <- c("NonDE", "DZI", "DE&DZI", "DE")
    result <- swarnAdjLrt(sce=testData,norm.method="log1p",RNAspike.use=FALSE)
    SwarnClass <- swarnClassDe(results = result, alpha = 0.01)
    expect_true(all(unique(SwarnClass$SwarnClassDE) %in% allowed_values))
})
test_that("The two matrices in the list have equal gene counts", {
    result <- swarnAdjLrt(sce=testData,norm.method="log1p",RNAspike.use=FALSE)
    top_genes <- swarnTopTags(results = result, m = 10)
    expect_length(top_genes, 2)
    expect_equal(nrow(top_genes[[1]]), nrow(top_genes[[2]]))
})




