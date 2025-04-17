library(TrajConserve)
exampledata = TrajConserve::load_example_data()
exampledataSce <- as.SingleCellExperiment(exampledata)
zellkonverter::writeH5AD(exampledataSce,"processed_data/toy_data/20250412_example_trajconserve.h5ad")
