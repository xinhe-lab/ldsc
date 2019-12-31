source("renv/activate.R")
cat("starting!\n")
library(daprcpp)

saveRDS(torus_fdr(snakemake@input[["gwasf"]],snakemake@input[["annof"]]),snakemake@output[["outputf"]])
