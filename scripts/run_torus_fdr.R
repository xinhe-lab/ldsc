source("renv/activate.R")
library(daprcpp)

saveRDS(torus_fdr(snakemake@input[["gwasf"]],snakemake@input[["annof"]]),snakemake@output[["outputf"]])
