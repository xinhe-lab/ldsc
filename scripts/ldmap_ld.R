source("renv/activate.R")
library(bigsnpr)
library(ldmap)
library(dplyr)
library(purrr)

data(ldetect_EUR)
rdsf <- snakemake@input[["rdslist"]]
ld_id <- snakemake@params[["region"]]
stopifnot(!is.null(ld_id),length(ld_id)==1)
ld_id <- as.integer(ld_id)
ldmr <- ldetect_EUR[ld_id]
srds <- ldmap:::subset_rds(ldmr = ldmr,reference_files = rdsf,pattern = paste0(ld_id,"_"))
saveRDS(panel_ld(srds,FALSE),snakemake@output[["ldf"]])
