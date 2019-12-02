source("renv/activate.R")
## save.image("ip.RData")
## stop()
  library(daprcpp)
  library(dplyr)
  library(purrr)
  library(readr)
  library(ldmap)
  library(fs)

  prior_r <- scan(snakemake@input[["prior_r"]],what = character())
  od <- snakemake@output[["outputd"]]
  torus_ret <- daprcpp:::run_torus_cmd(gf=snakemake@input[["gwasf"]],af=snakemake@input[["annof"]],torus_p=prior_r)

  saveRDS(torus_ret$df,snakemake@output[["outputf"]])

  iwalk(torus_ret$prior,function(pr,region_id) {
  write_tsv(pr,fs::path(od,region_id,ext = "txt.gz"))
  })
