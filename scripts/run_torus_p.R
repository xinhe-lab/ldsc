#source("renv/activate.R")

save.image("ip.RData")
str(snakemake)
library(daprcpp)
library(dplyr)
library(purrr)
library(readr)
library(ldmap)
library(fs)
library(tidyr)
library(stringr)


gf <- snakemake@input[["gwasf"]]
af <- snakemake@input[["annof"]]
if(is.null(af)){
af <- tempfile()
write_tsv(tibble::tibble(SNP=character()),af)
}

prior_rf <- snakemake@input[["prior_r"]]

  prior_r <- scan(prior_rf,what = character())
  od <- snakemake@output[["outputd"]]
  torus_ret <- daprcpp:::run_torus_cmd(gf=gf,af=af,torus_p=prior_r)

  saveRDS(torus_ret$df,snakemake@output[["outputf"]])
if (!dir.exists(od)) {
  fs::dir_create(od,recurse = TRUE)
  }
  iwalk(torus_ret$priors,function(pr,region_id) {
  trid <- region_id
  separate(pr,SNP,c("chrom","pos","ref","alt"),convert=TRUE) %>%
    dplyr::mutate(chrom=as.integer(stringr::str_replace(chrom,"chr",""))) %>%
    compact_snp_struct(snp_struct="SNP") %>% 
    saveRDS(fs::path(od,trid,ext = "RDS"))
  })
