source("renv/activate.R")
library(bigsnpr)
library(fs)
library(EigenH5)
library(ldmap)
library(dplyr)
library(purrr)

iff <- snakemake@input[["bedlist"]]
off <- snakemake@output[["h5"]]
bed2h5 <- function(inputf, outputf){
  td <- tempdir()
  pext <- fs::path_file(inputf)
  pext <- fs::path_ext_remove(pext)
  tbf <- fs::path(td, pext)
  rds <- snp_readBed(inputf, backingfile = tbf)

  bs <- bigsnpr::snp_attach(rds)
  tibble::as_tibble(bs$map) %>%
    dplyr::mutate(index = 1:dplyr::n()) %>%
    compact_snp_struct(chrom =  "chromosome",
                       pos =  "physical.pos",
                       ref = "allele2",
                       alt = "allele1",
                       remove = TRUE) %>%
    write_df_h5(outputf,"snp")

  bsx <- bs$genotypes
  bsxm <- bsx[,]
  write_matrix_h5(bsxm, outputf, "dosage")
  return(outputf)
}

bed2h5(iff,off)
