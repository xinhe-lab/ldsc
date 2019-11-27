library(readr)
library(dplyr)
library(purrr)
library(forcats)
library(ldmap)
## save.image("ppo.RData")
## stop("?")
## load("ppo.RData")
data(ldetect_EUR)


sumstat_h5f <- snakemake@input[["gwas_hf"]]
annof <- snakemake@input[["annot_f"]]
index_f <- snakemake@input[["bimf"]]

anno_n <- snakemake@params[["annot"]]
chromlist <- snakemake@params[["chroms"]]
outputf <- snakemake@output[["outputf"]]

region_l <- purrr::map(annof, ~read_bed(.x)$ldmap_range)
names(region_l) <- paste0(anno_n,"_d")
bc <- bim_cols(chrom = col_chromosome(prefix_chr=FALSE))
index_l <- purrr::walk(index_f, function(x) {
  fe <- file.exists(outputf)
  input_b <- read_plink_bim(x, cols = bc)$snp_struct
  snp_in_ranges(input_b, region_l) %>%
    rename(SNP = ldmap_snp) %>%
    write_delim(outputf, delim = " ", append = !fe)
  })
