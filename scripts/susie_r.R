source("renv/activate.R")
  library(ldmap)
  library(dplyr)
  library(purrr)
  library(vroom)

  ld_regions <- readRDS(snakemake@input[["ldgf"]]) %>%
    filter(region_id == as.integer(snakemake@params[["region_id"]]))
  stopifnot( nrow(ld_regions) == 1L )
  gwas_df <- read_df_h5(snakemake@input[["inputf"]],
                        datapath = "snp",
                        subcols = c("snp_struct", "beta", "se","N"),
                        offset =  ld_regions$offset,
                        datasize = ld_regions$datasize) %>% 
    align_reference(snakemake@input[["rdsf"]]) %>% mutate(snp_struct = as.character(snp_struct))
  fmt <- gwas_df$snp_struct
  R <- readRDS(snakemake@input[["ldf"]])[fmt,fmt]


  gwas_df <- vroom(snakemake@input[["priorf"]],delim = " ",col_names = c("snp_struct","prior","region_id")) %>% inner_join(gwas_df)


  saveRDS(susie_bhat(
    bhat = gwas_df$beta,
    shat = gwas_df$se,
    R = R,
    n = max(gwas_df$N), prior_weights = gwas_df$prior,
    L = 1,
    estimate_residual_variance = TRUE,
    estimate_prior_variance = FALSE
  ), snakemake@output[["outputf"]])

source("renv/activate.R")
saveRDS(purrr::map(snakemake@input,readRDS),snakemake@output[["outputf"]])
