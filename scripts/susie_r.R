source("renv/activate.R")
  library(ldmap)
  library(EigenH5)
  library(susieR)
  library(dplyr)
  library(purrr)


  align_reference <- function(gwas_df, reference_file, remove_missing = TRUE, read_map_fn =  identity){
  bsmap <- read_map_fn(reference_file)
  
  sc_df <- ldmap:::snp_cols(gwas_df)
  ret_df <- match_ref_panel(gwas_df,bsmap$snp_struct)%>%
  dplyr::select(-{{sc_df}}) %>% dplyr::rename({{sc_df}} := match)
  if(remove_missing){
  ret_df <- dplyr::filter(ret_df,!is.na(index))
  }else{
  stopifnot(all(!is.na(ret_df$index)))
  }
  return(ret_df)
  }

  read_map_fn <- function(x) read_df_h5(x,"snp")
  ld_regions <- readRDS(snakemake@input[["ldgf"]]) %>%
    filter(region_id == as.integer(snakemake@params[["region_id"]]))
  stopifnot( nrow(ld_regions) == 1L )
  p <- dim_h5(snakemake@input[["inputf"]],"snp/snp_struct")
  gwas_df <- read_df_h5(snakemake@input[["inputf"]],
                        datapath = "snp",
                        subcols = c("snp_struct", "beta", "se","N"),
                        offset =  ld_regions$offset,
                        datasize = ld_regions$datasize)  %>% 
    align_reference(snakemake@input[["ldf"]],read_map_fn = read_map_fn) %>% 
arrange(index) 

  R <- read_matrix_h5v(snakemake@input[["ldf"]],"R",gwas_df$index, gwas_df$index)


  gwas_df <- readRDS(snakemake@input[["priorf"]]) %>% 
    select(snp_struct=SNP,prior) %>% 
    inner_join(gwas_df)
    h <- 0.1
    h_p <- h/p
    prior_v <- h_p*nrow(gwas_df)
    sres <- susie_suff_stat(
      bhat = gwas_df$beta,
      shat=gwas_df$se,
      R = R,
      n = max(gwas_df$N),
      scaled_prior_variance=prior_v,
      prior_weights = gwas_df$prior,
      L = 1,
      null_weight=NULL
    )
    saveRDS(sres, snakemake@output[["outputf"]])
