library(tidyverse)
library(EigenH5)
library(readr)
library(ldmap)
## load("~/Dropbox/Repos/ldsc/workflow/tf.RData")

input_f <- snakemake@input[["inputf"]]
index_f <-  snakemake@input[["indexf"]]
chrom <- snakemake@params[["chrom"]]
stopifnot(!is.null(chrom))
schrom <- as.integer(chrom)
output_f <- snakemake@output[["outputf"]]


ind_spec <- cols(
  chrom = col_integer(),
  pos = col_double(),
  rsid = col_integer(),
  alleles_as_ambig = col_character()
)

gwas_type <- if_else(
  is.null(snakemake@params[["gwas_t"]]),
  "",
  paste0(".", snakemake@params[["gwas_t"]])
)


beta_col <- glue::glue("beta{gwas_type}")
se_col <- glue::glue("se{gwas_type}")
N_col <- glue::glue("N{gwas_type}")
P_col <- glue::glue("pval{gwas_type}")

sel_cols <- c("snp_struct",
              beta_col,
              "A1",
              "A2",
              se_col,
              N_col,
              P_col)

sel_cols <- stringr::str_replace(
                       sel_cols,
                       "\\.$",
                       "")

index_df <- vroom::vroom(
                     index_f,
                     delim = "\t",
                     col_names = names(ind_spec$cols),
                     col_types = ind_spec,
                     skip = 1L
                   )

chrom_df <- read_tibble_h5(input_f, "chrom_offset", list()) %>%
  filter(chrom == schrom) %>% mutate(offset = as.integer(offset), datasize = as.integer(datasize)) %>% 
  arrange(offset)

jdf <- pmap_dfr(chrom_df, function(chrom, datasize, offset) {
  subset_l <- seq(offset + 1, length.out = datasize)
  input_i <- EigenH5::read_df_h5(filename = input_f,
                          datapath = "snp",
                            subcols = sel_cols,
                            subset = subset_l) %>%
    mutate(subset = (1:n()) + offset)

    inner_join(index_df,  bind_cols(input_i,ldmap::ldmap_snp_2_dataframe(input_i$snp_struct)))
})

                                        #%>% mutate(snp_struct = as_ldmap_snp(snp_struct))  %>%
stopifnot(all(jdf$chrom == schrom))

jdf  %>% rename(beta =  {{beta_col}},
                se =  {{se_col}},
                N =  {{N_col}}) %>%
  dplyr::distinct(rsid, .keep_all = TRUE) %>% 
  dplyr::transmute(SNP = paste0("rs",rsid), N = N, Z = beta / se, A1 = A1, A2 = A2,P=pval) %>%
  vroom::vroom_write(output_f,delim = "\t")
