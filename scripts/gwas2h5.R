library(dplyr)
library(purrr)
library(readr)


  library(EigenH5)
  library(readr)
  library(ldmap)


  ## mc <- cols(
  ##     rsid = col_character(),
  ##     chrom = col_integer(),
  ##     pos = col_double(),
  ##     A1 = col_character(),
  ##     A2 = col_character(),
  ##     N = col_double(),
  ##     freq = col_double(),
  ##     beta = col_double(),
  ##     se = col_double(),
  ##     pval = col_double(),
  ##     Q = col_double(),
  ##     het = col_double(),
  ##     N.local = col_double(),
  ##     freq.local = col_double(),
  ##     beta.local = col_double(),
  ##     se.local = col_double(),
  ##     pval.local = col_double(),
  ##     N.23andMe = col_double(),
  ##     freq.23andMe = col_double(),
  ##     beta.23andMe = col_double(),
  ##     se.23andMe = col_double(),
  ##     pval.23andMe = col_double()
  ## )


  input_f <- snakemake@input[["inputf"]]
  output_f <- snakemake@output[["outputf"]]
  paramf <- snakemake@input[["paramf"]]
  stopifnot(!is.null(paramf))
  source(paramf)


  callback_fun <- function(df, filename, datapath, ...){
    write_df_h5(
      df = dplyr::slice(
                    dplyr::mutate(df,
                                  ref = fast_str2ascii(A2),
                                  alt = fast_str2ascii(A1),
                                  snp_struct =
                                    new_ldmap_snp(chrom, pos, ref, alt),
                                  rsid = fast_str2int(rsid, prefix = "rs"),
                                  ),
                    rank.ldmap_snp(snp_struct)),
      filename = filename, datapath = datapath, ... = ...)
  }

  stopifnot(!is.null(input_f),
            !is.null(output_f),
            file.exists(input_f),
            !file.exists(output_f))

  delim2h5(input_f,
           output_file = output_f,
           h5_args = list(datapath = "snp"),
           delim = data_delim,
           col_names = names(mc$cols),
           skip = 1L,
           callback_fun = callback_fun,
           col_types = mc,
           progress = TRUE,
           chunk_size = 150000)

  chrom_vec <- read_vector_h5v(output_f, "snp/chrom", i = integer())
  chrom_df <- rle2offset(chrom_vec) %>%
      dplyr::rename(chrom = value)
  write_df_h5(chrom_df,output_f,"chrom_offset")
