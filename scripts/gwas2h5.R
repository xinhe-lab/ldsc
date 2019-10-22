library(EigenH5)
library(readr)


mc <- cols(
    id = col_skip(),
    chr = col_integer(),
    pos = col_double(),
    A1 = col_character(),
    A2 = col_character(),
    N = col_double(),
    freq = col_double(),
    beta = col_double(),
    se = col_double(),
    pval = col_double(),
    Q = col_double(),
    het = col_double(),
    N.local = col_double(),
    freq.local = col_double(),
    beta.local = col_double(),
    se.local = col_double(),
    pval.local = col_double(),
    N.23andMe = col_double(),
    freq.23andMe = col_double(),
    beta.23andMe = col_double(),
    se.23andMe = col_double(),
    pval.23andMe = col_double()
)

input_f <- snakemake@input[["input"]]
output_f <- snakemake@output[["output"]]



callback_fun <- function(df, filename, datapath, ...){
    write_df_h5(df = dplyr::mutate(df,
                                   ref = fast_str2ascii(A2),
                                   alt = fast_str2ascii(A1),
                                   snp_struct = fast_snp_pos_struct(chr, pos, ref, alt),
                                   rsid = NA_integer_),
                filename = filename, datapath = datapath, ... = ...)
}

stopifnot(!is.null(input_f),
          !is.null(output_f),
          file.exists(input_f),
          !file.exists(output_f))

delim2h5(input_f,
         output_file = output_f,
         h5_args = list(datapath = "snp"),
         delim = "\t",
         col_names = names(mc$cols),
         callback_fun = callback_fun,
         col_types = mc,
         progress = TRUE,
         chunk_size = 15000)
