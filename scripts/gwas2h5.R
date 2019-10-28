library(EigenH5)
library(readr)
library(ldmap)


mc <- cols(
    id = col_skip(),
    chrom = col_integer(),
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

snp_f <- snakemake@input[["snp_f"]]
input_f <- snakemake@input[["input"]]
output_f <- snakemake@output[["output"]]


callback_fun <- function(df, filename, datapath, ...){
    write_df_h5(df = dplyr::slice(dplyr::mutate(df,
                                   ref = fast_str2ascii(A2),
                                   alt = fast_str2ascii(A1),
                                   snp_struct = new_ldmap_snp(chr, pos, ref, alt),
                                   rsid = NA_integer_), rank.ldmap_snp(snp_struct)),
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
         skip = 1L,
         callback_fun = callback_fun,
         col_types = mc,
         progress = TRUE,
         chunk_size = 15000)

chrom_vec <- read_vector_h5v(output_f, "snp/chrom", i = integer())
chrom_df <- rle2offset(chrom_vec) %>%
    rename(chrom = value)
write_df_h5(chrom_df,output_f,"chrom_offset")


gwas_offset <- filter(gwas_offset, !is.na(chrom)) %>%
    rename(offset_gwas = offset,
           datasize_gwas = datasize)


snpdb_offset <- read_df_h5(snp_f, "chrom_offset") %>%
    filter(!is.na(chrom)) %>%
    rename(offset_db = offset,
           datasize_db = datasize)

gwas_db_df <- inner_join(gwas_offset, snpdb_offset) %>%
    pwalk(
        function(chrom,offset_db,datasize_db,offset_gwas,datasize_gwas) {
            rseq_gwas <- seq(offset_gwas + 1,length.out = datasize_gwas)
            rseq_db <- seq(offset_db + 1,length.out = datasize_db)
            db_df <- read_df_h5(snp_f,
                                datapath = "snp",
                                subcols = c("snp_struct", "rsid"),
                                subset =  rseq_db)
            gw_df <- read_df_h5(output_f,
                                datapath = "snp",
                                subcols = c("snp_struct"),
                                subset =  rseq_gwas)
            ret_gl <- left_join(gw_df,db_df)
            stopifnot(ret_gl$snp_struct ==gw_df$snp_struct)
            write_df_h5(ret_gl, datapath = "snp", subset = rseq_gwas)})






}
