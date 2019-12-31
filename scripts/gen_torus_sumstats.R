## stop()
      library(readr)
      library(dplyr)
      library(purrr)
      library(forcats)
      library(ldmap)
      library(EigenH5)

## stop()
      data(ldetect_EUR)


      sumstat_h5f <- snakemake@input[["inputf"]]
      snplist <- snakemake@input[["snplist"]]
      chromlist <- snakemake@params[["chroms"]]
      outputf <- snakemake@output[["outputf"]]

      chrom_df <- read_df_h5(sumstat_h5f, "chrom_offset") %>% 
        dplyr::slice(1:22) %>% 
        dplyr::mutate(offset = as.integer(offset),
                      datasize = as.integer(datasize)) %>%
        dplyr::arrange(offset)

    bc <- bim_cols(chrom=col_chromosome(prefix_chr=FALSE))
    index_l <- purrr::map(snplist, ~read_plink_bim(.x,cols = bc)$snp_struct)
  mutate(chrom_df,snplist_l = index_l) %>%
    pwalk(
      function(chrom, offset, datasize, snplist_l, ...) {
        fe <- file.exists(outputf)
        idf <- EigenH5::read_df_h5(
                          filename = sumstat_h5f,
                          datapath = "snp",
                          subcols = c("snp_struct", "beta", "se"),
                          offset = offset,
                          datasize = datasize) %>%
          match_ref_panel(snplist_l) %>%
          filter(!is.na(index)) %>% 
          dplyr::transmute(SNP = match,
                           locus = snp_in_range(SNP, ldetect_EUR),
                           `z-vals` =  beta / se )
        stopifnot(all(!is.na(idf$locus)))


        write_delim(idf,outputf, delim = " ",append = fe)
      })
