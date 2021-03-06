library(dplyr)
library(purrr)
library(readr)


  library(ldmap)
  library(EigenH5)

  input_down <- snakemake@input[["input_down"]]
  input_up <- snakemake@input[["input_up"]]

  outputf <- snakemake@output[["bedf"]]

  dcols <- cols(
    chrom = col_factor(paste0("chr", c(as.character(1:22), "X"))),
    start = col_integer(),
    end = col_integer())

  diff_df <- vroom::vroom(c(input_up, input_down),
                          delim = "\t",
                          col_names = c("chrom", "start", "end"),
                          col_types = dcols)
  new_ldmap_range(diff_df$chrom,
                  diff_df$start,
                  diff_df$end) %>%
    split_ldmap_range_overlap() %>%
    ldmap_range_2_data_frame() %>%
    vroom::vroom_write(outputf, delim = "\t", col_names = FALSE)
