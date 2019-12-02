library(stringr)
library(Homo.sapiens)
library(ggbio)
library(RColorBrewer)
library(drake)
library(purrr)
library(dplyr)
library(biovizBase)
library(Homo.sapiens)

txdb <- TxDb(Homo.sapiens)

chrp <- snakemake@params[["chrom"]]
start <- snakemake@params[["start"]]
end <- snakemake@params[["end"]]

gr <- GenomicRanges::GRanges(seqnames = chrp,start = start,end = end)

suppressMessages(tg_df <-
                   OrganismDbi::selectByRanges(
                                  x = Homo.sapiens,
                                  ranges = gr,
                                  columns = c("SYMBOL", "TXNAME")) %>%
                   as_tibble() %>%
                   dplyr::select(tx_name = TXNAME, symbol = SYMBOL) %>%
                   tidyr::unnest(cols = c(tx_name, symbol)))

suppressMessages(gr.txdb <-
                   crunch(txdb, which = gr) %>%
                   as_tibble() %>%
                   mutate(tx_name = as.character(tx_name)) %>%
                   inner_join(tg_df) %>%
                   plyranges::as_granges() %>%
                   split(.$symbol))


saveRDS(gr.txdb,snakemake@output[["outf"]])
