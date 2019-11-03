library(vroom)

 input_f <- snakemake@input[["inputf"]]
 output <- snakemake@output[["outputf"]]

 vroom::vroom(input_f,delim="\t") %>% vroom_write(output,delim="\t")
