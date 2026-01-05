library(tidyverse)
library(Lamian)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

res <-
  input_file %>%
  readRDS()

set.seed(42)
result <- evaluate_uncertainty(res, n.permute = 1000)
saveRDS(result, output_file)
