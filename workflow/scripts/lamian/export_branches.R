library(tidyverse)

input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

x <- readRDS(input_file)

branches <- names(x[["order"]])
df <- lapply(branches, function(branch) {
  cells <- x[["order"]][[branch]]
  tibble(cell_id = cells, pseudotime = seq_along(cells), branch = branch)
}) %>%
  bind_rows()

write_csv(df, output_file)
