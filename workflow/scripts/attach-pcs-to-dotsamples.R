



# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)

samples_file <- snakemake@input$samples
cov_file <- snakemake@input$cov
dot_infile <- snakemake@input$dots
npc <- snakemake@params$num_pcs
dot_outfile <- snakemake@output$dots


#samples_file <- ".test/sample_subsets/all-fish.txt"
#cov_file <- "~/Documents/projects/yukon-chinookomes/BigAssocResults/pcangsd/out.cov"
#dot_infile <- ".test/data/dot-samples-all.tsv"

samples <- read_tsv(samples_file, col_names = "ID")
cov_mat <- matrix(scan(cov_file, what = numeric()), ncol = nrow(samples), byrow=TRUE)
indots <- read_tsv(dot_infile)

eig <- eigen(cov_mat)
colnames(eig$vectors) <- sprintf("PC-%d", 1:ncol(eig$vectors))

pc_tib <- as_tibble(eig$vectors[,1:npc]) %>%
  mutate(ID = samples$ID, .before = `PC-1`)

# now, add that first row to it with the variable types
pc_tib2 <- rbind(c("0", rep("C", npc)), pc_tib)

# then join things
outdots <- left_join(indots, pc_tib2, by = "ID")


# and finally, write that out
write_tsv(outdots, file = dot_outfile)

