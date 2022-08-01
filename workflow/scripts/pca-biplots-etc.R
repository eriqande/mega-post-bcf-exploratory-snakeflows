

#save.image(file = "/tmp/Rdata")

# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")


#results/pcangsd/Mondo-454-fish-maf-0.05/thin_6_1/maf_0.04/out.cov
#results/bcf/Mondo-454-fish-maf-0.05/thin_6_1/info/samples.txt
samps_file <- snakemake@input$samps
cov_file <- snakemake@input$cov
npp <- 4  # number of PCs to retain for plotting



library(tidyverse)

samples <- read_tsv(samps_file, col_names = "sample")
cov_mat <- matrix(scan(cov_file, what = numeric()), ncol = nrow(samples), byrow=TRUE)

eig <- eigen(cov_mat)
colnames(eig$vectors) <- sprintf("PC-%04d", 1:ncol(eig$vectors))


pca_tib <- tibble::as_tibble(eig$vectors[,1:npp]) %>%
  mutate(sample = samples$sample) %>%
  select(sample, everything())

pca_long <- pca_tib %>%
  tidyr::pivot_longer(cols = -sample, names_to = "PC", values_to = "val")

pca_pairs <- expand_grid(
  sample = pca_tib$sample,
  PCx = sprintf("PC-%04d", 1:npp),
  PCy = sprintf("PC-%04d", 1:npp)
) %>%
  dplyr::left_join(., pca_long, by = c("sample", "PCx" = "PC")) %>%
  dplyr::rename(val_x = val) %>%
  dplyr::left_join(pca_long, by = c("sample", "PCy" = "PC")) %>%
  dplyr::rename(val_y = val) 


# make the biplot
bp <- ggplot(
  data = pca_pairs, 
  mapping = aes(
    x = val_x, 
    y = val_y, 
  )
) + 
  geom_point(shape=21, fill=NA) +
  facet_grid(PCy ~ PCx) +
  theme_bw()
