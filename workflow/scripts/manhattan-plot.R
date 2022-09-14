


# redirect output and messages/errors to the log
log <- file(snakemake@log[[1]], open="wt")
sink(log, type = "output")
sink(log, type = "message")

library(tidyverse)
library(vroom)

sgfile <- snakemake@input$sg
lrt_file <- snakemake@input$lrt
outfile <- snakemake@output$mh_plot

#sgfile <- ".test/resources/scaffold_groups.tsv"
#lrt_file <- "results/bcf_yukon/filt_snps05/hasSDY/thin_0_0/do_asso/maf_0.05/age_cohort_4PCs/all-scaff-groups.lrt0.gz"
neg_log10_cutoff <- 0


# first, get the cumulative position of the start of each chromosome
sg <- read_tsv(sgfile) %>%
  mutate(
    cstart = 1 + lag(cumsum(stop), default = 0),
    cend = cstart + stop - 1
  )

# read in the LRTs and process them
r <- vroom(lrt_file) %>%
  filter(LRT != -999) %>% # remove missing ones
  mutate(neg_log10_p = -log10(pchisq(LRT, df = 1, lower.tail = FALSE))) %>%
  filter(neg_log10_p > neg_log10_cutoff) %>% 
  left_join(
    .,
    sg %>% select(chrom, mh_label, angsd_chrom, cstart),
    by = c("Chromosome" = "angsd_chrom")
  ) %>%
  mutate(xpos = Position + cstart - 1) %>%
  mutate(  # put each successive chromosome/scaff group into a separate group (a or b)
    cnf = factor(mh_label, levels = unique(mh_label)),
    cgroup = c("a", "b")[1 + as.integer(cnf) %% 2]
  )
  

# finally, we are going to want to get a tibble that tells
# us where the chromosome labels should go
cmids <- sg %>%
  group_by(mh_label) %>%
  summarise(xmid = (min(cstart) + max(cend)) / 2 )


# here we make the plot
g <- ggplot(r) +
  geom_point(aes(x = xpos / 1e6, y = neg_log10_p, colour = cgroup), size = 0.5) +
  theme_bw() +
  scale_colour_manual(values = c(a = "skyblue", b = "gray80")) +
  xlab("Megabases") +
  ylab("Negative log10 of p-value") +
  geom_text(
    data = cmids,
    mapping = aes(x = xmid / 1e6, y = 1, label = mh_label),
    size = 2.8,
    angle = 90,
    hjust = 0,
    vjust = 0.5
  ) +
  theme(legend.position = "none")


ggsave(g, filename = outfile, width = 10, height = 3)
