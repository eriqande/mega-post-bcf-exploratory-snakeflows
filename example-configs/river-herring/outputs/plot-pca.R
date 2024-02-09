
library(tidyverse)


cov <- read_table("example-configs/river-herring/outputs/out.cov", col_names = FALSE)


eig <- eigen(as.matrix(cov))


plot(eig$vectors[,1], eig$vectors[,2])
