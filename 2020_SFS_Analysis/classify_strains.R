library(tidyverse)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

sweeps <- data.table::fread("sweep_summary.tsv") %>%
  dplyr::select(-(I:X), -swept_chroms, -(II_hapshare:III_hapshare))

ce328_strains <- data.table::fread("ce328_strains.txt", header = F)

swept_isotypes <- tidyr::gather(sweeps, chrom, frac_sweep, -isotype) %>%
  dplyr::filter(frac_sweep >= 0.3, isotype %in% ce328_strains$V1) %>%
  dplyr::distinct(isotype) %>%
  dplyr::arrange(isotype) %>%
  dplyr::pull(isotype)
paste(swept_isotypes, collapse = ",")

divergent_isotypes <- tidyr::gather(sweeps, chrom, frac_sweep, -isotype) %>%
  # dplyr::filter(frac_sweep >=0.3) %>%
  dplyr::distinct(isotype) %>%
  dplyr::filter(!isotype%in%swept_isotypes, isotype %in% ce328_strains$V1) %>%
  dplyr::arrange(isotype) %>%
  dplyr::pull(isotype)

paste(divergent_isotypes, collapse = ",")


tidyr::gather(sweeps, chrom, frac_sweep, -isotype) %>%
  ggplot()+
  aes(x = frac_sweep)+
  geom_histogram(binwidth = 0.05)+
  geom_vline(aes(xintercept=0.3), color = "red")+
  facet_grid(chrom~.)
