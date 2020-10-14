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

classified_strains <- ce328_strains %>%
  dplyr::mutate(strain_set = ifelse(V1 == "XZ1516", "ancestor", 
                                    ifelse(V1 %in% swept_isotypes, "swept", "divergent"))) %>%
  dplyr::rename(strain = V1)

write.table(classified_strains, "classified_strains.tsv", quote=F,row.names = F, col.names = T, sep = "\t")

# tidyr::gather(sweeps, chrom, frac_sweep, -isotype) %>%
#   ggplot()+
#   aes(x = frac_sweep)+
#   geom_histogram(binwidth = 0.05)+
#   geom_vline(aes(xintercept=0.3), color = "red")+
#   facet_grid(chrom~.)
