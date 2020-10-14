#!/usr/bin/env Rscript

# args:
# 1 - SFS_INPUT.tsv file
# 2 - "no_indel" to exclude indels
# 3 - invariant_site_counts.tsv
# 4 - number of samples

# args <- c("SFS_INPUT.tsv", "no_indel", "invariant_site_by_region.tsv", "52")

args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(glue)

multi_dfe_out <- function(df, fname) {
  writeLines(
    c(
      nrow(df) - 1,
      paste(df$Selected, collapse =" "),
      paste(df$Neutral, collapse =" ")
    ),
    con = glue::glue(fname), 
    sep = "\n"
  )
}


if(args[2] == "no_indel"){
  print("Excluding Indels")
  sfs_df <- data.table::fread(args[1]) %>%
    dplyr::mutate(DERIVED_AF = ifelse(AA==REF, round(AF, digits = 5), 
                                      ifelse(AA == ALT, round(1-AF, digits = 5), NA))) %>%
    dplyr::filter(nchar(REF) == nchar(ALT)) %>%
    dplyr::filter(DERIVED_AF > 0)
} else {
  sfs_df <- data.table::fread(args[1]) %>%
    dplyr::mutate(DERIVED_AF = ifelse(AA==REF, round(AF, digits = 5), 
                                      ifelse(AA == ALT, round(1-AF, digits = 5), NA)))
}

# calculate theta for 4fold sites across genome
# extract 4FOLD sites
theta_df <- dplyr::filter(sfs_df, `4FOLD` == "4FOLD") %>%
  dplyr::group_by(GENOMIC_REGION, CHROM) %>%
  dplyr::summarise(theta = n()/(sum(1/1:(as.numeric(args[4])-1)))) %>%
  dplyr::arrange(desc(theta))
# segregating 4fold sites
seg_sites <- nrow(theta_df)
# denominator of formula
an <- (sum(1/1:(as.numeric(args[4])-1)))
# calculate genome-wide theta
seg_sites/an
# 18948.34

af_df <- data.frame(DERIVED_AF = round(seq(0,1,by = 1/as.numeric(args[4])), digits = 5))


sfs_out <- function(nclass, sclass, region, chrom){
  
  gl_sclass <- paste(sclass,collapse="_")
  
  save_name <- glue::glue("{nclass}_{gl_sclass}_{region}_CHROM-{chrom}")
  
  neutral <- sfs_df %>%
  {if (nclass == "4FOLD") dplyr::filter(.,`4FOLD` == nclass) else .} %>%
  {if (nclass == "intergenic_region") dplyr::filter(.,`SNPEFF_REGION` == nclass) else .} %>%
  {if (region != "GENOME") dplyr::filter(., GENOMIC_REGION == region) else .} %>%
  {if (chrom != "GENOME") dplyr::filter(., CHROM %in% chrom) else .} %>%
    dplyr::group_by(DERIVED_AF) %>%
    dplyr::summarise(Neutral = n()) %>%
    dplyr::left_join(af_df, .)
  
  if(length(sclass) == 1){
    
    selected <- sfs_df %>%
    {if (region != "GENOME") dplyr::filter(., GENOMIC_REGION == region) else .} %>%
    {if (sclass == "0FOLD") dplyr::filter(., `0FOLD` == sclass) else .} %>%
    {if (sclass != "0FOLD" & !("intergenic_region" %in% sclass)) dplyr::filter(., EFFECT == sclass) else .} %>%
    {if (chrom != "GENOME") dplyr::filter(., CHROM %in% chrom) else .} %>%
    {if ("intergenic_region" %in% sclass) dplyr::filter(., `SNPEFF_REGION` %in% sclass) else .} %>%
      dplyr::group_by(DERIVED_AF) %>%
      dplyr::summarise(Selected = n()) %>%
      dplyr::left_join(neutral, .)
    
  } else if (length(sclass) > 1) {
    
    if(grepl("0FOLD", sort(sclass)[1])){
      sclass_1 <- sort(sclass)[1]
      sclass_2 <- sort(sclass)[2:length(sclass)]
      
      selected <- sfs_df %>%
      {if (region != "GENOME") dplyr::filter(., GENOMIC_REGION == region) else .} %>%
      {if (chrom != "GENOME") dplyr::filter(., CHROM %in% chrom) else .} %>%
      {if ("intergenic_region" %in% sclass_2) dplyr::filter(., `0FOLD` == sclass_1 | EFFECT %in% sclass_2 | `SNPEFF_REGION` %in% sclass_2) else .} %>%
      {if (!("intergenic_region" %in% sclass_2)) dplyr::filter(., `0FOLD` == sclass_1 | EFFECT %in% sclass_2) else .} %>%
        dplyr::group_by(DERIVED_AF) %>%
        dplyr::summarise(Selected = n()) %>%
        dplyr::left_join(neutral, .)
    } else {
      selected <- sfs_df %>%
      {if (region != "GENOME") dplyr::filter(., GENOMIC_REGION == region) else .} %>%
      {if (chrom != "GENOME") dplyr::filter(., CHROM %in% chrom) else .} %>%
      {if ("intergenic_region" %in% sclass) dplyr::filter(., EFFECT %in% sclass | `SNPEFF_REGION` %in% sclass) else .} %>%
      {if (!("intergenic_region" %in% sclass)) dplyr::filter(., EFFECT %in% sort(sclass)) else .} %>%
        dplyr::group_by(DERIVED_AF) %>%
        dplyr::summarise(Selected = n()) %>%
        dplyr::left_join(neutral, .)
    }
    
  }
  
  selected[is.na(selected)] <- 0
  
  return(list(save_name, selected))
}

invariant_sites <- data.table::fread(args[3])

sel_class <- c("HIGH","MODERATE","LOW","0FOLD", "intergenic_region")

for(cr in c("GENOME")) {
  for(grab_s in 1:length(sel_class)){
    combos <- combn(sel_class, m = grab_s)
    for(sc in 1:ncol(combos)){
      s_groups <- as.vector(combos[,sc])
      for(rg in c("GENOME","ARM","CENTER")){
        print(s_groups)
        spectra <- sfs_out(nclass = "4FOLD", sclass = s_groups, region = rg, chrom = cr)
        
        if(grepl("0FOLD", spectra[[1]])) {
          invariant_spectra <- spectra[[2]]
          
          temp_invar_nsites <- invariant_sites %>%
            dplyr::filter(GENOMIC_REGION == rg, variant_fold == 4) %>%
            dplyr::pull(ct)
          
          temp_invar_ssites <- invariant_sites %>%
            dplyr::filter(GENOMIC_REGION == rg, variant_fold == 0) %>%
            dplyr::pull(ct)
          
          invariant_spectra$Neutral[1] <- temp_invar_nsites
          invariant_spectra$Selected[1] <- temp_invar_ssites
          
          multi_dfe_out(df = invariant_spectra,
                        fname = glue::glue("{spectra[[1]]}_INVARIANT.sfs"))
        }
        
        # save file
        multi_dfe_out(df = spectra[[2]],
                      fname = glue::glue("{spectra[[1]]}.sfs"))
        
        sfs_plot <- spectra[[2]] %>%
          dplyr::filter(DERIVED_AF != 0) %>%
          tidyr::gather(CLASS, COUNTS, -DERIVED_AF) %>%
          dplyr::group_by(CLASS) %>%
          dplyr::mutate(frq = COUNTS/sum(COUNTS)) %>%
          ggplot()+
          aes(x = DERIVED_AF, y = frq, color = CLASS)+
          geom_line()+
          theme_bw(15)+
          labs(x = "Derived AF", y = "Frequency")
        
        ggsave(sfs_plot, filename = glue::glue("{spectra[[1]]}.pdf"), height = 4, width = 10)
      }
    }
    
  }
}
