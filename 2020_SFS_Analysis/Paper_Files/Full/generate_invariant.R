library("tidyverse")
library("seqinr")

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

VARIANT_FILE="SFS_INPUT.tsv"
SAVE_FILE="invariant_site_by_region.tsv"
# HTA phenotyped strains
DIVERGED_MODIFIER=""
# Non-swept strains
DIVERGED_MODIFIER="_diverged"

VARIANT_FILE=glue::glue("SFS_INPUT{DIVERGED_MODIFIER}.tsv")
SAVE_FILE=glue::glue("invariant_site_by_region{DIVERGED_MODIFIER}.tsv")
DEGENERACY_FILE=glue::glue("gene_degeneracy{DIVERGED_MODIFIER}.Rda")

# load fasta files
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))
nt_fa <- seqinr::read.fasta(file = "WS261_spliced_cds.fa")

# load variant information
snv <- data.table::fread(VARIANT_FILE) %>%
  dplyr::filter(grepl("WBGene", WBGeneID), grepl("p.", aa_CHANGE)) %>%
  dplyr::mutate(DERIVED_AF = ifelse(AA == REF, round(AF, digits = 5), 
                                    ifelse(AA == ALT, round(1-AF, digits = 5), NA))) %>%
  dplyr::filter(nchar(REF) == 1, nchar(ALT) == 1, DERIVED_AF > 0) %>% # snps only and varies in the population
  tidyr::separate(CDS_POS, into = c("nt_POS", "gene_LENGTH"), sep = "/", convert = T)

# define codon degeneracy
degen_df <- data.frame(deg1 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                       deg2 = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                       deg3 = c(0, 4, 4, 0, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 0, 4, 0, 4, 0, 4),
                       aa   = c("M", "S", "R", "I", "H", "N", "G", "E", "D", "V", "Y", "K", "A", "W", "Q", "L", "F", "P", "C", "T"))

# function to chunk nucleotide vector into codons
chunk2 <- function(x,n) data.frame(t(data.frame(split(x, cut(seq_along(x), n, labels = FALSE)))))

# loop through all genes
degeneracy_list <- list()
for(gene in 1:length(nt_fa)) {
  
  # extract gene sequence
  tsequ <- nt_fa[[gene]]
  
  #  extract gene information
  tname <- strsplit(strsplit(getAnnot(tsequ), split = " ")[[1]][1], split = ">")[[1]][2]
  chrom <- strsplit(strsplit(strsplit(getAnnot(tsequ), split = " ")[[1]][2], split = "\\(")[[1]], split = ":")[[1]][2]
  sgene <- as.numeric(strsplit(strsplit(strsplit(getAnnot(tsequ), split = " ")[[1]][2], split = "\\)")[[1]][2], split = "-")[[1]][1])
  egene <- as.numeric(strsplit(strsplit(strsplit(getAnnot(tsequ), split = " ")[[1]][2], split = "\\)")[[1]][2], split = "-")[[1]][2])
  gdire <- strsplit(strsplit(strsplit(getAnnot(tsequ), split = " ")[[1]][2], split = "\\)")[[1]][1], split = "\\(")[[1]][2]
  
  # extract gene variant information
  gsnv <- snv %>%
    dplyr::filter(CHROM == chrom, POS < egene & POS > sgene) %>%
    dplyr::arrange(nt_POS)
  
  # deals with genes with no variants
  if(nrow(gsnv) == 0) {
    gene_length <- length(getSequence(tsequ))
  } else {
    gene_length <- unique(gsnv$gene_LENGTH)
    # deals with overlapping genes; grabs gene length corresponding to gene of interest: tname
    if(length(gene_length) > 1) {
      gene_length <- gene_length[which(gene_length == length(getSequence(tsequ)))]
      # deals with overlapping genes that are not the same length as the snpeff annotated gene
      if(length(gene_length) == 0){
        gene_length <- unique(gsnv$gene_LENGTH)[1]
      } else {
        gsnv <- gsnv %>%
          dplyr::filter(gene_LENGTH == gene_length)
      }
    }
  }
  
  # conditional deals with genes with multiple transcripts, only extract transcript corresponding to length in SNPEFF output
  if(length(getSequence(tsequ)) == gene_length){
    # replace variable sequence in gene with a "-"
    modified_sequence <- replace(getSequence(tsequ), gsnv$nt_POS, rep("-", times = length(gsnv$nt_POS)))
    
    # make gene data.frame
    gene_df <- data.frame(chunk2(modified_sequence, length(getTrans(modified_sequence))),
                          aa = getTrans(tsequ),
                          aa_mod = getTrans(modified_sequence)) %>%
      dplyr::rename(codon1 = X1, codon2 = X2, codon3 = X3) %>%
      dplyr::mutate(aa_pos = 1:n())
    
    # combine gene sequence with degeneracy table
    degeneracy_list[[gene]] <- dplyr::left_join(gene_df, degen_df, by = "aa") %>%
      tidyr::gather(codon, nt, -(aa:deg3)) %>%
      tidyr::gather(deg_pos, fold, -aa, -aa_mod, -aa_pos, -codon, -nt) %>%
      dplyr::mutate(codon = as.numeric(gsub("codon","", codon)),
                    deg_pos = as.numeric(gsub("deg","", deg_pos))) %>%
      dplyr::filter(codon == deg_pos) %>% # filter data frame to only keep degeneracy information for corresponding codon position
      dplyr::arrange(aa_pos) %>%
      dplyr::mutate(variant_fold = ifelse(nt == "-", NA, fold),
                    transcript_name = tname,
                    CHROM = chrom,
                    gene_start = sgene,
                    gene_end = egene,
                    gene_strand = gdire) # append gene information 
  } else {
    degeneracy_list[[gene]] <- NULL
  }
  
}

all_gene_degen_df <- dplyr::bind_rows(degeneracy_list[!sapply(degeneracy_list, is.null)] ) %>%
  tidyr::separate(transcript_name, into = c("cosmid", "transcript", "altsplice"), sep = "\\.") %>%
  dplyr::distinct(aa, aa_pos, codon, cosmid, transcript, .keep_all = T)

save(all_gene_degen_df, file = DEGENERACY_FILE)

load(DEGENERACY_FILE)

chr_regions <- readr::read_tsv("ARMS_CENTERS.bed.gz", col_names = F) %>%
  dplyr::rename(CHROM = X1, Spos = X2, Epos = X3, CLASS = X4)

append_region <- list()

for(region in 1:nrow(chr_regions)){
  temp_region <- chr_regions[region,]
  
  append_region[[region]] <- all_gene_degen_df %>%
    dplyr::filter(CHROM == temp_region$CHROM, gene_start < temp_region$Epos, gene_start > temp_region$Spos) %>%
    dplyr::mutate(GENOMIC_REGION = temp_region$CLASS) 
}

region_variant_fold <- dplyr::bind_rows(append_region) %>%
  dplyr::filter(nt != "-", aa != "*") %>%
  dplyr::group_by(GENOMIC_REGION, variant_fold) %>%
  dplyr::summarise(ct = n())

genome_invariant_fold <- all_gene_degen_df %>%
  dplyr::filter(nt != "-", aa != "*") %>%
  dplyr::group_by(variant_fold) %>%
  dplyr::summarise(ct = n()) %>%
  dplyr::mutate(GENOMIC_REGION = "GENOME") %>%
  dplyr::select(GENOMIC_REGION, variant_fold, ct)

invariant_sites <- dplyr::bind_rows(genome_invariant_fold, region_variant_fold)

write.table(invariant_sites, SAVE_FILE, col.names = T, row.names = F, quote = F, sep = "\t")
