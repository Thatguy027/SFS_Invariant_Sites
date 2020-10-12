############################################################################### GENERATE ANNOTATION FILES
# put yourself in annotation directory


bcftools filter -i N_MISSING=0 Ce328_From_2020GATK.vcf.gz |\
vcffixup - |\
bcftools filter -i 'AF>0' -Oz -o Ce328_From_2020GATK_noMissing.vcf.gz

tabix -p vcf Ce328_From_2020GATK_noMissing.vcf.gz


# make ancestor bed
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t[%TGT]\n' Ce328_From_2020GATK_noMissing.vcf.gz |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > ANC.bed.gz

tabix ANC.bed.gz

# parse snpeff
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' Ce328_From_2020GATK_noMissing.vcf.gz |\
awk -F"|" '$1=$1' OFS="\t" |\
awk '{print $1, $2, $3, $5}' OFS="\t" |\
awk '{if($4 == "") print $1, $2 = $2 - 1, $3, "intergenic_region";
	  else print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > SNPEFF_REGION.bed.gz

tabix SNPEFF_REGION.bed.gz

# make 4fold site annotation - 112045 sites
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' Ce328_From_2020GATK_noMissing.vcf.gz |\
grep "synonymous_variant" |\
grep "protein_coding" |\
grep -v "splice_" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' |\
awk '{if(($4 ~ "p.Ser" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Val" || $4 ~ "p.Leu" || $4 ~ "p.Gly" || $4 ~ "p.Arg") && $7 == 0) print $1, $2=$2-1, $3, "4FOLD"}' OFS="\t" |\
bgzip > 4FOLD_SITES.bed.gz

tabix 4FOLD_SITES.bed.gz


# make 0fold site annotation - 138537 sites
bcftools query --samples XZ1516 -f '%CHROM\t%POS\t%END\t%ANN\n' Ce328_From_2020GATK_noMissing.vcf.gz |\
grep "missense" |\
awk -F"||," '{print $1,$2,$3,$4}' OFS="\t" |\
cut -f-4 |\
cut -f1,11,13 -d"|" |\
sed 's/[A-Z]|//g' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk -F"/" '$1=$1' OFS="\t" |\
awk '{print $0, $5 % 3}' OFS="\t" |\
awk '$4 !~ "p.Ser" {print}' |\
awk '{if ($4 ~ "p.Met" || $4 ~ "p.Trp")  print $0, "0FOLD";
	  else if(($4 ~ "p.Phe" || $4 ~ "p.Ile" || $4 ~ "p.Val" || $4 ~ "p.Pro" || $4 ~ "p.Thr" || $4 ~ "p.Ala" || $4 ~ "p.Tyr" || $4 ~ "p.His" || $4 ~ "p.Gln" || $4 ~ "p.Asn" || $4 ~ "p.Lys" || $4 ~ "p.Asp" || $4 ~ "p.Glu" || $4 ~ "p.Cys" || $4 ~ "p.Gly") && ($7 == 1 || $7 == 2)) print $0, "0FOLD";
	  else if(($4 ~ "p.Leu" || $4 ~ "p.Arg") && $7 == 1) print $0, "0FOLD"}' OFS="\t" |\
cut -f1,2,3,8 |\
awk '{print $1, $2 = $2 - 1, $3, $4}' OFS="\t" |\
bgzip > 0FOLD_SITES.bed.gz

tabix 0FOLD_SITES.bed.gz


############################################################################### ANNOTATE VCF
# annotate vcf
vcfanno ANNOTATION_conf.toml Ce328_From_2020GATK_noMissing.vcf.gz |\
bcftools view -Oz -o Ce328_From_2020GATK_noMissing_ANNOTATED.vcf.gz

tabix -p vcf Ce328_From_2020GATK_noMissing_ANNOTATED.vcf.gz

# generate SFS input - remove MASKED regions from WS266
bcftools view --samples ^XZ1516 Ce328_From_2020GATK_noMissing_ANNOTATED.vcf.gz |\
vcffixup - |\
bcftools query  -f '%CHROM\t%POS\t%REF\t%ALT\t%AA\t%AC\t%AF\t%AN\t%SNPEFF_REGION\t%MASKED\t%4FOLD\t%0FOLD\t%GENOMIC_REGION\t%ANN\n' |\
awk -F"|" '$1=$1' OFS="\t" |\
awk '$10 != "Masked" {print $0}' OFS="\t" |\
awk 'BEGIN{OFS="\t"; print "CHROM", "POS", "REF", "ALT", "AA", "AC", "AF", "AN", "SNPEFF_REGION", "4FOLD" ,"0FOLD", "GENOMIC_REGION", "EFFECT", "WBGeneID", "aa_CHANGE", "CDS_POS", "aa_POS"};
	 {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $11, $12, $13, $16, $18, $24, $26, $27}' OFS="\t" > SFS_INPUT.tsv

############################################################################### load file in R and run processing script


# options
# 3_prime_UTR_variant
# 5_prime_UTR_premature_start_codon_gain_variant
# 5_prime_UTR_variant
# downstream_gene_variant
# initiator_codon_variant
# intergenic_region
# intron_variant
# missense_variant
# missense_variant&splice_region_variant
# non_coding_transcript_exon_variant
# splice_acceptor_variant&intron_variant
# splice_acceptor_variant&splice_donor_variant&intron_variant
# splice_acceptor_variant&splice_region_variant&intron_variant
# splice_donor_variant&intron_variant
# splice_region_variant
# splice_region_variant&intron_variant
# splice_region_variant&non_coding_transcript_exon_variant
# splice_region_variant&stop_retained_variant
# splice_region_variant&synonymous_variant
# start_lost
# start_lost&splice_region_variant
# stop_gained
# stop_gained&splice_region_variant
# stop_lost
# stop_lost&splice_region_variant
# stop_retained_variant
# synonymous_variant
# upstream_gene_variant

