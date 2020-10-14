############################################################################### GENERATE ANNOTATION FILES
# put yourself in annotation directory

bcftools view -s AB1,BRC20067,BRC20263,CB4851,CB4852,CB4853,CB4854,CB4855,CB4857,CB4858,CB4932,CX11254,CX11262,CX11264,CX11271,CX11276,CX11285,CX11292,CX11307,CX11314,CX11315,DL200,DL226,ECA348,ECA349,ECA592,ECA594,ECA640,ECA923,ECA928,ECA930,ED3005,ED3011,ED3012,ED3017,ED3040,ED3046,ED3048,ED3049,ED3052,ED3073,ED3077,EG4347,EG4349,EG4724,EG4725,EG4946,GXW1,JT11398,JU1088,JU1172,JU1200,JU1212,JU1213,JU1242,JU1246,JU1249,JU1395,JU1400,JU1409,JU1440,JU1491,JU1530,JU1543,JU1568,JU1581,JU1586,JU1652,JU1666,JU1792,JU1793,JU1808,JU1896,JU1934,JU2001,JU2007,JU2016,JU2017,JU2106,JU2131,JU2141,JU2234,JU2250,JU2257,JU2464,JU2466,JU2478,JU2513,JU2519,JU2522,JU2534,JU2565,JU2566,JU2570,JU2572,JU2575,JU2576,JU2578,JU258,JU2581,JU2586,JU2587,JU2592,JU2593,JU2600,JU2610,JU2619,JU2800,JU2811,JU2825,JU2829,JU2838,JU2841,JU2853,JU2862,JU2866,JU2878,JU2879,JU2906,JU2907,JU310,JU311,JU3125,JU3127,JU3128,JU3132,JU3134,JU3135,JU3137,JU3140,JU3141,JU3144,JU3166,JU3169,JU3224,JU3225,JU3228,JU323,JU3280,JU3282,JU3291,JU346,JU360,JU367,JU393,JU394,JU397,JU406,JU440,JU561,JU642,JU751,JU774,JU775,JU778,JU782,JU792,JU830,JU847,KR314,LKC34,MY1,MY10,MY18,MY2147,MY2212,MY2453,MY2530,MY2535,MY2573,MY2585,MY2693,MY2713,MY2741,MY518,MY679,MY772,MY795,MY920,N2,NIC1,NIC1049,NIC1107,NIC1119,NIC166,NIC195,NIC199,NIC2,NIC207,NIC231,NIC236,NIC242,NIC252,NIC255,NIC256,NIC258,NIC259,NIC260,NIC261,NIC262,NIC266,NIC267,NIC269,NIC271,NIC272,NIC274,NIC275,NIC276,NIC277,NIC3,NIC501,NIC511,NIC513,NIC514,NIC515,NIC522,NIC523,NIC526,NIC527,NIC528,NIC529,PB303,PB306,PS2025,PX179,QG2075,QG2811,QG2813,QG2818,QG2823,QG2824,QG2825,QG2827,QG2828,QG2832,QG2835,QG2836,QG2837,QG2838,QG2841,QG2843,QG2846,QG2850,QG2854,QG2855,QG2857,QG2873,QG2874,QG2875,QG2877,QG2932,QG536,QG556,QG557,QW947,QX1212,QX1233,QX1792,RC301,WN2001,WN2002,WN2033,WN2050,WN2063,WN2064,WN2066,XZ1515,XZ1516,XZ1672,XZ1734,XZ1735,XZ1756,XZ2018,XZ2020 Ce328_From_2020GATK.vcf.gz |\
bcftools filter -i N_MISSING=0 |\
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

