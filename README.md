# SFS_Invariant_Sites
Calculate Site Frequency Spectra With and Without Invariant Sites

All scripts used to generate the spectra are currently located in the repo, but more work needs to be done to make the scripts more extensible.

## Current Workflow

1. All analyses are run on three strain sets: swept, divergent, and the whole population. 
2. Scripts to generate spectra are in the `2020_SFS_Analysis` folder
3. Run `classify_strains.R` to generate the different strain populations
4. For each population do the following:
  * Paste strain names output by `classify_strains.R` on line 5 of `make_files.sh` in the population subdirectory 
  * Run `make_files.sh` using the pruned VCF: VCF pruning was performed by only including sites with no missing data from the 20200815 CeNDR release.
  * `make_files.sh` generates a file called `SFS_INPUT.tsv`
  * Run `generate_invariant.R`. This script takes a spliced CDS fasta file and the output from `make_files.sh` to generate counts of 0- and 4-fold sites that do are invariant across the tested population. This script needs to be updated to work within the repository file structure. 
  * Run `generate_spectra.R` to generate spectra for DFE analysis
5. See `Readme` in Scripts directory to perform DFE analysis
  

## Previous Workflow

1. Run `GENERATE_SFS_FILES.sh`. This takes a VCF and generates a processed data set that is used by `Generate_Spectra.R` to generate spectra files. This script needs to be updated in the following ways (currently doesn't run as a script, but all the commands are there):
  * Modify paths to work in the github repository file structure.
  * Modify script to take VCF, sample names file, and ancestor name as an input.
2. Run `Invariant_SFS.R`. This script takes a spliced CDS fasta file and the output from `GENERATE_SFS_FILES.sh` to generate counts of 0- and 4-fold sites that do are invariant across the tested population. This script needs to be updated to work within the repository file structure. 
3. Run `Generate_Spectra.R`. This generates `.sfs` files used for DFE analysis. The parameters for this script are:
  * The output from `GENERATE_SFS_FILES.sh`: `SFS_INPUT.tsv`
  * `no_indel` or `indel` to include indels in the spectra
  * The output from `Invariant_SFS.R`: `invariant_site_by_region.tsv` 
  
