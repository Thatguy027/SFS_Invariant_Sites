#!/bin/bash

# run DFE-alpha on Slim outputs

# for 1 individual slim output, need to:
#	calculate neutral and selected SFS
#	make config file
#	run est_dfe class 0
#	run est_dfe class 1 with output from class 0
#	run prop_muts_in_s_range on class 1 output

# if beneficials, run alpha_omega on est_dfe output class 0 and 1


# get all outputs in dir to run analyses on
echo "Provide the path to the directory containing the SLiM outputs to be analysed."
read dir

# put all those file names in a list
ls $dir | grep Fixed | sed "s/FixedOutput_//g" > base_InputNames.txt

# loop through the list and do each one at a time
for base_name in `cat base_InputNames.txt`
do
if echo $base_name | grep 20mbp
then
	genosize=20000000
fi
if echo $base_name | grep 24mbp
then
       	genosize=24000000
fi
if echo $base_name | grep 25mbp
then
	genosize=25000000
fi
if echo $base_name | grep 26mbp
then
	genosize=26000000
fi
if echo $base_name | grep 30mbp
then
	genosize=30000000
fi
if echo $base_name | grep 1e5
then
        popsize=100000
fi
if echo $base_name | grep 5e5
then
           popsize=500000
fi
if echo $base_name | grep 1e6
then
        popsize=1000000
fi
# turn the sample files into a format suitable for a "FullOutput..." file

## do this one if doing all generations sampled:
Rscript CommandLine_MakeModifiedFullSampleFile_ForDFEinput_throughTime.R $base_name $dir
done

# because there are multiple generations for the same file, make new list of base names:
ls $dir | grep Modified | sed "s/ModifiedSampleOutput_//g" > base_InputNames4.txt

# loop through the list and do each one at a time
for base_name in `cat base_InputNames4.txt`
do
# then analyze like normal:
Rscript CommandLine_RunSlimToDFEconversion_Celegans_ArmsVsCenters_UNFOLDED.R $base_name subsample $dir $genosize $popsize
# create divergence file (for all of them because easier to do in this loop anyway
Rscript CommandLine_RunSlimToAlphaOmega_Celegans_ArmsVsCenters.R $base_name subsample $dir $genosize $popsize
done

