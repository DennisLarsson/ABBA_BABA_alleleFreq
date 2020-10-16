# ABBA_BABA_alleleFreq version 1.0
calculate D statistics using allele frequencies

This script is a derivate of the tutorial written by Simon Martin here: http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/

Follow the instructions written there to understand what the script is doing in general. I have extended the script to do multiple analyses for several different populations and combinations of populations, as well integrating a bootstrap procedure to calculate a p-value instead of the jack-knife approached used in the original tutorial (Jack-knife makes sense for long continous sequences where loci might be linked (whole genome sequencing), while the bootstrap approach is more appropriate for the mostly unlinked loci in a RADseq dataset).

I have then also intregrated printing of the results into a single file for ease of overview.

Download the package by Simon martin using this command: git clone https://github.com/simonhmartin/genomics_general

Then transfer the scripts: parseVCF.py (in the folder VCF_processing) and freq.py to the workdirectory where you have your input vcf file.
