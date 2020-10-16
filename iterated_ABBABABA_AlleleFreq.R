#Genome wide ABBA BABA analysis
#http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/abba-baba-statistics/

#I am unsure of the package 'boot' is in R by default or if it needs to be downloaded. if it requires downloading use: install.packages("boot")
library("boot")
library("parallel")

#set the work directory where your input files are and where the outputfile will be produced
setwd("/path/to/workDirectory/")

#this runs a python script to parse the VCF file and produce a geno file from which allele frequencies will be calculated.
system("python2 parseVCF.py -i populations.snps.vcf | gzip > spicAll.geno.gz")

#This command runs a python script that will calculate allele frequencies from the geno file. It needs to be modifed for your purpose. 
# Each population that you want to analyze needs to be entered following a "-p". you then include a popmap after --popsFile in which you have defined which
# individuals you want in each population. It is possible to calculate allele frequencies in two way, 1) assigning all states held by the outgroup as ancestral or 
# 2) choose the ancestral state to always be the state of the minor allele. 

#this first command (and the recommended one) makes the last population (the last population following a "-p", in my case Ptet) the outgroup.
# We make sure the script understands that we want to do this by giving the argument --target the input 'derived'.
system("python2 freq.py -g spicAll.geno.gz -p PgalGMS128 -p PgalGMS137 -p PgalGMS141 -p PspGMS130 -p PspGMS138 -p PspD56 -p PspD47 -p PspD1 -p Ptet --popsFile /path/to/popmap --target derived -o spicAll.tsv.gz")

#In case you don't have an outgroup, you can still generate allele frequencies. I have not noticed significant biases in this but it is not the recommend approach. 
# To set this up make sure --target is set to 'minor'. 
#system("python2 freq.py -g spicAll.geno.gz -p PgalGMS128 -p PgalGMS137 -p PgalGMS141 -p PspGMS130 -p PspGMS138 -p PspD56 -p PspD47 -p PspD1 --popsFile /path/to/popmap --target minor -o spicAll.tsv.gz")

#Read in our frequencies.
freq_table = read.table("spicAll.tsv.gz", header=T, as.is=T)

#set which populations you want to be the "pure type" that has not experienced introgression
P1_species=c("PspD56", "PspD47", "PspD1")
#set this to the populations you want to test for introgression. it is recommened to also test what happens in of you put pure types in both P1 and P2. 
# They will act as a standard for what unintrogressed types look like.
P2_species=c("PspGMS130", "PspGMS138")
# this should be set to the populations that might have introgressed into P2
P3_species=c("PgalGMS128","PgalGMS137","PgalGMS141") 

#It might be wise to change the name of the output file for every test you make. the current is set for a test if gallicum introgresses into spicatum, 
# where all spicatum are being tested for introgression against every other spicatum
outFile="ABBABABA_result_gal_spi_spi.tab"

# NULL can be replaced with a number of cores if you don't want to use all. otherwise the if statement will try and detect how many you have and use them all
NrCores=NULL
if (is.null(NrCores)) {NrCores=detectCores()}

#How many bootstrap replicates to make. This won't change the results of the D-statistics, but it might impact the significance level if it is too low. 
# 1000 is usually fine.
nrReps = 1000

nrow(freq_table)

#Now, to compute D, we need to define populations P1, P2 and P3. W

# explanation of header: P1, P2, P3: are which populations are being analyzed. D-value is the value denoting how high the introgression is. F1 hybrids generally have 
# t0_Z-score, t0_p-value: are the Z-score and P-value for how much (and how likely) the mean of the bootstrapped D replicates varies from the original.
# Zero_Z-score, Zero_p-value: are how much (and how likely) mean of the bootstrapped D replicates are varying from 0 (no significant difference between ABBA and BABA)
# T-test_p-value: results of a T-test to see if mean of bootstrapped D is sigificantly different from 0. usually not very useful.
# I generally use a sigificance level of 0.001 rather than the standard scientific 0.05 because 5% probabilty of erronous estimate is quite high (every 20th test will be wrong).
# number_of_BABA_sites, number_of_ABBA_sites: self explanatory
#We set these populations and then compute D by extracting the derived allele frequencies for all SNPs for the three populations.

resultsOutput <- paste("P1","P2","P3","D-value","t0_Z-score","t0_p-value","Zero_Z-score","Zero_p-value","T-test_p-value","number_of_BABA_sites","number_of_ABBA_sites",sep = "\t")
cat(resultsOutput, file=outFile, sep="\n", append=FALSE)
for (P3 in P3_species) {
  for (P2 in P2_species) {
    for (P1 in P1_species) {
      if (P1 != P2 & P1 != P3 & P2 != P3) {
          #We then set up the population test where we calculate the allele frequences
        ABBA = (1 - freq_table[,P1]) * freq_table[,P2] * freq_table[,P3]
        BABA = freq_table[,P1] * (1 - freq_table[,P2]) * freq_table[,P3]
      
          #load the result into a dataframe
        ABBA_BABA_df = as.data.frame(cbind(ABBA,BABA))
        
          #We then define a function for D.
        D.stat = function(dataframe, i) (sum(dataframe$ABBA[i], na.rm = T) - sum(dataframe$BABA[i], na.rm = T)) / (sum(dataframe$ABBA[i], na.rm = T) + sum(dataframe$BABA[i], na.rm = T))
        
          #This function is then bootstrapped
        D_bootstrapped <- boot(ABBA_BABA_df, D.stat, R = nrReps, stype = "i", parallel = "multicore", ncpus = NrCores)
        print("P1 P2 P3")
        print(paste(P1,P2,P3))
        print(D_bootstrapped)
        #print(mean(D_bootstrapped$t))
          #Finally a p-value is calculated via z-score, does the mean of t deviate from 0?
        zero.z = (0 - mean(D_bootstrapped$t)) / sd(D_bootstrapped$t)
        #print(paste("Z-Score:", zero.z))
        zero.p.value <- (2*pnorm(-abs(zero.z))) #what is the probability that mean of t is 0?
        #print(paste("P-value:", zero.p.value))
        t0.z = (D_bootstrapped$t0 - mean(D_bootstrapped$t)) / sd(D_bootstrapped$t) # does the mean of the bootstrapped ts deviate from the original t?
        
        t0.p.value <- 1 - (2*pnorm(-abs(t0.z))) # what is the probability that the mean of the bootstrapped ts is the same as the original t?
        
        t.pvalue <- t.test(D_bootstrapped$t)
        
        print(paste("T-test p-value: ",t.pvalue$p.value))
          #prints the number of loci that are informative for the ABBA and the BABA configuration
        print(paste("number of BABA sites:", sum(BABA>0, na.rm = T)))
        print(paste("number of ABBA sites:", sum(ABBA>0, na.rm = T)))
        out <- paste(P1,P2,P3,D_bootstrapped$t0,t0.z, t0.p.value,zero.z,zero.p.value,t.pvalue$p.value,sum(BABA>0, na.rm = T),sum(ABBA>0, na.rm = T),sep = "\t")
        cat(out, file=outFile, sep="\n", append=TRUE)
      } #if statement
    } # P1
  } # P2
} # P3

