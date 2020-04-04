#Script for analysing the summary data with two-sample MR bi-directionally
#UPDATED for the newest version of MRBase R package
#For more information on commands/additional options see: https://mrcieu.github.io/TwoSampleMR/
#Modified from Robyn Wootton September 2018

######################################################################################################################################
# Load packages
######################################################################################################################################
rm(list=ls(all=TRUE)) #empties your R environment

#Load packages - you will have to install the packages first time you run the script
install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR") #re-run this to update MR BASE
library(TwoSampleMR)
install.packages("ggplot2")
library(ggplot2)
install.packages("knitr")
library(knitr)
install_github('qingyuanzhao/mr.raps')
library(mr.raps)
install.packages("xlsx")
library(xlsx)
install.packages("psych")
library(psych)


library(MRInstruments)

######################################################################################################################################
# Read exposure data
######################################################################################################################################
#The following code is for reading in your own data
#If instead, you want to use data from MR Base, see https://mrcieu.github.io/TwoSampleMR/

# Make sure headings are formatted like so: 
# "Phenotype 	SNP			CHR 	BP 			effect_allele 	other_allele 		eaf   beta 	  se 		  pval"
# "Neuroticism 	rs4653651 	1 		225862060 	A 				G 				0.32	-0.091 	0.0202 	6.443e-06"

#can also define headings yourself

#read in the data file

exp <- read.table("PATH_TO_FILE/file.txt",
                 header=T)


#create alternative SNP column if necessary

#exp$RSID <- paste(exp$CHR, exp$POS, "SNP", sep=":")

exp$phenotype <- rep("phenotype")
exp$sample_size <- rep(enter value)


#************************************************************************************************************************

#Read in this data file, formatted for MR Base
#Have told the package the name of the columns in the exposure data (in "") and what these correspond to
exp_dat <- format_data(exp, type="exposure",
                  snp_col = "SNP",
                  beta_col = "BETA",
                  se_col = "SE",
                  effect_allele_col = "EFFECT_ALLELE",
                  other_allele_col = "NON_EFFECT_ALLELE",
                  eaf_col = "EFFECT_ALLELE_FREQ",
                  phenotype_col = "phenotype",
                  pval_col = "PVAL",
                  samplesize_col = "N"
)

head(exp_dat) #N=#enter value



######################################################################################################################################
# Read in outcome data
######################################################################################################################################
#The following code is for using an outcome saved in MR Base

#get the library of available outcomes in MR Base
#ao <- available_outcomes()
#head(ao)

#Use the following commands to search for outcome of interest and save the ID and column numbers
#lung<-subset(ao, ao$consortium=="ILCCO") #id=966, c=1051
#heart<-subset(ao, ao$consortium=="CARDIoGRAMplusC4D") #id=7, c=756

#Alternatively, read in your own outcome GWAS summary dataset

outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "outcome_data.txt",
  sep="\t",
  #sep=" ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "Tested_Allele",
  other_allele_col = "Other_Allele",
  eaf_col = "Freq",
  pval_col = "P",
  samplesize_col = "N"
)

                     
outcome_dat$outcome <- rep("outcome_trait")

#define sample size if "N" was not present in the dataset
#outcome_dat$samplesize.outcome <- rep(enter value)

######################################################################################################################################
# Harmonise data
######################################################################################################################################
# harminising data - action 2 tries to align palindromic SNPs 
dat <- harmonise_data( 
	exposure_dat = exp_dat,
	outcome_dat = outcome_dat,
	action = 1
)

######################################################################################################################################
# Run 2 sample MR
######################################################################################################################################
# heterogeneity measures - significant heterogeneity in the SNP-exposure effects might be suggestive of pleiotropy
mr_het <- mr_heterogeneity(dat) #Cochran's Q
mr_ruck <- mr_rucker(dat) #Rucker's Q

# primary methods - this is the main MR analysis. See ?mr to see alternatives methods to include

#this is the default command:
res <- mr(dat)

#can use the command below to perform specific analyses:
res <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))


# Egger intercept - test of directional pleiotropy. A significant intercept suggests signifciant pleiotropy
mr_egger_int <- mr_pleiotropy_test(dat)

# Main results
res

#Turn into odds ratio (for binary outome)
or<-generate_odds_ratios(res)
or

#Save results
setwd("output_folder/") #change to your location
write.csv(or, "./Results_file.csv", row.names=F, quote=F)


#perform single SNP analysis:

# default single SNP analyses gives the wald ratio:
res_single <- mr_singlesnp(dat)

# leave one out analyses - by defalut uses IVW
res_loo <- mr_leaveoneout(dat)


#####################################################################################################################################

#PLOT RESULTS

p1 <- mr_scatter_plot(res, dat)
p1[[1]]

#save plot

ggsave(p1[[1]], file="./Results_file_MR_Egger_plot.pdf", width=7, height=7)

#create forest plot

p2 <- mr_forest_plot(res_single)
p2[[1]]

ggsave(p2[[1]], file="./Results_file_forest_plot.pdf", width=7, height=7)

#create plot for leave-one-out analysis
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]

ggsave(p3[[1]], file="./Results_file_loo_analysis.pdf", width=7, height=7)

#create funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]

ggsave(p4[[1]], file="./Results_file_funnel_plot.pdf", width=7, height=7)


#Save the full results report
mr_report(dat, output_path="./Results_file/", output_type = "html",
          author = "Author Name", study = "Two-sample  Results")
