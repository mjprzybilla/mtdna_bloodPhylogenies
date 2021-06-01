################################################################################################################################################
##                                                                                                                      
##  CHECK MEDIAN COVERAGE FOR THE MITOCHONDRIA IN FOETAL AND ADULT BLOOD
##                                                                                                                      
##  Date: 17 APRIL 2021                                                                                                                    
##  
##  Author: Moritz Przybilla                                                                                                                    
##                                                                                                                    
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
set.seed(1001)
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "reshape2", "stringr", "readr", "tidyverse", "ggpubr", "ComplexHeatmap", "BuenColors",
                      "bbmle", "emdbook", "extraDistr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

#####################################################################################
# READ IN THE COPY NUMBER DATA
#####################################################################################

# read in the copy number information
coverage.file <- read.csv("/Users/mp34/team268_martincorena/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv", header = T, sep = ",")
copynumber.file <- read.csv("/Users/mp34/team268_martincorena/data/whole_genome_coverage_pileup_and_bedtools_annotated.csv", header = T, sep = ",")
coverage.file <- coverage.file[, c("Sample", "Tissue", "pileup_mtDNA_coverage")]
copynumber.file <- copynumber.file[, c("Sample", "Tissue", "pileup_mtDNA_genomes")]
colnames(coverage.file) <- c("SampleID", "Tissue", "Coverage")
colnames(copynumber.file) <- c("SampleID", "Tissue", "copynumber")

# only keep the foetal blood infromation
foetal.blood.coverage.file <- coverage.file[coverage.file$Tissue == "blood_foetal",]
foetal.blood.coverage.file$PatientID <- substring(foetal.blood.coverage.file$SampleID, 1, 7)

# only keep the foetal blood infromation
foetal.blood.copynumber.file <- copynumber.file[copynumber.file$Tissue == "blood_foetal",]
foetal.blood.copynumber.file$PatientID <- substring(foetal.blood.copynumber.file$SampleID, 1, 7)

# calculate mean coverage
mean.blood.foetal.cov <- foetal.blood.coverage.file %>% group_by(PatientID) %>% summarize("Mean" = mean(Coverage))
mean.blood.foetal.copynumber <- foetal.blood.copynumber.file %>% group_by(PatientID) %>% summarize("Median" = median(copynumber))

## FOR ADULT BLOOD
# only keep the foetal blood infromation
adult.blood.coverage.file <- coverage.file[coverage.file$Tissue == "blood_emily",]
adult.blood.coverage.file$PatientID <- substring(adult.blood.coverage.file$SampleID, 1, 7)

# calculate mean coverage
mean.blood.adult.cov <- adult.blood.coverage.file %>% group_by(PatientID) %>% summarize("Mean" = mean(Coverage))

