################################################################################################################################################
##                                                                                                                      
##  ESTIMATE BETA-BINOMIAL DISTRIBUTION PARAMETER FROM FOETAL BLOOD MTDNA COPY NUMBER DISTRIBUTION
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
copynumber.file <- fread("/Users/mp34/team154_campbell/data/blood_foetal/mcmc_model_mtDNA_CopyNumber_tissue_Celltype.txt", header = T)

# visualize the copy number information
pdf("/Users/mp34/team154_campbell/plots/blood_foetal_copyNumber_hist.pdf")
hist(copynumber.file$copy_number, breaks = 200, col = "darkgreen", main = paste0("Mean copy number value, n = ", round(mean(copynumber.file$copy_number), 4)))
dev.off()

# read in the copy number information
copynumber.file <- fread("/Users/mp34/team154_campbell/data/blood_emily/mcmc_model_mtDNA_CopyNumber.txt", header = T)

# visualize the copy number information
pdf("/Users/mp34/team154_campbell/plots/blood_emily_copyNumber_hist.pdf")
hist(copynumber.file$copy_number, breaks = 200, col = "darkgreen", main = paste0("Mean copy number value, n = ", round(mean(copynumber.file$copy_number), 4)))
dev.off()

#####################################################################################
# EXPLORE THE BBMLE PACKAGE HERE 
#####################################################################################

# get the copy number vector here
cn.vector <- as.integer(round(copynumber.file$copy_number, 2))
foetal.blood.median <- median(cn.vector)

# get test data from a binomial distribution
# sim.vector <- rbetabinom(n=length(cn.vector),prob=0.3,size=7000,theta=1)
sim.vector <- rbetabinom(n=length(cn.vector),prob=0.1,size=5000,theta=25)

# visualize the observed and simulated 
pdf("/Users/mp34/team154_campbell/plots/blood_emily_simulated_copyNumber_hist.pdf", width = 12, height = 5)
par(mfrow=c(1, 2))
hist(cn.vector, breaks = 200, col = "darkgreen", main = paste0("\nMean observed copy number value,\n n = ", round(mean(cn.vector), 4)))
hist(sim.vector, breaks = 200, col = "darkblue", xlim = c(0, 5000), main = paste0("\nMean simulated copy number value,\n n = ", round(mean(sim.vector), 4)))
dev.off()

# check whether the parameters fit the observed histogram distribution
# i.e. they do not produce NAs
dbetabinom(cn.vector, prob = 0.2, size = 7000, theta = 1.25, log = T)

# provide the log-likelihood function from bbmle vignette
mtmp <- function(prob,size,theta){
  -sum(dbetabinom(cn.vector,prob,size,theta,log=TRUE))
}

# apply it to the data with best guess parameters
m0 <- mle2(mtmp, start=list(prob=0.2,theta=15),data=list(size=5000))

# check the model summary
summary(m0)

## OUTPUT
# Maximum likelihood estimation
# 
# Call:
#   mle2(minuslogl = mtmp, start = list(prob = 0.2, theta = 1), data = list(size = 7000))
# 
# Coefficients:
#   Estimate Std. Error z value     Pr(z)    
# prob  0.1152751  0.0065778  17.525 < 2.2e-16 ***
#   theta 5.2303402  0.4185043  12.498 < 2.2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# -2 log L: 5667.459

#####################################################################################
# ASSESS THE LOGLIKELIHOOD PROFILE
#####################################################################################

# construct the likelihood profile
p0 <- profile(m0)

# CHECK THE CI ESTIMATION
# quadratic approximation at the maximum likelihood estimate
confint(m0,method="quad")

# plot the each profile
par(mfrow=c(1,2))
plot(p0,plot.confstr=TRUE)

beta <- rbeta(1, 0.11, 21)

hist(test)
test <- rbinom(10, 5000, beta)


