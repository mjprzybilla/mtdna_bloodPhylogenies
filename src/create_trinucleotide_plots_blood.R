################################################################################################################################################
##                                                                                                                      
##  GENERATE OBSERVED VS EXPECTED MUTATION SPECTRA FOR BLOOD DATA
##                                                                                                                      
##  Date: 24 MARCH 2021                                                                                                                    
##  
##  Author: Andrew Lawson, adapted by Moritz Przybilla                                                                                                                    
##                                                                                                                        
##                                                                                                                      
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("BiocManager", "reshape2", "stringr", "readr", "tidyverse", "ggpubr", "SummarizedExperiment", "Rsamtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))

"%ni%" <- Negate("%in%")

#####################################################################################
# FUNCTION
#####################################################################################

trinucleotide_plot = function (mutations, file_name, analysis_type, analysis_region) {
  
  # subset the mutations to the respective columns
  mutations <- unique(mutations[,c("chr","pos","ref","mut","donor")])
  mutations <- mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")),]
  
  if(analysis_region == "coding"){
    mutations <- mutations[which(mutations$pos %in% coding_region),]
  }else if(analysis_region == "d_loop"){
    mutations <- mutations[which(mutations$pos %in% d_loop_region),]
  }else if(analysis_region != "all_mtDNA"){
    mutations <- NULL
  }
  
  mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
  
  # 2. Annotating the mutation from the pyrimidine base
  ntcomp = c(T="A",G="C",C="G",A="T")
  mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
  mutations$trinuc_ref_py = mutations$trinuc_ref
  for (j in 1:nrow(mutations)) {
    if (mutations$ref[j] %in% c("A","G")) { # Purine base
      mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
      mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
    }
  }
  
  # 3. Counting subs
  freqs_heavy = table(paste(mutations$sub[which(mutations$ref %in% c("A","G"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],3,3),sep="-"),sep=","))
  freqs_light = table(paste(mutations$sub[which(mutations$ref %in% c("C","T"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],3,3),sep="-"),sep=","))
  
  sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
  ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
  full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
  freqs_heavy_full = freqs_heavy[full_vec]; freqs_heavy_full[is.na(freqs_heavy_full)] = 0; names(freqs_heavy_full) = full_vec
  freqs_light_full = freqs_light[full_vec]; freqs_light_full[is.na(freqs_light_full)] = 0; names(freqs_light_full) = full_vec
  
  
  if(analysis_type == "obs_exp"){
    if(analysis_region == "coding"){
      heavy_base_freqs <- (mtdna_trinuc_freq["coding_heavy",] / sum(mtdna_trinuc_freq["coding_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["coding_light",] / sum(mtdna_trinuc_freq["coding_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else if(analysis_region == "d_loop"){
      heavy_base_freqs <- (mtdna_trinuc_freq["d_loop_heavy",] / sum(mtdna_trinuc_freq["d_loop_heavy",]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- (mtdna_trinuc_freq["d_loop_light",] / sum(mtdna_trinuc_freq["d_loop_light",]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else if(analysis_region == "all_mtDNA"){
      heavy_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_heavy","d_loop_heavy"),]))
      exp_heavy_counts <- sum(freqs_heavy) * as.numeric(c(rep(heavy_base_freqs[1:16]/3,times = 3),rep(heavy_base_freqs[17:32]/3,times = 3)))
      freqs_heavy_full <- freqs_heavy_full / exp_heavy_counts
      
      light_base_freqs <- colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]) / sum(colSums(mtdna_trinuc_freq[c("coding_light","d_loop_light"),]))
      exp_light_counts <- sum(freqs_light) * as.numeric(c(rep(light_base_freqs[1:16]/3,times = 3),rep(light_base_freqs[17:32]/3,times = 3)))
      freqs_light_full <- freqs_light_full / exp_light_counts
      
    }else{
      freqs_heavy_full = NULL
      freqs_light_full = NULL
    }
  }
  
  xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
  
  dev.new(width=10,height=4)
  colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
  y_heavy = freqs_heavy_full; y_light = freqs_light_full; maxy = max(c(y_heavy,y_light))
  
  if(analysis_type == "obs_exp"){
    ylab = "Mutation frequency (Obs/Exp)"
  }else{
    ylab = "Mutation count"
  }
  
  h_heavy = barplot(y_heavy, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab)
  h_light = barplot(-y_light, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, add = T)
  
  segments(y0 = maxy*1.5, y1 = maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = -maxy*1.5, y1 = -maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
  segments(y0 = 0, y1 = 0, x0 = 0.5, x1 = 192.5,  col = "black")
  abline(v = 0.5, col = "black")
  abline(v = 32.5, col = "black")
  abline(v = 64.5, col = "black")
  abline(v = 96.5, col = "black")
  abline(v = 128.5, col = "black")
  abline(v = 160.5, col = "black")
  abline(v = 192.5, col = "black")
  
  
  for (j in 1:length(sub_vec)) {
    xpos = h_heavy[c((j-1)*16+1,j*16)]
    rect(xpos[1]-0.5, maxy*1.25, xpos[2]+0.5, maxy*1.15, border=NA, col=colvec[j*16])
    text(x=mean(xpos), y=maxy*1.15, pos=3, labels=sub_vec[j])
  }    
  dev.copy(pdf,file_name,width=12,height=5)
  dev.off()
  dev.off()
}

#####################################################################################
# READ IN DATA AND SET THE PARAMETERS
#####################################################################################

# and output directory
input.dir <- "/Users/mp34/sanger/team268im/analysis/mtdna"
output.dir <- "/Users/mp34/team154_campbell/plots"
dir.create(paste0(output.dir, "/mut_spectra"))
dir.create(paste0(output.dir, "/mut_spectra/counts"))
dir.create(paste0(output.dir, "/mut_spectra/obs_exp"))

# list patient.txt files for each study that is finished
mtdna.variant.files <- list.files(input.dir, pattern = "all_mtDNA_variants.txt", recursive = T, full.names = T, all.files = T)
mtdna.variant.files <- mtdna.variant.files[grep("blood", mtdna.variant.files)]
mtdna.variant.files <- mtdna.variant.files[grep("/all", mtdna.variant.files)]
mtdna.variant.files <- mtdna.variant.files[-grep("transplant", mtdna.variant.files)]

# file path to mtdna fasta file
genomeFile <- "/Users/mp34/sanger/lustre/bed_ref/hs37d5.fa"

# read in the trinucleotide context
mtdna_trinuc_freq <- readRDS("/Users/mp34/sanger/team268im/metadata/mtDNA_trinuc_freqs_coding_dloop_heavy_light.Rds")

# set regions for coding and dloop regions
coding_region <- 577:16023
d_loop_region <- c(1:576,16024:16569)

# get the tissue ids
tissues <- unique(str_split_fixed(mtdna.variant.files, "/", 10)[,8])

#####################################################################################
# CREATE PLOTS FOR DIFFERENT TISSUES
#####################################################################################
i <- 1
for(i in 1:length(tissues)){
  
  # which tissue
  tissue.id <- tissues[i]
  print(tissue.id)
  
  # grep the tissue files
  tissue.mtdna.variant.files <- mtdna.variant.files[grep(paste0( "/", tissue.id, "/"), mtdna.variant.files)] 
  
  # read in the data from the variant calls
  tissue.mtdna.variant.files <- lapply(tissue.mtdna.variant.files, read.table, header = T, sep = "\t")
  
  for (i in 1:length(tissue.mtdna.variant.files)){
    
    tissue.mutations <- tissue.mtdna.variant.files[[i]]
    tissue.mutations$chr <- as.character(tissue.mutations$chr)
    tissue.mutations$ref <- as.character(tissue.mutations$ref)
    tissue.mutations$mut <- as.character(tissue.mutations$mut)
    tissue.mtdna.variant.files[[i]] <- tissue.mutations
  }
  
  tissue.mtdna.variant.files <- bind_rows(tissue.mtdna.variant.files)
  
  # get the important columns
  mutations <- tissue.mtdna.variant.files[,c("sampleID", "chr", "pos", "ref", "mut", "PID")]
  colnames(mutations)[ncol(mutations)] <- "donor"
  
  # remove sampleID column and deduplicate
  mutations$sampleID <- NULL
  mutations <- mutations[!duplicated(mutations),]
  
  # plot the raw counts
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/counts/", tissue.id,"_coding_region_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"counts","coding")
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/counts/", tissue.id,"_d_loop_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"counts","d_loop")
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/counts/", tissue.id,"_all_mtDNA_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"counts","all_mtDNA")
  
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/obs_exp/", tissue.id,"_coding_region_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"obs_exp","coding")
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/obs_exp/", tissue.id,"_d_loop_obs_exp_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"obs_exp","d_loop")
  trinucleotide_plot(mutations,paste0(output.dir, "/mut_spectra/obs_exp/", tissue.id,"_all_mtDNA_obs_below_exp_0.01_pct_cutoff_trinucleotide_plot.pdf"),"obs_exp","all_mtDNA")
  
  # plot trinucleotide context for each patient
  patient.ids <- unique(mutations$donor)
  
  for (j in 1:length(patient.ids)){
    
    patient <- patient.ids[j]
    pt.mutations <- mutations[mutations$donor == patient,]
    
    if (nrow(pt.mutations) > 3){
      
      # plot the raw counts
      trinucleotide_plot(pt.mutations,paste0(output.dir, "/mut_spectra/counts/", tissue.id, "_", patient, "_all_mtDNA_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"counts","all_mtDNA")
      trinucleotide_plot(pt.mutations,paste0(output.dir, "/mut_spectra/obs_exp/", tissue.id, "_", patient, "_all_mtDNA_obs_below_exp_0.01_pct_cutoff_trinucleotide_plot.pdf"),"obs_exp","all_mtDNA")
      
    } else {
      
      trinucleotide_plot(pt.mutations,paste0(output.dir, "/mut_spectra/counts/", tissue.id, "_", patient, "_all_mtDNA_counts_below_0.01_pct_cutoff_trinucleotide_plot.pdf"),"counts","all_mtDNA")
      cat("Not enough mutations detected in patient ", patient, "\n")
      next
    }
    
    
  }
  
}
