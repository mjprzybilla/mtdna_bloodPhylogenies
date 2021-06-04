################################################################################################################################################
##                                                                                                                      
##  VISUALIZE MITOCHONDRIAL MUTATION ON MAXIMUM PARSINOMY TREES FOR EMILYS BLOOD DATA
##                                                                                                                      
##  Date: 02 JULY 2021                                                                                                                    
##  
##  Author: Tim Coorens adapted by Moritz Przybilla                
##
################################################################################################################################################

# clear workspace beforehand
rm(list = ls())

# package dependencies, which have to be installed are checked and installed if not available
list.of.packages <- c("ape", "reshape2", "stringr", "readr", "tidyverse", "ggpubr", "ggtree", "BuenColors",
                      "ComplexHeatmap", "phyloch")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos = "http://cran.us.r-project.org")
if(length(new.packages)) BiocManager::install(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# ignore all "simple" diagnostic messages (warnings or errors)
suppressMessages(invisible(lapply(list.of.packages, require, character.only = TRUE)))
options(stringsAsFactors = F)

#####################################################################################
# FUNCTION
#####################################################################################

add.color.bar<-function(leg,cols, title= "VAF",lims=c(0, 0.1),digits=2, prompt=TRUE){
  digits=2
  y<-0.65
  x<-115
  X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
  Y<-cbind(rep(y,length(cols)),rep(y,length(cols))) 		
  lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=4+2,lend=2) 
  for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=4,lend=2)
  text(x=x,y=y,round(lims[1],digits),pos=3)
  text(x=x+leg,y=y,round(lims[2],digits),pos=3)
  if(is.null(title)) title<-"P(state=1)"
  text(x=(2*x+leg)/2,y=y,title,pos=3)
}

is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

#####################################################################################
# READ IN DATA FROM TREES & MITOCHONDRIA
#####################################################################################

# determine the input and output
input.dir <- "/Users/mp34/sanger/team268im/analysis/mtdna"
output.dir <- "/Users/mp34/team154_campbell/plots"
dir.create(paste0(output.dir, "/blood_emily"))

# read in sample metadata and everything else which is needed
sampleInfo=read.csv(paste0(input.dir, "Panbody_key_10082020.csv"))
tissue_colours=read.csv(paste0(input.dir, "tissuetype_colour.csv"), header=F)
tree=read.tree(paste0(input.dir, "PD43851_snv_tree_with_branch_length.tree"))

# get the tree files and the mutation files
blood_trees <- list.files("/Users/mp34/team154_campbell/data/blood_emily", pattern = ".tree", full.names = T)
blood_trees <- blood_trees[grep("_adj", blood_trees)]

# get the sample ids
blood_ids <- str_split_fixed(basename(blood_trees), "_", 2)[,1] 

# read in metadata
blood_metadata <- fread("/Users/mp34/team268_martincorena/data/blood_metadata.txt", header = T)
converted.ids <- blood_metadata[match(blood_ids, blood_metadata$Sample), "PDID"]
blood_ids <- converted.ids$PDID

# read in the trees and name them
blood_trees <- lapply(blood_trees, read.tree)
names(blood_trees) <- blood_ids

# read the variant matrices
NR_mtdna.files <- list.files("/Users/mp34/sanger/team268im/analysis/mtdna/blood_emily", pattern = "snp_NR_mtdna_all.txt", all.files = T, recursive = T, full.names = T)
NV_mtdna.files <- list.files("/Users/mp34/sanger/team268im/analysis/mtdna/blood_emily", pattern = "snp_NV_mtdna_all.txt", all.files = T, recursive = T, full.names = T)
NR_mtdna.files <- NR_mtdna.files[grep("/snp", NR_mtdna.files)]
NV_mtdna.files <- NV_mtdna.files[grep("/snp", NV_mtdna.files)]

#####################################################################################
# ITERATE OVER EACH PATIENT AND VISUALIZE MUTATIONS ON THE TREE
#####################################################################################
i <- 1
for (i in 1:length(blood_ids)){
  
  # which patient
  patient <- blood_ids[i]
  
  # create dir
  dir.create(paste0(output.dir, "/blood_emily"))
  dir.create(paste0(output.dir, "/blood_emily/", patient))
  
  # collapse the tree and make the edge length equal
  tree <- blood_trees[[i]]
  tree_collapsed=blood_trees[[i]]
  tree_collapsed$edge.length=rep(1,nrow(tree_collapsed$edge))
  tree_collapsed$edge.length[tree$edge.length==0]=0
  tree_collapsed$edge.length[tree$edge[,2]%in%(1:length(tree_collapsed$tip.label))]=1
  
  # tree dataframe
  tree_df=as.data.frame(fortify(tree_collapsed))
  
  pdf(paste0(output.dir, "/blood_emily/", patient,"/phylogenetic_tree.pdf"), width=14, height=5,useDingbats = F)
  plot(tree,label.offset=0.01*max(tree_df$x), show.tip.label=F, direction="downwards")
  axisPhylo(side = 2,backward=F)
  dev.off()
  
  pdf(paste0(output.dir, "/blood_emily/", patient, "/phylogenetic_tree_collapsed.pdf"),width=9,height=3,useDingbats = F)
  plot(tree_collapsed,label.offset=0.01*max(tree_df$x),show.tip.label=F,direction="downwards")
  dev.off()
  
  #####################################################################################
  # PROCESS THE MITOCHONDRIAL CALLS TO BRING THEM IN SHAPE FOR TREE PLOTTING
  #####################################################################################
  
  # check if mutation data are present for patient
  if (is.integer0(grep(patient, NR_mtdna.files))){
    
    cat("Patient ", patient, " does not have mutation data!\n")
    next()
  }
  
  # subset to patient only & read in the data
  pt_NR.mtdna.file <- read.table(NR_mtdna.files[grep(patient, NR_mtdna.files)], header = T, sep = "\t")
  pt_NV.mtdna.file <- read.table(NV_mtdna.files[grep(patient, NV_mtdna.files)], header = T, sep = "\t")
  
  # only keep the samples which are in the tree
  pt_NR.mtdna.file <- pt_NR.mtdna.file[, colnames(pt_NR.mtdna.file) %in% tree$tip.label]
  pt_NV.mtdna.file <- pt_NV.mtdna.file[, colnames(pt_NV.mtdna.file) %in% tree$tip.label]
  
  # calculate VAF
  pt_mtdna.vaf <- pt_NV.mtdna.file/pt_NR.mtdna.file
  pt_mtdna.vaf[round(pt_mtdna.vaf,2) < 0.01] <- 0
  
  # remove non-complete cases
  pt_mtdna.vaf <- pt_mtdna.vaf[complete.cases(pt_mtdna.vaf),]
  
  # get variants which are shared between at least two samples
  bin.pt.mdna <- pt_mtdna.vaf
  bin.pt.mdna[bin.pt.mdna >= 0.01] <- 1
  present.variants <- rownames(bin.pt.mdna[rowSums(bin.pt.mdna) > 1,])
  
  # subset the vaf matrix according to that
  pt_mtdna.vaf.fil <- pt_mtdna.vaf[present.variants,]
  
  # recurrent variants
  bin.pt.mdna <- pt_mtdna.vaf.fil
  bin.pt.mdna[bin.pt.mdna >= 0.01] <- 1
  recurrent.variants <- rownames(bin.pt.mdna[rowSums(bin.pt.mdna) < (ncol(bin.pt.mdna)-1),])
  
  # subset the vaf matrix according to that
  pt_mtdna.vaf.fil <- pt_mtdna.vaf[recurrent.variants,]
  
  # check the variant heatmap
  mtdna.vaf.heatmap <- as.matrix(pt_mtdna.vaf.fil)
  
  # limit to VAF of up to 50%
  mtdna.vaf.heatmap[mtdna.vaf.heatmap > 0.1] <- 0.1
  
  hm <- Heatmap(mtdna.vaf.heatmap,
                col=as.character(jdb_palette("solar_rojos",type="continuous")),
                show_row_names = TRUE, 
                cluster_columns = TRUE,
                name = "VAF", use_raster = FALSE,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 4),
                cluster_rows = TRUE, 
                show_column_names = TRUE, border = T)
  
  pdf(paste0(output.dir, "/blood_emily/", patient, "/filtered_mtdna_mutation_heatmap.pdf"), width = 20, height = 20)
  print(hm)
  dev.off()
  
  # make different matrices based on VAF thresholds
  low.pt_mtdna.bin <- med.pt_mtdna.bin <- high.pt_mtdna.bin <- pt_mtdna.vaf.fil
  
  # low = everything higher 1% VAF
  low.pt_mtdna.bin[low.pt_mtdna.bin < 0.01] <- 0
  low.pt_mtdna.bin[low.pt_mtdna.bin >= 0.1] <- 0.1
  
  # remove non recurrent variants
  bin.pt.mdna <- low.pt_mtdna.bin
  bin.pt.mdna[bin.pt.mdna >= 0.01] <- 1
  recurrent.variants <- rownames(bin.pt.mdna[rowSums(bin.pt.mdna) > 1,])
  
  # transpose and make rownmaes as columns
  low.pt_mtdna.bin <- as.data.frame(t(low.pt_mtdna.bin))
  low.pt_mtdna.bin <- low.pt_mtdna.bin[,recurrent.variants]
  low.pt_mtdna.bin$label <- rownames(low.pt_mtdna.bin)
  
  #####################################################################################
  # PLOT INDIVIDUAL VAF VARIANTS
  #####################################################################################
  
  # set the parameters for the plotting
  cols <- colorRampPalette(c("black", "yellow", "red"))(20)
  
  # VAF limits 
  xlims <- c(0, 0.1)
  digits <- 2
  
  # create the colors 
  breaks<-0:20/20*(xlims[2]-xlims[1])+xlims[1]
  whichColor<-function(p,cols,breaks){
    i<-1
    while(p>=breaks[i]&&p>breaks[i+1]) i<-i+1
    cols[i]
  }
  
  # create output.dir
  dir.create(paste0(output.dir, "/blood_emily/", patient, "/all"))
  
  i <- 1
  for (i in 1:(ncol(low.pt_mtdna.bin)-1)){
    
    # get the mutation of interest
    print(i)
    mutation <- colnames(low.pt_mtdna.bin)[i]
    group <- as.list(low.pt_mtdna.bin[low.pt_mtdna.bin[,mutation] > 0, "label"])
    colors<-sapply(as.numeric(low.pt_mtdna.bin[,i]), whichColor, cols=cols, breaks=breaks)
    names(colors) <- rownames(low.pt_mtdna.bin)
    
    # 
    use.colors <- colors[names(colors) %in% group]
    ecol <- edge.color(tree, group, col = use.colors)
    
    # COLLAPSED
    pdf(paste0(output.dir, "/blood_emily/", patient, "/all/", mutation, "_colour.pdf"),width=13,height=4,useDingbats = F)
    plot(tree, label.offset=0.01*max(tree_df$x), 
         edge.color = ecol, show.tip.label=F,
         no.margin = TRUE,  direction="downwards", edge.width = 1.5)
    add.color.bar(50,cols,title = "VAF",xlims,digits,prompt=FALSE)
    dev.off()
    
  }
}
