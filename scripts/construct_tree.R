#!/usr/bin/env Rscript
### Given mash distance construct phylogenetic tree from it

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0){
	stop("Use as follow ./construct_tree.R distance_file fasta_dir output_file")
}


disfile=args[1]
fastadir=args[2]
outfile=args[3]

fasta_suffix=ifelse( length(args)>3, args[4], "fa")


suffix_fasta=args[4]

#loading library 
library("tidyverse")
library("ape")

datdis  <- read.table(disfile,header=FALSE, stringsAsFactors =FALSE)


#rename the header
colnames(datdis)  <- c("anim1","anim2","distr","comp4","comp5")


# give correct assembly name

datdis$anim1c  <- str_extract(pattern ="([A-Z][A-Za-z]+)", string=datdis$anim1)
datdis$anim2c  <- str_extract(pattern ="([A-Z][A-Za-z]+)", string=datdis$anim2)

# make distance into a wide matrix 
datsel  <- datdis  %>% select(anim1c,anim2c, distr)
datwide  <- datsel  %>% pivot_wider(names_from = anim2c, values_from = distr)

datmat  <- as.matrix(datwide  %>% select(-anim1c))
rownames(datmat)  <- datwide$anim1c

# apply neighbor joining 

tr  <- nj(datmat)

# create seq file for cactus

write.tree(tr,file=outfile)

uniq_anim <- unique(datdis$anim1c)

for (anim in uniq_anim){
	cat(anim,paste0(fastadir,"/",anim,".",suffix_fasta),"\n",file=outfile,append=TRUE)
}


