library(maftools)

setwd("documents/Test Data/LGG_maf_files")
files <- list.files(pattern = "TCGA-*")
list_syn <- list()
list_nonsyn <- list()



samples <- sort((substr(files,1,13))) #duplicates ? 
genes <- sort(colnames(exp_table))
syn_muts <- matrix(0,length(samples),length(genes))
dimnames(syn_muts) <- list(samples,genes)
nonsyn_muts <- syn_muts

for(maf in files){
  update <- read.delim(maf,comment.char="#")
  nonsyn <- (read.maf(update))@data
  syn <- (read.maf(update))@maf.silent
  
  
  
  # list_nonsyn <- c(list_nonsyn,maftool@data)
  # list_syn <- c(list_syn,maftool@maf.silent)
}

nonsyn_matrix <- do.call(cbind,list_nonsyn)
syn_matrix <- do.call(rbind,list_syn)


nonsyn <- test@data
gene_names <- paste0(nonsyn$Hugo_Symbol,".",nonsyn$Entrez_Gene_Id)
syn <- test@maf.silent
tnon <- with(nonsyn,table(nonsyn$Tumor_Sample_Barcode,gene_names))
tnon <- tnon[,colnames(tnon) %in% colnames(nonsyn_muts)]



nonsyn_muts[substr(rownames(tnon),1,12),colnames(tnon)] <- tnon
