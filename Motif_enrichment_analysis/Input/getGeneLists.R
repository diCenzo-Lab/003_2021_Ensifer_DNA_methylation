# To consider only first gene of operon and non-operon genes (i.e., analysis part two)
operon_first <- read.delim("Downloads/Operons_First_Gene.txt", sep = "\t", header = TRUE)
operon_all <- read.delim("Downloads/Operons_Reorganized.txt", sep = "\t", header = TRUE)
carbon <- read.delim("Downloads/Carbon_source_regulated_genes.txt", sep = "\t", header = TRUE)
motif <- read.delim("Downloads/GANTC_250.txt", sep = "\t", header = TRUE)

i <- 0
carbon_new <- data.frame(matrix(,nrow=1, ncol=1))
names(carbon_new) <- c("Locus_tag")
for(n in 1:nrow(carbon)) {
  pos1 <- match(carbon[n,1], operon_all[,1])
  pos2 <- match(carbon[n,1], operon_first[,1])
  if(is.na(pos1)) {
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]]
  }
  if(!is.na(pos2)) {
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]]
  }
}

# To consider non-operon genes and first gene of operon for all genes in an operon (i.e., analysis part three)
operon_first <- read.delim("Downloads/Operons_First_Gene.txt", sep = "\t", header = TRUE)
operon_all <- read.delim("Downloads/Operons_Reorganized.txt", sep = "\t", header = TRUE)
carbon <- read.delim("Downloads/Carbon_source_regulated_genes.txt", sep = "\t", header = TRUE)
motif <- read.delim("Downloads/GANTC_250.txt", sep = "\t", header = TRUE)

i <- 0
carbon_new <- data.frame(matrix(,nrow=1, ncol=1))
names(carbon_new) <- c("Locus_tag")
for(n in 1:nrow(carbon)) {
  pos1 <- match(carbon[n,1], operon_all[,1])
  if(is.na(pos1)) {
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]]
  } else {
    i <- i + 1
    pos2 <- match(operon_all[pos1,2], operon_first[,2])
    carbon_new[i,1] <- levels(operon_first$Gene)[operon_first[pos2,1]]
  }
}
carbon_new <- unique(carbon_new)





