setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
operon_first <- read.delim("Operons_First_Gene.txt", sep = "\t", header = TRUE) #contains only the first operon
operon_all <- read.delim("Operons_Reorganized.txt", sep = "\t", header = TRUE) #contains all the operons, column 1= gene name, column 2= operon groups
sym <- read.delim("Symbiosis_regulated_genes.txt", sep = "\t", header = TRUE) #symbiosis
genome <- read.delim("CGCANNNNNGTG_250.txt", sep = "\t", header = TRUE) #one of the eight motifs
genome$Motif_count= as.numeric(as.character(genome$Motif_count)) #it was read as integer, but we need it to be numeric

#Symbiosis list
#making a new dataframe
sym_new <- data.frame(matrix(nrow=1, ncol=1)) 
names(sym_new) <- c("Locus_tag")
i = 0
for(n in 1:nrow(sym)) {
  pos1 <- match(sym[n,1], operon_all[,1]) #match with all the genes that s part of an operon
  if(is.na(pos1)) {
    i <- i + 1
    sym_new[i,1] <- levels(sym$Locus_tag)[sym[n,1]] #put it in the list
  } else { #if in the operon list
    i <- i + 1
    pos2 <- match(operon_all[pos1,2], operon_first[,2])  #find the first gene of that operon
    sym_new[i,1] <- levels(operon_first$Gene)[operon_first[pos2,1]] #put the first gene of that operon in the list
  }
}
sym_new = unique(sym_new) #remove duplicates
#merging the data to match the motif count back
x = merge(genome, sym_new, by.x="Gene_name", by.y="Locus_tag", all=FALSE)
sym_new = data.frame("Locus_tag"=x$Gene_name, "Motif_count"=x$Motif_count)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: Motif count - M.C
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(sym_new)){ 
  if (sym_new[r,2]>0) { #row 1, column 2 is the motif number
    sym_new[r,3] = 1
    M.C=M.C+1
  } else {
    sym_new[r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
colnames(sym_new) = c("Locus_tag", "Motif_count", "Motif_Promoter")

#Genome list
genome <- genome[!(genome$Gene_name %in% sym$Locus_tag),]
#new dataframe
genome_new<- data.frame(matrix(nrow=1, ncol=1)) 
names(genome_new) <- c("Locus_tag") #naming the column
I <- 0
for(n in 1:nrow(genome)) { 
  pos1 <- match(genome[n,1], operon_all[,1])
  if(is.na(pos1)) {
    I <- I + 1
    genome_new[I,1] <- levels(genome$Gene_name)[genome[n,1]] #put it in the list
  } else { #if in the operon list
    I <- I + 1
    pos2 <- match(operon_all[pos1,2], operon_first[,2])  #find the first gene of that operon
    genome_new[I,1] <- levels(operon_first$Gene)[operon_first[pos2,1]] #put the first gene of that operon in the list
  }
}
genome_new <- unique(genome_new) #remove duplicated
#finding the motif_count for these genes
y = merge(genome_new, genome, by.x="Locus_tag",by.y ="Gene_name", all=FALSE)
#this "y" is 85 less than genome_new, checked, it is all good 
genome_new =data.frame('Locus_tag'=y$Locus_tag,'Motif_count'=y$Motif_count)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
R=1 #row number
#counting the genes: M(motif) G(genome)
M.G=0 #number of genes with a motif in the promoter 
nM.G=0 #number of genes without a motif in the promoter (ie: no motif)
for (i in 1:nrow(genome_new)){
  
  if (genome_new[R,2]>0) {
    genome_new[R,3] = 1
    M.G=M.G+1
  } else {
    genome_new[R,3] = 0
    nM.G=nM.G+1
  }
  R=R+1
}
#naming the column
colnames(genome_new)=c("Locus_tag","Motif_count", "Motif_Promoter")

#Fisher test overall
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
Fisher.data = data.frame("Motif"=Motif, "No.Motif"=No.Motif)
row.names(Fisher.data) = c("Symbiosis", "Genome")
Fisher.data
fisher.test(Fisher.data)

#T-test
n.motif=c(sym_new$Motif_count, genome_new$Motif_count) #need to line them up
repS = rep("S", nrow(sym_new))
repG = rep ("G", nrow(genome_new))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#dealing with the Down and Up clusters
#merging the cluster Down and Up back
XX = merge(sym_new, sym, by.x = "Locus_tag", by.y = "Locus_tag")

sub1 = subset(XX, Cluster == "Down")
c1 = data.frame("Locus_tag"=sub1$Locus_tag, "Motif_count"=sub1$Motif_count, "Motif_Promoter"=sub1$Motif_Promoter)

sub2 = subset(XX, Cluster =="Up")
c2 = data.frame("Locus_tag"=sub2$Locus_tag, "Motif_count"=sub2$Motif_count, "Motif_Promoter"=sub2$Motif_Promoter)

#Down (sub 1)
#counting the number of genes with motif
subM.C = subset(c1, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c1)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c1$Motif_count, genome_new$Motif_count)
repS = rep("S", nrow(c1))
repG = rep ("G", nrow(genome_new))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#Up (sub 2)
#counting the number of genes with motif
subM.C = subset(c2, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c2)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c2$Motif_count, genome_new$Motif_count)
repS = rep("S", nrow(c2))
repG = rep ("G", nrow(genome_new))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)









