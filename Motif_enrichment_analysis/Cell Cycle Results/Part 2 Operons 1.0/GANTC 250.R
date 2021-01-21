setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
operon_first <- read.delim("Operons_First_Gene.txt", sep = "\t", header = TRUE) #contains only the first operon
operon_all <- read.delim("Operons_Reorganized.txt", sep = "\t", header = TRUE) #contains all the operons, column 1= gene name, column 2= operon groups
cell <- read.delim("Cell_cycle_regulated_genes.txt", sep = "\t", header = TRUE) #cell regulation
genome <- read.delim("GANTC_250.txt", sep = "\t", header = TRUE) #one of the eight motifs
genome$Motif_count= as.numeric(as.character(genome$Motif_count)) #it was read as integer, but we need it to be numeric

#For each gene in cell regulation list, determine if the gene is part of an operon
#only keep those are not and those that are the first gene of an operon
i <- 0
#making a new dataframe
cell_new <- data.frame(matrix(nrow=1, ncol=1)) 
names(cell_new) <- c("Locus_tag")

for(n in 1:nrow(cell)) {    #runs for as many rows there are in the cell regulation gene list 
  pos1 <- match(cell[n,1], operon_all[,1])  #matching the first gene in cell regulation list to any row in operon list column 1
  pos2 <- match(cell[n,1], operon_first[,1]) 
  if(is.na(pos1)) { #if it is NA = not in position 1 (ie any operon)
    i <- i + 1
    cell_new[i,1] <- levels(cell$Locus_tag)[cell[n,1]] #put this gene in the new dataframe
  }
  if(!is.na(pos2)) { #if it is not NA for match in position 2 (ie. not not in first operon list = yes it is in first operon list)
    i <- i + 1
    cell_new[i,1] <- levels(cell$Locus_tag)[cell[n,1]] #put this gene in the new dataframe
  }
}

#merging the data to match the motif count etc back to the cell cycle gene list
x = merge(genome, cell_new, by.x="Gene_name", by.y="Locus_tag", all=FALSE)
cell_new = data.frame("Locus_tag"=x$Gene_name, "Motif_count"=x$Motif_count)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon)
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(cell_new)){ 
  if (cell_new[r,2]>0) { #row 1, column 2 is the motif number
    cell_new[r,3] = 1
    M.C=M.C+1
  } else {
    cell_new[r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
colnames(cell_new) = c("Locus_tag", "Motif_count", "Motif_Promoter")

#considering operons in GANTC 
genome <- genome[!(genome$Gene_name %in% cell$Locus_tag),]

#For each gene in genome, determine if the gene is part of an operon
#only keep those are not and those that are the first gene of an operon
I <- 0
genome_new<- data.frame(matrix(nrow=1, ncol=1)) 
names(genome_new) <- c("Locus_tag") #naming the column

for(n in 1:nrow(genome)) {    #runs for as many rows there are in the carbon gene list (ie 40)
  pos1 <- match(genome[n,1], operon_all[,1])  #matching the first gene in carbon list row n column 1 to any row in operon list column 1
  pos2 <- match(genome[n,1], operon_first[,1]) 
  if(is.na(pos1)) { #if it is NA = not in position 1 (ie any operon)
    I <- I + 1
    genome_new[I,1] <- levels(genome$Gene_name)[genome[n,1]] #put this gene in the new matrix
  }
  if(!is.na(pos2)) { #if it is not NA for match in position 2 (ie. not not in first operon list = yes it is in first operon list)
    I <- I + 1
    genome_new[I,1] <- levels(genome$Gene_name)[genome[n,1]] #put this gene in the new matrix as well
  }
}

#finding the motif_count for these genes
y = merge(genome_new, genome, by.x="Locus_tag",by.y ="Gene_name", all=FALSE)
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

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
Fisher.data = data.frame("Motif"=Motif, "No.Motif"=No.Motif)
row.names(Fisher.data) = c("Cell_cycle", "Genome")
Fisher.data
fisher.test(Fisher.data)

#T-test
n.motif=c(cell_new$Motif_count, genome_new$Motif_count) #need to line them up
repC = rep("C", nrow(cell_new))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#merging the original cell regulation document with cell_new to put the clusters back
XX= merge(cell, cell_new, by.x="Locus_tag", by.y="Locus_tag", all=FALSE)
#creating subset for the cluster groups
sub1 = subset(XX, Cluster == 1)
c1 = data.frame("Locus_tag"=sub1$Locus_tag, "Motif_count"=sub1$Motif_count, "Motif_Promoter"=sub1$Motif_Promoter)

sub2 = subset(XX, Cluster == 2)
c2 = data.frame("Locus_tag"=sub2$Locus_tag, "Motif_count"=sub2$Motif_count,"Motif_Promoter"=sub2$Motif_Promoter)

sub3 = subset(XX, Cluster == 3)
c3 = data.frame("Locus_tag"=sub3$Locus_tag, "Motif_count"=sub3$Motif_count, "Motif_Promoter"=sub3$Motif_Promoter)

sub4 = subset(XX, Cluster == 4)
c4 = data.frame("Locus_tag"=sub4$Locus_tag, "Motif_count"=sub4$Motif_count, "Motif_Promoter"=sub4$Motif_Promoter)

sub5 = subset(XX, Cluster == 5)
c5 = data.frame("Locus_tag"=sub5$Locus_tag, "Motif_count"=sub5$Motif_count, "Motif_Promoter"=sub5$Motif_Promoter)

sub6 = subset(XX, Cluster == 6)
c6 = data.frame("Locus_tag"=sub6$Locus_tag, "Motif_count"=sub6$Motif_count, "Motif_Promoter"=sub6$Motif_Promoter)
#all of these cluster groups already have the Motif_Promoter coloumn with 1 and 0 allocated to presence of motif

#cluster 1
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
repC = rep("C", nrow(c1))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 2
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
repC = rep("C", nrow(c2))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 3
#counting the number of genes with motif
subM.C = subset(c3, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c3)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c3$Motif_count, genome_new$Motif_count)
repC = rep("C", nrow(c3))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 4
#counting the number of genes with motif
subM.C = subset(c4, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c4)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c4$Motif_count, genome_new$Motif_count)
repC = rep("C", nrow(c4))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 5
#counting the number of genes with motif
subM.C = subset(c5, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c5)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c5$Motif_count, genome_new$Motif_count)
repC = rep("C", nrow(c5))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 6
#counting the number of genes with motif
subM.C = subset(c6, Motif_Promoter == 1)
M.C = as.numeric(nrow(subM.C))
nM.C = as.numeric(nrow(c6)) - M.C

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c6$Motif_count, genome_new$Motif_count)
repC = rep("C", nrow(c6))
repG = rep ("G", nrow(genome_new))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)
