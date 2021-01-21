setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
cell <- read.delim("Cell_cycle_regulated_genes.txt", sep = "\t", header = TRUE) #cell cycle genes
GANTC <- read.delim("GANTC_100.txt", sep = "\t", header = TRUE)

#the cell cycle dataframe
#merging the data to match the motif count etc back to the cell cycle gene list
x = merge(GANTC, cell, by.x="Gene_name", by.y="Locus_tag", all=FALSE)
#taking out the gene name and motif count column to make a dataframe
dfcell= data.frame("Locus_tag"=x$Gene_name, "Motif_count"=x$Motif_count)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(cell)
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(dfcell)){ 
  if (dfcell[r,2]>0) { #row 1, column 2 is the motif number
    dfcell[r,3] = 1
    M.C=M.C+1
  } else {
    dfcell [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
#naming the column
colnames(dfcell)= c("Locus_tag", "Motif_count", "Motif_Promoter")

#GANTC dataframe
#taking out gene name and motif count
dfGANTC = data.frame("Gene_name"=GANTC$Gene_name, "Motif_count"=GANTC$Motif_count)
#remove all the genes in the cell cycle list
dfGANTC <- dfGANTC[!(dfGANTC$Gene_name %in% dfcell$Locus_tag),]

#allocating 1 and 0 to third column for presence of GANTC in promoter (if motif count >0)
R = 1
M.G = 0
nM.G = 0
for (i in 1:nrow(dfGANTC)) {
  if (dfGANTC[R,2]>0) {
    dfGANTC [R,3] = 1
    M.G = M.G +1
  } else {
    dfGANTC [R,3] = 0
    nM.G = nM.G +1
  }
  R = R + 1
}
#naming the columns
colnames(dfGANTC)= c("Locus_tag", "Motif_count", "Motif_Promoter")

#Fisher's test overall  
Motif = c(M.C, M.G)
No.Motif = c(nM.C, nM.G)
Fisher.data = data.frame("Motif"=Motif, "No.Motif"=No.Motif)
row.names(Fisher.data) = c("Cell_cycle", "Genome")
Fisher.data
fisher.test(Fisher.data)

#T-test
n.motif=c(dfcell$Motif_count, dfGANTC$Motif_count) #need to line them up
repC = rep("C", nrow(dfcell))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#create subset for clusters
sub1 = subset(x, Cluster == 1)
c1 = data.frame("Locus_tag"=sub1$Gene_name, "Motif_count"=sub1$Motif_count)

sub2 = subset(x, Cluster == 2)
c2 = data.frame("Locus_tag"=sub2$Gene_name, "Motif_count"=sub2$Motif_count)

sub3 = subset(x, Cluster == 3)
c3 = data.frame("Locus_tag"=sub3$Gene_name, "Motif_count"=sub3$Motif_count)

sub4 = subset(x, Cluster == 4)
c4 = data.frame("Locus_tag"=sub4$Gene_name, "Motif_count"=sub4$Motif_count)

sub5 = subset(x, Cluster == 5)
c5 = data.frame("Locus_tag"=sub5$Gene_name, "Motif_count"=sub5$Motif_count)

sub6 = subset(x, Cluster == 6)
c6 = data.frame("Locus_tag"=sub6$Gene_name, "Motif_count"=sub6$Motif_count)

#cluster 1 analysis 
#allocating 1 and 0 to the third coloumn
r= 1
M.C=0 
nM.C=0 
for (i in 1:nrow(c1)){ 
  if (c1[r,2]>0) { #row 1, column 2 is the motif number
    c1 [r,3] = 1
    M.C=M.C+1
  } else {
    c1 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}

#naming the column
colnames(c1)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c1$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c1))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 2
#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(cell)
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(c2)){ 
  if (c2[r,2]>0) { #row 1, column 2 is the motif number
    c2[r,3] = 1
    M.C=M.C+1
  } else {
    c2 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
#naming the column
colnames(c2)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c2$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c2))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 3
#allocating 1 and 0 to the third coloumn
r= 1
M.C=0 
nM.C=0 
for (i in 1:nrow(c3)){ 
  if (c3[r,2]>0) { #row 1, column 2 is the motif number
    c3 [r,3] = 1
    M.C=M.C+1
  } else {
    c3 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}

#naming the column
colnames(c3)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c3$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c3))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 4
#allocating 1 and 0 to the third coloumn
r= 1
M.C=0 
nM.C=0 
for (i in 1:nrow(c4)){ 
  if (c4[r,2]>0) { #row 1, column 2 is the motif number
    c4 [r,3] = 1
    M.C=M.C+1
  } else {
    c4 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}

#naming the column
colnames(c4)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c4$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c4))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 5
#allocating 1 and 0 to the third coloumn
r= 1
M.C=0 
nM.C=0 
for (i in 1:nrow(c5)){ 
  if (c5[r,2]>0) { #row 1, column 2 is the motif number
    c5 [r,3] = 1
    M.C=M.C+1
  } else {
    c5 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}

#naming the column
colnames(c5)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c5$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c5))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#cluster 6
#allocating 1 and 0 to the third coloumn
r= 1
M.C=0 
nM.C=0 
for (i in 1:nrow(c6)){ 
  if (c6[r,2]>0) { #row 1, column 2 is the motif number
    c6 [r,3] = 1
    M.C=M.C+1
  } else {
    c6 [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}

#naming the column
colnames(c6)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c6$Motif_count, dfGANTC$Motif_count)
repC = rep("C", nrow(c6))
repG = rep ("G", nrow(dfGANTC))
source = c(repC, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)



