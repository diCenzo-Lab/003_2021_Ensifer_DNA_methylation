setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
cell <- read.delim("Cell_cycle_regulated_genes.txt", sep = "\t", header = TRUE) #cell cycle genes
GANTC <- read.delim("GANTC_250.txt", sep = "\t", header = TRUE)

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
M.C1=0 
nM.C1=0 
for (i in 1:nrow(c1)){ 
  if (c1[r,2]>0) { #row 1, column 2 is the motif number
    c1 [r,3] = 1
    M.C1=M.C1+1
  } else {
    c1 [r,3] = 0
    nM.C1=nM.C1+1
  }
  r=r+1
}

#naming the column
colnames(c1)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.1=c(M.C1, M.G)
No.Motif.1=c(nM.C1, nM.G)
f1=data.frame('Motif'=Motif.1, 'No.Motif'=No.Motif.1)
row.names(f1) = c("Cell_cycle", "Genome")
f1
fisher.test(f1) 

#T-test 
c1.motif=c(c1$Motif_count, dfGANTC$Motif_count)
repC1 = rep("C", nrow(c1))
repG = rep ("G", nrow(dfGANTC))
source = c(repC1, repG)
t1 = data.frame("Motif"=c1.motif, "Source"=source)
t.test(Motif ~ Source, data=t1)

#cluster 2 analysis
#allocating 1 and 0 to the third coloumn
r= 1
M.C2=0 
nM.C2=0 
for (i in 1:nrow(c2)){ 
  if (c2[r,2]>0) { #row 1, column 2 is the motif number
    c2 [r,3] = 1
    M.C2=M.C2+1
  } else {
    c2 [r,3] = 0
    nM.C2=nM.C2+1
  }
  r=r+1
}
#naming the column
colnames(c2)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.2=c(M.C2, M.G)
No.Motif.2=c(nM.C2, nM.G)
f2=data.frame('Motif'=Motif.2, 'No.Motif'=No.Motif.2)
row.names(f2) = c("Cell_cycle", "Genome")
f2
fisher.test(f2) 

#T-test 
c2.motif=c(c2$Motif_count, dfGANTC$Motif_count)
repC2 = rep("C", nrow(c2))
repG = rep ("G", nrow(dfGANTC))
source = c(repC2, repG)
t2 = data.frame("Motif"=c2.motif, "Source"=source)
t.test(Motif ~ Source, data=t2)

#cluster 3
r= 1
M.C3=0 
nM.C3=0 
for (i in 1:nrow(c3)){ 
  if (c3[r,2]>0) { #row 1, column 2 is the motif number
    c3 [r,3] = 1
    M.C3=M.C3+1
  } else {
    c3 [r,3] = 0
    nM.C3=nM.C3+1
  }
  r=r+1
}
#naming the column
colnames(c3)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.3=c(M.C3, M.G)
No.Motif.3=c(nM.C3, nM.G)
f3=data.frame('Motif'=Motif.3, 'No.Motif'=No.Motif.3)
row.names(f3) = c("Cell_cycle", "Genome")
f3
fisher.test(f3) 

#T-test 
c3.motif=c(c3$Motif_count, dfGANTC$Motif_count)
repC3 = rep("C", nrow(c3))
repG = rep ("G", nrow(dfGANTC))
source = c(repC3, repG)
t3 = data.frame("Motif"=c3.motif, "Source"=source)
t.test(Motif ~ Source, data=t3)

#cluster 4
#allocating 1 and 0 to the third coloumn
r= 1
M.C4=0 
nM.C4=0 
for (i in 1:nrow(c4)){ 
  if (c4[r,2]>0) { #row 1, column 2 is the motif number
    c4 [r,3] = 1
    M.C4=M.C4+1
  } else {
    c4 [r,3] = 0
    nM.C4=nM.C4+1
  }
  r=r+1
}

#naming the column
colnames(c4)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.4=c(M.C4, M.G)
No.Motif.4=c(nM.C4, nM.G)
f4=data.frame('Motif'=Motif.4, 'No.Motif'=No.Motif.4)
row.names(f4) = c("Cell_cycle", "Genome")
f4
fisher.test(f4) 

#T-test 
c4.motif=c(c4$Motif_count, dfGANTC$Motif_count)
repC4 = rep("C", nrow(c4))
repG = rep ("G", nrow(dfGANTC))
source = c(repC4, repG)
t4 = data.frame("Motif"=c4.motif, "Source"=source)
t.test(Motif ~ Source, data=t4)

#cluster 5
r= 1
M.C5=0 
nM.C5=0 
for (i in 1:nrow(c5)){ 
  if (c5[r,2]>0) { #row 1, column 2 is the motif number
    c5 [r,3] = 1
    M.C5=M.C5+1
  } else {
    c5 [r,3] = 0
    nM.C5=nM.C5+1
  }
  r=r+1
}
#naming the column
colnames(c5)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.5=c(M.C5, M.G)
No.Motif.5=c(nM.C5, nM.G)
f5=data.frame('Motif'=Motif.5, 'No.Motif'=No.Motif.5)
row.names(f5) = c("Cell_cycle", "Genome")
f5
fisher.test(f5) 

#T-test 
c5.motif=c(c5$Motif_count, dfGANTC$Motif_count)
repC5 = rep("C", nrow(c5))
repG = rep ("G", nrow(dfGANTC))
source = c(repC5, repG)
t5 = data.frame("Motif"=c5.motif, "Source"=source)
t.test(Motif ~ Source, data=t5)

#cluster 6
r= 1
M.C6=0 
nM.C6=0 
for (i in 1:nrow(c6)){ 
  if (c6[r,2]>0) { #row 1, column 2 is the motif number
    c6 [r,3] = 1
    M.C6=M.C6+1
  } else {
    c6 [r,3] = 0
    nM.C6=nM.C6+1
  }
  r=r+1
}
#naming the column
colnames(c6)=c("Gene_name","Motif_count", "Motif_Promoter")

#Fisher's exact test
Motif.6=c(M.C6, M.G)
No.Motif.6=c(nM.C6, nM.G)
f6=data.frame('Motif'=Motif.6, 'No.Motif'=No.Motif.6)
row.names(f6) = c("Cell_cycle", "Genome")
f6
fisher.test(f6) 

#T-test 
c6.motif=c(c6$Motif_count, dfGANTC$Motif_count)
repC6 = rep("C", nrow(c6))
repG = rep ("G", nrow(dfGANTC))
source = c(repC6, repG)
t6 = data.frame("Motif"=c6.motif, "Source"=source)
t.test(Motif ~ Source, data=t6)



