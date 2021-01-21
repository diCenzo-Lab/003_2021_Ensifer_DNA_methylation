setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
sym <- read.delim("Symbiosis_regulated_genes.txt", sep = "\t", header = TRUE) #symbiosis
genome <- read.delim("GANTC_100.txt", sep = "\t", header = TRUE)
sym = unique(sym)

#Symbiosis list
#merging the data to match the motif count 
x = merge(genome, sym, by.x="Gene_name", by.y="Locus_tag", all=FALSE)
#XX <- sym[!(sym$Locus_tag %in% x$Gene_name),] #looking for the unmatched ones
dfsym = data.frame("Locus_tag"=x$Gene_name, "Motif_count"=x$Motif_count)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) S(symbiosis)
M.S=0 #number of genes with a motif in promoter
nM.S=0 #number of genes without a motif in promoter
for (i in 1:nrow(dfsym)){ 
  if (dfsym[r,2]>0) { #row 1, column 2 is the motif number
    dfsym[r,3] = 1
    M.S=M.S+1
  } else {
    dfsym [r,3] = 0
    nM.S=nM.S+1
  }
  r=r+1
}
#naming the column
colnames(dfsym)= c("Locus_tag", "Motif_count", "Motif_Promoter")

#genome
#taking out gene name and motif count
dfgenome = data.frame("Gene_name"=genome$Gene_name, "Motif_count"=genome$Motif_count)
#remove all the genes in the symbiosis list 
dfgenome <- dfgenome[!(dfgenome$Gene_name %in% sym$Locus_tag),]

#allocating 1 and 0 to third column for presence of GANTC in promoter (if motif count >0)
R = 1
M.G = 0
nM.G = 0
for (i in 1:nrow(dfgenome)) {
  if (dfgenome[R,2]>0) {
    dfgenome [R,3] = 1
    M.G = M.G +1
  } else {
    dfgenome [R,3] = 0
    nM.G = nM.G +1
  }
  R = R + 1
}
#naming the columns
colnames(dfgenome)= c("Locus_tag", "Motif_count", "Motif_Promoter")

#Fisher's test overall  
Motif = c(M.S, M.G)
No.Motif = c(nM.S, nM.G)
Fisher.data = data.frame("Motif"=Motif, "No.Motif"=No.Motif)
row.names(Fisher.data) = c("Symbiosis", "Genome")
Fisher.data
fisher.test(Fisher.data)

#T-test
n.motif=c(dfsym$Motif_count, dfgenome$Motif_count) #need to line them up
repS = rep("S", nrow(dfsym))
repG = rep ("G", nrow(dfgenome))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#dealing with the Down and Up clusters
#merging the dfsym back to original sym with cluster Down and Up
XX = merge(dfsym, sym, by.x = "Locus_tag", by.y = "Locus_tag")

sub1 = subset(XX, Cluster == "Down")
c1 = data.frame("Locus_tag"=sub1$Locus_tag, "Motif_count"=sub1$Motif_count, "Motif_Promoter"=sub1$Motif_Promoter)

sub2 = subset(XX, Cluster =="Up")
c2 = data.frame("Locus_tag"=sub2$Locus_tag, "Motif_count"=sub2$Motif_count, "Motif_Promoter"=sub2$Motif_Promoter)

#Down (sub 1)
#counting the number of genes with motif
subM.S = subset(c1, Motif_Promoter == 1)
M.S = as.numeric(nrow(subM.S))
nM.S = as.numeric(nrow(c1)) - M.S

#Fisher's exact test
Motif=c(M.S, M.G)
No.Motif=c(nM.S, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c1$Motif_count, dfgenome$Motif_count)
repS = rep("S", nrow(c1))
repG = rep ("G", nrow(dfgenome))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

#Up (sub 2)
#counting the number of genes with motif
subM.S = subset(c2, Motif_Promoter == 1)
M.S = as.numeric(nrow(subM.S))
nM.S = as.numeric(nrow(c2)) - M.S

#Fisher's exact test
Motif=c(M.S, M.G)
No.Motif=c(nM.S, nM.G)
f=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(f) = c("Cell_cycle", "Genome")
f
fisher.test(f) 

#T-test 
n.motif=c(c2$Motif_count, dfgenome$Motif_count)
repS = rep("S", nrow(c2))
repG = rep ("G", nrow(dfgenome))
source = c(repS, repG)
t = data.frame("Motif"=n.motif, "Source"=source)
t.test(Motif ~ Source, data=t)

