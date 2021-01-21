setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
operon_first <- read.delim("Operons_First_Gene.txt", sep = "\t", header = TRUE) #contains only the first operon
operon_all <- read.delim("Operons_Reorganized.txt", sep = "\t", header = TRUE) #contains all the operons, column 1= gene name, column 2= operon groups
carbon <- read.delim("Carbon_source_regulated_genes.txt", sep = "\t", header = TRUE) #carbon source gene list - 40 genes
genome <- read.delim("GANTC_100.txt", sep = "\t", header = TRUE) #one of the eight motifs
genome$Motif_count= as.numeric(as.character(genome$Motif_count)) #it was read as integer, but we need it to be numeric

#For each gene in carbon source list, determine if the gene is part of an operon
#only keep those are not and those that are the first gene of an operon
i <- 0
#making a new dataframe, this will have all the genes that are either not in operon or is the first gene in operons
carbon_new <- data.frame(matrix(nrow=1, ncol=1)) 
names(carbon_new) <- c("Locus_tag") #naming the column

for(n in 1:nrow(carbon)) {    #runs for as many rows there are in the carbon gene list (ie 40)
  pos1 <- match(carbon[n,1], operon_all[,1])  #matching the first gene in carbon list row n column 1 to any row in operon list column 1
  pos2 <- match(carbon[n,1], operon_first[,1]) 
  if(is.na(pos1)) { #if it is NA = not in position 1 (ie any operon)
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]] #put this gene in the new matrix
  }
  if(!is.na(pos2)) { #if it is not NA for match in position 2 (ie. not not in first operon list = yes it is in first operon list)
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]] #put this gene in the new matrix as well
  }
}
carbon_new #checking 
dim(carbon_new) #check - should be less than 40 

#finding the motif_count for these genes
#merging, only showing the interacting variables ("by.x" and "by.y"- the gene names) that are present in both
x = merge(carbon_new, genome, by.x="Locus_tag", by.y="Gene_name", all=FALSE)
x
#taking out the gene name and motif count column to make a dataframe
carbon_new = data.frame('Locus_tag'=x$Locus_tag, 'Motif_count'=x$Motif_count)
head(carbon_new) #checking

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon)
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(carbon_new)){ 
  if (carbon_new[r,2]>0) { #row 1, column 2 is the motif number
    carbon_new[r,3] = 1
    M.C=M.C+1
  } else {
    carbon_new[r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
M.C
nM.C
#naming the column
colnames(carbon_new)=c("Locus_tag","Motif_count", "Motif_Promoter")
carbon_new #the final 

#working on the genome list
#remove all the genes that are in the carbon source list
genome<- genome[!(genome$Gene_name %in% carbon$Locus_tag),]
dim(genome)
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
dim(genome_new) #check - should be less than total
tail(genome_new)

#finding the motif_count for these genes
y = merge(genome_new, genome, by.x="Locus_tag",by.y ="Gene_name", all=FALSE)
dim(y)
genome_new =data.frame('Locus_tag'=y$Locus_tag,'Motif_count'=y$Motif_count)
str(genome_new)

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
head(genome_new)
dim(genome_new)
genome_new #the final

#Fisher's exact test
#gather the data needed:#GANTC w motif, #GANTC w/out motif, #carbon w motif, #carbon w/out motif
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)

Fisher_Data=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(Fisher_Data) = c("Carbon", "Genome")
Fisher_Data
fisher.test(Fisher_Data)

#T-test 
#gether and arrange the data
n.motif=c(carbon_new$Motif_count,genome_new$Motif_count) #need to line up number of motifs
t=data.frame ('Motif'=n.motif)
dim(t)
#assigning C (carbon) and G(Genome)
t[1:22,2]= c('C')
t[23:3847,2]= c('G')
tail(t) #check
colnames(t)=c("Motif", "Source")
head(t)
t.test(Motif ~ Source, data=t)

##
#Considering only the genes upregulated in Glucose
SubGlucose=subset(carbon,Cluster=='Glucose') #Picking out the list of genes from carbon that are Cluster=glucose
SubGlucose #there should be 29
#merge it with the list of carbon that are selected to be non-operon or first in operon
g=merge(carbon_new, SubGlucose, by.x = "Locus_tag", by.y = "Locus_tag") 
g[,4]=NULL #removing the 4th column (which is the cluster column indicating glucose or succinate)
g #dataframe: glucose upregulated + non-operon/first in operon

#counting number of genes with motifs in promoter
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon)
M.CG=0 #number of genes with a motif in promoter
nM.CG=0 #number of genes without a motif in promoter
for (i in 1:nrow(g)){ 
  if (g[r,3]>0) { #row 1, column 2 is the motif number
    M.CG=M.CG+1
  } else {
    nM.CG=nM.CG+1
  }
  r=r+1
}

#Fisher's exact test
G_Motif=c(M.CG, M.G)
G_No.Motif=c(nM.CG, nM.G)
G.Fisher_Data=data.frame('Motif'=G_Motif, 'No.Motif'=G_No.Motif)
row.names(G.Fisher_Data) = c("Carbon", "Genome")
G.Fisher_Data
fisher.test(G.Fisher_Data) 

#T-test 
#gether and arrange the data
G.n.motif=c(g$Motif_count,genome_new$Motif_count) #need to line them up
Gt=data.frame ('Motif'=G.n.motif)
dim(Gt)
#assigning C (carbon) and G(Genome)
dim(g)
Gt[1:15,2]= c('C')
Gt[16:3840,2]= c('G')
tail(Gt) #check
colnames(Gt)=c("Motif", "Source")
head(Gt)
t.test(Motif ~ Source, data=Gt)

##
#Considering only the genes upregulated in Succinate
SubSucc=subset(carbon,Cluster=='Succinate')
dim(SubSucc) #should be 11
#merge it with the list of carbon that are selected to be non-operon or first in operon
s=merge(carbon_new, SubSucc, by.x = "Locus_tag", by.y = "Locus_tag") 
s[,4]=NULL
s #final succ dataframe

#counting number of genes with motifs in promoter
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon)
M.CS=0 #number of genes with a motif in promoter
nM.CS=0 #number of genes without a motif in promoter
for (i in 1:nrow(s)){ 
  if (s[r,3]>0) { #row 1, column 2 is the motif number
    M.CS=M.CS+1
  } else {
    nM.CS=nM.CS+1
  }
  r=r+1
}

#Fisher's exact test
S_Motif=c(M.CS, M.G)
S_No.Motif=c(nM.CS, nM.G)
S.Fisher_Data=data.frame('Motif'=S_Motif, 'No.Motif'=S_No.Motif)
row.names(S.Fisher_Data) = c("Carbon", "Genome")
S.Fisher_Data
fisher.test(S.Fisher_Data) 

#T-test 
#gether and arrange the data
S.n.motif=c(s$Motif_count,genome_new$Motif_count) #need to line them up
St=data.frame ('Motif'=S.n.motif)
dim(St) #check to make sure it has the total number of genes 
#assigning C (carbon) and G(Genome)
dim(s)
St[1:7,2]= c('C')
St[8:3832,2]= c('G')
tail(St) #check
colnames(St)=c("Motif", "Source")
head(St)
t.test(Motif ~ Source, data=St)


