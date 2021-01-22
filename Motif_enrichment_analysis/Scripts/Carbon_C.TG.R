#Load the files
# To consider non-operon genes and first gene of operon for all genes in an operon (i.e., analysis part three)
operon_first <- read.delim("Operons_First_Gene.txt", sep = "\t", header = TRUE)
operon_all <- read.delim("Operons_Reorganized.txt", sep = "\t", header = TRUE)
carbon <- read.delim("Carbon_source_regulated_genes.txt", sep = "\t", header = TRUE)
genome <- read.delim("CGCANNNNNGTG_125.txt", sep = "\t", header = TRUE)
genome$Motif_count= as.numeric(as.character(genome$Motif_count)) #it was read as integer, but we need it to be numeric

i <- 0 
carbon_new <- data.frame(matrix(nrow=1, ncol=1))
names(carbon_new) <- c("Locus_tag")
for(n in 1:nrow(carbon)) {
  pos1 <- match(carbon[n,1], operon_all[,1]) #match with all the genes that s part of an operon
  if(is.na(pos1)) {
    i <- i + 1
    carbon_new[i,1] <- levels(carbon$Locus_tag)[carbon[n,1]] #put it in the list
  } else { #if in the operon list
    i <- i + 1
    pos2 <- match(operon_all[pos1,2], operon_first[,2])  #find the first gene of that operon
    carbon_new[i,1] <- levels(operon_first$Gene)[operon_first[pos2,1]] #put the first gene of that operon in the list
  }
}
carbon_new <- unique(carbon_new) #remove duplicated

#finding the motif_count for these genes
#merging, only showing the interacting variables ("by.x" and "by.y"- the gene names) that are present in both
x = merge(carbon_new, genome, by.x="Locus_tag", by.y="Gene_name", all=FALSE)
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

#naming the column
colnames(carbon_new)=c("Locus_tag","Motif_count", "Motif_Promoter")

#working on the genome list
#remove all the genes that are in the carbon source list
genome<- genome[!(genome$Gene_name %in% carbon$Locus_tag),]

I<-0
genome_new<-data.frame(matrix(nrow = 1, ncol = 1))
names(genome_new)=c("Locus_tag")

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

#merging to look for all the motif count
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
head(genome_new)
genome_new #the final

#Fisher's exact test
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
t[1:25,2]= c('C')
t[26:3851,2]= c('G')
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
Gt[16:3841,2]= c('G')
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
St[8:3833,2]= c('G')
tail(St) #check
colnames(St)=c("Motif", "Source")
head(St)
t.test(Motif ~ Source, data=St)

