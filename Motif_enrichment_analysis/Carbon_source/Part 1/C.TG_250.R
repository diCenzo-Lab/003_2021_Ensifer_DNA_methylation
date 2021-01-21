setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Loading Carbon source file
CarbonData=read.delim("Carbon_source_regulated_genes.txt", header=TRUE, sep="\t")
#Loading CGCANNNNNGTG_250.txt
C.TG_250=read.delim("CGCANNNNNGTG_250.txt", header=TRUE, sep="\t")
C.TG_250$Motif_count = as.numeric(as.character(C.TG_250$Motif_count))
head(C.TG_250)
str(C.TG_250)

#Making carbon Data frame
#merging the two datasets, only showinng the interacting variables (by.x and by.y) that are present in both
x = merge(C.TG_250, CarbonData, by.x="Gene_name", by.y="Locus_tag", all=FALSE)
x
dfCarbon = data.frame('Gene_name'=x$Gene_name, 'Motif_Count'=x$Motif_count) #taking out the two columns we want
dfCarbon 
str(dfCarbon)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon)
M.C=0 #number of genes with a motif in promoter
nM.C=0 #number of genes without a motif in promoter
for (i in 1:nrow(dfCarbon)){ 
  
  if (dfCarbon[r,2]>0) { #row 1, column 2 is the motif number
    dfCarbon[r,3] = 1
    M.C=M.C+1
  } else {
    dfCarbon [r,3] = 0
    nM.C=nM.C+1
  }
  r=r+1
}
M.C
nM.C

#naming the column
colnames(dfCarbon)=c("Gene_name","Motif_count", "Motif_Promoter")
dfCarbon #the final dfCarbon

#making the dfC.TG_250 dataframe
dfC.TG_250= data.frame('Gene_name'=C.TG_250$Gene_name, 'Motif_count'= C.TG_250$Motif_count)
dim(dfC.TG_250) #check number of lines
head(dfC.TG_250)

#remove all the genes that are in the carbon source list
dfC.TG_250<- dfC.TG_250[!(dfC.TG_250$Gene_name %in% dfCarbon$Gene_name),]
dim(dfC.TG_250) #checking- works - then continue

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
r=1
#counting the genes: M(motif) G(genome)
M.G=0 #number of genes with a motif in the promoter 
nM.G=0 #number of genes without a motif in the promoter (ie: no motif)
for (i in 1:nrow(dfC.TG_250)){
  
  if (dfC.TG_250[r,2]>0) {
    dfC.TG_250[r,3] = 1
    M.G=M.G+1
  } else {
    dfC.TG_250[r,3] = 0
    nM.G=nM.G+1
  }
  r=r+1
}
M.G
nM.G 
#naming the column
colnames(dfC.TG_250)=c("Gene_name","Motif_count", "Motif_Promoter")
head(dfC.TG_250)
dfC.TG_250 #the final

#Fisher's exact test
Motif=c(M.C, M.G)
No.Motif=c(nM.C, nM.G)

Fisher_Data=data.frame('Motif'=Motif, 'No.Motif'=No.Motif)
row.names(Fisher_Data) = c("Carbon", "Genome")
Fisher_Data
fisher.test(Fisher_Data)

#T-test 
#gether and arrange the data
n.motif=c(dfCarbon$Motif_count,dfC.TG_250$Motif_count) #need to line them up
t=data.frame ('Motif'=n.motif)
dim(t) #check to make sure it has the totall number of genes 
#assigning C (carbon) and G(Genome)
t[1:40,2]= c('C')
t[41:6256,2]= c('G')
tail(t) #check
colnames(t)=c("Motif", "Source")
head(t)
t.test(Motif ~ Source, data=t)

##
#Considering only the genes upregulated in Glucose
#Picking out the list of genes that are upregulated in Glucose
SubGlucose=subset(x,Cluster=='Glucose')
head(SubGlucose)
dfGlucose=data.frame('Gene_name'=SubGlucose$Gene_name, 'Motif_count'=SubGlucose$Motif_count)
head(dfGlucose) #just checking

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
rG=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon) G(Glucose)
M.CG=0 #number of genes with a motif in promoter
nM.CG=0 #number of genes without a motif in promoter
for (i in 1:nrow(dfGlucose)){ 
  if (dfGlucose[rG,2]>0) { #row 1, column 2 is the motif number
    dfGlucose[rG,3] = 1
    M.CG=M.CG+1
  } else {
    dfGlucose [rG,3] = 0
    nM.CG=nM.CG+1
  }
  rG=rG+1
}
M.CG
nM.CG
#naming the column
colnames(dfGlucose)=c("Gene_name","Motif_count", "Motif_Promoter")
dfGlucose

#Fisher's exact test
#gather the data
G_Motif=c(M.CG, M.G)
G_No.Motif=c(nM.CG, nM.G)
G.Fisher_Data=data.frame('Motif'=G_Motif, 'No.Motif'=G_No.Motif)
row.names(G.Fisher_Data) = c("Carbon", "Genome")
G.Fisher_Data
fisher.test(G.Fisher_Data) 

#T-test 
#gether and arrange the data
G.n.motif=c(dfGlucose$Motif_count,dfC.TG_250$Motif_count) #need to line them up
Gt=data.frame ('Motif'=G.n.motif)
dim(Gt) #check to make sure it has the totall number of genes 
#assigning C (carbon) and G(Genome)
dim(dfGlucose)
Gt[1:29,2]= c('C')
Gt[30:6245,2]= c('G')
tail(Gt) #check
colnames(Gt)=c("Motif", "Source")
head(Gt)
t.test(Motif ~ Source, data=Gt)

##
#Considering only the genes upregulated in Succinate
#Picking out the list of genes that are upregulated in Succinate
SubSucc=subset(x,Cluster=='Succinate')
head(SubSucc)
dfSucc=data.frame('Gene_name'=SubSucc$Gene_name, 'Motif_count'=SubSucc$Motif_count)
head(dfSucc)

#allocating 1 and 0 to the third coloumn, indicating the presence of motif in promoter region
rS=1 #starting at row 1 of the dataframe
#counting the genes: M(motif) C(Carbon) G(Glucose)
M.CS=0 #number of genes with a motif in promoter
nM.CS=0 #number of genes without a motif in promoter

for (i in 1:nrow(dfSucc)){ 
  if (dfSucc[rS,2]>0) { #row 1, column 2 is the motif number
    dfSucc[rS,3] = 1
    M.CS=M.CS+1
  } else {
    dfSucc [rS,3] = 0
    nM.CS=nM.CS+1
  }
  rS=rS+1
}
M.CS
nM.CS
#naming the column
colnames(dfSucc)=c("Gene_name","Motif_count", "Motif_Promoter")
dfSucc #final

#Fisher's exact test
#gather the data
S_Motif=c(M.CS, M.G)
S_No.Motif=c(nM.CS, nM.G)
S.Fisher_Data=data.frame('Motif'=S_Motif, 'No.Motif'=S_No.Motif)
row.names(S.Fisher_Data) = c("Carbon", "Genome")
S.Fisher_Data
fisher.test(S.Fisher_Data)

#T-test 
#gether and arrange the data
S.n.motif=c(dfSucc$Motif_count,dfC.TG_250$Motif_count) #need to line them up
St=data.frame ('Motif'=S.n.motif)
dim(St) #check to make sure it has the totall number of genes
#assigning C (carbon) and G(Genome)
dim(dfSucc)
St[1:11,2]= c('C')
St[12:6227,2]= c('G')
tail(St) #check
colnames(St)=c("Motif", "Source")
head(St)
t.test(Motif ~Source, data=St)
