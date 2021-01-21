#Making a new table
setwd("/Users/janischeng/Desktop/Methylation Project/Analysis")
#Load the files
operon_first <- read.delim("Operons_First_Gene.txt", sep = "\t", header = TRUE) #contains only the first operon
operon_all <- read.delim("Operons_Reorganized.txt", sep = "\t", header = TRUE) #contains all the operons, column 1= gene name, column 2= operon groups
cell <- read.delim("Cell_cycle_regulated_genes.txt", sep = "\t", header = TRUE) #cell cycle genes
GANTC_250 <- read.delim("GANTC_250.txt", sep = "\t", header = TRUE)
GANTC_100 <- read.delim("GANTC_100.txt", sep = "\t", header = TRUE)
sub.100= subset(GANTC_100,Motif_count > 0)
sub.250 = subset(GANTC_250, Motif_count > 0)

cell_new=data.frame("Locus_tag"=cell$Locus_tag, "Cluster"=cell$Cluster)

#adding a column to indicate whether the genes are: not in operon(0), first in operon(1), others(2)
i = 0
for(n in 1:nrow(cell)) {   
  pos1 <- match(cell[n,1], operon_all[,1])  #matching with column 1 to any row in operon list column 1
  pos2 <- match(cell[n,1], operon_first[,1]) #matching with list of first in operons
  if(is.na(pos1)) { #if it is NA = not in position 1 (ie not in any operon)
    i <- i + 1
    cell_new[i,3] = 0 #labeled as "0"
  }
  else if(!is.na(pos2)) { #if it is not NA for match in position 2 (ie. not not in first operon list = yes it is in first operon list)
    i <- i + 1
    cell_new[i,3] = 1 #labeled as "1"
  } else { #everything else
    i <- i + 1
    cell_new[i,3] = 2 #labeled as "2"
  }
}
head(cell_new)
#cell_new = na.omit(cell_new) #no need for this now

#Column 4: whether it is in GANTC 100 with a motif count of < 0
I = 0
n1 = 0
n2 = 0
for(n in 1:nrow(cell)) {   
  pos3 <- match(cell[n,1], sub.100[,1])  #matching with column 1 to any row in operon list column 1
  if(!is.na(pos3)) { #not NA = in the list
    I <- I + 1
    n1 = n1 + 1
    cell_new[I,4] = 1 #labeled as "1"
  } 
  if(is.na(pos3)) { 
    I <- I + 1
    n2 = n2 +1
    cell_new[I,4] = 0 #labeled as "0"
  }
}

#Column 5: whether it is in GANTC 250 with a motif count of < 0
I = 0
n3 = 0
n4 = 0
for(n in 1:nrow(cell)) {   
  pos4 <- match(cell[n,1], sub.250[,1])  #matching with column 1 to any row in operon list column 1
  if(!is.na(pos4)) { #not NA = in the list
    I <- I + 1
    n3 = n3 + 1
    cell_new[I,5] = 1 #labeled as "1"
  } 
  if(is.na(pos4)) { 
    I <- I + 1
    n4 = n4 +1
    cell_new[I,5] = 0 #labeled as "0"
  }
}

#naming the new column
colnames(cell_new)=c("Locus_tag", "Cluster", "Operon_status", "GANTC_100", "GANTC_250")
head (cell_new)

write.table(cell_new, file = "New_cell_cycle", append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE, qmethod = c("escape", "double"))


