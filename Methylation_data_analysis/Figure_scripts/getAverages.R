davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt),
                   V6 = length(d1$Start_nt), V7 = length(d1$Start_nt), V8 = length(d1$Start_nt), V9 = length(d1$Start_nt), V10 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
  davg[i,7] <- (d1$Motif_count[i] + d2$Motif_count[i] + d3$Motif_count[i]) / 3
  davg[i,8] <- d4$Motif_count[i]
  davg[i,9] <- d5$Motif_count[i]
  davg[i,10] <- (d6$Motif_count[i] + d7$Motif_count[i] + d8$Motif_count[i]) / 3
  davg[i,11] <- (d9$Motif_count[i] + d10$Motif_count[i] + d11$Motif_count[i]) / 3
  
  davg[i,12] <- (d12$Extent_methylation[i] + d13$Extent_methylation[i] + d14$Extent_methylation[i]) / 3
  davg[i,13] <- (d12$Motif_count[i] + d13$Motif_count[i] + d14$Motif_count[i]) / 3
}

sum(davg[,2] * davg[,7], na.rm=TRUE) / sum(davg[,7], na.rm=TRUE)
sum(davg[,3] * davg[,8], na.rm=TRUE) / sum(davg[,8], na.rm=TRUE)
sum(davg[,4] * davg[,9], na.rm=TRUE) / sum(davg[,9], na.rm=TRUE)
sum(davg[,5] * davg[,10], na.rm=TRUE) / sum(davg[,10], na.rm=TRUE)
sum(davg[,6] * davg[,11], na.rm=TRUE) / sum(davg[,11], na.rm=TRUE)

sum(davg[,12] * davg[,13], na.rm=TRUE) / sum(davg[,13], na.rm=TRUE)
