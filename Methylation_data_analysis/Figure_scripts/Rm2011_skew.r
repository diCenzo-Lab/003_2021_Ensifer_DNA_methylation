require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# 2011 chromosome succinate
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d5$Extent_methylation[i] + d6$Extent_methylation[i] + d7$Extent_methylation[i]) / 3
  davg[i,4] <- (d8$Extent_methylation[i] + d9$Extent_methylation[i] + d10$Extent_methylation[i]) / 3
  davg[i,5] <- (d11$Extent_methylation[i] + d12$Extent_methylation[i] + d13$Extent_methylation[i]) / 3
}
plot1 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)
d4b <- data.frame(V1 = length(d4$Start_nt), V2 = length(d4$End_nt))
for(i in 1:length(d4$Start_nt)) {
  d4b[i,1] <- d4$End_nt[i] / 1000000
  if(i == 1) {
    d4b[i,2] <- d4$GC_skew[i]
  }
  else {
    d4b[i,2] <- d4b$V2[i-1] + d4$GC_skew[i]
  }
}
plot2 <- ggplot(d4b, aes(V1, V2)) +
  geom_point(x = d4b$V1, y = d4b$V2, size=0.1, color="black") +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  scale_y_continuous(limits = c(-2,3),
                     labels=scaleFUN,
                     "Cumulative GC skew"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 0.4)
plot7 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)

# 2011 pSymB succinate
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_GANTC_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_GANTC_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_GANTC_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_CGCANNNNNGTG_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_CGCANNNNNGTG_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_CGCANNNNNGTG_CP004139_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d5$Extent_methylation[i] + d6$Extent_methylation[i] + d7$Extent_methylation[i]) / 3
  davg[i,4] <- (d8$Extent_methylation[i] + d9$Extent_methylation[i] + d10$Extent_methylation[i]) / 3
  davg[i,5] <- (d11$Extent_methylation[i] + d12$Extent_methylation[i] + d13$Extent_methylation[i]) / 3
}
plot3 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)
d4b <- data.frame(V1 = length(d4$Start_nt), V2 = length(d4$End_nt))
for(i in 1:length(d4$Start_nt)) {
  d4b[i,1] <- d4$End_nt[i] / 1000000
  if(i == 1) {
    d4b[i,2] <- d4$GC_skew[i]
  }
  else {
    d4b[i,2] <- d4b$V2[i-1] + d4$GC_skew[i]
  }
}
plot4 <- ggplot(d4b, aes(V1, V2)) +
  geom_point(x = d4b$V1, y = d4b$V2, size=0.1, color="black") +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  scale_y_continuous(limits = c(-2,3),
                     labels=scaleFUN,
                     "Cumulative GC skew"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 0.4)
plot8 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)

# 2011 pSymA succinate
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_GANTC_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_GANTC_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_GANTC_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_CGCANNNNNGTG_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_CGCANNNNNGTG_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_CGCANNNNNGTG_CP004138_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d5$Extent_methylation[i] + d6$Extent_methylation[i] + d7$Extent_methylation[i]) / 3
  davg[i,4] <- (d8$Extent_methylation[i] + d9$Extent_methylation[i] + d10$Extent_methylation[i]) / 3
  davg[i,5] <- (d11$Extent_methylation[i] + d12$Extent_methylation[i] + d13$Extent_methylation[i]) / 3
}
plot5 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)
d4b <- data.frame(V1 = length(d4$Start_nt), V2 = length(d4$End_nt))
for(i in 1:length(d4$Start_nt)) {
  d4b[i,1] <- d4$End_nt[i] / 1000000
  if(i == 1) {
    d4b[i,2] <- d4$GC_skew[i]
  }
  else {
    d4b[i,2] <- d4b$V2[i-1] + d4$GC_skew[i]
  }
}
plot6 <- ggplot(d4b, aes(V1, V2)) +
  geom_point(x = d4b$V1, y = d4b$V2, size=0.1, color="black") +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  scale_y_continuous(limits = c(-2,3),
                     labels=scaleFUN,
                     "Cumulative GC skew"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 0.4)
plot9 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="grey") +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 5)

combinedPlot <- ggarrange(plot1, plot3, plot5, plot7, plot8, plot9, plot2, plot4, plot6, 
                          ncol = 3, nrow = 3,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                          font.label = list(size = 10, color = "black"))

svg(filename = "../Figures/Rm2011_skew.svg")
combinedPlot
dev.off()
