require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# GANTC - chromosome
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sativa_distal_1/Rm2011_sativa_distal_motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sativa_proximal_1/Rm2011_sativa_proximal_motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sativa_whole_1/10734879-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_sativa_whole_2/10734885-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_sativa_whole_3/10734888-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_GANTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
}
plot1 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot2 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot3 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot4 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot5 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V6, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)

# RCGCCTC - chromosome
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sativa_distal_1/Rm2011_sativa_distal_motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sativa_proximal_1/Rm2011_sativa_proximal_motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sativa_whole_1/10734879-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_sativa_whole_2/10734885-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_sativa_whole_3/10734888-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_RCGCCTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_RCGCCTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_RCGCCTC_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
}
plot6 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot7 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot8 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot9 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot10 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V6, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)

# CGCANNNNNGTG - chromosome
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sativa_distal_1/Rm2011_sativa_distal_motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sativa_proximal_1/Rm2011_sativa_proximal_motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sativa_whole_1/10734879-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/Rm2011_sativa_whole_2/10734885-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/Rm2011_sativa_whole_3/10734888-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/Rm2011_stationary_1/analysis-2011_rep1-180-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/Rm2011_stationary_2/analysis-2011_rep2-179-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/Rm2011_stationary_3/analysis-2011_rep3-178-motifs_CGCANNNNNGTG_CP004140_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
}
plot11 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V2, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot12 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V3, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot13 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V4, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot14 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V5, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)
plot15 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V6, size=0.1, color="grey") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)

# Plot
combinedPlot <- ggarrange(plot1, plot11, plot6, plot5, plot10, plot15, plot2, plot12, plot7, plot3, plot13, plot8, plot4, plot14, plot9,
                          ncol = 3, nrow = 5,
                          labels = c("A", "F", "K", "B", "G", "L", "C", "H", "M", "D", "I", "N", "E", "J", "O"),
                          font.label = list(size = 10, color = "black"))
svg(filename = "../Figures/Rm2011_chromosome_nodule_skew.svg")
combinedPlot
dev.off()

