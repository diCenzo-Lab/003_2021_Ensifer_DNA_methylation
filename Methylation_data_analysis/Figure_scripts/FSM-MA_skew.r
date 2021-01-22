require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# 2011 chromosome succinate
d1 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_GANTC_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_GANTC_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_GANTC_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_TCGANNNNNNNNTCGA_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_TCGANNNNNNNNTCGA_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_TCGANNNNNNNNTCGA_CP019584_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_GANTC_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_GANTC_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_GANTC_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_TCGANNNNNNNNTCGA_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_TCGANNNNNNNNTCGA_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_TCGANNNNNNNNTCGA_CP019584_1_10kb.txt", header = TRUE, sep = "\t")
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
d1 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_GANTC_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_GANTC_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_GANTC_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_TCGANNNNNNNNTCGA_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_TCGANNNNNNNNTCGA_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_TCGANNNNNNNNTCGA_CP019586_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_GANTC_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_GANTC_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_GANTC_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_TCGANNNNNNNNTCGA_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_TCGANNNNNNNNTCGA_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_TCGANNNNNNNNTCGA_CP019586_1_10kb.txt", header = TRUE, sep = "\t")
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
d1 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_GCskew.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
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

svg(filename = "../Figures/FSM-MA_skew.svg")
combinedPlot
dev.off()
