require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# GANTC - chromosome
d1 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_sativa_distal_1/FSM-MA_sativa_distal_motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_sativa_proximal_1/FSM-MA_sativa_proximal_motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_sativa_whole_1/10635471-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/FSM-MA_sativa_whole_2/10635472-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/FSM-MA_sativa_whole_3/11388042-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/FSM-MA_truncatula_whole_1/10538740-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/FSM-MA_truncatula_whole_2/10538749-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/FSM-MA_truncatula_whole_3/10538758-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d14 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt), V7 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
  davg[i,7] <- (d12$Extent_methylation[i] + d13$Extent_methylation[i] + d14$Extent_methylation[i]) / 3
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
plot11 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V7, size=0.1, color="black") +
  scale_y_continuous(limits = c(0.6,1),
                     "Extent of methylation"
  ) +
  scale_x_continuous(limits = c(0,4),
                     "Nucleotide position (Mb)"
  ) +
  theme_classic() +
  theme(axis.text = black.10.text, axis.title = black.12.text) +
  coord_equal(ratio = 6)

# TCGANNNNNNNNTCGA - chromosome
d1 <- read.table("../Output2/FSM-MA_succinate_1/10735010-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_succinate_2/10735001-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_succinate_3/10738420-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_sativa_distal_1/FSM-MA_sativa_distal_motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_sativa_proximal_1/FSM-MA_sativa_proximal_motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_sativa_whole_1/10635471-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d7 <- read.table("../Output2/FSM-MA_sativa_whole_2/10635472-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d8 <- read.table("../Output2/FSM-MA_sativa_whole_3/11388042-motifs_GANTC_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d9 <- read.table("../Output2/FSM-MA_truncatula_whole_1/10538740-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d10 <- read.table("../Output2/FSM-MA_truncatula_whole_2/10538749-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d11 <- read.table("../Output2/FSM-MA_truncatula_whole_3/10538758-motifs_TCGANNNNNNNNTCGA_CP019585_1_Sinorhizobium_meliloti_strain_CCMM_B554_FSM_MA_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d12 <- read.table("../Output2/FSM-MA_stationary_1/analysis-FSM_rep1-184-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d13 <- read.table("../Output2/FSM-MA_stationary_2/analysis-FSM_rep2-183-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d14 <- read.table("../Output2/FSM-MA_stationary_3/analysis-FSM_rep3-182-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt), V7 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- d4$Extent_methylation[i]
  davg[i,4] <- d5$Extent_methylation[i]
  davg[i,5] <- (d6$Extent_methylation[i] + d7$Extent_methylation[i] + d8$Extent_methylation[i]) / 3
  davg[i,6] <- (d9$Extent_methylation[i] + d10$Extent_methylation[i] + d11$Extent_methylation[i]) / 3
  davg[i,7] <- (d12$Extent_methylation[i] + d13$Extent_methylation[i] + d14$Extent_methylation[i]) / 3
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
plot12 <- ggplot(davg, aes(V1, V2)) +
  geom_point(x = davg$V1, y = davg$V7, size=0.1, color="grey") +
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
combinedPlot <- ggarrange(plot1, plot6, plot11, plot12, plot2, plot7, plot3, plot8, plot4, plot9, plot5, plot10,
                          ncol = 2, nrow = 6,
                          labels = c("A", "G", "B", "H", "C", "I", "D", "J", "E", "K", "F", "L"),
                          font.label = list(size = 10, color = "black"))
svg(filename = "../Figures/FSM-MA_pSymA_nodule_skew.svg")
combinedPlot
dev.off()


