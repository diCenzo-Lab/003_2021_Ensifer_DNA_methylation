require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# GANTC - chromosome
d1 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-A17_wt_b-534-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf1_b-533-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf2_b-532-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf4_b-531-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf5_b-530-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf7_b-529-motifs_GANTC_CP019585_1_10kb.txt", header = TRUE, sep = "\t")

davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt), V7 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- d1$Extent_methylation[i]
  davg[i,3] <- d2$Extent_methylation[i]
  davg[i,4] <- d3$Extent_methylation[i]
  davg[i,5] <- d4$Extent_methylation[i]
  davg[i,6] <- d5$Extent_methylation[i]
  davg[i,7] <- d6$Extent_methylation[i]
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
plot6 <- ggplot(davg, aes(V1, V2)) +
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


# CGANNNNNNNNTCGA - chromosome
d1 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-A17_wt_b-534-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf1_b-533-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf2_b-532-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf4_b-531-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf5_b-530-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/FSM-MA_dnf_nodules_1/analysis-dnf7_b-529-motifs_TCGANNNNNNNNTCGA_CP019585_1_10kb.txt", header = TRUE, sep = "\t")

davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt), V4 = length(d1$Start_nt), V5 = length(d1$Start_nt), V6 = length(d1$Start_nt), V7 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- d1$Extent_methylation[i]
  davg[i,3] <- d2$Extent_methylation[i]
  davg[i,4] <- d3$Extent_methylation[i]
  davg[i,5] <- d4$Extent_methylation[i]
  davg[i,6] <- d5$Extent_methylation[i]
  davg[i,7] <- d6$Extent_methylation[i]
}
plot7 <- ggplot(davg, aes(V1, V2)) +
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
plot8 <- ggplot(davg, aes(V1, V2)) +
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
plot9 <- ggplot(davg, aes(V1, V2)) +
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
plot10 <- ggplot(davg, aes(V1, V2)) +
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
plot11 <- ggplot(davg, aes(V1, V2)) +
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
plot12 <- ggplot(davg, aes(V1, V2)) +
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

# Plot
combinedPlot <- ggarrange(plot2, plot8, plot5, plot11, plot3, plot9, plot6, plot12, plot4, plot10, plot1, plot7,
                          ncol = 2, nrow = 6,
                          labels = c("A", "G", "B", "H", "C", "I", "D", "J", "E", "K", "F", "L"),
                          font.label = list(size = 10, color = "black"))
svg(filename = "../Figures/FSM-MA_dnf_pSymA_skew.svg")
combinedPlot
dev.off()


