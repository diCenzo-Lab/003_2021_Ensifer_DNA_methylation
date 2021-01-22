require(ggplot2)
require(ggpubr)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
black.10.text <- element_text(color = "black", size = 8)
black.12.text <- element_text(color = "black", size = 8)
scaleFUN <- function(x) sprintf("%.1f", x)

# 2011 chromosome GANTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_GANTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
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

# 2011 pSymB GANTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_GANTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot2 <- ggplot(davg, aes(V1, V2)) +
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

# 2011 pSymA GANTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_GANTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
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

# 2011 chromosome CGCANNNNNGTG
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_CGCANNNNNGTG_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot4 <- ggplot(davg, aes(V1, V2)) +
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

# 2011 pSymB CGCANNNNNGTG
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_CGCANNNNNGTG_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
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

# 2011 pSymA CGCANNNNNGTG
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_CGCANNNNNGTG_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot6 <- ggplot(davg, aes(V1, V2)) +
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

# 2011 chromosome RCGCCTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_RCGCCTC_CP004140_1_Sinorhizobium_meliloti_2011_complete_genome_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot7 <- ggplot(davg, aes(V1, V2)) +
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

# 2011 pSymB RCGCCTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_RCGCCTC_CP004139_1_Sinorhizobium_meliloti_2011_plasmid_pSymB_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot8 <- ggplot(davg, aes(V1, V2)) +
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

# 2011 pSymA RCGCCTC
d1 <- read.table("../Output2/Rm2011_succinate_1/10738447-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d2 <- read.table("../Output2/Rm2011_succinate_2/10738438-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d3 <- read.table("../Output2/Rm2011_succinate_3/10738429-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d4 <- read.table("../Output2/Rm2011_sucrose_1/10734983-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d5 <- read.table("../Output2/Rm2011_sucrose_2/10734986-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
d6 <- read.table("../Output2/Rm2011_sucrose_3/10734974-motifs_RCGCCTC_CP004138_1_Sinorhizobium_meliloti_2011_plasmid_pSymA_complete_sequence_10kb.txt", header = TRUE, sep = "\t")
davg <- data.frame(V1 = length(d1$Start_nt), V2 = length(d1$Start_nt), V3 = length(d1$Start_nt))
for(i in 1:length(d1$Start_nt)) {
  davg[i,1] <- d1$End_nt[i] / 1000000
  davg[i,2] <- (d1$Extent_methylation[i] + d2$Extent_methylation[i] + d3$Extent_methylation[i]) / 3
  davg[i,3] <- (d4$Extent_methylation[i] + d5$Extent_methylation[i] + d6$Extent_methylation[i]) / 3
}
plot9 <- ggplot(davg, aes(V1, V2)) +
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

# Make combined plot
combinedPlot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, 
                          ncol = 3, nrow = 3,
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                          font.label = list(size = 10, color = "black"))

svg(filename = "../Figures/Rm2011_carbon_skew.svg")
combinedPlot
dev.off()
