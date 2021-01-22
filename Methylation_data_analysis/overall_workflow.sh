# Make output directories
mkdir Output/
mkdir Output2/

# Rm2011 succinate
mkdir Output2/Rm2011_succinate_1/
mkdir Output2/Rm2011_succinate_2/
mkdir Output2/Rm2011_succinate_3/
perl Scripts/modify_2011.pl Data/10738447-motifs.gff > Data/10738447-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10738438-motifs.gff > Data/10738438-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10738429-motifs.gff > Data/10738429-motifs_modified.gff
sh initAnalysis.sh Data/10738437-motifs.csv Data/10738447-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_succinate_1/
sh initAnalysis.sh Data/10738437-motifs.csv Data/10738438-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_succinate_2/
sh initAnalysis.sh Data/10738437-motifs.csv Data/10738429-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_succinate_3/

# Rm2011 sucrose
mkdir Output2/Rm2011_sucrose_1/
mkdir Output2/Rm2011_sucrose_2/
mkdir Output2/Rm2011_sucrose_3/
perl Scripts/modify_2011.pl Data/10734983-motifs.gff > Data/10734983-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10734986-motifs.gff > Data/10734986-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10734974-motifs.gff > Data/10734974-motifs_modified.gff
sh initAnalysis.sh Data/10734973-motifs.csv Data/10734983-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sucrose_1/
sh initAnalysis.sh Data/10734973-motifs.csv Data/10734986-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sucrose_2/
sh initAnalysis.sh Data/10734973-motifs.csv Data/10734974-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sucrose_3/

# Rm2011 sativa whole nodules
mkdir Output2/Rm2011_sativa_whole_1/
mkdir Output2/Rm2011_sativa_whole_2/
mkdir Output2/Rm2011_sativa_whole_3/
perl Scripts/modify_2011.pl Data/10734879-motifs.gff > Data/10734879-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10734885-motifs.gff > Data/10734885-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10734888-motifs.gff > Data/10734888-motifs_modified.gff
sh initAnalysis.sh Data/10734878-motifs.csv Data/10734879-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sativa_whole_1/
sh initAnalysis.sh Data/10734878-motifs.csv Data/10734885-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sativa_whole_2/
sh initAnalysis.sh Data/10734878-motifs.csv Data/10734888-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sativa_whole_3/

# Rm2011 sativa distal nodules
mkdir Output2/Rm2011_sativa_distal_1/
perl Scripts/modify_2011.pl Data/Rm2011_sativa_distal_motifs.gff > Data/Rm2011_sativa_distal_motifs_modified.gff
sh initAnalysis.sh Data/Rm2011_sativa_distal_motifs.csv Data/Rm2011_sativa_distal_motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sativa_distal_1/

# Rm2011 sativa proximal nodules
mkdir Output2/Rm2011_sativa_proximal_1/
perl Scripts/modify_2011.pl Data/Rm2011_sativa_proximal_motifs.gff > Data/Rm2011_sativa_proximal_motifs_modified.gff
sh initAnalysis.sh Data/Rm2011_sativa_proximal_motifs.csv Data/Rm2011_sativa_proximal_motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_sativa_proximal_1/

# Rm2011 ∆pSymAB sucrose
mkdir Output2/pSymAB_sucrose_1/
mkdir Output2/pSymAB_sucrose_2/
mkdir Output2/pSymAB_sucrose_3/
perl Scripts/modify_2011.pl Data/10635469-motifs.gff > Data/10635469-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10635470-motifs.gff > Data/10635470-motifs_modified.gff
perl Scripts/modify_2011.pl Data/10635502-motifs.gff > Data/10635502-motifs_modified.gff
sh initAnalysis.sh Data/10635465-motifs.csv Data/10635469-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/pSymAB_sucrose_1/
sh initAnalysis.sh Data/10635465-motifs.csv Data/10635470-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/pSymAB_sucrose_2/
sh initAnalysis.sh Data/10635465-motifs.csv Data/10635502-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/pSymAB_sucrose_3/

# OV14 succinate
mkdir Output2/OV14_succinate_1/
mkdir Output2/OV14_succinate_2/
mkdir Output2/OV14_succinate_3/
perl Scripts/modify_OV14.pl Data/10734925-motifs.gff > Data/10734925-motifs_modified.gff
perl Scripts/modify_OV14.pl Data/10734927-motifs.gff > Data/10734927-motifs_modified.gff
perl Scripts/modify_OV14.pl Data/10740678-motifs.gff > Data/10740678-motifs_modified.gff
sh initAnalysis.sh Data/10734924-motifs.csv Data/10734925-motifs_modified.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.fna 1
mv Output/* Output2/OV14_succinate_1/
sh initAnalysis.sh Data/10734924-motifs.csv Data/10734927-motifs_modified.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.fna 1
mv Output/* Output2/OV14_succinate_2/
sh initAnalysis.sh Data/10734924-motifs.csv Data/10740678-motifs_modified.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.gff Genomes/GCA_000583045.1_ASM58304v1_genomic.fna 1
mv Output/* Output2/OV14_succinate_3/

# NGR234 succinate
mkdir Output2/NGR234_succinate_1/
mkdir Output2/NGR234_succinate_2/
mkdir Output2/NGR234_succinate_3/
perl Scripts/modify_NGR234.pl Data/10734882-motifs.gff > Data/10734882-motifs_modified.gff
perl Scripts/modify_NGR234.pl Data/10734884-motifs.gff > Data/10734884-motifs_modified.gff
perl Scripts/modify_NGR234.pl Data/10734920-motifs.gff > Data/10734920-motifs_modified.gff
sh initAnalysis.sh Data/10734874-motifs.csv Data/10734882-motifs_modified.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.fna 1
mv Output/* Output2/NGR234_succinate_1/
sh initAnalysis.sh Data/10734874-motifs.csv Data/10734884-motifs_modified.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.fna 1
mv Output/* Output2/NGR234_succinate_2/
sh initAnalysis.sh Data/10734874-motifs.csv Data/10734920-motifs_modified.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.gff Genomes/GCA_000018545.1_ASM1854v1_genomic.fna 1
mv Output/* Output2/NGR234_succinate_3/

# FSM-MA succinate
mkdir Output2/FSM-MA_succinate_1/
mkdir Output2/FSM-MA_succinate_2/
mkdir Output2/FSM-MA_succinate_3/
perl Scripts/modify_FSMMA.pl Data/10735010-motifs.gff > Data/10735010-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/10735001-motifs.gff > Data/10735001-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/10738420-motifs.gff > Data/10738420-motifs_modified.gff
sh initAnalysis.sh Data/10735009-motifs.csv Data/10735010-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_succinate_1/
sh initAnalysis.sh Data/10735009-motifs.csv Data/10735001-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_succinate_2/
sh initAnalysis.sh Data/10735009-motifs.csv Data/10738420-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_succinate_3/

# FSM-MA sativa whole nodules
mkdir Output2/FSM-MA_sativa_whole_1/
mkdir Output2/FSM-MA_sativa_whole_2/
mkdir Output2/FSM-MA_sativa_whole_3/
perl Scripts/modify_FSMMA.pl Data/10635471-motifs.gff > Data/10635471-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/10635472-motifs.gff > Data/10635472-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/11388042-motifs.gff > Data/11388042-motifs_modified.gff
sh initAnalysis.sh Data/10635467-motifs.csv Data/10635471-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_sativa_whole_1/
sh initAnalysis.sh Data/10635467-motifs.csv Data/10635472-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_sativa_whole_2/
sh initAnalysis.sh Data/10635467-motifs.csv Data/11388042-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_sativa_whole_3/

# FSM-MA truncatula whole nodules
mkdir Output2/FSM-MA_truncatula_whole_1/
mkdir Output2/FSM-MA_truncatula_whole_2/
mkdir Output2/FSM-MA_truncatula_whole_3/
perl Scripts/modify_FSMMA.pl Data/10538740-motifs.gff > Data/10538740-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/10538749-motifs.gff > Data/10538749-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/10538758-motifs.gff > Data/10538758-motifs_modified.gff
sh initAnalysis.sh Data/10538739-motifs.csv Data/10538740-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_truncatula_whole_1/
sh initAnalysis.sh Data/10538739-motifs.csv Data/10538749-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_truncatula_whole_2/
sh initAnalysis.sh Data/10538739-motifs.csv Data/10538758-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_truncatula_whole_3/

# FSM-MA sativa distal nodules
mkdir Output2/FSM-MA_sativa_distal_1/
perl Scripts/modify_FSMMA.pl Data/FSM-MA_sativa_distal_motifs.gff > Data/FSM-MA_sativa_distal_motifs_modified.gff
sh initAnalysis.sh Data/FSM-MA_sativa_distal_motifs.csv Data/FSM-MA_sativa_distal_motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_sativa_distal_1/

# FSM-MA sativa proximal nodules
mkdir Output2/FSM-MA_sativa_proximal_1/
perl Scripts/modify_FSMMA.pl Data/FSM-MA_sativa_proximal_motifs.gff > Data/FSM-MA_sativa_proximal_motifs_modified.gff
sh initAnalysis.sh Data/FSM-MA_sativa_proximal_motifs.csv Data/FSM-MA_sativa_proximal_motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_sativa_proximal_1/

# Rm2011 stationary
mkdir Output2/Rm2011_stationary_1/
mkdir Output2/Rm2011_stationary_2/
mkdir Output2/Rm2011_stationary_3/
perl Scripts/modify_2011.pl Data/analysis-2011_rep1-180-motifs.gff > Data/analysis-2011_rep1-180-motifs_modified.gff
perl Scripts/modify_2011.pl Data/analysis-2011_rep2-179-motifs.gff > Data/analysis-2011_rep2-179-motifs_modified.gff
perl Scripts/modify_2011.pl Data/analysis-2011_rep3-178-motifs.gff > Data/analysis-2011_rep3-178-motifs_modified.gff
sh initAnalysis.sh Data/analysis-2011_rep1-180-motifs.csv Data/analysis-2011_rep1-180-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_stationary_1/
sh initAnalysis.sh Data/analysis-2011_rep2-179-motifs.csv Data/analysis-2011_rep2-179-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_stationary_2/
sh initAnalysis.sh Data/analysis-2011_rep3-178-motifs.csv Data/analysis-2011_rep3-178-motifs_modified.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.gff Genomes/GCA_000346065.1_ASM34606v1_genomic.fna 1
mv Output/* Output2/Rm2011_stationary_3/

# FSM-MA succinate
mkdir Output2/FSM-MA_stationary_1/
mkdir Output2/FSM-MA_stationary_2/
mkdir Output2/FSM-MA_stationary_3/
perl Scripts/modify_FSMMA.pl Data/analysis-FSM_rep1-184-motifs.gff > Data/analysis-FSM_rep1-184-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/analysis-FSM_rep2-183-motifs.gff > Data/analysis-FSM_rep2-183-motifs_modified.gff
perl Scripts/modify_FSMMA.pl Data/analysis-FSM_rep3-182-motifs.gff > Data/analysis-FSM_rep3-182-motifs_modified.gff
sh initAnalysis.sh Data/analysis-FSM_rep1-184-motifs.csv Data/analysis-FSM_rep1-184-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_stationary_1/
sh initAnalysis.sh Data/analysis-FSM_rep2-183-motifs.csv Data/analysis-FSM_rep2-183-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_stationary_2/
sh initAnalysis.sh Data/analysis-FSM_rep3-182-motifs.csv Data/analysis-FSM_rep3-182-motifs_modified.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.gff Genomes/GCA_002215195.1_ASM221519v1_genomic.fna 1
mv Output/* Output2/FSM-MA_stationary_3/
