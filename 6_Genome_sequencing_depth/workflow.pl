#Download BAM files
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_22/download/_JAMO/5dd347f80a3792142a000d76/10734900-mapped.alignmentset.bam' -b cookies > Rm2011_alfalfa_1.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_23/download/_JAMO/5dd495020a3792142a001b7d/10738440-mapped.alignmentset.bam' -b cookies > Rm2011_succinate_1.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_20/download/_JAMO/5ebd719eddecee1bc444ac9c/mapped.bam' -b cookies > FSM-MA_alfalfa_1.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_21/download/_JAMO/5d8417ea95f4dcd30aeabe0d/10538760-mapped.alignmentset.bam' -b cookies > FSM-MA_truncatula_1.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_24/download/_JAMO/5dd347fd0a3792142a000df6/10735012-mapped.alignmentset.bam' -b cookies > FSM-MA_succinate_1.bam

curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_22/download/_JAMO/5dd347f80a3792142a000d71/10734895-mapped.alignmentset.bam' -b cookies > Rm2011_alfalfa_2.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_23/download/_JAMO/5dd494f80a3792142a001b73/10738431-mapped.alignmentset.bam' -b cookies > Rm2011_succinate_2.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_20/download/_JAMO/5da9fbcbaa74fab996decb7e/10635480-mapped.alignmentset.bam' -b cookies > FSM-MA_alfalfa_2.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_21/download/_JAMO/5d8417dd95f4dcd30aeabe03/10538751-mapped.alignmentset.bam' -b cookies > FSM-MA_truncatula_2.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_24/download/_JAMO/5dd4947f0a3792142a001b69/10738422-mapped.alignmentset.bam' -b cookies > FSM-MA_succinate_2.bam

curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_22/download/_JAMO/5dd347f70a3792142a000d61/10734881-mapped.alignmentset.bam' -b cookies > Rm2011_alfalfa_3.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_23/download/_JAMO/5dd4950b0a3792142a001b87/10738449-mapped.alignmentset.bam' -b cookies > Rm2011_succinate_3.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_20/download/_JAMO/5da9fbcbaa74fab996decb7d/10635479-mapped.alignmentset.bam' -b cookies > FSM-MA_alfalfa_3.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_21/download/_JAMO/5d8417d195f4dcd30aeabdf9/10538742-mapped.alignmentset.bam' -b cookies > FSM-MA_truncatula_3.bam
curl 'https://genome.jgi.doe.gov/portal/ext-api/downloads/get_tape_file?blocking=true&url=/Sinmelmdetection_24/download/_JAMO/5dd347fc0a3792142a000dec/10735003-mapped.alignmentset.bam' -b cookies > FSM-MA_succinate_3.bam

mv *.bam BAM_files/

#Determine coverage
samtools coverage BAM_files/Rm2011_alfalfa_1.bam > Output/Rm2011_alfalfa_1.txt
samtools coverage BAM_files/Rm2011_succinate_1.bam > Output/Rm2011_succinate_1.txt
samtools coverage BAM_files/FSM-MA_alfalfa_1.bam > Output/FSM-MA_alfalfa_1.txt
samtools coverage BAM_files/FSM-MA_truncatula_1.bam > Output/FSM-MA_truncatula_1.txt
samtools coverage BAM_files/FSM-MA_succinate_1.bam > Output/FSM-MA_succinate_1.txt

samtools coverage BAM_files/Rm2011_alfalfa_2.bam > Output/Rm2011_alfalfa_2.txt
samtools coverage BAM_files/Rm2011_succinate_2.bam > Output/Rm2011_succinate_2.txt
samtools coverage BAM_files/FSM-MA_alfalfa_2.bam > Output/FSM-MA_alfalfa_2.txt
samtools coverage BAM_files/FSM-MA_truncatula_2.bam > Output/FSM-MA_truncatula_2.txt
samtools coverage BAM_files/FSM-MA_succinate_2.bam > Output/FSM-MA_succinate_2.txt

samtools coverage BAM_files/Rm2011_alfalfa_3.bam > Output/Rm2011_alfalfa_3.txt
samtools coverage BAM_files/Rm2011_succinate_3.bam > Output/Rm2011_succinate_3.txt
samtools coverage BAM_files/FSM-MA_alfalfa_3.bam > Output/FSM-MA_alfalfa_3.txt
samtools coverage BAM_files/FSM-MA_truncatula_3.bam > Output/FSM-MA_truncatula_3.txt
samtools coverage BAM_files/FSM-MA_succinate_3.bam > Output/FSM-MA_succinate_3.txt

samtools coverage BAM_files/analysis-A17_wt_b-446-mapped.bam > Output/A17.txt
samtools coverage BAM_files/analysis-dnf1_b-445-mapped-002.bam > Output/dnf1.txt
samtools coverage BAM_files/analysis-dnf2_b-444-mapped.bam > Output/dnf2.txt
samtools coverage BAM_files/analysis-dnf4_b-450-mapped-002.bam > Output/dnf4.txt
samtools coverage BAM_files/analysis-dnf5_b-449-mapped-002.bam > Output/dnf5.txt
samtools coverage BAM_files/analysis-dnf7_b-448-mapped-002.bam > Output/dnf7.txt
