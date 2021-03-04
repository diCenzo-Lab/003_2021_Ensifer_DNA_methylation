# Set up directories
mkdir Proteomes
mkdir ProteomesHMM
mkdir HMMsearch
mkdir HMMsearchParsed
mkdir HMMsearchHits
mkdir HMMscan
mkdir HMMscanParsed
mkdir HMMscanTop
mkdir OutputFiles
mkdir HMMscanProteins
mkdir intermediaryFiles
mkdir hmmDatabaseFiles
mkdir HMMscanTopLists

# Prepare protein files
cp ../genome_reannotation/intermediaryFiles/genomeList.txt intermediaryFiles # Get the list of genomes
mv hmm_list.txt intermediaryFiles # Move this file
cp ../genome_reannotation/Output/*/*.faa Proteomes # Get the proteomes
perl Scripts/switchNames.pl # Switch the names of the proteins in the faa files
cat ProteomesHMM/*.faa > intermediaryFiles/combined_proteomes_HMM.faa # Combine the faa files into one file
perl Scripts/modifyFasta.pl intermediaryFiles/combined_proteomes_HMM.faa > intermediaryFiles/combined_proteomes_HMM_modified.faa # Modify the fasta file for easy extraction

# Download and prepare HMM libraries
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz # get the Pfam HMM files
wget ftp://ftp.jcvi.org/pub/data/TIGRFAMs//TIGRFAMs_15.0_HMM.LIB.gz # get the TIGRFAM HMM files
gunzip Pfam-A.hmm.gz # unzip the Pfam files
gunzip TIGRFAMs_15.0_HMM.LIB.gz # unzip the TIGRFAM files
mv Pfam-A.hmm hmmDatabaseFiles/Pfam-A.hmm # move the Pfam files
mv TIGRFAMs_15.0_HMM.LIB hmmDatabaseFiles/TIGRFAMs_15.0_HMM.LIB # move the TIGRFAM files
hmmconvert hmmDatabaseFiles/Pfam-A.hmm > hmmDatabaseFiles/Pfam-A_converted.hmm # convert the database to the necessary format
hmmconvert hmmDatabaseFiles/TIGRFAMs_15.0_HMM.LIB > hmmDatabaseFiles/TIGRFAM_converted.hmm # convert the database to the necessary format
cat hmmDatabaseFiles/Pfam-A_converted.hmm hmmDatabaseFiles/TIGRFAM_converted.hmm > hmmDatabaseFiles/converted_combined.hmm # combined all hidden Markov models into a single file
hmmpress hmmDatabaseFiles/converted_combined.hmm # prepare files for hmmscan searches

# Perform the HMMsearch screens
perl Scripts/performHMMsearch.pl # A short script to repeat for all HMM files, the build, hmmsearch, parsing, and hit extraction

# Perform the HMM scan screens
cp /datadisk1/georged/Other/HMM_Database_Files/* hmmDatabaseFiles/ # get the HMM database
perl Scripts/performHMMscan.pl # a short script to repeat for all the HMM search output files, to perform hmmscan, parse, and hit extraction

# Determine strains with each protein
perl Scripts/determineProteinPresence.pl > OutputFiles/proteinPresenceAbsence.txt # determine which of the six proteins are in each of the strains
