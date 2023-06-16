# Genome Base Call Summary programm

This program takes in an aligned multi-FASTA file and annotates the number of calls at each position in the genome.
It then categorises these into 'Nt', 'Gaps', 'Mixed', and 'Ambiguous'. It calculates the percentage of samples with a given group call each position and generates a summary area plot for those categories. A dataframe is saved as a csv.
