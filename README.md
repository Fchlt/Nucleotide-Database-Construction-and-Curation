# Nucleotide-Database-Construction-and-Curation
R scripts to extract genes from genomes and generate a reference taxonomy database

The first R script (GetGenefromGenome.R) is used to extract gene sequences (nucleotide) from a genome coding sequence genomic file (cds_from_genomic.fna) based on user-defined search terms.


**How it works:** 

**_1. Extract gene sequences from genomes:_**
  
First, download genomic files from the NCBI datasets website (https://www.ncbi.nlm.nih.gov/datasets/genome/). Then, in the GetGenefromGenome.R, indicate the name of the gene you want to extract from your genomic files (using Target=; e.g. Target="rpoB"). If the search mode is set to strict (strict=TRUE), the search term will have an added ']' at the end. This allows closely related matched to be excluded (e.g. rpoB'). Optionally, if the initial search has not returned a hit (use Use.Alternative.Queries=TRUE), a list of alternative search terms stored as a separate .txt file (no header) where each line is a single search term (see the rpoB_search_terms.txt file for examples. 
The loop will then search each genomic file stored in the user defined PATH directory (assuming that the PATH directory contains genomic files, with each genome being a separate folder containing a single 'cds_from_genomic.fna' file). 
For the loop to work, a data_summary.tsv should also be present in the PATH directory (downloaded as part of the package download on NCBI datasets). 

The loop should output a '.fasta' file

_Summary of the arguments:_

PATH (character): the path to the directory containing the genome folders

genome.info (data.frame): data.frame containing genome info, downloaded as part of the package from NCBI datasets

Target (character): The abbreviation name of the gene of interest (e.g. 'amoA'; 'rpoB')

strict (logical): Whether an additional ']' character s added at the end of the 'Target' name to avoid hits with closely related genes

Use.Alternative.Queries (Logical): whether alternative gene names should be used to query the genomic files

Queries (file): text file containing the alternative gene names to be used to query the genomic files. (use one name per line)

**_2. Full ORF check_**

The next step is to make sure only full length genes (Open Reading Frames; ORF) are retained. To do so, the ORFchecker.R will look at each  sequence provided in the 'sequences.fasta' file and 1) check for the presence of a START codon at the beginning of the sequence (the list of accepted STRAT codons can be manually changed by setting the START.list argument) 2) check the presence of a STOP codon at the end of the sequence (by translating to an amino acid sequence and checking for a STOP codon * at the end of the sequence) and 3) check for the absence of a STOP codon within the sequence. Additionally, the sequences can be filtered to retain only those superior to the Length_threshold value (use Length_threshold=0 if no length filtering is necessary).

The loop should output a "full_ORF.fasta" file

_Summary of the arguments:_

Target (character): the name of the target gene

sequences (DNAStringSet): DNA sequences in fasta format (uploaded via the Biostrings package).

Length_threshold (numeric): the minimum length for the nucleotide sequences

START.list (character) any combination of the following: 'ATG','GTG','TTG','CTG','ATT','ATC'

**_3. Find species name_**

This quick intermediate step is used to extract the Genus and species names of each nucleotide sequence from their header. The loop can currently extract the species names from fasta files downloaded from NCBI datasets, NCBI (nucleotide), Microbial Genome Data Base (MGDB) or European Nucleotide Archive (ENA). If using the MGDB, the genomelist.txt file needs to be used. 

The loop should output a 'species_list.csv' file containing the Genus and Species names for each sequence in the database

_Summary of the arguments:_

Target (character): the name of the target gene

sequences (DNAStringSet or sequences): DNA sequences in fasta format (uploaded via the Biostrings or seqinr package, respectively). 

genome.list (file): file containing the species abbreviations and corresponding full names

**_4. Generate full taxonomy_**

The last step is to generate a full taxonomy file for each sequence. This is achieved by searching the Species names generated at the previous step against the NCBI Taxonomy database using the "classification" function from the taxize package. The loop will first search for the species name provided for each sequence. If no hit is found, the loop will try modifying the species name (e.g. changing "Nitrosomonas sp." to "uncultured Nitrosomonas"). If still no hit are found, the loop will use the genus name instead. 

The loop should output three files: a taxonomy file in cdv format and two taxonomy text files, one compatible with NBC classifier and one for BLCA. 





