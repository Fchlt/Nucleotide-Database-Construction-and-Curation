**Nucleotide-Database-Construction-and-Curation (NDCC.R)**

R function to extract genes from annotated genomes, separate between full and partial ORF and generate a reference taxonomy database.
This function works using genomic files downloaded from the NCBI datasets website (https://www.ncbi.nlm.nih.gov/datasets/genome/) and the original annotation for with genomic files re-annotated using prokka.

**How it works:**

The function contains three main loops:

**_Loop 1: Extract gene sequences from genomes:_**

This first loop extracts genes from your genomic files based on the user-defined "Target" argument (using Target=; e.g. Target="rpoB"). If the search mode is set to strict (strict=TRUE), the search term will have an added "]" at the end. This allows closely related matched to be excluded (e.g. rpoB'). Optionally, if the initial search has not returned a hit (use Use.Alternative.Queries=TRUE), a list of alternative search terms stored as a separate .txt file (no header) where each line is a single search term (see the rpoB_search_terms.txt file for examples. The loop will then search each genomic file stored in the user defined PATH directory, assuming that the PATH directory contains genomic files, with each genome being a separate folder containing a single 'cds_from_genomic.fna' file (when using the original NCBI datasets annotations) or a single .ffn file (when using Prokka re-annotated files). For the loop to work, a genome summary file (.tsv) should also be present in the PATH directory (downloaded as part of the package download on NCBI datasets).

The loop should output a "[Target].fasta" file as well as a "[Target] Summary search.csv" file which recaps which genome returned hits.

**_Loop 2. Full ORF check_**

The next step is to make sure only full length genes (Open Reading Frames; ORF) are retained. To do so, the ORFchecker.R will look at each sequence provided in the "[Target].fasta" file and 1) check for the presence of a START codon at the beginning of the sequence (the list of accepted STRAT codons can be manually changed by setting the START.list argument) 2) check the presence of a STOP codon at the end of the sequence (by translating to an amino acid sequence and checking for a STOP codon * at the end of the sequence) and 3) check for the absence of a STOP codon within the sequence. Additionally, the sequences can be filtered to retain only those superior to the Length_threshold value (use Length_threshold=0 if no length filtering is necessary).

The loop should output a "[Target]_NCBI_datasets_full_ORF.fasta" file (and "[Target]_NCBI_datasets_partial_ORF.fasta" / "[Target]_NCBI_datasets_NON_ORF.fasta" if partial or non-ORF genes are found)

**_Loop 3. Find species name_**

This quick intermediate step is used to extract the Genus and species names of each nucleotide sequence from their header.

The loop should output a "[Target]_NCBI_datasets_full_ORF_species_list.csv" file containing the Genus and Species names for each sequence in the database

**_Loop 4. Generate full taxonomy_**

The last step is to generate a full taxonomy file for each sequence. This is achieved by searching the Species names generated at the previous step against the NCBI Taxonomy database using the "classification" function from the taxize package. The loop will first search for the species name provided for each sequence. If no hit is found, the loop will try modifying the species name (e.g. changing "Nitrosomonas sp." to "uncultured Nitrosomonas"). If still no hit are found, the loop will use the genus name instead.

The loop should output three files: a "[Target]_NCBI_datasets_full_ORF_Taxonomy.csv" file and two taxonomy text files, ("[Target]_NCBI_datasets_full_ORF_Taxonomy_NBC.txt and [Target]_NCBI_datasets_full_ORF_Taxonomy_BLCA.txt") one compatible with NBC classifier and one for BLCA.


_Summary of the arguments:_

PATH (character): the path to the directory containing the genome folders

annotations: whether the original NCBI_datasets files are used (annotation='ncbi') or the Prokka re-annotated files are used (annotation='prokka') [DEFAULT='ncbi']

Target (character): the name of the target gene 

Length_threshold (numeric): the minimum length for the nucleotide sequences [DEFAULT=0]

Use.Alternative.Queries (Logical): whether alternative gene names should be used to query the genomic files [DEFAULT=FALSE]

Queries (file): text file containing the alternative gene names to be used to query the genomic files. (use one name per line) [DEFAULT=NULL]

strict (logical): Whether an additional ']' character s added at the end of the 'Target' name to avoid hits with closely related genes [DEFAULT=FALSE]

START.list (character) any combination of the following: 'ATG','GTG','TTG','CTG','ATT','ATC'
