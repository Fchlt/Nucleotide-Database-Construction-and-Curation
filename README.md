# Nucleotide-Database-Construction-and-Curation
R scripts to extract genes from genomes and generate a reference taxonomy database

The first R script (GetGenefromGenome.R) is used to extract gene sequences (nucleotide) from a genome coding sequence genomic file (cds_from_genomic.fna) based on user-defined search terms.

_How it works_: First, download genomic files from the NCBI datasets website (https://www.ncbi.nlm.nih.gov/datasets/genome/). Then, in the GetGenefromGenome.R, indicate the name of the gene you want to extract from your genomic files (using Target=; e.g. Target="rpoB") and, optionally (use Use.Alternative.Queries=TRUE), a list of alternative search terms stored as a separate .txt file (no header) where each line is a single search term (see the rpoB_search_terms.txt file for examples. 
The loop will then search each genomic file stored in the user defined PATH directory (assuming that the PATH directory contains genomic files, with each genome being a separate folder containing a single 'cds_from_genomic.fna' file). 
For the loop to work, a data_summary.tsv should also be present in the PATH directory (downloaded as part of the package download on NCBI datasets). 

