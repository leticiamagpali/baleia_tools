"""
This script calls each step of the get_sequences pipeline.
Each step is clearly marked and can be commented out with triple-quotations if it is not
needed.

Authors: Letícia Magpali, Rafael Copstein, Yuri Kulish
"""

import os
from typing import List, Tuple
from Bio import SeqIO
from get_sequences import extract_seqs
from get_sequences import process_isoforms

# Modify this variable to point to the root directory where you will run the analyses
# or leave as it is and the analysis will run on your working directory
# (i.e. the same directory where this file is)
WORKING_DIR = "/Users/leticiamagpali/Google Drive/My Drive/phd_leticia/research_project/scripts/sandbox/get-seqs-test-1"
# os.path.dirname(os.path.realpath(__file__))

# Modify this path to point to the folder containing the genomic data you wish to analyse
GENOMICS_DIR = f"{WORKING_DIR}/test-sequences"

# Modify the paths below to point to your results folders
# you can also leave it as it is and they will be created for you inside the working directory
# in this case the results will be saved on the same directory where the code is running

RESULTS_DIR = f"{WORKING_DIR}/results"
os.makedirs(RESULTS_DIR, exist_ok=True)

RESULTS_ISOFORMS = f"{WORKING_DIR}/results/largest-isoforms"
os.makedirs(RESULTS_ISOFORMS, exist_ok=True)

RESULTS_EXTRACTED = f"{WORKING_DIR}/results/extracted-genes"
os.makedirs(RESULTS_EXTRACTED, exist_ok=True)

RESULTS_FILTERED = f"{WORKING_DIR}/results/filtered-lowq"
os.makedirs(RESULTS_FILTERED, exist_ok=True)

RESULTS_NAME_CHANGED_GENES = f"{WORKING_DIR}/results/genes-name-changed"
os.makedirs(RESULTS_NAME_CHANGED_GENES, exist_ok=True)

RESULTS_MULTISPECIES_GENES = f"{WORKING_DIR}/results/genes-multispecies"
os.makedirs(RESULTS_MULTISPECIES_GENES, exist_ok=True)

RESULTS_PROTEINS = f"{WORKING_DIR}/results/proteins"
os.makedirs(RESULTS_PROTEINS, exist_ok=True)

RESULTS_MULTISPECIES_PROT = f"{WORKING_DIR}/results/proteins-multispecies"
os.makedirs(RESULTS_MULTISPECIES_PROT, exist_ok=True)

# modify this path to point to file containing the list of genes you want to extract
LIST_PATH = f"{WORKING_DIR}/cognition_genes_test.txt"


##########################
# Process nucleotide file
##########################

def process_nucleotide_file() -> Tuple[str, List[str]]:
    '''Makes a list of file paths for the DNA files'''

    # creates empty list to store paths to each cds_from_genomic.fna file
    raw_data_file_nuc: List[str] = []

    # finds all cds_from_genomic.fna.fna files in the end of each path inside "genomics_directory"
    # writes a path to those files by joining root and the file name
    # adds path to empty list
    for current_dir, _, files in os.walk(GENOMICS_DIR):
        for file in files:
            if file == "cds_from_genomic.fna":
                raw_data_file_nuc.append(os.path.join(current_dir, file))

    # Uncomment the part below and modify to run on a single file instead of whole dataset
    # raw_data_file =
    # [genomics_directory + "Delphinapterus_leucas/GCF_002288925.2/cds_from_genomic.fna"]

    return LIST_PATH, raw_data_file_nuc


#############################
# Extracting longest isoform
#############################

def extract_longest_isoform(raw_data_file_nuc: List[str]):
    '''Runs the module to process isoforms (and extract the longest ones)
    on all cds_from_genomic.fna files '''

    for filepath in raw_data_file_nuc:
        dataset_i = process_isoforms.keep_longest_isoform(filepath)

        # Writes the resulting dataset to a FASTA file named after the organism
        species_name = filepath.split('/')[-3]
        with open(f"{RESULTS_ISOFORMS}/{species_name}_isoforms.fna", "w") as output_file:
            SeqIO.write(dataset_i, output_file, "fasta")


##############################
# Extracting genes of interest
##############################

def extract_genes_of_interest(list_path: str):
    '''Runs the extract_seqs.extract_genes() function on the results from the step above
    (extracting longest isoform)'''

    # locates directory storing the results from previous step (remove isoforms)
    # Then, for each file on this directory:
    # runs the function extract_genes
    # and stores the results on the variable dataset_IE
    for filename in os.listdir(RESULTS_ISOFORMS):
        if filename.startswith("."):
            continue

        path_result = os.path.join(RESULTS_ISOFORMS, filename)
        with open(path_result, "r", encoding="ascii", errors="surrogateescape") as file:
            dataset_ie = extract_seqs.extract_genes(file, list_path)

        # Writes the results stored on dataset_IE on a fasta file
        # Names file accordingly to each species (takes name of species from filepath)
        species_name = filename.replace("isoforms", "extracted")
        with open(f"{RESULTS_EXTRACTED}/{species_name}", "w") as output_file:
            SeqIO.write(dataset_ie, output_file, "fasta")


################################
# Removing low quality sequences
################################

def remove_low_quality_sequences():
    '''Runs the extract_seqs.fiter_low_quality() function on the results from the step above
    (extracting genes of interest)'''

    # locates directory storing the results from previous step (extracting genes of interest)
    # Then, for each file on this directory:
    # runs the function filter_low_quality
    # and stores the results on the variable dataset_IEF

    for filename in os.listdir(RESULTS_EXTRACTED):
        if filename.startswith("."):
            continue

        path_result = os.path.join(RESULTS_EXTRACTED, filename)
        with open(path_result, "r", encoding="ascii", errors="surrogateescape") as file:
            dataset_ief = extract_seqs.fiter_low_quality(file)

        species_name = filename.replace("extracted", "filtered")
        with open(f"{RESULTS_FILTERED}/{species_name}", "w") as output_file:
            SeqIO.write(dataset_ief, output_file, "fasta")


#############################################
# Edit sequence names (description) on files
#############################################

def edit_sequence_names():
    '''Runs the extract_seqs.change_seq_description() function on the results from the step above
(removing low quality sequences)'''
    for filename in os.listdir(RESULTS_FILTERED):
        if filename.startswith("."):
            continue

        path_result = os.path.join(RESULTS_FILTERED, filename)
        species_name = '_'.join(filename.split('_')[:2])
        extract_seqs.change_seq_description(
            path_result, species_name, RESULTS_NAME_CHANGED_GENES)


################################################
# Make files for one gene and multiple species
################################################

def make_multispecies_gene_file():
    '''Runs the extract_seqs.multispecies_gene() function on the results from the step above
    (edit sequence names)'''

    for filename in os.listdir(RESULTS_NAME_CHANGED_GENES):
        if filename.startswith("."):
            continue

        path_result = os.path.join(RESULTS_NAME_CHANGED_GENES, filename)
        extract_seqs.multispecies_gene(
            path_result, LIST_PATH, RESULTS_MULTISPECIES_GENES)


########################
# Process protein files
########################

def process_protein_file() -> List[str]:
    '''Makes a list of file paths for the protein files'''

    # all paths to protein.faa files inside the genomics_directory
    raw_data_file_prot: List[str] = []

    # finds all files in the end of each path inside "genomics_directory"
    # writes a path to those files by joining root and the file name
    # adds path to empty list
    for current_dir, _, files in os.walk(GENOMICS_DIR):
        for file in files:
            if file == "protein.faa":
                raw_data_file_prot.append(os.path.join(current_dir, file))

    return raw_data_file_prot


#################################
# Extracting proteins of interest
#################################

def extract_proteins_of_interest(raw_data_file_prot: List[str]):
    '''Extract proteins from .faa files based on the genes that have been extracted
    by the pipeline above'''

    # creates empty list to store protein ID numbers
    proteins_list: List[str] = []

    # locates directory storing the results from previous step (results_filtered)
    # ignores hidden files (such as .DS_store)
    for filename in os.listdir(RESULTS_FILTERED):
        if filename.startswith("."):
            continue

        # for each file on "results_filtered":
        # runs the function create_protein_list and adds results to the list
        path_result = os.path.join(RESULTS_FILTERED, filename)
        with open(path_result, "r", encoding="ascii", errors="surrogateescape") as file:
            proteins_list.extend(extract_seqs.create_protein_list(file))

    # Runs "extract_proteins" for each protein.faa file on the path list (aka "raw_data_file"),
    # stores the results on dataset_P and writes them on a fasta file
    # Names file accordingly to each species (takes name of species from filepath)
    for filepath in raw_data_file_prot:
        dataset_p = extract_seqs.extract_proteins(filepath, proteins_list)
        species_name = filepath.split('/')[-3]
        with open(f"{RESULTS_PROTEINS}/{species_name}_prot.fna", "w") as output_file:
            SeqIO.write(dataset_p, output_file, "fasta")


################################################
# Make files for one protein and multiple species
################################################

def make_multispecies_prot_file():
    '''Runs the extract_seqs.multispecies_gene() function on the results from the step above
    (removing low quality sequences)'''

    for filename in os.listdir(RESULTS_PROTEINS):
        if filename.startswith("."):
            continue

        filename_altered = filename.replace("_prot.fna", "")
        gene_file = os.path.join(
            f"{RESULTS_FILTERED}/{filename_altered}_filtered.fna")

        path_result = os.path.join(RESULTS_PROTEINS, filename)
        extract_seqs.multispecies_prot(gene_file,
                                       path_result, RESULTS_MULTISPECIES_PROT)


########################
# Main Program Execution
########################

# list_path, raw_data_file = process_nucleotide_file()
# extract_longest_isoform(raw_data_file)
# extract_genes_of_interest(list_path)
# remove_low_quality_sequences()
# edit_sequence_names()
make_multispecies_gene_file()
# raw_data_file = process_protein_file()
# extract_proteins_of_interest(raw_data_file)
# make_multispecies_prot_file()
