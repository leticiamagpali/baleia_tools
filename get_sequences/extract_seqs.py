'''
04-03-2022

Authors: Leticia Magpali & Rafael Copstein

This module extracts individual sequences from a multi-sequence file 
according to functional terms and removes low quality sequences.
To utilize these tools, access through run_get_seqs.py
'''
from typing import List
from Bio import SeqIO
from file_editing import TextToList


###############################
# Extract genes from .fna files
###############################

def extract_genes(file, genes_list):
    '''Insert function description'''

    genes_list = TextToList.create_list(genes_list)

    # Transforms genes_list into a set (data structure that
    # is optimized to check whether an element is whithin it or not)
    genes_set = set(genes_list)

    # creates empty list to store the sequence records of the extracted genes
    extracted_genes = []

    # reads sequence file
    for record in SeqIO.parse(file, "fasta"):

        # returns 1st field of sequence description split by whitespaces (aka [gene=GENE_NAME])
        # edits [gene=GENE_NAME] to display only GENE_NAME
        # if GENE_NAME matches the gene names on a list (user-provided):
        # appends sequence record to the extracted_genes list
        # OBS: the result of seq_record.description.split()[1] is a string
        if record.description.split()[1].lstrip("[gene=").rstrip("]") in genes_set:
            extracted_genes.append(record)
    return (extracted_genes)


def fiter_low_quality(filepath):
    '''Insert function description'''

    filtered_genes = []

# parses cds.fna sequence file for each species, from filepath in genomics directory
    # write error if parse function doesn't find the sequence
    for record in SeqIO.parse(filepath, "fasta"):

        # prints only sequences which do not have the string protein=LOW QUALITY PROTEIN" on their descriptions
        if "protein=LOW QUALITY PROTEIN" not in record.description:
            filtered_genes.append(record)
    return (filtered_genes)


#############################################
# Edit sequence names (description) on files
#############################################

def change_seq_description(gene_file, species_name, output_folder):
    '''Parses a fasta file and chages the description to include only gene name and species name'''
    for record in SeqIO.parse(gene_file, "fasta"):
        gene_name = record.description.split()[1].lstrip("[gene=").rstrip("]")
        record.id = f"{gene_name}_{species_name}"
        record.description = ""
        with open(f"{output_folder}/{species_name}_nuc.fasta", 'a+', encoding="utf-8") as output_gene_file:
            SeqIO.write(record, output_gene_file, "fasta")


##############################################
# Make files for one gene and multiple species
##############################################

def multispecies_gene(seq_file, genes_list, output_folder):
    '''Takes a file with multiple genes for a single species 
    and turns into one file for each gene for multiple species'''

    genes = TextToList.create_list(genes_list)

    for record in SeqIO.parse(seq_file, "fasta"):
        gene_name = record.description.split("_")[0]
        if gene_name in genes:
            with open(f"{output_folder}/{gene_name}_nuc.fasta", 'a+', encoding="utf-8") as gene_file:
                SeqIO.write(record, gene_file, "fasta")


##################################
# Extract proteins from .faa files
##################################

def create_protein_list(file) -> List[str]:
    '''Insert function description'''

    # creates empty list to store protein names
    prot_list = []

    # parses fasta file to extract protein name from header
    # by splitting the header with blank space as separator
    # extracting the element that contains "protein_id="
    # and then removing unecessary terms ([protein= and ])
    # finally, adds protein name on a list
    for record in SeqIO.parse(file, "fasta"):
        for element in record.description.split():
            if "protein_id=" in element:
                prot_list.append(element.lstrip("[protein_id=").rstrip("]"))

    return prot_list


def extract_proteins(file, prot_list):
    '''Insert function description'''

    # creates a list to store the sequence records of the extracted proteins
    extracted_proteins = []

    # Transforms prot_list into a set (data structure that
    # is optimized to check whether an element is whithin it or not)
    prot_set = set(prot_list)

    # reads sequence file
    for record in SeqIO.parse(file, "fasta"):
        # splits sequence description/header by whitespaces and returns 1st item (aka protein ID)
        # if protein name matches the IDs on the list (created above):
        # appends sequence record to the extracted_proteins list
        if record.description.split()[0].lstrip(">") in prot_set:
            extracted_proteins.append(record)

    return (extracted_proteins)


##############################################
# Make files for one gene and multiple species
##############################################

def make_gene_prot_dict(gene_file):
    '''Makes a dictionary relating gene name to protein id'''
    gene_prot_dict = {}

    for record in SeqIO.parse(gene_file, "fasta"):
        gene_name = record.description.split()[1].lstrip("[gene=").rstrip("]")
        prot_id = record.description.split("[")[4].lstrip(
            "protein_id=").rstrip("] ")
        gene_prot_dict.setdefault(gene_name, prot_id)

    return (gene_prot_dict)


def multispecies_prot(gene_file, prot_file, output_folder):
    '''Takes a file with multiple genes for a single species 
    and turns into one file for each gene for multiple species'''

    gene_prot = make_gene_prot_dict(gene_file)
    proteins = create_protein_list(gene_file)

    for record in SeqIO.parse(prot_file, "fasta"):
        prot_record = record.description.split()[0]
        if prot_record in proteins:
            for key, value in gene_prot.items():
                if prot_record == value:
                    gene_name = key
                    break
            with open(f"{output_folder}/{gene_name}_prot.fasta", 'a+', encoding="utf-8") as gene_file:
                SeqIO.write(record, gene_file, "fasta")

        # makes sure functions will work when called from this file AND also when called from run.py
        # if __name__ == '__main__':
        # extract_genes()
        # fiter_low_quality()
        # multispecies_gene()
        # multispecies_prot()
