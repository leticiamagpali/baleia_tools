#########################
# MANIPULATING ALIGNMENTS
#########################

import os
import extract_seqs_id
from Bio import SeqIO
from Bio import AlignIO

# Run folders for alignment extraction
RUN_FOLDER = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/genes-multispecies"
RESULTS_FOLDER = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/alignments"
LIST_PATH = "/home/chordata/Desktop/Bielawski_Lab/MSc/results/maleXfemale/alignments/maleXfemale.list"


######################################
# Extracting sequences from fasta file
######################################

for sequence_file in os.listdir(RUN_FOLDER):
    if sequence_file.endswith(".fasta"):
        file_path = os.path.join(RUN_FOLDER, sequence_file)

        extracted_seqs = extract_seqs_id.extract_seqs_id(
            file_path, LIST_PATH)
        with open(f"{RESULTS_FOLDER}/{sequence_file}", "w", encoding="utf-8") as output_file:
            SeqIO.write(extracted_seqs, output_file, "fasta")


######################################
# Extracting sequences from alignment
######################################

# for alignment in os.listdir(RUN_FOLDER):
    # if alignment.endswith(".fasta"):
        # alignment_path = os.path.join(RUN_FOLDER, alignment)

        # extracted_seqs = ExtractSeqsID.extract_seqs_id(
            # alignment_path, LIST_PATH)
        # with open(f"{RESULTS_FOLDER}/{alignment}", "w", encoding="utf-8") as output_file:
            # AlignIO.write(extracted_seqs, output_file, "fasta")
