#### Editing alignments ####
############################

'''This module has functions to edit multispecies alignments'''

# Importing necessary packages
import os
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse
from filedit import TextToList

def is_alignment(file_path, format="fasta"):
    '''This function checks if a file is an alignment.'''
    # Check if file is an alignment by trying to read it with AlignIO
    # If it raises a ValueError, it's not an alignment
    # If it raises an IOError, it's not a valid file
    try:
        AlignIO.read(file_path, format)
        print(f"Your sequence {file_path} is an alignment.")
        return True
    except ValueError:
        print(f"Your sequence {file_path} is not an alignment.")
        return False
    except IOError as e:
        print(f"Error reading file {file_path}: {e}")
        return False

def write_sequences(output_path, sequences, format="fasta"):
    ''' This function writes sequences to a file in the specified format.'''
    with open(output_path, "w", encoding="utf-8") as output_file:
        if isinstance(sequences, MultipleSeqAlignment):
            AlignIO.write(sequences, output_file, format)
        else:
            SeqIO.write(sequences, output_file, format)

def extract_seqs_id(list_path, seq_folder):
    '''This function extracts sequences from a file based on a list of the sequence IDs'''
    seqs_list = TextToList.create_list(list_path)
    seqs_set = set(seqs_list)

    # For each sequence in the folder, get the path
    # and the gene name from the sequence name
    # and create the output file path
    for sequence in os.listdir(seq_folder):
        if sequence.endswith(".fasta"):
            sequence_path = os.path.join(seq_folder, sequence)
            gene_name = sequence.split("_")[0]
            output_file_path = os.path.join(seq_folder, f"{gene_name}_codon_aligned_d2.fasta")

        # Check if the sequence is an alignment or not
        # If it is an alignment, use AlignIO to parse it
        if is_alignment(sequence_path):
            # Read the single alignment into a MultipleSeqAlignment object
            # Filter the records based on IDs in seqs_set
            # Write the filtered sequences to the output file
            alignment = AlignIO.read(sequence_path, "fasta")
            extracted_seqs = [record for record in alignment if record.id in seqs_set]
            write_sequences(output_file_path, extracted_seqs, "fasta")
        else:
            sequence_file = SeqIO.parse(sequence_path, "fasta")
            extracted_seqs = [record for record in sequence_file if record.id in seqs_set]
            write_sequences(output_file_path, extracted_seqs, "fasta")

# ensures code is only run when the script is executed directly
# and not when imported as a module
if __name__ == "__main__":
    # Create the argument parser
    # and add arguments for the list path and sequence folder
    # This allows the script to be run from the command line
    parser = argparse.ArgumentParser(description="Extract sequences by ID.")
    parser.add_argument("list_path", help="Path to the list of sequence IDs.")
    parser.add_argument("sequence_folder", help="Folder containing sequence files.")
    args = parser.parse_args()

    extract_seqs_id(args.list_path, args.sequence_folder)
