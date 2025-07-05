###################################
# Running IQ-TREE in parallel
###################################

# This script runs IQ-TREE in parallel for multiple alignments.
# It takes a folder of FASTA files as input and generates trees for each alignment.

# Author: Letícia Magpali

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
import argparse

# Define a function to run IQ-TREE
def run_iqtree(sequence, seq_type, model, prefix, threads):
    command = f"iqtree -s {sequence} -st {seq_type} -m {model} -pre {prefix} -nt {threads}"
    print(f"Running: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Completed: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error running {command}: {e}")

if __name__ == "__main__":
    # Create the argument parser
    parser = argparse.ArgumentParser(description="Runs IQTree in parallel. Usage: python run_iqtree.py <run_folder> <out_folder> <seq_type> <model> <threads>")
    # Add arguments
    parser.add_argument("run_folder", help="Path to folder containing fasta alignments used as input for IQTree.")
    parser.add_argument("out_folder", help="Path to folder where the output files will be saved.")
    parser.add_argument("seq_type", help="Type of sequence (DNA, CODON, AA, or AUTO).")
    parser.add_argument("model", help="Model from IQTree used for the analysis.")
    parser.add_argument("threads", type=int, help="Number of threads to use for each IQTree run.")
    parser.add_argument("--max_parallel_runs", type=int, default=10, help="Maximum number of parallel IQTree runs.")
    args = parser.parse_args()

    # Store arguments as variables
    run_folder = args.run_folder
    out_folder = args.out_folder
    seq_type = args.seq_type
    model = args.model
    threads = args.threads
    max_parallel_runs = args.max_parallel_runs

    # Use ThreadPoolExecutor to parallelize the runs
    with ThreadPoolExecutor(max_workers=max_parallel_runs) as executor:
        # Creates an empty list to story future objects (asynchronous execution of tasks submitted to the executor)
        # Allows to keep track of the tasks that have been submitted
        # and their completion status
        futures = []
        for alignment in os.listdir(run_folder):
            if alignment.endswith(".fasta"):
                gene_name = alignment.split("_")[0]
                alignment_path = os.path.join(run_folder, alignment)
                prefix_path = f"{out_folder}/{gene_name}_{model}_tree"
                futures.append(executor.submit(run_iqtree, alignment_path, seq_type, model, prefix_path, threads))

        # Wait for all futures to complete
        for future in futures:
            future.result()

    print("All IQ-TREE runs completed.")