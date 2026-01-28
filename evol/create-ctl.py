###################################
# Generates control file for codeml
###################################

# Authors: Letícia Magpali & Rafael Copstein

'''This script generates a control file for codeml by editing a 
file template and making multiple copies of it for each gene.
Each control file for each gene/analysis is stored in its own folder.
'''

import os
from pathlib import Path
import argparse


if __name__ == "__main__":
    # Adding main arguments
    parser = argparse.ArgumentParser(description="Generate codeml control files for multiple genes. Requires a run folder, a paths file with paths to alignments and trees, an analysis name, a hypothesis number, and a control file template.")
    parser.add_argument("run_folder")
    parser.add_argument("paths_file")
    parser.add_argument("analysis")
    parser.add_argument("hypothesis")
    parser.add_argument("control_file_template")

    # Adding codeml options as optional arguments  
    parser.add_argument("--model", type=str, help="model option for codeml")
    parser.add_argument("--NSsites", type=str, help="NSsites option for codeml")
    parser.add_argument("--omega", type=str, help="omega option for codeml")
    parser.add_argument("--fix_omega", type=str, help="fix_omega option for codeml")
    parser.add_argument("--CodonFreq", type=str, help="CodonFreq option for codeml")
    parser.add_argument("--ncatg", type=str, help="ncatG option for codeml")
    parser.add_argument("--noisy", type=str, help="noisy option for codeml")
    parser.add_argument("--verbose", type=str, help="verbose option for codeml")
    parser.add_argument("--seqtype", type=str, help="seqtype option for codeml")
    parser.add_argument("--ndata", type=str, help="ndata option for codeml")
    parser.add_argument("--icode", type=str, help="icode option for codeml")
    parser.add_argument("--cleandata", type=str, help="cleandata option for codeml")
    parser.add_argument("--estFreq", type=str, help="estFreq option for codeml")
    parser.add_argument("--clock", type=str, help="clock option for codeml")
    args = parser.parse_args()

# Assigning main arguments to variables
RUN_FOLDER = args.run_folder # folder where the control files will be stored
PATHS_FILE = args.paths_file # file with paths to alignments folder and tree folder one per line
HYPOTHESIS = args.hypothesis # number/code of your hypothesis to be written in folders and outfiles
CONTROL_FILE_TEMPLATE = args.control_file_template # path or name of control file template
ANALYSIS = args.analysis # analysis name to be written in folders and outfiles


# Reading input files and assigning each line of the file to a path variable
with open(PATHS_FILE, "r") as paths_file:
    paths = [line.strip() for line in paths_file]

ALIGNMENTS_FOLDER, TREES_FOLDER = paths[:2]

# Ensure the paths are valid
RUN_FOLDER = Path(RUN_FOLDER)
ALIGNMENTS_FOLDER = Path(ALIGNMENTS_FOLDER)
TREES_FOLDER = Path(TREES_FOLDER)

# Print to verify
print(f"run folder: {RUN_FOLDER}")
print(f"alignments folder: {ALIGNMENTS_FOLDER}")
print(f"treees folder: {TREES_FOLDER}")

# Summarized way to print codeml options and their defaults
defaults = {
    "model": "0",
    "NSsites": "0",
    "CodonFreq": "7",
    "omega": "0.5",
    "fix_omega": "0",
    "ncatg": "None",
}

for opt, default in defaults.items():
    val = getattr(args, opt, None)
    print(f"{opt} = {val if val is not None else default} {'(default)' if val is None else ''}")

# Creates a dictionary of the codeml template file
codeml_dictionary = {}
with open(CONTROL_FILE_TEMPLATE) as control_file_template:
    for line in control_file_template:
        option, comment = line.split("*")
        (key, val) = option.split("=")
        key = key.strip()
        val = val.strip()
        codeml_dictionary[key] = val

# List of supported codeml options and their values from args
codeml_options = {
    "model": args.model,
    "NSsites": args.NSsites,
    "omega": args.omega,
    "fix_omega": args.fix_omega,
    "CodonFreq": args.CodonFreq,
    "ncatG": args.ncatg,
    "noisy": args.noisy,
    "verbose": args.verbose,
    "seqtype": args.seqtype,
    "ndata": args.ndata,
    "icode": args.icode,
    "cleandata": args.cleandata,
    "estFreq": args.estFreq,
    "clock": args.clock,
}

# This loop creates a control file from the dictionary
# and stores it on a run folder
for alignment in os.listdir(ALIGNMENTS_FOLDER):
    if alignment.endswith(".phy"):
        alignment_path = ALIGNMENTS_FOLDER/alignment
        # Extracts gene name from alignment file name
        gene_name = alignment.split("_")[0]
        # Take tree path from trees folder
        trees = os.listdir(TREES_FOLDER)
        found_tree = None
        for tree_name in trees:
            if tree_name.startswith(gene_name):
                found_tree = tree_name
                break  # Found the first matching tree, no need to search further
        if found_tree:
            # A tree starting with "gene name" was found
            tree_path = TREES_FOLDER/found_tree
        else:
            # No tree starting with "id1" was found in the list.
            tree_path = TREES_FOLDER/trees[0]

        # Modifies dictionary adding the values you specified above
        codeml_dictionary["seqfile"] = alignment_path
        codeml_dictionary["treefile"] = tree_path
        codeml_dictionary["outfile"] = f"out_{gene_name}_{ANALYSIS}-{HYPOTHESIS}.txt"

        # Update the dictionary with the codeml options
        for option, value in codeml_options.items():
            if value is not None:
                codeml_dictionary[option] = value

        # Creates a folder to store the ctl file
        run_subfolder = RUN_FOLDER/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}"

        run_subfolder.mkdir(parents=True, exist_ok=True)

        # and writes the dictionary to a text file to create a control file
        with open(run_subfolder/f"{ANALYSIS}-{HYPOTHESIS}_{gene_name}.ctl", "w", encoding="utf-8") as control_file:
            for key, val in codeml_dictionary.items():
                control_file.write(f"{key} = {val}\n")
