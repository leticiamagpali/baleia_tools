####### Label phylogenetic trees #########
###########################################

# Author: Letícia Magpali
''' This script labels tree nodes and branches.
To run it, type python label-trees.py <foreground branches> <label> <path to trees folder>
Requires trees in Newick format (.tre) and the ETE3 library.
'''

import os
import argparse
from re import sub
from ete3 import EvolTree


class EvolTree2(EvolTree):

    def write(self, features=None, outfile=None, format=13):
        if format == 13:
            nwk = super().write(format=11)
            nwk = sub(r'\[&&NHX:mark=([ #\$0-9.]*)\]', r'\1', nwk)
        else:
            nwk = super().write(features=features, format=format)

        if outfile is not None:
            with open(outfile, "w", encoding="UTF-8") as output_file:
                output_file.write(nwk)
            return nwk
        else:
            return nwk

if __name__ == "__main__":
    # Adding main arguments
    parser = argparse.ArgumentParser(description="Formats gene trees to codeml and labels foreground branches. Requires a run folder, a paths file with paths to alignments and trees, an analysis name, a hypothesis number, and a control file template.")
    parser.add_argument("hypothesis", type=str, help="Hypothesis number or code to be written in the output file name (e.g., 'M0', 'H1', 'H2')")
    parser.add_argument("tree_type", type=str, help="Type of tree, e.g., 'genetree', 'speciestree'")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing the trees")
    parser.add_argument("output_folder", type=str, help="Path to the output folder where labeled trees will be saved")

    # Adding optional arguments  
    parser.add_argument("--label", type=str, default="", help="Label to be added to the foreground branches (default: '')")
    parser.add_argument("--leaves", type=str, default="", help="Comma-separated list of leaves to label in the tree(default: '')")
    parser.add_argument("--clades", type=str, default="", help="Space-separated list of clades to label in the tree(default: '')")
    # Making a variable to store the arguments
    args = parser.parse_args()


LEAVES = args.leaves
CLADES = args.clades
LABEL = args.label  
HYPOTHESIS = args.hypothesis  
TREE_TYPE = args.tree_type  
TREES_FOLDER = args.input_folder  
OUTPUT_FOLDER = args.output_folder

# Takes user input, a string of with the names of foreground branches
# separated by spaces, and turns it into a list
foreground_leaves = [leaf.strip() for leaf in LEAVES.split(",") if leaf.strip()]
foreground_clades = [clade.strip() for clade in CLADES.split(" ") if clade.strip()]

# This loop goes through each tree in your trees folder,
# parses the tree and checks only for the leafs
# if a particular leaf is in the foreground list,
# it will add the user specified label
for file in os.listdir(TREES_FOLDER):
    if file.endswith(".tre"):
        gene_name = file.split("_")[0]
        tree_path = os.path.join(TREES_FOLDER, file)
        tree = EvolTree2(tree_path)

        if LABEL and LEAVES:
            # If a label and foreground is specified, mark the foreground nodes
            for leaf in foreground_leaves:
                leaf_to_label = tree.get_leaves_by_name(leaf)
                if leaf_to_label:
                    tree.mark_tree([leaf_to_label[0].node_id], marks=[f"#{LABEL}"])
                
        if CLADES:
            for clade in foreground_clades:
                clade_species = clade.split(",") # Now a list of species names
                monophyletic = tree.check_monophyly(clade_species, target_attr="name")[0] 
                if monophyletic:
                    # If the clade is monophyletic, label the common ancestor
                    clade_to_label = tree.get_common_ancestor(clade_species)
                    tree.mark_tree([clade_to_label.node_id], marks=[f"${LABEL}"])
                if not monophyletic:
                    subclades = tree.get_monophyletic(clade_species, target_attr="name")
                    for node in subclades:
                        if node.is_leaf():
                        # Species out of clade, label leaf
                            tree.mark_tree([node.node_id], marks=[f"#{LABEL}"])
                        else:
                        # Subclade, label common ancestor
                            tree.mark_tree([node.node_id], marks=[f"${LABEL}"])
                

        # Always write the tree, labeled or not
        tree.write(outfile=f"{OUTPUT_FOLDER}/{gene_name}_{TREE_TYPE}_{HYPOTHESIS}.tre", format=13)
