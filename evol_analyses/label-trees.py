####### Label phylohgenetic trees #########
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

    def get_descendants_by_name(self, name: str):
        for node in self.traverse():
            if node.name == name:
                return node
        return None

if __name__ == "__main__":
    # Adding main arguments
    parser = argparse.ArgumentParser(description="Formats gene trees to codeml and labels foreground branches. Requires a run folder, a paths file with paths to alignments and trees, an analysis name, a hypothesis number, and a control file template.")
    parser.add_argument("hypothesis", type=str, help="Hypothesis number or code to be written in the output file name (e.g., 'M0', 'H1', 'H2')")
    parser.add_argument("tree_type", type=str, help="Type of tree, e.g., 'genetree', 'speciestree'")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing the trees")
    parser.add_argument("output_folder", type=str, help="Path to the output folder where labeled trees will be saved")

    # Adding optional arguments  
    parser.add_argument("--label", type=str, default="", help="Label to be added to the foreground branches (default: '')")
    parser.add_argument("--foreground", type=str, default="", help="Comma-separated list of foreground branches or clades to label in the tree (e.g., 'leaf1,leaf2,clade_limit_1:clade_limit_2')(default: '')")

    # Making a variable to store the arguments
    args = parser.parse_args()


FOREGROUND = args.foreground  
LABEL = args.label  
HYPOTHESIS = args.hypothesis  
TREE_TYPE = args.tree_type  
TREES_FOLDER = args.input_folder  
OUTPUT_FOLDER = args.output_folder

# Takes user input, a string of with the names of foreground branches
# separated by spaces, and turns it into a list
foreground_nodes = FOREGROUND.split(",")

# This loop goes through each tree in your trees folder,
# parses the tree and checks only for the leafs
# if a particular leaf is in the foreground list,
# it will add the user specified label
for file in os.listdir(TREES_FOLDER):
    if file.endswith(".tre"):
        gene_name = file.split("_")[0]
        tree_path = os.path.join(TREES_FOLDER, file)
        tree = EvolTree2(tree_path)

        if LABEL and FOREGROUND:
            # If a label or foreground is specified, mark the foreground nodes
            for element in foreground_nodes:
                if ":" not in element and element != "":
                    # this is a leaf
                    leaf = tree.get_descendants_by_name(element)
                    if leaf is not None:
                        tree.mark_tree([leaf.node_id], marks=[f"#{LABEL}"])
                elif ":" in element:
                    # this is a clade
                    clade_limits = element.split(":")
                    lim_left = tree.get_descendants_by_name(clade_limits[0])
                    lim_right = tree.get_descendants_by_name(clade_limits[1])
                    if lim_left is not None and lim_right is not None:
                        common = lim_left.get_common_ancestor(lim_right)
                        tree.mark_tree([common.node_id], marks=[f"${LABEL}"])

        # Always write the tree, labeled or not
        tree.write(outfile=f"{OUTPUT_FOLDER}/{gene_name}_{TREE_TYPE}_{HYPOTHESIS}.tre", format=13)
