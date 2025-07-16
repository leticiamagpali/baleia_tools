####### Label phylogenetic trees #########
###########################################

# Authors: Letícia Magpali & Rafael Copstein
''' This script formats trees to codeml and labels foreground nodes (leaves and clades).

It handles discrepancies between gene and species trees, allowing to tailor the labeling to the topology of the tree.
For example, if you have a gene tree with species A, B, C, and D, but the species tree has A, B, and C as a clade,
you can label the branch leading to D as a foreground branch, while labeling the clade A, B, C as another foreground branch.

To run it, type python label-trees.py <hypothesis> <tree_type> <input_folder> <output_folder> --options
Requires trees in Newick format (.tre) and the ete3 library.

OBS: The EvolTree class of ete3 was slightly modified to allow writing trees in codeml format, 
without branch lengths and with labels.
'''


# Import necessary libraries
import os
import argparse
from re import sub
from itertools import combinations
from ete3 import EvolTree


# Define a subclass of EvolTree to handle specific writing formats
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

    ############# Argument parser #############
    parser = argparse.ArgumentParser(description="Formats trees to codeml and labels foreground branches. Requires a an input folder with tree files in Newick format (.tre) and the ETE3 library.")
    
    # Positional arguments
    parser.add_argument("hypothesis", type=str, help="Hypothesis number or code to be written in the output file name (e.g., 'H0', 'H1', 'H2')")
    parser.add_argument("tree_type", type=str, help="Type of tree, e.g., 'genetree', 'speciestree'")
    parser.add_argument("input_folder", type=str, help="Path to the folder containing the newick trees to be labeled")

    # Optional arguments
    parser.add_argument("--output_folder", type=str, help="Path to the output folder where labeled trees will be saved")
    parser.add_argument("--label", type=str, default="", help="Label to be added to the foreground branches (default: '')")
    parser.add_argument("--leaves", type=str, default="", help="Comma-separated list of leaves to label in the tree(default: ''), example: species1,species2,species3")
    parser.add_argument("--clades", type=str, default="", help="Space-separated list of clades to label in the tree(default: ''), example: species1,species2 species3,species4 (species1 and species2 are a clade, species3 and species4 are another clade)")
    
    # Making a variable to store the arguments
    args = parser.parse_args()

# Assigning arguments to variables
LEAVES = args.leaves
CLADES = args.clades
LABEL = args.label  
HYPOTHESIS = args.hypothesis  
TREE_TYPE = args.tree_type  
TREES_FOLDER = args.input_folder  
OUTPUT_FOLDER = args.output_folder

if OUTPUT_FOLDER is None:
    OUTPUT_FOLDER = os.path.join(TREES_FOLDER, "labeled_trees")


# Extracting foreground leaves and clades from the input arguments
# Leaves are expected to be a comma-separated string, and clades a space-separated string
foreground_leaves = [leaf.strip() for leaf in LEAVES.split(",") if leaf.strip()]
foreground_clades = [clade.strip() for clade in CLADES.split(" ") if clade.strip()]

# This loop goes through each tree in your trees folder,
# parses the tree and checks for the foreground leaves and clades to be labeled.
# if a particular branch/clade is in the foreground list,
# it will add the user specified label
for file in os.listdir(TREES_FOLDER):
    if file.endswith(".tre"):
        gene_name = file.split("_")[0]
        tree_path = os.path.join(TREES_FOLDER, file)
        tree = EvolTree2(tree_path)

        # Checking all nodes of the tree
        # this way, if the user inputs a foreground that's not present in the tree
        # it won't cause an error, just won't label anything
        present_node_ids = [node.node_id for node in tree.get_descendants()]

        # Labels to leaves and clades will only be added when specified by the user
        if LABEL and LEAVES:
            for leaf in foreground_leaves:
                leaf_to_label = tree.get_leaves_by_name(leaf)
                if leaf_to_label and leaf_to_label[0].node_id in present_node_ids:
                    tree.mark_tree([leaf_to_label[0].node_id], marks=[f"#{LABEL}"])
        
        # This section labels clades in the tree
        # It checks if the clade is monophyletic and labels the common ancestor
        # If the clade is not monophyletic, it finds the most inclusive subclade
        # and labels it with $1, while labeling the outsiders with #1      
        if CLADES:
            for clade in foreground_clades:
                clade_species = clade.split(",") 
                monophyletic = tree.check_monophyly(clade_species, target_attr="name", ignore_missing=True)[0] 
                if monophyletic:
                    clade_to_label = tree.get_common_ancestor(clade_species)
                    if clade_to_label and clade_to_label.node_id in present_node_ids:
                        tree.mark_tree([clade_to_label.node_id], marks=[f"${LABEL}"])
                if not monophyletic:
                    list_of_subclades = []
                    for r in range(2, len(clade_species)+1):
                        for subset in combinations(clade_species, r):
                            subclade_nodes = list(tree.get_monophyletic(list(subset), target_attr="name"))
                            if subclade_nodes:
                                list_of_subclades.append((subset, subclade_nodes))

                    if list_of_subclades:
                        # Sort by length of subset, descending
                        most_inclusive = max(list_of_subclades, key=lambda x: len(x[0]))
                        # Tag all nodes in the most inclusive subclade
                        for node in most_inclusive[1]:
                            if node.node_id in present_node_ids:
                                tree.mark_tree([node.node_id], marks=[f"${LABEL}"])

                        outsiders = [species for species in clade_species if species not in most_inclusive[0]]
                        for leaf in outsiders:
                            outsider_to_label = tree.get_leaves_by_name(leaf)
                            if outsider_to_label and outsider_to_label[0].node_id in present_node_ids:
                                tree.mark_tree([outsider_to_label[0].node_id], marks=[f"#{LABEL}"])

        # Always writes the tree, labeled or not -- useful for changing the tree to codeml format automatically
        tree.write(outfile=f"{OUTPUT_FOLDER}/{gene_name}_{TREE_TYPE}_{HYPOTHESIS}.tre", format=13)


# Future improvements:
# - Add section to process and label single ancestral nodes
# - Add a check for the existence of the output folder and create it if it doesn't exist
# - Add error handling for file reading/writing
# - Add more flexible input options for clades and leaves (e.g., from a file)
# - Add a summary of labeled nodes at the end of the script