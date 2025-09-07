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

## This is an algorithm that implements labling criterion for gene-tree species-tree disagreements.


# Import necessary libraries
import os
import argparse
from re import sub
from itertools import combinations
from ete3 import EvolTree


# Define a subclass of EvolTree to handle specific writing formats (codeml trees without branch lengths)
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
    parser = argparse.ArgumentParser(description="Formats trees to codeml and labels foreground branches. " \
    "Requires a an input folder with tree files in Newick format (.tre or .treefile) and the ETE3 library." \
    "Input trees can be gene trees or species trees. " \
    "Input trees should be named as <gene_name>_something.tre or <gene_name>_something.treefile. ")

    # Positional arguments
    parser.add_argument("input", type=str, help="Path to the folder containing the newick trees to be labeled")
    
    # Optional arguments
    parser.add_argument("--out", type=str, default="", help="Path to the output folder where labeled trees will be saved")
    parser.add_argument("--h", type=str, default="", help="Hypothesis number or code to be written in the output file name (e.g., 'H0', 'H1', 'H2'). If no hypothesis is specified, field will be left blank")
    parser.add_argument("--t", type=str, default="tree", help="Type of tree, e.g., 'genetree', 'speciestree'. If no type is specified, 'tree' will be used")
    parser.add_argument("--tag", type=str, default="", help="Label to be added to the foreground branches (default: '')")
    parser.add_argument("--leaves", type=str, default="", help="Comma-separated list of leaves to label in the tree (default: ''), example: species1,species2,species3")
    parser.add_argument("--ancestor", type=str, default="", help="Space-separated list of clades whose common ancestors will be labeled in the tree (default: ''), example: species1,species2 species3,species4 (species1 and species2 are a clade, species3 and species4 are another clade)")
    parser.add_argument("--clades", type=str, default="", help="Space-separated list of clades to label in the tree (default: ''), example: species1,species2 species3,species4 (species1 and species2 are a clade, species3 and species4 are another clade)")
    
    # Making a variable to store the arguments
    args = parser.parse_args()

# Assigning arguments to variables
LEAVES = args.leaves
CLADES = args.clades
ANCESTOR = args.ancestor
LABEL = args.label
HYPOTHESIS = args.hypothesis
TREE_TYPE = args.tree_type  
TREES_FOLDER = args.input_folder  
OUTPUT_FOLDER = args.output_folder

# If no output folder is specified, create one inside the input folder
if OUTPUT_FOLDER is None:
    OUTPUT_FOLDER = os.path.join(TREES_FOLDER, "labeled_trees")

# Extracting foreground leaves and clades from the input arguments
# Leaves are expected to be a comma-separated string, and clades a space-separated string
foreground_leaves = [leaf.strip() for leaf in LEAVES.split(",") if leaf.strip()]
foreground_clades = [clade.strip() for clade in CLADES.split(" ") if clade.strip()]
foreground_ancestors = [clade.strip() for clade in ANCESTOR.split(" ") if clade.strip()]

# Ensure the output folder exists before writing files
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# This loop goes through each tree in your trees folder,
# parses the tree and checks for the foreground leaves and clades to be labeled.
# Trees should be in Newick format with .tre or .treefile extensions.
# if a particular branch/clade is in the foreground list,
# it will add the user specified label
for file in os.listdir(TREES_FOLDER):
    if file.endswith(".tre") or file.endswith(".treefile"):
        gene_name = file.split("_")[0]
        tree_path = os.path.join(TREES_FOLDER, file)
        tree = EvolTree2(tree_path)

        # This line reads all nodes of the tree
        # this way, if the user inputs a foreground that's not present in the tree
        # it won't cause an error, just won't label anything
        present_node_ids = [node.node_id for node in tree.get_descendants()]

        # Labels to leaves and clades will only be added when specified by the user
        if LABEL and LEAVES:
            for leaf in foreground_leaves:
                leaf_to_label = tree.get_leaves_by_name(leaf)
                if leaf_to_label:
                    tree.mark_tree([leaf_to_label[0].node_id], marks=[f"#{LABEL}"])
                else:
                    print(f"Warning: Leaf {leaf} not found in tree {gene_name}. No leaves labelled.")
        
        # This section labels clades in the tree
        # It checks if the clade is monophyletic and labels the common ancestor
        # If the clade is not monophyletic, it finds the most inclusive subclade
        # and labels it with $1, while labeling the outsiders with #1      
        if LABEL and CLADES:
            for clade in foreground_clades:
                clade_species = [species for species in clade.split(",") if species in tree.get_leaf_names()]
                # Only proceed if at least one species is present
                if clade_species:  
                    monophyletic = tree.check_monophyly(clade_species, target_attr="name", ignore_missing=True)[0]
                    if monophyletic:
                        ancestor_to_label = tree.get_common_ancestor(clade_species)
                        if ancestor_to_label:
                            tree.mark_tree([ancestor_to_label.node_id], marks=[f"${LABEL}"])
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

                        if not list_of_subclades:
                            # If no subclades were found, label the entire clade
                            for leaf in clade_species:
                                leaf_to_label = tree.get_leaves_by_name(leaf)
                                if leaf_to_label and leaf_to_label[0].node_id in present_node_ids:
                                    tree.mark_tree([leaf_to_label[0].node_id], marks=[f"#{LABEL}"])

        if LABEL and ANCESTOR:
            for clade in foreground_ancestors:
                ancestor_species = [species for species in clade.split(",") if species in tree.get_leaf_names()]
                if ancestor_species: 
                    monophyletic = tree.check_monophyly(ancestor_species, target_attr="name", ignore_missing=True)[0]
                    if monophyletic:
                        ancestor_to_label = tree.get_common_ancestor(ancestor_species)
                        if ancestor_to_label:
                            tree.mark_tree([ancestor_to_label.node_id], marks=[f"#{LABEL}"])
                    else:
                        list_of_subclades = []
                        for r in range(2, len(ancestor_species)+1):
                            for subset in combinations(ancestor_species, r):
                                subclade_nodes = list(tree.get_monophyletic(list(subset), target_attr="name"))
                                if subclade_nodes:
                                    list_of_subclades.append((subset, subclade_nodes))

                        if list_of_subclades:
                            # Sort by length of subset, descending
                            most_inclusive = max(list_of_subclades, key=lambda x: len(x[0]))
                            # Tag all nodes in the most inclusive subclade
                            for node in most_inclusive[1]:
                                if node.node_id in present_node_ids:
                                    tree.mark_tree([node.node_id], marks=[f"#{LABEL}"])
                        
                        if not list_of_subclades:
                            print(f"Warning: No monophyletic clade found for ancestor {clade} in tree {gene_name}. No ancestors labelled.")
        
        # Always writes the tree, labeled or not -- useful for changing the tree to codeml format automatically
        tree.write(outfile=f"{OUTPUT_FOLDER}/{gene_name}_{TREE_TYPE}_{HYPOTHESIS}.tre", format=13)

# Future improvements:
# - Add error handling for file reading/writing
# - Add more flexible input options for clades and leaves (e.g., from a file)
# - Add a summary of labeled nodes at the end of the script