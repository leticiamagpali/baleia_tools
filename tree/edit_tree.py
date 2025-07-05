import re
import os

TREE_NAME = "NLGN3_codon.treefile"
TREE_PATH = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/phd_leticia/research_project/chapter2_genes/trees/gene_trees/small_dataset/NLGN3_tree/NLGN3_codon.treefile"
RESULTS_FOLDER = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/phd_leticia/research_project/chapter2_genes/evol_models/codeml/small_dataset/cml_trees/gene_trees/M0_gene_trees"

###### Remove branch lengths ######
###################################


def remove_branch_lengths():
    # Output file
    output_file = TREE_NAME.replace(".treefile", ".tre")

    # Read input file
    with open(TREE_PATH) as file:
        tree_with_branch = file.read()

    # Remove branch lengths
    tree_without_branch = re.sub(r":[0-9.]+", '', tree_with_branch)

    output_file = os.path.join(RESULTS_FOLDER, output_file)
    with open(output_file, "w", encoding="utf-8") as file:
        file.write(tree_without_branch)


remove_branch_lengths()
