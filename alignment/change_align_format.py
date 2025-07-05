import os
from Bio import AlignIO

# Run folders for changing alignment format
ALIGNMENT_FOLDER = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/PhD_Letícia/research_project/chapter2_genes/alignments/codon_align_phylip_codeml/codon_alignments_v4_Jan24"
ALIGNMENT_FOLDER_CHANGED = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/PhD_Letícia/research_project/chapter2_genes/alignments/codon_align_fasta"

#################
# Phylip to fasta
#################

for alignment in os.listdir(ALIGNMENT_FOLDER):
    if alignment.endswith(".phy"):
        alignment_path = os.path.join(ALIGNMENT_FOLDER, alignment)
        original_alignment = AlignIO.parse(alignment_path, "phylip-relaxed")
        with open(f"{ALIGNMENT_FOLDER_CHANGED}/{alignment}.fasta", "w") as reformated_alignment:
            AlignIO.write(original_alignment, reformated_alignment, "fasta")
