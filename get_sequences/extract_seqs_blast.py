import sys
from Bio import SeqIO

# Print error message if user inputs wrong number of arguments
if len(sys.argv) != 4:
    print("Use: extract_seqs_blast.py <file.fasta> <GENE_IDs_n_Coords.list> <Orgasnism.list>")
    sys.exit(1)

# This script has 3 arguments:
fasta_file = sys.argv[1]
geneIDs_coords = sys.argv[2]
organisms_file = sys.argv[3]

geneIDs_coords_dict = dict()
organisms_list = list()

# Make a list of organisms from which sequences will be downloaded
with open(organisms_file) as file:
    for line in file:
        if line not in organisms_list:
            organisms_list.append(line.strip())

# Add gene IDs and coords to a dictionary:
with open(geneIDs_coords) as file:
    for line in file:
        if line.split()[0] not in geneIDs_coords_dict:
            geneIDs_coords_dict[line.split()[0].strip()] = []
        geneIDs_coords_dict[line.split()[0]].append(
            (line.split()[1].strip(), line.split()[2].strip()))

# For each organism in the list, and for each gene ID and coordinate in the list,
# this loop retrieves the fasta file from the alignment you inputed
for org in organisms_list:
    for key, value in geneIDs_coords_dict.items():
        for rec in SeqIO.parse(fasta_file, "fasta"):
            if org == rec.id:
                with open(key + '.fasta', 'a') as out_file:
                    out_file.write(
                        f">{org}_{key}\n{rec.seq[int(value[0][0]):int(value[0][1])+1]}\n")
