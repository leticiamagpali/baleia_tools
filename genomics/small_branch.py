treefile = '/Users/leticiamagpali/Google Drive/My Drive/PhD_Letícia/research_project/chapter2_genes/trees/gene_trees/iqtree_v4_Jan24/MF_codon_tree/CNTNAP2_MF.treefile'

f = open(treefile, 'r')
tree = f.read()


output = ""

for entry in tree.split(','):
    species = str(entry.split(':')).split(',')[0]
    species = species.replace('[', '')
    species = species.replace('(', '')
    species = species.replace("'", '')

    distance = str(entry.split(':')).split(',')[1]
    distance = distance.replace(']', '')
    distance = distance.replace(')', '')
    distance = distance.replace("'", '')
    distance = distance.replace("\\n", '')
    distance = float(distance.replace(";", ''))
    if (distance < 0.0002):
        output += species + "\t" + str(distance) + "\n"

newfile = "CNTNAP2_short_species.txt"

f = open(newfile, 'w')

f.write(output)
