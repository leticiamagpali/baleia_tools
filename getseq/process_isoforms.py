#### Remove Isoforms #####
##########################

# Author: Rafael Copstein
'''Fasta file is parsed, sequences are ordered from longest to shortest, 
and then for each sequence record the loop will check if a record with that gene name 
(aka distinguishing attribute) has already been added to the list attributes_added 
- if not, it will add it, and also add the seq record to the new_dataset list.'''


from os import PathLike
from typing import Optional
from Bio import SeqIO, SeqRecord


def sortseq_descending_length(sequences):
    '''Sorts sequences according to length'''
    return sorted(sequences, key=lambda r: -len(r))


def is_gene(path: PathLike) -> bool:
    '''Return true if path ends with .fna - i.e., if file is a dna file'''
    return path.endswith("cds_from_genomic.fna")


def is_protein(path: PathLike) -> bool:
    '''Return true if path ends with .faa - i.e., if file is a protein file'''
    return path.endswith("protein.faa")


def get_gene_name(seq_record: SeqRecord) -> Optional[str]:
    '''Parses the sequence record description (like the one below) and extracts the gene name.
    >lcl|NW_011888782.1_cds_XP_023384417.1_10 [gene=AMOTL2] [db_xref=GeneID:105306059] 
    [protein=angiomotin-like protein 2 isoform X3] [protein_id=XP_023384417.1] 
    [location=join(89845..90578,93231..93537,94580..94724,95044..95136,99130..99425,
100197..100504,101254..101471,102998..103059)] [gbkey=CDS]'''
    gene_description = seq_record.description.split(" ")
    for element in gene_description:
        if element.startswith("[gene="):
            return element[6:-1]
        if element.startswith("[locus_tag="):
            return element[11:-1]
    return None


def get_protein_name(seq_record: SeqRecord) -> Optional[str]:
    '''Parses the sequence record description (like the one below) and returns the protein name.
    OBS: sometimes the protein will have an isoform, which we are not interested in.
    >NP_001292124.1 cyclin-dependent kinase inhibitor 1B [Pteropus vampyrus] 
    >NP_001292124.1 cyclin-dependent kinase inhibitor 1B isoform X1 [Pteropus vampyrus]'''
    protein_description = seq_record.description.split(" ")

    # Updates the protein description to remove the first element (aka the ID number)
    protein_description = protein_description[1:]

    for element, index in protein_description:
        # accomodates for proteins that have isoform or not
        if element == "isoform" or element.startswith("["):
            # takes everything in the protein description until that index
            protein_description = protein_description[:index]
            break  # stops the loop

    # protein_description is a list, and the elements are the words that make up the protein's name
    # we need to join all these elements into a string, separated by spaces, using the line below:
    return " ".join(protein_description)


def keep_longest_isoform(filepath):
    '''Sorts sequences according to length from longest to shortest.
    Parses the sequence file and checks if 
    the gene/protein name is on the sequence_names set. 
    If not, it adds the name to the set and adds the sequence record to the dataset list'''
    sequence_names = set()
    dataset = []

    sequences = list(SeqIO.parse(filepath, "fasta"))
    sequences_ordered = sortseq_descending_length(sequences)

    for seq_record in sequences_ordered:
        # guarantees that the loop starts with the variable empty
        seq_name: Optional[str] = None

        if is_gene(filepath):
            seq_name = get_gene_name(seq_record)
        if is_protein(filepath):
            seq_name = get_protein_name(seq_record)

        # Prints an error message if the gene or protein name is not found
        if seq_name is None:
            print(f"[!] SeqRecord attribute not found: {seq_record}")
            continue

        if seq_name not in sequence_names:
            sequence_names.add(seq_name)
            dataset.append(seq_record)

    return dataset
