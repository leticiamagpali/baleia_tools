from Bio import AlignIO
from Bio import SeqIO

INPUT_FILE = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/phd_leticia/phd_project/chapter1_acoustic_genes/alignments/big_dataset/aligned/codon_align/codon_align_fasta/RYR2_codon_aligned.fasta"
OUTPUT_FILE = "/Users/leticiamagpali/Library/CloudStorage/GoogleDrive-leticiamagpali@gmail.com/My Drive/phd_leticia/phd_project/chapter1_acoustic_genes/alignments/big_dataset/aligned/codon_align/codon_align_fasta/trimmed/RYR2_codon_aligned_filtered.fasta"
MAX_GAP_PERCENTAGE = 50.0  # maximum allowed gap percentage

# Read the alignment
alignment = AlignIO.read(INPUT_FILE, "fasta")
alignment_length = alignment.get_alignment_length()
max_gaps = (MAX_GAP_PERCENTAGE / 100.0) * alignment_length

# Filter sequences
filtered_sequences = []
for record in alignment:
    if record.seq.count("-") <= max_gaps:
        filtered_sequences.append(record)

# Write the filtered alignment
SeqIO.write(filtered_sequences, OUTPUT_FILE, "fasta")
