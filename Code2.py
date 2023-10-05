from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Read in gene sequences and their corresponding names from file
with open('output.fasta', 'r') as f:
    gene_data = [line.strip() for line in f if line.strip()]

# Split gene data into sequence names and sequences
sequence_names = [line[1:] for i, line in enumerate(gene_data) if i % 2 == 0]
gene_sequences = [line for i, line in enumerate(gene_data) if i % 2 == 1]

# Calculate pairwise similarities between all sequences
similarities = np.zeros((len(gene_sequences), len(gene_sequences)))
for i, j in combinations(range(len(gene_sequences)), 2):
    seq1, seq2 = gene_sequences[i], gene_sequences[j]
    if len(seq1) > 0:
        similarity = sum(1 for a, b in zip(seq1, seq2) if a == b) / len(seq1)
    else:
        similarity = 0.0
    similarities[i, j] = similarity
    similarities[j, i] = similarity

# Perform hierarchical clustering on the similarities
Z = linkage(similarities, method='complete')

# Plot dendrogram with sequence names
plt.figure(figsize=(10, 5))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Gene Index')
plt.ylabel('Distance')
dendrogram(Z, leaf_rotation=90., leaf_font_size=8., labels=sequence_names)
plt.show()
