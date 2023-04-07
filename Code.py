# The output of the code would be a dendrogram plot showing
# the hierarchical clustering of the gene sequences based on their pairwise similarities.
# The plot will have the x-axis as the gene index and the y-axis as the distance between the clusters.
# The labels for the gene sequences will be shown on the leaves of the dendrogram.
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage

# Read in gene sequences from file
with open('seq.txt', 'r') as f:
    gene_sequences = [line.strip() for line in f if line.strip()]

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

# Plot dendrogram
plt.figure(figsize=(10, 5))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Gene Index')
plt.ylabel('Distance')
dendrogram(Z, leaf_rotation=90., leaf_font_size=8., labels=range(len(gene_sequences)))
plt.show()



