# Persistent homology on Lioness single sample networks
Pipeline of my work on the usage of persistence homology to cluster single sample Lioness networks.

1. Download data from TCGA using tcga_download.R
2. Expression matrix is subsetted for the L1000 genes.
3. Expression counts are log transformed.
4. For each sample we run Lioness. Transform the score into probabilities. Artificially set the diagonal to 1. We consider this sort of a "correlation similarity" matrix.
5. This is transformed into a "distance"/similarity network.
6. Run Persistence homology using Gudhi
7. Export all figures and diagonal.
8. Create distance matrices using the bottleneck distance. (Maybe this can be improved to use the Wasserstein distance)
9. Cluster using k-means, generally with k=2.
10. Create design file to be used with Limma in R.
11. In R run Limma and export the t-scores.
12. Pre-ranked gene set enrichment test (gsea).
