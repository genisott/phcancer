import numpy as np
import pandas as pd
from functions import *


# Read the expression file directly exported form the R downloader.
# Auxiliar file attached.

ctype = "ACC"
exp = pd.read_csv("expression_"+ctype+".txt").transpose() # Samples x genes

#Subset genes from the L1000 genes.
file = open("1000genes.txt")
genes = [gene[:-1] for gene in file.readlines()]

genes = sorted(set(genes).intersection(exp.columns))
samples = exp.index
nsamples = len(samples)
ngenes = len(genes)

exp = exp[genes] #Samples x genes

# Log transform
exp = np.log2(exp+1)
exp.to_csv("expression_"+ctype+"_1000log.txt")

# Main loop. 
# Computes lioness for each sample and then runs PH.
# Then generates the output files.

agg = exp.corr() # genes x genes

for i in range(0,nsamples):
    net = lioness_sample(exp,agg,samples,i) #Run lioness for the sample i.
    net = norm.cdf(net) #Transform into probabilities.
    for j in range(0,ngenes):
       net[j,j] = 1  #Set the diagonal to 1. (artificial tho).
    diag = persistence_coexp(net) #Computes the persistence of the network.
    plot_persistence(diag,str(i),str(ctype))
    save_diag(diag,str(i),str(ctype))

# Create distance matrices:
# The distance matrices are created from the files since I don't accumulate in memory the different diags.
d0,d1 = distance_matrix(nsamples, samples, ctype)

# Cluster with K-means

k02 = KMeans(n_clusters=2).fit(d0)
k12 = KMeans(n_clusters=2).fit(d1)



cl = pd.DataFrame(k02.labels_)
cl.columns=["k02"]
cl["k12"]=k12.labels_
cl.index=d0.index
cl.to_csv("clusterings_"+ctype+".txt")

# Create design files to be used in Limma.
for i in cl.columns:
    design(cl,i,ctype)


