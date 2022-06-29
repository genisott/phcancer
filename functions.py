import pandas as pd
import numpy as np
import gudhi
import matplotlib.pyplot as plt
from scipy.stats import norm
import os
from sklearn.cluster import KMeans

#Functions

def persistence_coexp(matrix):
   """Function that computes the 0- and 1-dimensial Persistent Homology for a network given by "kinda" correlation matrix. """
   distance_matrix = np.array(1 - matrix) # Transform correlation matrix into distance matrix.
   rips_complex = gudhi.RipsComplex(distance_matrix=distance_matrix) # Create Rips Complex from the distance matrix.
   simplex_tree = rips_complex.create_simplex_tree(max_dimension=2) # Create the Simplex Tree from the Rips complex. This takes A LOT of memory. For our networks max_dim>2 crashes the server.
   diag = simplex_tree.persistence() # Compute the persistence homology of the simplex tree.
   return diag # Return the persistence homology information.

def plot_persistence(diag,label,ctype):
    """Plot all the persistence homology plots for a given network."""

    # Makes sure the output folders exist.
    cwd = os.getcwd()
    dir_density = os.path.join(cwd,"densities_"+ctype)
    dir_barcode = os.path.join(cwd,"barcodes_"+ctype)
    dir_diagram = os.path.join(cwd,"diagrams_"+ctype)

    if not os.path.exists(dir_density):
        os.mkdir(dir_density)
    if not os.path.exists(dir_barcode):
        os.mkdir(dir_barcode)
    if not os.path.exists(dir_diagram):
        os.mkdir(dir_diagram)
    
    # Plotting.
    gudhi.plot_persistence_density(diag)
    plt.savefig("densities_"+ctype+"/s"+label+"_density.png")
    plt.clf()

    gudhi.plot_persistence_barcode(diag)
    plt.savefig("barcodes_"+ctype+"/s"+label+"_barcode.png")
    plt.clf()

    gudhi.plot_persistence_diagram(diag)
    plt.savefig("diagrams_"+ctype+"/s"+label+"_diagram.png")
    plt.clf()

def save_diag(diag,label,ctype):
    """ Save the diagonal file so that PH has not to be run again..."""

    # Making sure the folder exists.
    cwd = os.getcwd()
    dir_diags = os.path.join(cwd,"diags_"+ctype)
    if not os.path.exists(dir_diags):
        os.mkdir(dir_diags)
    
    # Saves persistence data (birth/death).


    DF = pd.DataFrame(diag)
    pv = list(zip(*list(pd.DataFrame(diag)[1])))
    DF[1] = pv[0]
    DF[2] = pv[1]
    DF.to_csv("diags_"+ctype+"/s"+label+"_diag.txt")

def lioness_sample(exp,agg,samples,i):
    """ Computes lioness for a specific sample from expression data.

    exp <- Expression matrix, samples x genes.
    agg <- Aggregated correlation matrix. Passed as argument so that is is not computed again and again. And can run PH individually.
    samples <- name of the samples.
    i <- index of the current sample in the samples list.
    """

    ss = exp.iloc[exp.index!=samples[i],:].corr() # Compute the single sample correlation.
    ls = len(samples)*(agg-ss)+ss # Apply Lioness equation.
    return ls

def design(cluster,column,ctype):
    """ From the cluster dataframe generate design file to be used by the Limma script. """
    des = cl[[column]]
    des["m"] = np.abs(1-des)
    des.to_csv("design_"+ctype+"_"+column+".txt",header=None)

def distance_matrix(nsamples,samples,ctype): # HAS TO BE REVISED!!!
    dfs = list()
    for i in range(0,nsamples):
        dfs.append(pd.read_csv("diags_"+ctype+"/s"+str(i)+"_diag.txt",index_col=0))

    D0 = np.zeros((nsamples,nsamples)) #Takes the 0dim part of the diagonal file.
    for i in range (0,nsamples):
        for j in range(0,i+1):
            dfx = dfs[i][dfs[i]["0"]==0][["1","2"]]
            dfy = dfs[j][dfs[j]["0"]==0][["1","2"]]
            D0[i,j] = gudhi.bottleneck_distance(dfx,dfy)
            D0[j,i] = D0[i,j]

    D0=pd.DataFrame(D0)
    D0.columns,D0.index=samples,samples
    D0.to_csv(ctype+"_0dim_distance.txt")

    D1 = np.zeros((nsamples,nsamples)) #Takes the 1dim part of the diagonal file.
    for i in range (0,nsamples):
        for j in range(0,i+1):
            dfx = dfs[i][dfs[i]["0"]==1][["1","2"]]
            dfy = dfs[j][dfs[j]["0"]==1][["1","2"]]
            D1[i,j] = gudhi.bottleneck_distance(dfx,dfy)
            D1[j,i] = D1[i,j]
            
    D1=pd.DataFrame(D1)
    D1.columns,D1.index=samples,samples
    D1.to_csv(ctype+"_1dim_distance.txt")

    return D0,D1

