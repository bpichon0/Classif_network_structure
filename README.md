# Code for running analyses of inference of ecological networks from their structure

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform the analyses for both main text and supplementary informations.
All the code was made on R (*v4.1.0*).

Here is the different steps to reproduce the figures:


## `Replicating the figures`

The file `Make_figs.R` contains the different steps to replicate all analysis in the paper. It is organized into different chunks of code, each doing a separate figure. 
**Make sure** that you have all the packages needed to perform the analyses (see `Functions_network_inference.R` for the list).



## `Computing the metrics on a single network`


Below, we provide an example on how to compute the different metrics of network structure, from micro-scale metrics to macro-scale ones (Fig. 1).  

<p align="center">
    <img src="https://github.com/bpichon0/Classif_network_structure/blob/master/Example/Scale_network_structure.jpg" width="1000">
</p>

<p align="center">
    <b>Figure 1: Scale dependant analysis of ecological network structure  .</b>
</p>


We take as an example a plant-pollinator network:

```R

include("Fonctions_network_inference.R") #loading the functions

network=read.table("./Example/A_HP_001.csv",sep=",")>0 #bipartite unweighted network

g = graph_from_incidence_matrix(network)
E(g)$arrow.mode = "-"
V(g)$color = ifelse(V(g)$type == TRUE, "lightblue", "lightgreen")

plot(g,layout=layout_as_bipartite, arrow.mode=0,
     vertex.label=NA,vertex.size=6,asp=0.3)

```

<p align="center">
    <img src="https://github.com/bpichon0/Classif_network_structure/blob/master/Example/PL_network.svg" width="600">
</p>

<p align="center">
    <b>Figure 2: Structure of plant-pollinator network.</b>
</p>


We can compute modularity, nestedness as well as the spectral density of the Laplacian graph.
```R
computeModules(network,method = "Beckett") #modularity using Beckett method
NODFc(network) #nestedness using NODFc metric

band=0.025 #defining a bandwidth
spectre= spectR_network_abs(network, n = 200,bandwith=band) #we keep n=200 points along the density

plot_spectR_network(spectR = spectre,bw = band)


```



<p align="center">
    <img src="https://github.com/bpichon0/Classif_network_structure/blob/master/Example/spectral_density.svg" width="600">
</p>

<p align="center">
    <b>Figure 3: Spectral density of the Laplacian graph of the plant-pollinator network.</b>
</p>

In addition, we can extract the motif composition of this network using the *mcount* function from the **bmotif** package.



```R
motif_freq=mcount(six_node = T, normalisation=T, mean_weight=F, standard_dev=F)
Plot_motif_frequency(motif_freq) #we normalize each motif occurrence by the total count of motifs
```



<p align="center">
    <img src="https://github.com/bpichon0/Classif_network_structure/blob/master/Example/Motif_frequency.svg" width="600">
</p>

<p align="center">
    <b>Figure 3: Motif composition in the plant-pollinator network (frequency sqrt transformed). </b>
</p>



