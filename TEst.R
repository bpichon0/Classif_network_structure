
source("Fonctions_network_inference.R")

network= #bipartite unweighted network

network=matrix(as.numeric(read.table("./Example/A_HP_001.csv",sep=",")>0),
               10,
               18)
colnames(network)=paste0("P_",1:ncol(network))
rownames(network)=paste0("PL_",1:nrow(network))

g <- graph_from_incidence_matrix(network)

E(g)$arrow.mode = "-"
V(g)$color <- ifelse(V(g)$type == TRUE, "lightblue", "lightgreen")

plot(g, 
     layout=layout_as_bipartite, 
     arrow.mode=0,
     vertex.label=NA,
     vertex.size=6,
     asp=2)


#laplacian graph
band=0.025
spectre= spectR_network_abs(network, n = 200,bandwith=band)
plot_spectR_network(spectR = spectre,bw = band)
motif_freq=mcount(,
           six_node = T, normalisation=T, mean_weight=F, standard_dev=F)


#simulating the bipartite evol model



