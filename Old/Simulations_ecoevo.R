rm(list=ls())

library(reshape2)
library(ggplot2)
library(stringr)
library(FactoMineR)
library(factoextra)
library(neuralnet)
library(nnet)

setwd("~/APT 2/Stage/Interactions planctoniques/Scripts")
source("Fonctions.R")
dir.create("../Data/bipartite_sim")
dir.create("../Data/bipartite_sim/Data")


#######################################################################################
#  This script use the simulation from an eco-evolutionary model (Maliet et al., 2020).
#  It's extracted, then we create a data frame containing all the data to train the
#  neuralnetwork. 
#######################################################################################



#EXTRACTION OF THE DATA

list_simul=list.files(path = "../Data/bipartite_sim/network")

connectance_simul=NULL
network_simul=NULL
value_simul=NULL
sum_simul=NULL
diff_simul=NULL
type_simul=c()

u=1
p=1
for (i in list_simul){
  network=read.csv(paste0("../Data/bipartite_sim/network/",i),sep=";",header=T)
  network=as.matrix(network)
  type=na.omit(str_extract(i,c("mutualism","antagonism","neutral")))[1]
  network=weighted_sample(network,taux=70)[[1]]
  if (dim(network)[1]>9 & dim(network)[2]>9 & type!="neutral"){
    network[network>0]=1
    network=network[rowSums(network)>0,]
    network=network[,colSums(network)>0]
    spectre= spectR_network_abs(network, n = 200)
    if (spectre$bw>0.04){spectre= spectR_network_abs(network,bandwith = 0.025, n = 200)}
    connectance_simul[u]=sum(network)/(dim(network)[1]*dim(network)[2])
    network_simul[[u]]=network
    value_simul = rbind.data.frame(value_simul, spectre$y)
    diff_simul[u]=abs(dim(network)[1]-dim(network)[2])
    sum_simul[u]=abs(dim(network)[1]+dim(network)[2])
    type_simul=c(type_simul,as.character(type))
    u=u+1
    print(u)
    
  }
  p=p+1
}

value_simul=value_simul[1:100]

type_simul[type_simul=="mutualism"]="BipartiteEvol Mutualistic"
type_simul[type_simul=="antagonism"]="BipartiteEvol Antagonistic"
type_simul[type_simul=="neutral"]="BipartiteEvol Neutral"



#motif analysis

count_motif_ecoevo=NULL
for (i in network_simul){
  i=as.matrix(i)
  dat=mcount(i, six_node = T, normalisation=T, mean_weight=F, standard_dev=F)
  count_motif_ecoevo=rbind(count_motif_ecoevo,dat$normalise_sum)
}

c=c()
for (m in 1:dim(count_motif_ecoevo)[2]){
  c=c(c,paste0("Motif_",m))
}

count_motif_ecoevo=as.data.frame(count_motif_ecoevo)
colnames(count_motif_ecoevo)=c
write.table(count_motif_ecoevo,"../Data/bipartite_sim/Data/count_motif_BipartiteEvol.csv",sep=";")



#NEstedness modularity

Nestedness=NODF_network(network_simul)
Modularity=modularity_network(network_simul)



data_bipartiteevol=data.frame(category_DB=type_simul,
                          type_DB="BipartiteEvol Simulations",
                          diff_DB=diff_simul,
                          sum_DB=sum_simul,
                          connectance_DB=connectance_simul,
                          value_DB=value_simul,
                          nestedness=Nestedness,
                          modularity=Modularity,
                          motif=count_motif_ecoevo,
                          stringsAsFactors=FALSE)

colnames(data_bipartiteevol)=c("Category","Type","Diff","Sum","Connectance",paste0("Var",1:100),"Nestedness","Modularity",
                           paste0("Motif_",1:44))


write.table(data_bipartiteevol,file="../Data/bipartite_sim/Data/data_BipartiteEvol.csv",sep=";")
save(network_simul,file="../Data/bipartite_sim/Data/network_simul.RData")
