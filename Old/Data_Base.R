rm(list = ls())

library(ggpubr)
library(ggrepel)
library(readr)
library(bipartite)
library(maxnodf)
library(igraph)
library(bmotif)
setwd("~/APT 2/Stage/Interactions planctoniques/Scripts")
source("Fonctions.R")


#______________________________________________________________________________#
#  This script create a data frame containing the spectral characteristics of all
#  the network present in Michalska-Smith et al., 2019 + the trophic network in 
#  Web of life and some mycorhize networks. At the end, we compile XXX networks 
#  depending on the size limit fixed for empirical networks
#______________________________________________________________________________#




#### MYCORHIZAL NETWORKS ------


value_myco=NULL #value of the spectral density
connectance_myco=NULL#connectance of empirical network 
empiric_myco=NULL#stock the empirical network  
diff_myco=NULL#dimention of each network
sum_myco=NULL#dimention of each network
type_myco=NULL
category_myco=NULL
name_myco=NULL
band=0.025
list_files <- list.files(pattern = c('network'),path = '../Data/Web_of_life/Myco' )
metadata=read_csv("../Data/Web_of_life/Myco/metadata_mycorrhizal.csv")

m=1
for(i in list_files){
  network <- read.table(paste0("../Data/Web_of_life/Myco/" ,i),sep=",")
  if (dim(network)[1]>9 & dim(network)[2]>9){
    network[is.na(network)]<-0 # Only numbers
    name_myco[m]=i
    #remove empty rows and columns
    network=network[rowSums(network)!=0,] 
    network=network[,colSums(network)!=0]
    network[network>0]=1 # make the networks binary
    empiric_myco[[m]]=network
    spectre= spectR_network_abs(network, n = 200,bandwith=band)
    if (spectre$bw>0.04){spectre= spectR_network_abs(network,bandwith = band, n = 200)}
    value_myco = rbind.data.frame(value_myco, spectre$y)
    connectance_myco[m] = sum(network) / (dim(network)[1] * dim(network)[2])
    sum_myco[m]=dim(network)[1]+dim(network)[2]
    diff_myco[m]=abs(dim(network)[1]-dim(network)[2])
    m=m+1
  }
}

name_myco=gsub('.csv','', name_myco)
name_myco=gsub('network_','', name_myco)

p=1
type_myco=rep("Mycorhizal",length(name_myco))
category_myco=rep("Mutualistic",length(type_myco))




#### NETWORKS FROM MICHASLKA AND ALLESINA (2019; Plos Computational biology) ------


network_DB=NULL
eigen_DB=NULL
value_DB=NULL
connectance_DB=NULL 
diff_DB=NULL
sum_DB=NULL
count=1
type_DB=NULL
category_DB=NULL
metadata <- read_csv("../Data/Ploscomp_2019_network/Metadata.csv")


list_DB=list.files(pattern = '.csv',path = '../Data/Ploscomp_2019_network/network')

list_DB <- list_DB[which(!list_DB %in% c("Orford_2016_Field-19.csv","Orford_2016_Field-55.csv"))] 
#we don't use this networks because they're causing error


for (i in list_DB){
  edge_list<-read_csv(paste0("../Data/Ploscomp_2019_network/network/" ,i), col_names=FALSE)
  edge=edge_list  
  g <- graph_from_data_frame(edge_list, directed=FALSE)
  V(g)$type <- V(g)$name %in% edge_list$X2
  network <- get.incidence(g)
  if (dim(network)[1]>9 & dim(network)[2]>9 & dim(network)[2]<400 & dim(network)[1]<400){
    network_DB[[count]]=network
    spectre= spectR_network_abs(network, n = 200,bandwith=band)
    if (spectre$bw>0.04){spectre= spectR_network_abs(network,bandwith = band, n = 200)}
    eigen_DB[[count]]=spectre$eigenvalues
    value_DB = rbind.data.frame(value_DB, spectre$y)
    connectance_DB[count] = sum(network) / (dim(network)[1] * dim(network)[2])
    sum_DB[count]=dim(network)[1]+dim(network)[2]
    diff_DB[count]=abs(dim(network)[1]-dim(network)[2])
    i=gsub('.csv','', i)
    type_DB[count]=metadata$feature_2[which(metadata$name==i)]
    category_DB[count]=metadata$feature_1[which(metadata$name==i)]
    
    count=count+1
  }
}


category_DB[which(category_DB=="mutualism")]="Mutualistic"
category_DB[which(category_DB=="antagonism")]="Antagonistic"


type_DB<- gsub('herbivory*',replacement = 'Herbivory',x = type_DB)
type_DB<- gsub('parasitism',replacement = 'Host-Parasite',x = type_DB)
type_DB<- gsub('anemone-fish',replacement = 'Anemone-Fish',x = type_DB)
type_DB<- gsub('ant-plant',replacement = 'Plant-Ant',x = type_DB)
type_DB<- gsub('pollination',replacement = 'Pollination',x = type_DB)
type_DB<- gsub('seed dispersal',replacement = 'Seed dispersal',x = type_DB)
type_DB<- gsub('bacteria-phage',replacement = 'Bacteria-Phage',x = type_DB)
type_DB<- gsub('host-parasitoid',replacement = 'Host-Parasitoid',x = type_DB)


c=c()
for (i in 1:(ncol(value_DB))) {
  c = c(c, paste0("Var", i))#nomme les 200 valeurs de la densité spectrale
}

colnames(value_DB)=c
colnames(value_myco)=c


networks=c(network_DB,empiric_myco)
value_DB=rbind(value_DB,value_myco)
connectance_DB=c(connectance_DB, connectance_myco)
diff_DB=c(diff_DB,diff_myco)
sum_DB=c(sum_DB,sum_myco)
type_DB=c(type_DB,type_myco)
category_DB=c(category_DB,category_myco)


value_DB=value_DB[1:100]

# modularity=modularity_network(networks)
# nestedness=NODF_network(networks)

count_motif=NULL
for (i in networks){
  i=as.matrix(i)
  dat=mcount(i, six_node = T, normalisation=T, mean_weight=F, standard_dev=F)
  count_motif=rbind(count_motif,dat$normalise_sum)
}

c=c()
for (m in 1:dim(count_motif)[2]){
  c=c(c,paste0("Motif_",m))
}

count_motif=as.data.frame(count_motif)
colnames(count_motif)=c

write.table(count_motif,"../Data/Ploscomp_2019_network/Data/count_motif_BD.csv",sep=";")


data_BD=data.frame(category_DB=category_DB,
                 type_DB=type_DB,
                 diff_DB=diff_DB,
                 sum_DB=sum_DB,
                 connectance_DB=connectance_DB,
                 value_DB=value_DB,
                 #nestedness=nestedness,
                 #modularity=modularity,
                 motif=count_motif,
                 stringsAsFactors=FALSE)



ind_myco=grep(pattern = "network",x=data_BD$category_DB)
data_BD$category_DB[ind_myco]="Mutualistic"
data_BD$type_DB[ind_myco]="Mycorhizal"
colnames(data_BD)=c("Category","Type","Diff","Sum","Connectance",paste0("Var",1:100),"Nestedness","Modularity",
                    paste0("Motif_",1:44))

dir.create("../Data/Ploscomp_2019_network/Data/")
write.table(data_BD,"../Data/Ploscomp_2019_network/Data/data_BD.csv",sep=";")

