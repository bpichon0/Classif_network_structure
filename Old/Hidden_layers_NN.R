rm(list=ls())
library(neuralnet)
library(viridis)
library(MLmetrics)
setwd("~/APT 2/Stage/Interactions planctoniques/Scripts")
source("Fonctions.R")

################################################################################
# This script aims to identify the optimal parameters of the intermediate layer 
# of the neural network for each of the data (motifs, global metrics). For 
# spectral density, the intermediate neuron size has been fixed to 60. For the 
# combination we will add the number of neurons between data. This works as an
# additive method.
################################################################################



################################################################################
#                      PARAMETRISATION On GLOBAL METRICS                          #
################################################################################
data_BD <- read.csv("../Data/Ploscomp_2019_network/Data/data_BD.csv",sep=";")
data_BD=data_BD[,c(1:5,106,107)]

AA=NULL
MM=NULL
taille=NULL
p=1
for (x in c(3,4,5,6,7,8,9,10,15,20)){
  aa=NULL
  mm=NULL
  g=1
  for (l in 1:10){
    data.test=NULL
    data.train=NULL
    set.seed(sample(1:10000,1))
    sample_test = sample(1:nrow(data_BD),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(data_BD),
                         replace = F)
    data.test <- data_BD[sample_test,]
    data.train <- data_BD[-sample_test,]
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = x , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn <- predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test <- prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) # If value>0.5, classify as nested
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    g=g+1
  }
  AA[p]=mean(aa)
  MM[p]=mean(mm)
  taille[p]=x
  p=p+1
}

plot(y=(AA+MM)/2,x=taille)






################################################################################
#                      PARAMETRISATION        ON Motifs                           #
################################################################################
data_BD <- read.csv("APT 2/Stage/Interactions planctoniques/Wgithub/data_BD.csv",sep=" ")
data_BD=data_BD[,c(1,2,108:151)]

AA=NULL
MM=NULL
taille=NULL
p=1
for (x in c(10,15,20,25,30,35,40,45,50,60,70)){
  aa=NULL
  mm=NULL
  g=1
  for (l in 1:20){
    data.test=NULL
    data.train=NULL
    set.seed(sample(1:10000,1))
    sample_test = sample(1:nrow(data_BD),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(data_BD),
                         replace = F)
    data.test <- data_BD[sample_test,]
    data.train <- data_BD[-sample_test,]
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = x , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn <- predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test <- prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) # If value>0.5, classify as nested
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    g=g+1
  }
  AA[p]=mean(aa)
  MM[p]=mean(mm)
  taille[p]=x
  p=p+1
}

plot(y=(AA+MM)/2,x=taille)





################################################################################
#                      PARAMETRISATION     ON    Motifs + spectral density        #
################################################################################
data_BD <- read.csv("APT 2/Stage/Interactions planctoniques/Wgithub/data_BD.csv",sep=" ")
data_BD=data_BD[,-c(3:5,106,107)]

AA=NULL
MM=NULL
taille=NULL
p=1
for (x in c(60,70,80,90,100,110,120,130)){
  aa=NULL
  mm=NULL
  g=1
  for (l in 1:20){
    data.test=NULL
    data.train=NULL
    set.seed(sample(1:10000,1))
    sample_test = sample(1:nrow(data_BD),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(data_BD),
                         replace = F)
    data.test <- data_BD[sample_test,]
    data.train <- data_BD[-sample_test,]
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = x , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn <- predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test <- prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) # If value>0.5, classify as nested
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    g=g+1
  }
  AA[p]=mean(aa)
  MM[p]=mean(mm)
  taille[p]=x
  p=p+1
}

plot(y=(AA+MM)/2,x=taille)







