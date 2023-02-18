rm(list=ls())


source("Fonctions_network_inference.R")

dir.create("./Figures/",showWarnings = F)
dir.create("./Figures/SI",showWarnings = F)


#empirical data 

load("./network_BD.RData")
data_empiric = read.csv("./Data/data_BD.csv",sep=";")
data_empiric$Category=as.character(data_empiric$Category)
data_empiric$Category[which(data_empiric$Category=="Empirical antagonistic")]="Antagonistic"
data_empiric$Category[which(data_empiric$Category=="Empirical mutualistic")]="Mutualistic"
save=data_empiric



#Data from simulations

save_bipartite=read.csv("../Data/bipartite_sim/Data/data_BipartiteEvol.csv",sep=";")

combination=list(c(1:3,4,5,106,107), # global metrics
                 c(1,2,108:151), # motifs
                 c(1,2,6:105), #spectral density
                 c(1:105), #global metrics + spectral density
                 c(1:5,108:151), #global metrics + motifs
                 c(1,2,5:106,108:151), #spectral density + motifs
                 c(1:151)) #spectral density + motifs + global metrics

name_combination=list("global metrics",
                      "motifs",
                      "spectral density",
                      "global metrics + spectral density",
                      "global metrics + motifs",
                      "spectral density + motifs",
                      "spectral density + motifs + global metrics")



#******************************************************************************#

# Main figures ----

#******************************************************************************#

## >> Fig metrics comparison mutualistic, antagonistic ----

#First global metrics

data_empiric=save[order(data_empiric$Category),c(1:2,5,106,107)] 
data_empiric$Number=as.character(1:343)
global_metric_empiric=melt(data_empiric)
colnames(global_metric_empiric)=c("Category","Type","Number","Variable","Value")

p1=ggplot(global_metric_empiric%>%
               mutate(., Category=recode_factor(Category,"Mutualistic"="Mutualistic networks","Antagonistic"="Antagonistic networks" )),
             aes(x=Variable, y=Value,fill=Variable)) + xlab("")+ylab("log(value)")+
  scale_y_log10()+facet_wrap( ~ Category, ncol = 2) +
  geom_line(aes(group=Number),alpha = 0.5, colour = "darkgrey") +
  geom_violin(alpha=0.8,trim=T)+theme_minimal()+scale_fill_manual(values=c("#ba4a00","#d68910","#229954"))+
  theme(legend.position="none")+
  geom_boxplot(width=0.1,outlier.shape = NA,color="black")+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        strip.text.x = element_text(colour = "black", size=12),
        axis.title.x = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12))




# secondly, laplacian graph


#Difference in eigen density between antagonistic and mutualistic networks
laplacian_data=save[,c(1,2,6:105)]

eigen_pos=NULL;freq=NULL;Nature=NULL;Dat=NULL


for (i in 3:ncol(laplacian_data)){
  data_mixed_model=laplacian_data[,c(1,2,i)]
  colnames(data_mixed_model)[3]="eigen"
  if (wilcox.test(data_mixed_model$eigen[which(data_mixed_model$Category=="Mutualistic")],
                  data_mixed_model$eigen[which(data_mixed_model$Category=="Antagonistic")])$p.value<0.01){ #we keep those that significantly differ using mixed effect models
    
    eigen_pos=rbind(as.numeric(i-2),as.numeric(i-2)) 
    freq=rbind(aggregate(x = laplacian_data[,i],by = list(laplacian_data$Category),FUN = mean)[2][2,],
               aggregate(x = laplacian_data[,i],by = list(laplacian_data$Category),FUN = mean)[2][1,])
    Nature=rbind("Mutualistic","Antagonistic")
    Dat=rbind(Dat,cbind(as.numeric(eigen_pos),as.numeric(freq),as.character(Nature)))
  }
}



Dat=data.frame(Eig=as.numeric(Dat[,1]),Freq=as.numeric(Dat[,2]),Nature=as.character(Dat[,3]))
data_2 = melt(as.matrix(laplacian_data[,3:ncol(laplacian_data)]))
data_2$value=sqrt(data_2$value)
data_2$Category = rep(laplacian_data$Category, ncol(laplacian_data)-2)

data_2 = data_2[which(data_2$Var2 %in%  paste0("Var",unique(Dat$Eig))),]

data_2$Var2=round(seq(0,1,length.out=100)[as.numeric(str_remove(data_2$Var2,pattern = "Var"))],3)
p2=ggplot(data_2, aes(x=as.factor(Var2), y=sqrt(value), fill=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic()+# theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+
  xlab("Eigenvalues") +ylab(expression(sqrt('Density of eigen values')))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+ylim(0,2)+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")+
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) 




#lastly: motifs

count_motif=save[,c(1,2,108:151)]

Motif_name=NULL;freq=NULL;Nature=NULL;Dat=NULL


for (i in 3:ncol(count_motif)){
  data_mixed_model=count_motif[,c(1,2,i)]
  colnames(data_mixed_model)[3]="Motif"
  if (wilcox.test(data_mixed_model$Motif[which(data_mixed_model$Category=="Mutualistic")],
                  data_mixed_model$Motif[which(data_mixed_model$Category=="Antagonistic")])$p.value<0.01){ #we keep those that significantly differ using mixed effect models
    
    Motif_name=rbind(as.numeric(i-2),as.numeric(i-2)) #motif number starts at 1 and index at 3
    freq=rbind(aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][2,],
               aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][1,])
    Nature=rbind("Mutualistic","Antagonistic")
    Dat=rbind(Dat,cbind(as.numeric(Motif_name),as.numeric(freq),as.character(Nature)))
  }
}


Dat=data.frame(Motif=as.numeric(Dat[,1]),Freq=as.numeric(Dat[,2]),Nature=as.character(Dat[,3]))

data_2 = melt(as.matrix(count_motif[,3:ncol(count_motif)]))
data_2$value=sqrt(data_2$value)
data_2$Category = rep(count_motif$Category, ncol(count_motif)-2)

data_2 = data_2[which(data_2$Var2 %in%  paste0("Motif_",unique(Dat$Motif))),]

p3=ggplot(data_2%>%
            mutate(., Category=recode_factor(Category,"Mutualistic"="Mutualistic networks","Antagonistic"="Antagonistic networks" )),
          aes(x=reorder(Var2, -value), y=value, fill=Category,col=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+
  scale_color_manual(values=c('black','black'))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")



Fig1=ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.1,1,.1)),
               p2+theme(legend.position = "none"),
               ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),ncol=3,widths=c(.05,1,.05)),
               nrow=3,labels = LETTERS[1:3],font.label = list(size=20))

ggsave("./Figures/Fig1.pdf",Fig1,width = 9,height = 12)
png("./Figures/Fig1.png",Fig1,width = 9,height = 12)


## >> Fig PCA for the 3 type of metrics ----



# PCA global metrics 

data_empiric=save[,c(1:3,4,5,106,107)]
res.pca=PCA(data_empiric[,3:dim(data_empiric)[2]], scale.unit = T, ncp = 2,  graph=F)

p1=pca_metrics=fviz_pca_ind(res.pca, geom.ind = "point", 
                               axes=c(1,2), 
                               pointsize = 2, 
                               col.ind = data_empiric$Category,
                               palette=c('#DECF3F','#8fd175'),
                               addEllipses = F,
                               label = "var",
                               repel = T,mean.point = FALSE)+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
       y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),col.ind="")+
  ggtitle("") +theme_classic()+theme(legend.position = "bottom",
                                     legend.text = element_text(size=14))

km.res = kmeans(res.pca$ind$coord, 2, nstart = 9)
Kmeans_table=table(km.res$cluster,data_empiric$Category)
write.table(Kmeans_table,"./Figures/Kmeans_globalmetrics.csv",sep=";")
chisq.test(Kmeans_table,simulate.p.value = T)




# Laplacian eigen values

data_empiric=save[,c(1,2,6:105)]
res.pca=PCA(data_empiric[,3:dim(data_empiric)[2]], scale.unit = T, ncp = 2,  graph=F)

p2=fviz_pca_ind(res.pca, geom.ind = "point", 
             axes=c(1,2), 
             pointsize = 2, 
             col.ind = data_empiric$Category,
             palette=c('#DECF3F','#8fd175'),
             #addEllipses = TRUE,ellipse.level = .95,
             label = "var",
             repel = T,mean.point = FALSE)+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
       y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),col.ind="")+
  ggtitle("") +theme_classic()+theme(legend.position = "bottom",
                                     legend.text = element_text(size=14))+
  xlim(-12,24)



km.res = kmeans(res.pca$ind$coord, 2, nstart = 9)
Kmeans_table=table(km.res$cluster,data_empiric$Category)
write.table(Kmeans_table,"./Figures/Kmeans_Laplacian.csv",sep=";")

chisq.test(Kmeans_table,simulate.p.value = T)



# Motifs frequency

data_empiric=save[,c(1,2,108:151)]%>%
  mutate(., Category=recode_factor(Category,"Mutualistic"="Mutualistic networks","Antagonistic"="Antagonistic networks" ))
res.pca=PCA(data_empiric[,3:dim(data_empiric)[2]], scale.unit = T, ncp = 2,  graph=F)

p3=fviz_pca_ind(res.pca, geom.ind = "point", 
                              axes=c(1,2), 
                              pointsize = 2, 
                              col.ind = data_empiric$Category,
                              palette=c('#DECF3F','#8fd175'),
                              addEllipses = F,
                              label = "var",
                              repel = T,mean.point = FALSE)+
  labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
       y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),col.ind="")+
  ggtitle("") +theme_classic()+theme(legend.position = "bottom",
                                     legend.text = element_text(size=14))
km.res = kmeans(res.pca$ind$coord, 2, nstart = 9)
Kmeans_table=table(km.res$cluster,data_empiric$Category)
write.table(Kmeans_table,"./Figures/Kmeans_Motifs.csv",sep=";")
chisq.test(Kmeans_table)

Fig2=ggarrange(p1+theme(legend.position = "none")+ggtitle("Global metrics")+theme(plot.title = element_text(size=18)),
               p2+theme(legend.position = "none")+ggtitle("Laplacian graph")+theme(plot.title = element_text(size=18)),
               p3+labs(color="",fill="",shape="")+ggtitle("Motifs frequency")+theme(plot.title = element_text(size=18))+
                 guides(color = guide_legend(override.aes = list(size = 4)),shape="none"),
               nrow=3,labels = LETTERS[1:3],font.label = list(size=20))
ggsave("./Figures/Fig2.pdf",Fig2,width = 7,height = 12)






## >> Table 1 : Classification of empirical networks using empirical networks ----

data_empiric=save
data_empiric$Category=as.character(data_empiric$Category)
data_empiric$Category[which(data_empiric$Category=="Antagonistic network")]="Antagonistic"
data_empiric$Category[which(data_empiric$Category=="Mutualistic network")]="Mutualistic"

percentage_classif=rep(0,nrow(data_empiric))
presence=rep(0,nrow(data_empiric))
g=1
nb_tours=10000
for (b in sample(1:10000,nb_tours)){
  data.test=NULL
  data.train=NULL
  set.seed(b)
  sample_test = sample(1:nrow(data_empiric),size = 0.2*nrow(data_empiric),replace = F)
  
  
  data.test = data_empiric[sample_test,]
  data.train = data_empiric[-sample_test,]
  data.train=rbind(data.train)
  
  
  formule = as.formula(c('Category ~ ', paste0(colnames(data.train[, 3:length(data.train)]), collapse = "+"))) #classify networks according to the type
  nn_1layer = neuralnet(formule,data.train,hidden = 60 , 
                        stepmax = 10 ** 8,lifesign = 'full', 
                        linear.output=F,rep=1) # number of max step to train
  
  prediction_nn = predict(nn_1layer, data.test[,3:length(data.test)])
  prediction_nn_test = prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
  
  u=1
  for (i in prediction_nn_test){
    
    if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
    if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
    u=u+1
  }
  prediction_nn_test=as.data.frame(prediction_nn_test[,2])
  prediction_nn_test$classif=data_empiric$Category[as.numeric(rownames(prediction_nn_test))]
  colnames(prediction_nn_test)=c("Prediction","Classification") #Classification = TRue category
  presence[as.numeric(rownames(prediction_nn_test))]=presence[as.numeric(rownames(prediction_nn_test))]+1
  for (k in 1:dim(prediction_nn_test)[1]){
    colnames(prediction_nn_test)=c("Prediction","Classification") 
    if (prediction_nn_test$Prediction[k]==prediction_nn_test$Classification[k]){
      percentage_classif[as.numeric(rownames(prediction_nn_test[k,]))]=percentage_classif[as.numeric(rownames(prediction_nn_test[k,]))]+1
    }
    
  }
  g=g+1
}

confiance=data.frame(percentage=percentage_classif/presence,category=data_empiric$Category,type=data_empiric$Type)

write.table(confiance,"./Figures/Data/data_Table_1.csv",sep=";")

confiance=read.table("./Figures/Data/data_Table_1.csv",sep=";")

confiance$percentage=confiance$percentage*100
confiance=confiance[order(confiance$category,-confiance$percentage),]
confiance$number=c(1:table(confiance$category)[1],1:table(confiance$category)[2])
confiance$classif_final=" " 

# We classify according to the mean of the results : 
# For example, if a mutualistic network, is classified less than 
# 50% of the time mutualistic, then, the final classification would
# be antagonistic. Furthermore, we distinguish between high (between 75 and 100%) and low confidence (between 50 and 75)

for (i in 1:nrow(confiance)){
  if (confiance$percentage[i]>=50 && confiance$percentage[i]<75){
    if (confiance$category[i]=="Antagonistic"){
      confiance$classif_final[i]="Antagonistic_low"
    }
    else {confiance$classif_final[i]="Mutualistic_low"}
  }
  if (confiance$percentage[i]>=75){
    if (confiance$category[i]=="Antagonistic"){
      confiance$classif_final[i]="Antagonistic_high"
    }
    else {confiance$classif_final[i]="Mutualistic_high"}
  }
  if (confiance$percentage[i]<50 && confiance$percentage[i]>25) {
    if (confiance$category[i]=="Antagonistic"){
      confiance$classif_final[i]="Mutualistic_low"
    }
    else {confiance$classif_final[i]="Antagonistic_low"}  
  }
  if (confiance$percentage[i]<=25) {
    if (confiance$category[i]=="Antagonistic"){
      confiance$classif_final[i]="Mutualistic_high"
    }
    else {confiance$classif_final[i]="Antagonistic_high"}  
  }
}

results = table(confiance$category,confiance$classif_final)

confiance$type=data_empiric$Type[as.numeric(rownames(confiance))]
Table_1=table(confiance$type,confiance$classif_final)
write.table(Table_1,"./Figures/Table_1.csv",sep=";")




















#******************************************************************************#

# SI Figs ----

#******************************************************************************#



## >> Tlobal metrics of empirical networks ----


data_empiric=save[,c(1:5,106,107)]
colnames(data_empiric)=c("Nature of the network",colnames(data_empiric)[2:7])


data_empiric$`Nature of the network`=as.character(data_empiric$`Nature of the network`)
data_empiric$`Nature of the network`[which(data_empiric$`Nature of the network`=="Empirical antagonistic")]="Antagonistic"
data_empiric$`Nature of the network`[which(data_empiric$`Nature of the network`=="Empirical mutualistic")]="Mutualistic"

panel_a=ggplot(data_empiric,aes(x=Type,y=Connectance,fill=`Nature of the network`,col=`Nature of the network`))+
  geom_boxplot()+
  theme_classic()+xlab("")+   ylab("Connectance")+theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+
  scale_color_manual(values=c('black','black'))+theme(legend.position="none",
                                                      axis.text.x = element_blank(),
                                                      axis.ticks = element_blank(),
                                                      axis.title.y = element_text(size=12))

panel_b=ggplot(data_empiric,aes(x=Type,y=Nestedness,fill=`Nature of the network`,col=`Nature of the network`))+
  geom_boxplot()+
  theme_classic()+xlab("")+  ylab("Nestedness")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+
  scale_color_manual(values=c('black','black'))+theme(legend.position="none",
                                                      axis.text.x = element_blank(),
                                                      axis.ticks = element_blank(),
                                                      axis.title.y = element_text(size=12))

panel_c=ggplot(data_empiric,aes(x=Type,y=Modularity,fill=`Nature of the network`,col=`Nature of the network`))+
  geom_boxplot()+
  theme_classic()+xlab("")+ ylab("Modularity")+labs(color="",fill="")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+
  scale_color_manual(values=c('black','black'))+theme(legend.position="none",
                                                      axis.text.x = element_blank(),
                                                      axis.ticks = element_blank(),
                                                      axis.title.y = element_text(size=12))

panel_d=ggplot(data_empiric,aes(x=Type,y=Diff,fill=`Nature of the network`,col=`Nature of the network`))+
  geom_boxplot()+
  theme_classic()+xlab("")+ ylab("Difference \nof guilds size")+labs(color="",fill="")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+ylim(0,100)+
  scale_color_manual(values=c('black','black'))+theme(legend.position="none",
                                                      axis.text.x = element_blank(),
                                                      axis.ticks = element_blank(),
                                                      axis.title.y = element_text(size=12))


panel_e=ggplot(data_empiric,aes(x=Type,y=Sum,fill=`Nature of the network`,col=`Nature of the network`))+
  geom_boxplot()+
  theme_classic()+xlab("")+ ylab("Sum of \nguilds size")+labs(color="",fill="")+ theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_fill_manual(values=c('#DECF3F','olivedrab3'))+ylim(0,200)+
  scale_color_manual(values=c('black','black'))+theme(legend.position="bottom",
                                                      axis.text.x = element_text(size=13),
                                                      legend.text = element_text(size=14),
                                                      axis.title.y = element_text(size=12))


panel_abc=ggarrange(panel_a,panel_b,panel_c,panel_d,nrow =4,labels=LETTERS[1:4])
p=ggarrange(panel_abc,panel_e,nrow =2,heights = c(1.85, 1),labels=c("",LETTERS[5]))
ggsave("./Figures/SI/Global_metric_empirical_N.pdf", p, height = 29.7/1.3, width = 17, units = "cm")


## >> Txample of Laplacian density ----

pdf("./Figures/SI/Example_Laplacian_density.pdf",width = 5,height = 5)

#A : example network


load("./Data/network_BD.RData")
network=network_bdd[[277]]
matrix_network = as.matrix(network)
A = rbind(cbind(matrix(
  0, nrow = nrow(network), ncol = nrow(network)
), matrix_network),
cbind(t(matrix_network), matrix(
  0, nrow = ncol(network), ncol = ncol(network)
)))
# degree matrix
D = matrix(0, nrow = nrow(A), ncol = ncol(A))
diag(D) = c(rowSums(A))
# Laplacian graph

MGL = D - A #laplacian matrix
sqrt_D = eigen(D)$vectors %*% diag(1 / sqrt(eigen(D)$values)) %*% t(eigen(D)$vectors) #used to normalised the laplacian matrix
nMGL = diag(nrow(D)) - sqrt_D %*% A %*% sqrt_D #normalized laplacian matrix

print(image(t(network_bdd[[277]][nrow(network_bdd[[277]]):1,]),col=c("white","forestgreen"),axes=F))
print(image(t(nMGL[nrow(nMGL):1,]),col=c(pal(5)[1:2],"white",pal(5)[3:5])))
spectre=spectR_network_abs(network,bandwith = .025)
print(plot_spectR_network(spectre,bw=.025))

# B: Modular network

modular_simulation=module_matrix(20,20,5)
image(modular_simulation,col = c("white","#DECF3F"),axes=FALSE)
n=NODFc(modular_simulation)
m=computeModules(modular_simulation)@likelihood
spectre=spectR_network_abs(modular_simulation,bandwith = 0.025)
plot_spectR_network(spectre,bw=0.025)
text(x=1.5,y=12,round(n,2))
text(x=1.5,y=10,round(m,2))


# C: Nested network

nested_simulation=nested_matrix(20,20,3)
image(nested_simulation,col = c("white","olivedrab3"),axes=FALSE)
n=NODFc(nested_simulation)
m=computeModules(nested_simulation)@likelihood
spectre=spectR_network_abs(nested_simulation,bandwith = 0.025)
plot_spectR_network(spectre,bw=0.025)
text(x=1.5,y=12,round(n,2))
text(x=1.5,y=10,round(m,2))


# D : Antagonistic parasitism network : network 22

antagonistic_empirical=network_bdd[[22]]
image(antagonistic_empirical,col = c("white","#DECF3F"),axes=FALSE)
n=NODFc(antagonistic_empirical)
m=computeModules(antagonistic_empirical)@likelihood
spectre=spectR_network_abs(antagonistic_empirical,bandwith = 0.025)
plot_spectR_network(spectre,bw=0.025)
text(x=1.5,y=12,round(n,2))
text(x=1.5,y=10,round(m,2))


# E: Mutualistic pollination network : network 97

mutualisitc_empirical=network_bdd[[97]]
image(mutualisitc_empirical,col = c("white","olivedrab3"),axes=FALSE)
n=NODFc(mutualisitc_empirical)
m=computeModules(mutualisitc_empirical)@likelihood
spectre=spectR_network_abs(mutualisitc_empirical,bandwith = 0.025)
plot_spectR_network(spectre,bw=0.025)
text(x=1.5,y=12,round(n,2))
text(x=1.5,y=10,round(m,2))

dev.off()







## >> Tifference in spectral density in empirical networks for the different types of ecology ----


data_laplacian=save[,c(2,6:105)]

eigen_name=NULL;freq=NULL;Type=NULL;Dat=NULL
data_laplacian=filter(data_laplacian,Type %in% c("Pollination","Host-Parasite","Seed dispersal","Herbivory"))

for (i in 2:ncol(data_laplacian)){
  data_mixed_model=data_laplacian[,c(1,i)]
  colnames(data_mixed_model)[2]="Eigen"
  if (Anova(lm(data=mutate(data_mixed_model,Type=as.factor(Type)),
               formula=Eigen~Type))$`Pr(>F)`[1] < 0.01){ 
    
    eigen_name=c(rep(seq(0,1,length.out=100)[i-1],length(unique(data_laplacian$Type)))) #motif number starts at 1 and index at 3
    freq=aggregate(x = data_laplacian[,i],by = list(data_laplacian$Type),FUN = mean)$x
    Nature=aggregate(x = data_laplacian[,i],by = list(data_laplacian$Type),FUN = mean)$Group.1
    Dat=rbind(Dat,cbind(as.numeric(eigen_name),as.numeric(freq),as.character(Nature)))
  }
}


Dat=data.frame(Eigen=round(as.numeric(Dat[,1]),3),Freq=as.numeric(Dat[,2]),Type=as.character(Dat[,3]))
colnames(data_laplacian)[2:ncol(data_laplacian)]=as.character(round(seq(0,1,length.out=100),3))
data_2=melt(as.matrix(data_laplacian[,2:ncol(data_laplacian)]))
data_2$value=sqrt(data_2$value)
data_2$Type=rep(data_laplacian$Type, ncol(data_laplacian)-1)
data_2=data_2[which(data_2$Var2 %in%  unique(Dat$Eigen)),]

colnames(data_2)=c("Network","Eigen","Density","Type")

#We split in two graphs
p1=ggplot(data_2[which(as.numeric(data_2$Eigen)<0.25),], aes(x=as.character(Eigen), y=Density, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#C57700","#DECF3F","#0CA206","olivedrab3"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")

p2=ggplot(data_2[which(as.numeric(data_2$Eigen)>0.25 & as.numeric(data_2$Eigen)<.5),], aes(x=as.character(Eigen), y=Density, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#C57700","#DECF3F","#0CA206","olivedrab3"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")


p3=ggplot(data_2[which(as.numeric(data_2$Eigen)>0.5 & as.numeric(data_2$Eigen)<.75),], aes(x=as.character(Eigen), y=Density, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#C57700","#DECF3F","#0CA206","olivedrab3"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")

p4=ggplot(data_2[which(as.numeric(data_2$Eigen)>0.75),], aes(x=as.character(Eigen), y=Density, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#C57700","#DECF3F","#0CA206","olivedrab3"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")

p=ggarrange(p1,p2,p3,p4,nrow=4)

ggsave("./Figures/SI/Spectral_density_type_ecology.pdf",p,width = 10,height = 16)


## >> Totifs difference ecology and type of interaction ----

count_motif=save[,c(1,2,108:151)]


for (type_test in c("A","B","C")){
  
  Motif_name=NULL;freq=NULL;Nature=NULL;Dat=NULL
  
  if (type_test=="A"){# Testing the dichotomy (A/M), such as in Fig-3B
    for (i in 3:ncol(count_motif)){
      data_mixed_model=count_motif[,c(1,2,i)]
      colnames(data_mixed_model)[3]="Motif"
      if (wilcox.test(data_mixed_model$Motif[which(data_mixed_model$Category=="Mutualistic")],
                      data_mixed_model$Motif[which(data_mixed_model$Category=="Antagonistic")])$p.value<0.01){ #we keep those that significantly differ using mixed effect models
        
        Motif_name=rbind(as.numeric(i-2),as.numeric(i-2)) #motif number starts at 1 and index at 3
        freq=rbind(aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][2,],
                   aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][1,])
        Nature=rbind("Mutualistic","Antagonistic")
        Dat=rbind(Dat,cbind(as.numeric(Motif_name),as.numeric(freq),as.character(Nature)))
      }
    }
  } else if (type_test=="B"){ #Testing the dichotomy when controlling for the type of interaction
    
    for (i in 3:ncol(count_motif)){
      data_mixed_model=count_motif[,c(1,2,i)]
      colnames(data_mixed_model)[3]="Motif"
      if (Anova(lmer(data=mutate(data_mixed_model,Type=as.factor(Type),Category=as.factor(Category)),
                     formula=Motif~Category+(1|Type)))$`Pr(>Chisq)` < 0.01){ #we keep those that significantly differ using mixed effect models
        
        Motif_name=rbind(as.numeric(i-2),as.numeric(i-2)) #motif number starts at 1 and index at 3
        freq=rbind(aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][2,],
                   aggregate(x = count_motif[,i],by = list(count_motif$Category),FUN = mean)[2][1,])
        Nature=rbind("Mutualistic","Antagonistic")
        Dat=rbind(Dat,cbind(as.numeric(Motif_name),as.numeric(freq),as.character(Nature)))
      }
    }
  }
  
  
  Dat=data.frame(Motif=as.numeric(Dat[,1]),Freq=as.numeric(Dat[,2]),Nature=as.character(Dat[,3]))
  
  data_2=melt(as.matrix(count_motif[,3:ncol(count_motif)]))
  data_2$value=sqrt(data_2$value)
  data_2$Category=rep(count_motif$Category, ncol(count_motif)-2)
  
  data_2=data_2[which(data_2$Var2 %in%  paste0("Motif_",unique(Dat$Motif))),]
  
  assign(paste0("panel_",type_test), ggplot(data_2, aes(x=Var2, y=value, fill=Category,col=Category)) + 
           geom_boxplot(outlier.size = 0) + theme_classic() + 
           theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
           xlab("") + ylab(expression(sqrt('Motif frequency')))+
           scale_fill_manual(values=c('#DECF3F','olivedrab3'))+
           scale_color_manual(values=c('black','black'))+
           theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color=""))
  
  
}

count_motif=save[,c(1,2,108:151)]
Motif_name=NULL;freq=NULL;Nature=NULL;Dat=NULL;Type=NULL
count_motif=filter(count_motif,Type %in% c("Pollination","Host-Parasite","Seed dispersal","Herbivory"))

for (i in 3:ncol(count_motif)){
  data_mixed_model=count_motif[,c(1,2,i)]
  colnames(data_mixed_model)[3]="Motif"
  if (Anova(lm(data=mutate(data_mixed_model,Type=as.factor(Type),Category=as.factor(Category)),
               formula=Motif~Type))$`Pr(>F)`[1] < 0.01){ 
    
    Motif_name=c(rep(i-2,length(unique(count_motif$Type)))) #motif number starts at 1 and index at 3
    freq=aggregate(x = count_motif[,i],by = list(count_motif$Type),FUN = mean)$x
    Nature=aggregate(x = count_motif[,i],by = list(count_motif$Type),FUN = mean)$Group.1
    Dat=rbind(Dat,cbind(as.numeric(Motif_name),as.numeric(freq),as.character(Nature)))
  }
}


Dat=data.frame(Motif=as.numeric(Dat[,1]),Freq=as.numeric(Dat[,2]),Type=as.character(Dat[,3]))
data_2=melt(as.matrix(count_motif[,3:ncol(count_motif)]))
data_2$value=sqrt(data_2$value)
data_2$Type=rep(count_motif$Type, ncol(count_motif)-2)
data_2=data_2[which(data_2$Var2 %in%  paste0("Motif_",unique(Dat$Motif))),]

panel_C= ggplot(data_2, aes(x=Var2, y=value, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c("#C57700","#DECF3F","#0CA206","olivedrab3"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")




pdf("./Figures/SI/Significant_motifs_btw_networks_AB.pdf",width = 8,height = 5)
print(panel_A)
print(panel_B)
dev.off()
pdf("./Figures/SI/Significant_motifs_btw_networks_C.pdf",width = 12,height = 5)
print(panel_C)
dev.off()


## >> TCA & Kmeans on empirical alone + (empirical + simulated) network ----

set.seed(123)

combination=list(c(1:3,4,5,106,107), # global metrics
                 c(1,2,6:105), #Laplacian spectral density
                 c(1,2,108:151)) # motifs

fig_number=c(4,5,6)
u=1
for (i in c("Global_metrics","Laplacian","Motifs")){ #Fig number
  set.seed(123)
  
  pdf(paste0("./Figures/SI/Fig_PCA_Simu_empirical_",i,".pdf"),height = 4,width = 7)
  
  data_BipartiteEvol=save_bipartite[sample(1:nrow(save_bipartite),size=300),combination[[u]]]
  data_empiric=save[,combination[[u]]]
  data_empiric$Category=as.character(data_empiric$Category)
  data_empiric$Category[which(data_empiric$Category=="Antagonistic")]="Empirical antagonistic"
  data_empiric$Category[which(data_empiric$Category=="Mutualistic")]="Empirical mutualistic"
  
  
  
  
  data_all=rbind(data_empiric,data_BipartiteEvol) #data_nested,data_modular,data_BipartiteEvol)
  
  
  res.pca=PCA(data_all[,3:(dim(data_all)[2])], scale.unit = TRUE, ncp = 2,  graph=F)
  
  print(fviz_pca_ind(res.pca, geom.ind = "point", 
                     axes=c(1,2), pointshape = 21,
                     pointsize = 1.5, 
                     col.ind=data_all$Category,
                     fill.ind=data_all$Category,
                     palette = c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"),
                     addEllipses = F,
                     alpha.ind = 1,
                     label = "var",
                     repel = T,mean.point = FALSE,
                     legend.title = "")+
          labs(title="")+
          theme_classic()+
          xlab(paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"))+
          ylab(paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"))+
          theme(plot.title = element_text(hjust = 0.2)))
  dev.off()
  
  set.seed(1)
  km.res=kmeans(res.pca$ind$coord, 4, nstart = 50) #before : 6 clusters with topological simulations
  Kmeans_table=table(km.res$cluster,data_all$Category) #table of CAH
  print(chisq.test(Kmeans_table,simulate.p.value = T))
  
  pdf(paste0("./Figures/SI/Fig_PCA_ecology_",i,".pdf"),height = 4,width = 7)
  
  res.pca=PCA(data_empiric[,3:(dim(data_empiric)[2])], scale.unit = TRUE, ncp = 2,  graph=F)
  
  print(fviz_pca_ind(res.pca, geom.ind = "point", 
                     axes=c(1,2), 
                     pointsize = 1.5, 
                     col.ind=data_empiric$Type,
                     fill.ind=data_empiric$Type,
                     addEllipses = F,
                     alpha.ind = 1,
                     label = "var",
                     repel = T,mean.point = FALSE,
                     legend.title = "Network ecology")+labs(title="")+
          theme_classic()+scale_shape_manual(values=0:8)+
          xlab(paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"))+
          ylab(paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"))+
          theme(plot.title = element_text(hjust = 0.2)))
  
  dev.off()
  
  set.seed(1)
  
  km.res=kmeans(res.pca$ind$coord, 2, nstart = 50)
  Kmeans_table=table(km.res$cluster,data_empiric$Type) #table of CAH
  
  chisq.test(Kmeans_table,simulate.p.value = T)
  
  classif_out_cluster=c()
  for (k in 1:ncol(Kmeans_table)){
    max_clust=which(Kmeans_table[,k]==max(Kmeans_table[,k]))
    classif_out_cluster=c(classif_out_cluster,sum(Kmeans_table[-max_clust,k])/sum(Kmeans_table[,k]))
  }
  
  Kmeans_table=rbind(Kmeans_table,round(100*classif_out_cluster,1))
  
  print(mean(classif_out_cluster*apply(Kmeans_table, 2, sum))) #moyenne pondérée totale
  print(mean(classif_out_cluster[c(1,6,7,8,9)]*apply(Kmeans_table, 2, sum)[c(1,6,7,8,9)])) #moyenne pondérée mutualistic
  print(mean(classif_out_cluster[c(2,3,4,5)]*apply(Kmeans_table, 2, sum)[c(2:5)])) #moyenne pondérée antagonistic
  
  write.table(Kmeans_table,paste0("./Figures/SI/Table_PCA_kmeans_ecology_",i,".csv"),sep=";")
  
  
  u=u+1
}


## >> Testing influence of overfitting on our results ----




data_BD=save[,c(1,2,108:151)]

size_layer=c(0,1,5,9,13,17,21,25,29,33,37,40)
table_results=data.frame(AA_test=0,MM_test=0,AA_train=0,MM_train=1:length(size_layer))
count=1
for (k in size_layer){
  g=1
  aa=NULL
  mm=NULL
  AA=NULL
  MM=NULL
  for (x in 1:50){
    data.test=NULL
    data.train=NULL
    set.seed(sample(1:10000,1))
    sample_test = sample(1:nrow(data_BD),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(data_BD),
                         replace = F)
    data.test=data_BD[sample_test,]
    data.train=data_BD[-sample_test,]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = k , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2])  
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    
    prediction_nn=predict(nn_1layer, data.train[, 3:length(data.train)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.train$Category, prediction_nn_test[,2])  
    
    
    AA[g]=results[1]/(results[3]+results[1])*100
    MM[g]=results[4]/(results[2]+results[4])*100
    
    
    g=g+1
  }
  table_results$AA_test[count]=mean(aa)
  table_results$AA_train[count]=mean(AA)
  table_results$MM_test[count]=mean(mm)
  table_results$MM_train[count]=mean(MM)
  count=count+1
}
table_results=cbind(size_layer,table_results)

# figure of overfitting
figure_data=(190*table_results$MM_train+153*table_results$AA_train)/343  #number of networks
figure_data=cbind(figure_data,table_results$MM_test,table_results$AA_test)
figure_data=cbind(figure_data,size_layer)
colnames(figure_data)=c("Classification network training set","Classification mutualistic test set","Classification antagonistic test set","Hidden layer size")
figure_data=as.data.frame(figure_data)
figure_data=melt(data=figure_data,id.vars=4)

p=ggplot(figure_data)+
  geom_line(aes(x=`Hidden layer size`,y=value,col=variable))+ylim(c(0,100))+
  geom_point(aes(x=`Hidden layer size`,y=value))+theme_minimal()+
  scale_color_manual(values=c("black",'#8fd175','#DECF3F'))+labs(y="",color="")

ggsave(p,filename = "./Figures/SI/Overfitting.pdf",width = 9,height = 4)




## >> Tontrolling for the number of ecologies (A) and the asymetry antagonist/mutualist in the training set ----


#Part A : controlling for the number of ecology in the training set

data_empiric=save[,c(1,2,108:151)]

size_training=c()
MM_ecology_controled=NULL;MM=NULL
AA_ecology_controled=NULL;AA=NULL
p=1
N_max_training=NULL
for (n_max in seq(5,20,3)){
  
  aa_ecology_controled=NULL;aa=NULL
  mm_ecology_controled=NULL;mm=NULL
  g=1
  
  for (b in sample(1:10000,50)){
    set.seed(b)
    
    
    #We first create training set with the same number of network per type of ecology
    
    sample_train=c()
    for (j in unique(data_empiric$Type)){
      
      data_j=data_empiric[which(data_empiric$Type==j),]
      if(nrow(data_j)<n_max){
        sample_ecology= rownames(data_j)
      } else{
        sample_ecology=sample(rownames(data_j),n_max)#maximum n_max networks from the same type of ecology in the training set
      }
      sample_train=c(sample_train,sample_ecology)
    }
    
    sample_train=as.numeric(sample_train)
    
    data.test=data_empiric[-sample_train,]
    data.train=data_empiric[sample_train,]
    size_training=c(size_training,nrow(data.train))
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) 
    
    
    aa_ecology_controled[g]=100*(results[1]/(results[1]+results[3]))
    mm_ecology_controled[g]=100*(results[4]/(results[2]+results[4]))
    
    
    # We compare the latter classification with the normal one, when we don't control for the 
    # number of networks in each ecology
    
    sample_train=sample(1:343,length(sample_train)) #same length
    
    data.test=data_empiric[-sample_train,]
    data.train=data_empiric[sample_train,]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) 
    
    
    
    aa[g]=100*results[1]/(results[1]+results[3])
    mm[g]=100*results[4]/(results[2]+results[4])
    
    g=g+1
    
  }
  AA_ecology_controled[p]=mean(aa_ecology_controled)
  MM_ecology_controled[p]=mean(mm_ecology_controled)
  AA[p]=mean(aa)
  MM[p]=mean(mm)
  N_max_training[p]=n_max
  p=p+1
}

results_NN=tibble(Antagonist_controlled=AA_ecology_controled,Mutualist_controlled=MM_ecology_controled,N_max_training=N_max_training,
                  Antagonist=AA,Mutualist=MM)
results_NN=melt(results_NN,id.vars = 3)


pdf("./Figures/SI/Controlling_number_ecology.pdf",width = 7,height = 4)
print(ggplot(results_NN)+geom_point(aes(x=N_max_training,y=value,color=variable),size=2)+theme_classic()+
        geom_line(aes(x=N_max_training,y=value,color=variable,linetype=variable),lwd=1)+
        labs(x="Maximal number of network per type of ecology",y="Classification accuracy",color="",linetype="")+ylim(0,100)+
        scale_linetype_manual(values = c(2,2,1,1))+ theme(legend.position="bottom",panel.grid.major=element_line(color = "gray75",size=.5))+
        scale_color_manual(values=c("#DECF3F","olivedrab3","#DECF3F","olivedrab3"))+
        annotate(geom = "text",x=seq(5,20,3),y=90,label=paste0(( unique(size_training))))
      
      
      
)
dev.off()

write.table(results_NN,"./Figures/SI/Controlling_number_ecology.csv",sep=";")



#Part B : controlling for the asymetry antagonist/mutualist in the training set

data_empiric=save[,c(1,2,108:151)]


MM_ecology_controled=NULL;MM=NULL
AA_ecology_controled=NULL;AA=NULL
p=1
N_max_training=NULL
for (n_max in seq(5,100,15)){
  
  aa_ecology_controled=NULL;aa=NULL
  mm_ecology_controled=NULL;mm=NULL
  g=1
  
  for (b in sample(1:10000,50)){
    set.seed(b)
    
    
    #We first create training set with the same number of network per type of interaction (n_max for each)
    
    sample_train=c(sample(rownames(data_empiric[which(data_empiric$Category=="Mutualistic network"),]),n_max),
                   sample(rownames(data_empiric[which(data_empiric$Category=="Antagonistic network"),]),n_max))
    
    sample_train=as.numeric(sample_train)
    
    data.test=data_empiric[-sample_train,]
    data.train=data_empiric[sample_train,]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2]) 
    
    
    
    #aa_ecology_controled[g]=F1_Score(data.test$Category,prediction_nn_test[,2],positive = "Antagonistic")
    #mm_ecology_controled[g]=F1_Score(data.test$Category,prediction_nn_test[,2],positive = "Mutualistic")
    aa_ecology_controled[g]=100*results[1]/(results[1]+results[3])
    mm_ecology_controled[g]=100*results[4]/(results[2]+results[4])
    
    
    # We compare the latter classification with the normal one, when we don't control for the 
    # number of networks in each ecology
    
    sample_train=sample(1:343,length(sample_train)) #same length
    
    data.test=data_empiric[-sample_train,]
    data.train=data_empiric[sample_train,]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2])  
    
    
    
    #aa[g]=F1_Score(data.test$Category,prediction_nn_test[,2],positive = "Antagonistic")
    #mm[g]=F1_Score(data.test$Category,prediction_nn_test[,2],positive = "Mutualistic")
    aa[g]=100*results[1]/(results[1]+results[3])
    mm[g]=100*results[4]/(results[2]+results[4])
    
    g=g+1
    
  }
  AA_ecology_controled[p]=mean(aa_ecology_controled)
  MM_ecology_controled[p]=mean(mm_ecology_controled)
  AA[p]=mean(aa)
  MM[p]=mean(mm)
  N_max_training[p]=n_max
  p=p+1
}

results_NN=tibble(Antagonist_controlled=AA_ecology_controled,Mutualist_controlled=MM_ecology_controled,N_max_training=N_max_training,
                  Antagonist=AA,Mutualist=MM)
results_NN=melt(results_NN,id.vars = 3)


pdf("./Figures/SI/Asymetry_mutua_anta.pdf",width = 7,height = 4)
print(ggplot(results_NN)+geom_point(aes(x=N_max_training,y=value,color=variable),size=2)+theme_classic()+
        geom_line(aes(x=N_max_training,y=value,color=variable,linetype=variable),lwd=1)+
        labs(x="Number of networks in the training set",y="Classification accuracy",color="",linetype="")+ylim(0,100)+
        scale_linetype_manual(values = c(2,2,1,1))+ theme(legend.position="bottom",panel.grid.major=element_line(color = "gray75",size=.5))+
        scale_color_manual(values=c("#DECF3F","olivedrab3","#DECF3F","olivedrab3"))+
        scale_x_continuous(breaks= c(25,50,75),labels = c("50","100","150"))
)
dev.off()

write.table(results_NN,"./Figures/SI/Asymetry_mutua_anta_data.csv",sep=";")




## >> Tlobal metrics of empirical & simulated networks  ----


data_empiric=save[,c(1,3:5,106,107)]
col=c()
for (k in data_empiric$Category){
  if (k=="Mutualistic") col=c(col,'#8fd175')
  else col=c(col,'#DECF3F')
}

data_empiric = cbind("Empirical networks",data_empiric$Category,data_empiric[,2:6],col)
colnames(data_empiric)=c("Category","Nature","Diff","Sum","Connectance","Nestedness","Modularity","Color")


#BipartiteEVol networks
Nature=c();color=c()
for (k in save_bipartite$Category){
  if (k=="BipartiteEvol mutualistic") {Nature=c(Nature,'Mutualistic')
  color=c(color,"#bb8fce")}
  else {Nature=c(Nature,'Antagonistic') 
  color=c(color,"#5B2C6F")}
}
data_BipartiteEvol=data.frame("Category"="BipartiteEvol networks","Nature"=Nature,
                              "Diff"=save_bipartite[,3],"Sum"=save_bipartite[,4],
                              "Connectance"=save_bipartite[,5],"Nestedness"=save_bipartite[,106],
                              "Modularity"=save_bipartite[,107],"Color"=color) %>% filter(Nestedness<12) #we filter for graphical reason, so that other plots don't get flaten


#Merging the data and plotting the results
info_data=rbind(data_empiric,data_BipartiteEvol) #,data_nested,data_modular)

panel_a=qplot(Nature, Connectance, data = info_data, fill=Color,alpha=.8)+geom_violin(trim=T)+theme_classic()+geom_boxplot(alpha=1,width=0.1,outlier.shape = NA,color="black")+
  facet_wrap( ~ Category, ncol = 3) +  labs(x="",y="Connectance")+
  scale_fill_manual(values=c("#5B2C6F",'olivedrab3','#bb8fce','#DECF3F'))+
  theme(strip.text.x = element_text(colour = "black", face = "italic"),legend.position="none",axis.text.x = element_blank(),
        axis.ticks = element_blank())


panel_b=qplot(Nature, Nestedness, data = info_data, fill=Color,alpha=.8)+geom_violin(trim=T)+theme_classic()+geom_boxplot(alpha=1,width=0.1,outlier.shape = NA,color="black")+
  facet_wrap( ~ Category, ncol = 3) +  labs(x="",y="Nestedness")+
  scale_fill_manual(values=c("#5B2C6F",'olivedrab3','#bb8fce','#DECF3F'))+
  theme(strip.text.x = element_blank(),legend.position="none",axis.text.x = element_blank(),
        axis.ticks = element_blank(),strip.background = element_blank())


panel_c=qplot(Nature, Modularity, data = info_data, fill=Color,alpha=.8)+geom_violin(trim=T)+theme_classic()+geom_boxplot(alpha=1,width=0.1,outlier.shape = NA,color="black")+
  facet_wrap( ~ Category, ncol = 3) +  labs(x="",y="Modularity")+
  scale_fill_manual(values=c("#5B2C6F",'olivedrab3','#bb8fce','#DECF3F'))+
  theme(strip.text.x = element_blank(),legend.position="none",axis.text.x = element_blank(),
        axis.ticks = element_blank(),strip.background = element_blank())


panel_d=qplot(Nature, Diff, data = info_data, fill=Color,alpha=.8)+geom_violin(trim=T)+theme_classic()+geom_boxplot(alpha=1,width=0.1,outlier.shape = NA,color="black")+
  facet_wrap( ~ Category, ncol = 3) +  labs(x="",y="Difference of guilds size")+
  scale_fill_manual(values=c("#5B2C6F",'olivedrab3','#bb8fce','#DECF3F'))+
  theme(strip.text.x = element_blank(),legend.position="none",axis.text.x = element_blank(),
        axis.ticks = element_blank(),strip.background = element_blank())

panel_e=qplot(Nature, Sum, data = info_data, fill=Color,alpha=.8)+geom_violin(trim=T)+theme_classic()+geom_boxplot(alpha=1,width=0.1,outlier.shape = NA,color="black")+
  facet_wrap( ~ Category, ncol = 3) +  labs(x="",y="Sum of guilds size")+
  scale_fill_manual(values=c("#5B2C6F",'olivedrab3','#bb8fce','#DECF3F'))+
  theme(strip.text.x = element_blank(),legend.position="none",
        strip.background = element_blank())


panel_abcd=ggarrange(panel_a,panel_b,panel_c,panel_d,nrow =4,labels=LETTERS[1:4])
p=ggarrange(panel_abcd,panel_e,nrow =2,heights = c(1, .285),labels=c("",LETTERS[5]))

ggsave("./Figures/SI/Global_metrics_simu_empirical.pdf", p, height = 27.7, width = 19, units = "cm")


#Test associated with the figures: 
#Comparison for nestedness, connectance and modularity

wilcox.test(save$Connectance,save_bipartite$Connectance);mean(save_bipartite$Connectance);mean(save$Connectance)
wilcox.test(save$Nestedness,save_bipartite$Nestedness);mean(save_bipartite$Nestedness);mean(save$Nestedness)
wilcox.test(save$Modularity,save_bipartite$Modularity);mean(save_bipartite$Modularity);mean(save$Modularity)
wilcox.test(save$Sum,save_bipartite$Sum);mean(save_bipartite$Sum);mean(save$Sum)
wilcox.test(save$Diff,save_bipartite$Diff);mean(save_bipartite$Diff);mean(save$Diff)

#******************************************************************************#







#******************************************************************************#

## >> Totifs frequencies between simulations and empirical networks ----


data_empiric=save[,c(1,2,108:151)]
data_empiric$Category=as.character(data_empiric$Category)
data_empiric$Category[which(data_empiric$Category=="Antagonistic")]="Empirical antagonistic"
data_empiric$Category[which(data_empiric$Category=="Mutualistic")]="Empirical mutualistic"
Motif_data=rbind(data_empiric,save_bipartite[,c(1,2,108:151)]) #,save_topo[,c(1,2,108:151)])
colnames(Motif_data)[1]="Category"
Motif_data[,c(3:46)]=sqrt(Motif_data[,c(3:46)]);Motif_data$Simulation=c(rep("Empirical",343),rep("Simulated",1175))
Motif_data=Motif_data[,c(1,2,47,3:46)]
save_motif_data=Motif_data#saving the file to reuse it


Motif_name=NULL;freq=NULL;Nature=NULL;Dat=NULL;Type=NULL

for (i in 4:ncol(Motif_data)){
  data_mixed_model=Motif_data[,c(1,2,3,i)]
  colnames(data_mixed_model)[4]="Motif"
  if (Anova(lm(data=mutate(data_mixed_model,Simulation=as.factor(Simulation),Category=as.factor(Category)),
               formula=Motif~Simulation))$`Pr(>F)`[1] < 0.01){ 
    
    Motif_name=c(rep(i-2,length(unique(Motif_data$Category)))) #motif number starts at 1 and index at 3
    freq=aggregate(x = Motif_data[,i],by = list(Motif_data$Category),FUN = mean)$x
    Nature=aggregate(x = Motif_data[,i],by = list(Motif_data$Category),FUN = mean)$Group.1
    Dat=rbind(Dat,cbind(as.numeric(Motif_name),as.numeric(freq),as.character(Nature)))
  }
}


Dat=data.frame(Motif=as.numeric(Dat[,1]),Freq=as.numeric(Dat[,2]),Type=as.character(Dat[,3]))
data_2=melt(as.matrix(Motif_data[,4:ncol(Motif_data)]))
data_2$value=sqrt(data_2$value)
data_2$Type=rep(Motif_data$Category, ncol(Motif_data)-3)
data_2=data_2[which(data_2$Var2 %in%  paste0("Motif_",unique(Dat$Motif))),]

p1= ggplot(data_2[which(data_2$Var2 %in% as.factor(paste0("Motif_",1:10))),], aes(x=Var2, y=value, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p2= ggplot(data_2[which(data_2$Var2 %in% as.factor(paste0("Motif_",11:20))),], aes(x=Var2, y=value, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p3= ggplot(data_2[which(data_2$Var2 %in% as.factor(paste0("Motif_",21:30))),], aes(x=Var2, y=value, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p4= ggplot(data_2[which(data_2$Var2 %in% as.factor(paste0("Motif_",30:44))),], aes(x=Var2, y=value, fill=Type)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Motif frequency')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")


p=ggarrange(p1,p2,p3,p4,nrow = 4, labels = LETTERS[1:4])
ggsave("./Figures/SI/Motifs_frequency_empiric_simulations.pdf",p,width = 12,height = 14)


## >> Tifference in spectral density between simulated and empirical networks ----

data_empiric=save[,c(1,2,6:105)]
data_empiric$Category=as.character(data_empiric$Category)
data_empiric$Category[which(data_empiric$Category=="Antagonistic")]="Empirical antagonistic"
data_empiric$Category[which(data_empiric$Category=="Mutualistic")]="Empirical mutualistic"
Spectral_data=rbind(data_empiric,save_bipartite[,c(1,2,6:105)])
colnames(Spectral_data)[1]="Category";Spectral_data$Simulation=c(rep("Empirical",343),rep('Simulation',1175))
Spectral_data=Spectral_data[,c(1,2,103,3:102)]
save_spectral_data=Spectral_data#saving the file to reuse it

colnames(Spectral_data)[4:ncol(Spectral_data)]=as.character(round(seq(0,1,length.out=100),3))

eigen_name=NULL;freq=NULL;Type=NULL;Dat=NULL

for (i in 4:ncol(Spectral_data)){
  data_mixed_model=Spectral_data[,c(1,2,3,i)]
  colnames(data_mixed_model)[4]="Eigen"
  if (Anova(lm(data=mutate(data_mixed_model,Category=as.factor(Category),Simulation=as.factor(Simulation)),
               formula=Eigen~Simulation))$`Pr(>F)`[1] < 0.01){ 
    
    eigen_name=c(rep(seq(0,1,length.out=100)[i-3],length(unique(Spectral_data$Category)))) #motif number starts at 1 and index at 3
    freq=aggregate(x = Spectral_data[,i],by = list(Spectral_data$Category),FUN = mean)$x
    Nature=aggregate(x = Spectral_data[,i],by = list(Spectral_data$Category),FUN = mean)$Group.1
    Dat=rbind(Dat,cbind(as.numeric(eigen_name),as.numeric(freq),as.character(Nature)))
  }
}


Dat=data.frame(Eigen=round(as.numeric(Dat[,1]),3),Freq=as.numeric(Dat[,2]),Type=as.character(Dat[,3]))
colnames(Spectral_data)[4:ncol(Spectral_data)]=as.character(round(seq(0,1,length.out=100),3))
data_2=melt(as.matrix(Spectral_data[,4:ncol(Spectral_data)]))
data_2$value=sqrt(data_2$value)
data_2$Type=rep(Spectral_data$Category, ncol(Spectral_data)-3)
data_2=data_2[which(data_2$Var2 %in%  unique(Dat$Eigen)),]

colnames(data_2)=c("Network","Eigen","Density","Category")

#We split in two graphs

p1=ggplot(filter(data_2,as.numeric(as.character(Eigen))<0.25), aes(x=as.factor(Eigen), y=Density, fill=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p2=ggplot(filter(data_2,as.numeric(as.character(Eigen))>0.25 & as.numeric(as.character(Eigen))<0.5), aes(x=as.factor(Eigen), y=Density, fill=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p3=ggplot(filter(data_2,as.numeric(as.character(Eigen))<0.75 & as.numeric(as.character(Eigen))>0.5), aes(x=as.factor(Eigen), y=Density, fill=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "none",legend.text = element_text(size=14))+labs(fill="",color="")

p4=ggplot(filter(data_2,as.numeric(as.character(Eigen))>0.75), aes(x=as.factor(Eigen), y=Density, fill=Category)) + 
  geom_boxplot(outlier.size = 0) + theme_classic() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,size=10)) +
  xlab("") + ylab(expression(sqrt('Laplacian density')))+
  scale_fill_manual(values=c("#5b2c6f","#bb8fce","#DECF3F","olivedrab3","#1a5276","#5dade2"))+
  theme(legend.position = "bottom",legend.text = element_text(size=14))+labs(fill="",color="")

p=ggarrange(p1,p2,p3,p4,nrow=4)

ggsave("./Figures/SI/Spectral_density_simulations_empirical.pdf",p,width = 10,height = 16)



## >> Table: AIC mixed effect models type of ecology ----


# For connectance
mod_connectance_ecology=lm(data=data_empiric,
                           formula=Connectance~Type)
mod_connectance_dichotomy_corrected=lmer(data=data_empiric,
                                         formula=Connectance~Category+(1|Type))
mod_connectance_dichotomy=lm(data=data_empiric,
                             formula=Connectance~Category)

test_c=AIC(mod_connectance_ecology,mod_connectance_dichotomy_corrected,mod_connectance_dichotomy)$AIC


#For nestedness
mod_nestedness_ecology=lm(data=data_empiric,
                          formula=Nestedness~Type)
mod_nestedness_dichotomy_corrected=lmer(data=data_empiric,
                                        formula=Nestedness~Category+(1|Type))
mod_nestedness_dichotomy=lm(data=data_empiric,
                            formula=Nestedness~Category)

test_n=AIC(mod_nestedness_ecology,mod_nestedness_dichotomy_corrected,mod_nestedness_dichotomy)$AIC


#For modularity
mod_modularity_ecology=lm(data=data_empiric,
                          formula=Modularity~Type)
mod_modularity_dichotomy_corrected=lmer(data=data_empiric,
                                        formula=Modularity~Category+(1|Type))
mod_modularity_dichotomy=lm(data=data_empiric,
                            formula=Modularity~Category)

test_m=AIC(mod_modularity_ecology,mod_modularity_dichotomy_corrected,mod_modularity_dichotomy)$AIC


#For Difference in guild sizes
mod_diff_ecology=lm(data=data_empiric,
                    formula=Diff~Type)
mod_diff_dichotomy_corrected=lmer(data=data_empiric,
                                  formula=Diff~Category+(1|Type))
mod_diff_dichotomy=lm(data=data_empiric,
                      formula=Diff~Category)

test_d=AIC(mod_diff_ecology,mod_diff_dichotomy_corrected,mod_diff_dichotomy)$AIC

#For sum of guild sizes
mod_sum_ecology=lm(data=data_empiric,
                   formula=Sum~Type)
mod_sum_dichotomy_corrected=lmer(data=data_empiric,
                                 formula=Sum~Category+(1|Type))
mod_sum_dichotomy=lm(data=data_empiric,
                     formula=Sum~Category)

test_s=AIC(mod_sum_ecology,mod_sum_dichotomy_corrected,mod_sum_dichotomy)$AIC

table_=data.frame(Dichotomy_corrected=c(test_n[2],test_c[2],test_m[2],test_d[2],test_s[2]),
                    Interaction_type=c(test_n[1],test_c[1],test_m[1],test_d[1],test_s[1]),
                    Dichotomy_uncorrected=c(test_n[3],test_c[3],test_m[3],test_d[3],test_s[3]));rownames(Table_S7)=c("Nestedness","Connectance","Modularity","Diff","Sum")


write.table(table_,"./Figures/SI/AIC_global_metrics.csv",sep=";")



## >> Methods for classifying empirical network ----




# >>> Probit classification

MM=NULL ; AA=NULL;F1_score=NULL
p=1
method=NULL
for (m in combination){
  data_BD=save[,m]
  data_BD$Category=as.factor(data_BD$Category)
  aa=NULL ;  mm=NULL;F1=NULL
  g=1
  for (b in sample(1:10000,50)){
    
    set.seed(b)
    
    sample_train = sample(1:nrow(data_BD), 4*nrow(data_BD)/5)
    train=data_BD[sample_train,]
    test=data_BD[-sample_train,]
    
    formule = as.formula(c('Category ~ ',
                           paste0(colnames(data_BD[, 3:length(data_BD)]),
                                  collapse = "+"))) 
    
    res.glm=glm(formula=formule,family =binomial,data=train)
    
    
    pred=predict.glm(res.glm,newdata=test,type="response")
    pred=pred > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in pred){
      
      if (i=="FALSE"){pred[u]="Antagonistic"}
      if( i=="TRUE"){pred[u]="Mutualistic"}
      u=u+1
    }
    results=table(test$Category,pred)
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    F1[g]=F1_Score(test$Category,pred)
    g=g+1
  }
  AA[p]=paste0(round(mean(aa),0)," +/- ",round(sd(aa),1))
  MM[p]=paste0(round(mean(mm),0)," +/- ", round(sd(mm),1))
  F1_score[p]=mean(F1)
  method[p]=name_combination[[p]]
  p=p+1
}

results_probit=data.frame(Antagonist=AA,Mutualist=MM,F1=F1_score,method=method)
write.table(results_probit,"./Figures/SI/results_probit.csv",sep=";")


# >>> Lasso classification

MM=NULL ;AA=NULL;F1_score=NULL
p=1
method=NULL
for (m in combination){
  data_BD=save[,m]
  data_BD$Category=as.factor(data_BD$Category)
  
  aa=NULL;F1=NULL
  mm=NULL
  g=1
  for (b in sample(1:10000,50)){
    
    set.seed(b)
    
    formule = as.formula(c('Category ~ ',
                           paste0(colnames(data_BD[, 3:length(data_BD)]),
                                  collapse = "+"))) 
    
    
    mod_lasso=model.matrix(formule, data_BD)
    category=data_BD$Category
    lambda_seq=10^seq(2, -2, by = -.1)
    
    # Splitting the data into test and train. Only using matrix and no data frame
    train = sample(1:nrow(mod_lasso), 4*nrow(mod_lasso)/5)
    
    x_train=mod_lasso[train,]
    y_train=category[train]
    
    cv_output=cv.glmnet(x_train, y_train,
                           alpha = 1, lambda = lambda_seq, 
                           nfolds = 5,family="binomial")
    
    # identifying best lambda
    best_lam=cv_output$lambda.min
    
    
    
    x_test = mod_lasso[-train,]
    y_test = category[-train]
    
    # Rebuilding the model with best lambda value identified
    lasso_best=glmnet(mod_lasso[train,], category[train], 
                         alpha = 1,
                         lambda = best_lam,family="binomial")
    pred=predict(lasso_best, s = best_lam, newx = x_test)
    pred=pred>0
    pred[which(pred==F)]="Antagonistic"
    pred[which(pred==T)]="Mutualistic"
    
    results=table(y_test, pred)
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    F1[g]=F1_Score(y_test,pred)
    
    g=g+1
  }
  AA[p]=paste0(round(mean(aa),0)," +/- ",round(sd(aa),1))
  MM[p]=paste0(round(mean(mm),0)," +/- ", round(sd(mm),1))
  F1_score[p]=mean(F1)
  method[p]=name_combination[[p]]
  p=p+1
}

results_lasso=data.frame(Antagonist=AA,Mutualist=MM,F1=F1_score,method=method)
write.table(results_lasso,"./Figures/SI/results_lasso.csv",sep=";")



# >>> Neural network classification 

size_layer=c(5,20,60,65,25,80,85) # the size of the intermediate layer depends on the data

F1_score=NULL
MM=NULL
AA=NULL
p=1
method=NULL
for (m in combination[-1]){
  data_BD=save[,m]
  data_BD$Category=as.factor(data_BD$Category)
  
  aa=NULL;F1=NULL
  mm=NULL
  g=1
  for (b in sample(1:10000,50)){
    
    set.seed(b)
    sample_test = sample(1:nrow(data_BD),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(data_BD),
                         replace = F)
    data.test=data_BD[sample_test,]
    data.train=data_BD[-sample_test,]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = size_layer[p] , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2])  
    
    
    
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    F1[g]=F1_Score(data.test$Category,prediction_nn_test[,2])
    g=g+1
  }
  AA[p]=paste0(round(mean(aa),0)," +/- ",round(sd(aa),1))
  MM[p]=paste0(round(mean(mm),0)," +/- ", round(sd(mm),1))
  F1_score[p]=mean(F1)
  method[p]=name_combination[[p]]
  p=p+1
}

results_NN=data.frame(Antagonist=AA,Mutualist=MM,F1=F1_score,method=method)
write.table(results_NN,"./Figures/SI/results_NN.csv",sep=";")




## >> Classifying empirical networks using simulations ----


# >>> Probit classification

for (v in c("BipartiteEvol_simulations")){
  MM_emp=NULL
  AA_emp=NULL
  MM_topo=NULL
  AA_topo=NULL
  p=1
  method=NULL
  for (m in combination){
    data_topological = save_bipartite
    data_BD=save[,m]
    data_topological=data_topological[,m];data_topological$Category=as.factor(data_topological$Category)
    
    aa_emp=NULL
    mm_emp=NULL
    aa_topo=NULL
    mm_topo=NULL
    g=1
    for (b in sample(1:10000,100)){
      
      set.seed(b)
      
      sample_train = sample(1:nrow(data_topological), 4*nrow(data_topological)/5)
      train=data_topological[sample_train,]
      test=data_topological[-sample_train,]
      
      formule = as.formula(c('Category ~ ',
                             paste0(colnames(data_topological[, 3:length(data_topological)]),
                                    collapse = "+"))) 
      
      res.glm=glm(formula=formule,family =binomial,data=train)
      
      
      pred=predict.glm(res.glm,newdata=test,type="response")
      pred=pred > 0.5# the value of 0.5 is the threshold of the network to classify them
      
      u=1
      for (i in pred){
        if (i=="FALSE"){pred[u]="Antagonistic"}
        if( i=="TRUE"){pred[u]="Mutualistic"}
        u=u+1
      }
      results=table(test$Category,pred)
      
      
      aa_topo[g]=results[1]/(results[3]+results[1])*100
      mm_topo[g]=results[4]/(results[2]+results[4])*100
      
      #and predict on empirical data
      pred=predict.glm(res.glm,newdata=data_BD,type="response")
      pred=pred > 0.5# the value of 0.5 is the threshold of the network to classify them
      
      u=1
      for (i in pred){
        if (i=="FALSE"){pred[u]="Antagonistic"}
        if( i=="TRUE"){pred[u]="Mutualistic"}
        u=u+1
      }
      results=table(data_BD$Category,pred)
      
      
      aa_emp[g]=results[1]/(results[3]+results[1])*100
      mm_emp[g]=results[4]/(results[2]+results[4])*100
      f[g]=
        g=g+1
    }
    AA_emp[p]=paste0(round(mean(aa_emp),0)," +/- ",round(sd(aa_emp),1))
    MM_emp[p]=paste0(round(mean(mm_emp),0)," +/- ", round(sd(mm_emp),1))
    AA_topo[p]=paste0(round(mean(aa_topo),0)," +/- ", round(sd(aa_topo),1))
    MM_topo[p]=paste0(round(mean(mm_topo),0)," +/- ", round(sd(mm_topo),1))
    method[p]=name_combination[[p]]
    p=p+1
  }
  
  results_probit=data.frame(Antagonist_topo=AA_topo,Mutualist_topo=MM_topo,Antagonist_emp=AA_emp,Mutualist_emp=MM_emp,method=method)
  print(results_probit)
  write.table(results_probit,
              paste0("./Figures/SI/results_probit_",v,".csv"),sep=";")
  
}



# >>> Lasso classification

for (v in c("BipartiteEvol_simulations")){
  MM_emp=NULL
  AA_emp=NULL
  MM_topo=NULL
  AA_topo=NULL
  p=1
  method=NULL
  for (m in combination){
    data_BD=save
    
    data_topological=save_bipartite
    data_BD=data_BD[,m]
    data_topological=data_topological[,m]
    
    aa_emp=NULL
    mm_emp=NULL
    aa_topo=NULL
    mm_topo=NULL
    g=1
    for (b in sample(1:10000,100)){
      
      set.seed(b)
      formule = as.formula(c('Category ~ ',
                             paste0(colnames(data_topological[, 3:length(data_topological)]),
                                    collapse = "+"))) 
      
      
      mod_lasso=model.matrix(formule, data_topological)
      category=data_topological$Category
      lambda_seq=10^seq(2, -2, by = -.1)
      
      # Splitting the data into test and train. Only using matrix and no data frame
      train = sample(1:nrow(mod_lasso), 4*nrow(mod_lasso)/5)
      
      x_train=mod_lasso[train,]
      y_train=category[train]
      
      cv_output=cv.glmnet(x_train, y_train,
                             alpha = 1, lambda = lambda_seq, 
                             nfolds = 5,family="binomial")
      
      # identifying best lambda
      best_lam=cv_output$lambda.min
      
      
      
      x_test = mod_lasso[-train,]
      y_test = category[-train]
      
      # Rebuilding the model with best lambda value identified
      lasso_best=glmnet(mod_lasso[train,], category[train], 
                           alpha = 1,
                           lambda = best_lam,family="binomial")
      pred=predict(lasso_best, s = best_lam, newx = x_test)
      pred=pred>0
      pred[which(pred==F)]="Antagonistic"
      pred[which(pred==T)]="Mutualistic"
      
      results=table(y_test, pred)
      
      
      aa_topo[g]=results[1]/(results[3]+results[1])*100
      mm_topo[g]=results[4]/(results[2]+results[4])*100
      
      #and predict on empirical data
      mod_lasso=model.matrix(formule, data_BD)
      
      pred=predict(lasso_best, s = best_lam, newx = mod_lasso)
      pred=pred>0
      pred[which(pred==F)]="Antagonistic"
      pred[which(pred==T)]="Mutualistic"
      
      results=table(data_BD$Category, pred)
      
      aa_emp[g]=results[1]/(results[3]+results[1])*100
      mm_emp[g]=results[4]/(results[2]+results[4])*100
      
      
      
      
      g=g+1
    }
    AA_emp[p]=paste0(round(mean(aa_emp),0)," +/- ",round(sd(aa_emp),1))
    MM_emp[p]=paste0(round(mean(mm_emp),0)," +/- ", round(sd(mm_emp),1))
    AA_topo[p]=paste0(round(mean(aa_topo),0)," +/- ", round(sd(aa_topo),1))
    MM_topo[p]=paste0(round(mean(mm_topo),0)," +/- ", round(sd(mm_topo),1))
    method[p]=name_combination[[p]]
    p=p+1
  }
  
  results_lasso=data.frame(Antagonist_topo=AA_topo,Mutualist_topo=MM_topo,Antagonist_emp=AA_emp,Mutualist_emp=MM_emp,method=method)
  write.table(results_lasso,
              paste0("./Figures/SI/results_lasso_",v,".csv"),sep=";")
  
}




# >>> Neural network classification



size_layer=c(5,20,60,65,25,80,85) # the size of the intermediate layer depends on the data

for (v in c("BipartiteEvol_simulations")){
  MM_emp=NULL
  AA_emp=NULL
  MM_topo=NULL
  AA_topo=NULL
  p=1
  method=NULL
  for (m in combination){
    data_BD=save
    data_topological=save_bipartite
    data_BD=data_BD[,m]
    data_topological=data_topological[,m]
    
    
    
    aa_emp=NULL
    mm_emp=NULL
    aa_topo=NULL
    mm_topo=NULL
    g=1
    for (b in sample(1:10000,10)){
      
      set.seed(b)
      set.seed(b)
      sample_test = sample(1:nrow(data_topological),
                           size = 0.2*nrow(data_topological),
                           replace = F)
      data.test=data_topological[sample_test,]
      data.train=data_topological[-sample_test,]
      
      
      
      #Training the neural network
      formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
      nn_1layer = neuralnet(formule,
                            data.train,
                            hidden = size_layer[p] , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                            stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
      
      #Test the neural network
      prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
      prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
      
      u=1
      for (i in prediction_nn_test){
        
        if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
        if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
        u=u+1
      }
      
      results = table(data.test$Category, prediction_nn_test[,2])  
      
      
      aa_topo[g]=results[1]/(results[3]+results[1])*100
      mm_topo[g]=results[4]/(results[2]+results[4])*100
      
      #Test the neural network
      prediction_nn=predict(nn_1layer, data_BD[, 3:length(data.test)])
      prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
      
      u=1
      for (i in prediction_nn_test){
        
        if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
        if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
        u=u+1
      }
      
      results = table(data_BD$Category, prediction_nn_test[,2])  
      
      aa_emp[g]=results[1]/(results[3]+results[1])*100
      mm_emp[g]=results[4]/(results[2]+results[4])*100
      
      g=g+1
    }
    AA_emp[p]=paste0(round(mean(aa_emp),0)," +/- ",round(sd(aa_emp),1))
    MM_emp[p]=paste0(round(mean(mm_emp),0)," +/- ", round(sd(mm_emp),1))
    AA_topo[p]=paste0(round(mean(aa_topo),0)," +/- ", round(sd(aa_topo),1))
    MM_topo[p]=paste0(round(mean(mm_topo),0)," +/- ", round(sd(mm_topo),1))
    method[p]=name_combination[[p]]
    p=p+1
  }
  
  results_NN=data.frame(Antagonist_topo=AA_topo,Mutualist_topo=MM_topo,Antagonist_emp=AA_emp,Mutualist_emp=MM_emp,method=method)
  print(results_NN)
  write.table(results_NN,
              paste0("./Figures/SI/results_NN_",v,".csv"),sep=";")
  
}



## >> Table : Global classification without one ecology ----

data_empiric=save

MM=NULL 
AA=NULL
k=1
Name=NULL

for (p in names(table(value2$Type))){  
  g=1
  aa=NULL
  mm=NULL
  value2=save[,c(1,2,108:151)]
  value2$Type=as.character(value2$Type)
  value2=value2[-which(value2$Type==p),]
  for (i in 1:50){
    set.seed(sample(1:10000,1))
    sample_test = sample(1:nrow(value2),#nrow = nombre total de réseaux simulés
                         size = 0.2*nrow(value2),
                         replace = F)
    data.test=value2[sample_test,]
    data.train=value2[-sample_test,]
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Category, prediction_nn_test[,2])  
    aa[g]=results[1]/(results[3]+results[1])*100
    mm[g]=results[4]/(results[2]+results[4])*100
    
    g=g+1
  }
  AA[k]=paste0(mean(na.omit(aa))," +/- ", sd(na.omit(aa)))
  MM[k]=paste0(mean(na.omit(mm))," +/- ", sd(na.omit(mm)))
  Name[k]=p
  k=k+1
}
results_S5A=data.frame(aa=AA,mm=MM,name=Name)

write.table(results_S5A,"./Figures/SI/Table_without_ecology.csv",sep=";")




## >> Classification of the excluded network ----



#Part 1 : Classification accuracy of the excluded networks
classif_accuracy=NULL
n=1
name=NULL
for (p in names(table(save$Type))[-1]){  
  classif=NULL
  g=1
  for (k in 1:50){
    value2=save[,c(1,2,108:151)]
    value2$Category=as.character(value2$Category)
    value2$Type=as.character(value2$Type)
    
    set.seed(k)
    
    data.test=value2[which(value2$Type==p),]
    data.train=value2[-which(value2$Type==p),]
    data.train=data.train[sample(1:dim(data.train)[1],0.8*dim(data.train)[1]),]
    
    
    
    #Training the neural network
    formule = as.formula(c('Category ~ ', paste0(colnames(data.test[, 3:length(data.test)]), collapse = "+"))) #classify networks according to the type
    nn_1layer = neuralnet(formule,
                          data.train,
                          hidden = 20 , 
                          stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) 
    
    #Test the neural network
    prediction_nn=predict(nn_1layer, data.test[, 3:length(data.test)])
    
    prediction_nn_test=prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
    
    u=1
    for (i in prediction_nn_test){
      
      if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
      if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
      u=u+1
    }
    
    results = table(data.test$Type, prediction_nn_test[,2])  
    
    
    classif[k]  = ifelse(dim(results)[2]==1,1,results[,unique(save$Category[which(save$Type==p)])]/sum(results))
    g=g+1
  }
  name[n]=p
  classif_accuracy[n]=paste0(mean(na.omit(classif))," +/- ", sd(na.omit(classif)))
  n=n+1
}
results=data.frame('Classif_accuracy'=classif_accuracy,name=name)
write.table(results,"./Figures/SI/Table_excluded_network.csv",sep=";")

#Part 2: Classification accuracy with all the networks in the training

original_classif=read.table("./Figures/data_Table_1.csv",sep=";")
complementary=unlist(lapply(original_classif$category,function(x){
  ifelse(x=="Mutualistic","Antagonistic","Mutualistic")
}))

original_classif$final_classif=unlist(lapply(1:nrow(original_classif),function(x){
  ifelse(original_classif$percentage[x]>.5,original_classif$category[x],complementary[x])
}))


tab=table(original_classif$type,original_classif$final_classif)
original_result=c()
for (i in 1:nrow(tab)){
  original_result=c(original_result,
                    ifelse(unique(original_classif$category[which(original_classif$type==rownames(tab)[i])])=="Mutualistic",
                           tab[i,][2]/sum(tab[i,]),tab[i,][1]/sum(tab[i,])))
}

results=data.frame('Classification_accuracy'=original_result,name=rownames(tab))

write.table(results,"./Figures/SI/Table_total_robustness_fig.csv",sep=";")




## >> Overfitting in Lasso and Probit regression ----

#******************************************************************************#

# testing overfitting with Lasso regression

data_bdd=save[,c(1,2,108:151)]


MM=NULL
AA=NULL
g=1
for (b in sample(1:10000,50)){
  
  set.seed(b)
  
  formule = as.formula(c('Category ~ ',
                         paste0(colnames(data_bdd[, 3:length(data_bdd)]),
                                collapse = "+"))) 
  
  
  mod_lasso=model.matrix(formule, data_bdd)
  category=data_bdd$Category
  lambda_seq=10^seq(2, -2, by = -.1)
  
  # Splitting the data into test and train. Only using matrix and no data frame
  train = sample(1:nrow(mod_lasso), 4*nrow(mod_lasso)/5)
  
  x_train=mod_lasso[train,]
  y_train=category[train]
  
  cv_output=cv.glmnet(x_train, y_train,
                         alpha = 1, lambda = lambda_seq, 
                         nfolds = 5,family="binomial")
  
  # identifying best lambda
  best_lam=cv_output$lambda.min
  
  
  
  
  x_test = mod_lasso[train,]
  y_test = category[train]
  
  # Rebuilding the model with best lambda value identified
  lasso_best=glmnet(mod_lasso[train,], category[train], 
                       alpha = 1,
                       lambda = best_lam,family="binomial")
  pred=predict(lasso_best, s = best_lam, newx = x_test)
  pred=pred>0
  pred[which(pred==F)]="Antagonistic"
  pred[which(pred==T)]="Mutualistic"
  
  results=table(y_test, pred)
  
  AA[g]=results[1]/(results[3]+results[1])*100
  MM[g]=results[4]/(results[2]+results[4])*100
  
  g=g+1
}
lasso_result=data.frame(MM=paste0(mean(MM),"+/-",sd(MM)),
                        AA=paste0(mean(AA),"+/-",sd(AA)),
                        MA=paste0(100-mean(MM),"+/-",sd(MM)),
                        AM=paste0(100-mean(AA),"+/-",sd(AA)))

write.table(lasso_result,"./Figures/SI/Overfitting_Lasso.csv",sep=";")

# testing overfitting with logit regression

data_bdd=save


MM=NULL
AA=NULL
aa=NULL
mm=NULL
g=1
for (b in sample(1:10000,50)){
  
  set.seed(b)
  
  
  sample_train = sample(1:nrow(data_bdd), 4*nrow(data_bdd)/5)
  train=data_bdd[sample_train,]
  test=data_bdd[-sample_train,]
  
  formule = as.formula(c('Category ~ ',
                         paste0(colnames(data_bdd[, 3:length(data_bdd)]),
                                collapse = "+"))) 
  
  res.glm=glm(formula=formule,family =binomial,data=train)
  
  
  pred=predict.glm(res.glm,newdata=train,type="response")
  pred=pred > 0.5# the value of 0.5 is the threshold of the network to classify them
  
  u=1
  for (i in pred){
    
    if (i=="FALSE"){pred[u]="Antagonistic"}
    if( i=="TRUE"){pred[u]="Mutualistic"}
    u=u+1
  }
  results=table(train$Category,pred)
  
  
  AA[g]=results[1]/(results[3]+results[1])*100
  MM[g]=results[4]/(results[2]+results[4])*100
  
  g=g+1
}
probit_result=data.frame(MM=paste0(mean(MM),"+/-",sd(MM)),
                         AA=paste0(mean(AA),"+/-",sd(AA)),
                         MA=paste0(100-mean(MM),"+/-",sd(MM)),
                         AM=paste0(100-mean(AA),"+/-",sd(AA)))

write.table(probit_result,"./Figures/SI/Overfitting_probit.csv",sep=";")






# adding the post hoc tests


data_empiric=data_empiric[-which(data_empiric$Type %in% c("Anemone-Fish","Plant-Ant")),]
write.table(as.data.frame(TukeyHSD(aov(Connectance ~ Type, data=data_empiric))$Type),"./Figures/SI/Post_hoc_Connectance_type_ecology.csv",sep=";")
write.table(as.data.frame(TukeyHSD(aov(Nestedness ~ Type, data=data_empiric))$Type),"./Figures/SI/Post_hoc_Nestedness_type_ecology.csv",sep=";")
write.table(as.data.frame(TukeyHSD(aov(Modularity ~ Type, data=data_empiric))$Type),"./Figures/SI/Post_hoc_Modularity_type_ecology.csv",sep=";")




## >> Classification of empirical networks using simulated networks ----


#BIPARTITEEVOL SIMULATIONS

data_BipartiteEvol=read.table("../Data/bipartite_sim/Data/data_BipartiteEvol.csv",sep=";")[,c(1,105:148)]
data_BipartiteEvol$Category=as.character(data_BipartiteEvol$Category)
data_BipartiteEvol$Category[which(data_BipartiteEvol$Category=="BipartiteEvol Antagonistic")]="Antagonistic"
data_BipartiteEvol$Category[which(data_BipartiteEvol$Category=="BipartiteEvol Antagonistic")]="Mutualistic"

aa=NULL
mm=NULL
AA=NULL
MM=NULL
g=1


for (k in 1:50){
  
  #Create 2 groups, one for training neural networks and one for test
  data.test=NULL
  data.train=NULL
  set.seed(sample(1:100000,1))
  
  sample_test = sample(1:nrow(data_BipartiteEvol),#nrow = nombre total de réseaux simulés
                       size = 0.2*nrow(data_BipartiteEvol),
                       replace = F)
  data.test = data_BipartiteEvol[sample_test,]
  data.train = data_BipartiteEvol[-sample_test,]
  
  
  #Training the neural network
  formule = as.formula(c('Category ~ ', paste0(colnames(data.train[, 2:length(data.train)]), collapse = "+"))) #classify networks according to the type
  nn_1layer = neuralnet(formule,
                        data.train,
                        hidden = 40 , # hidden = (output+input)/2 max= Nsample/(2(output+input))
                        stepmax = 10 ** 8,lifesign = 'full', linear.output=F,rep=1) # number of max step to train
  
  #Test the neural network
  prediction_nn = predict(nn_1layer, data.test[, 2:length(data.test)])
  prediction_nn_test = prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
  
  u=1
  for (i in prediction_nn_test){
    
    if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
    if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
    u=u+1
  }
  
  results1 = table(data.test$Category, prediction_nn_test[,2])  
  
  aa[g]=results1[1]/(results1[3]+results1[1])*100
  mm[g]=results1[4]/(results1[2]+results1[4])*100
  
  data.test=data_empiric
  
  prediction_nn = predict(nn_1layer, data.test[, 2:length(data.test)])
  prediction_nn_test = prediction_nn > 0.5# the value of 0.5 is the threshold of the network to classify them
  
  u=1
  for (i in prediction_nn_test){
    
    if (i=="FALSE"){prediction_nn_test[u]="Antagonistic"}
    if( i=="TRUE"){prediction_nn_test[u]="Mutualistic"}
    u=u+1
  }
  
  results2 = table(data.test$Category, prediction_nn_test[,2])  
  AA[g]=results2[1]/(results2[3]+results2[1])*100
  MM[g]=results2[4]/(results2[2]+results2[4])*100
  g=g+1
}

final_results_BipartiteEvol=data.frame(Simulations=rbind(paste0("Classification of Mutualistic : ",mean(mm)," +/- ",sd(mm)), 
                                           paste0('Classification of Antagonistic : ',mean(aa)," +/- ",sd(aa))),
                         Empiric=rbind(paste0("Classification of Mutualistic : ",mean(MM)," +/- ",sd(MM)), 
                                       paste0('Classification of Antagonistic : ',mean(AA)," +/- ",sd(AA))))

write.table(final_results_BipartiteEvol,"./Figures/SI/Classif_with_bipartiteevol.csv",sep=";")


#*********************************************************************#

# Tests used in the main text ----

#*********************************************************************#


data_empiric=save

## >> Wilcox test ----

# Wilcox test for nestedness
wilcox.test(data_empiric$Nestedness~data_empiric$Category, data=data_empiric) #test
ggplot(data=data_empiric,aes(x=Category,y=Nestedness))+geom_boxplot()+theme_classic() #plot
aggregate(data_empiric$Nestedness,by=list(data_empiric$Category),mean) #mean values
aggregate(data_empiric$Nestedness,by=list(data_empiric$Category),sd) #sd values

#Wilcox test for connectance
wilcox.test(data_empiric$Connectance~data_empiric$Category, data=data_empiric)
ggplot(data=data_empiric,aes(x=Category,y=Connectance))+geom_boxplot()+theme_classic()
aggregate(data_empiric$Connectance,by=list(data_empiric$Category),mean)
aggregate(data_empiric$Connectance,by=list(data_empiric$Category),sd)

#Wilcox test for modularity
wilcox.test(data_empiric$Modularity~data_empiric$Category, data=data_empiric)
ggplot(data=data_empiric,aes(x=Category,y=Modularity))+geom_boxplot()+theme_classic()
aggregate(data_empiric$Modularity,by=list(data_empiric$Category),mean)
aggregate(data_empiric$Modularity,by=list(data_empiric$Category),sd)

#Wilcox test for diff network size
wilcox.test(data_empiric$Diff~data_empiric$Category, data=data_empiric)
ggplot(data=data_empiric,aes(x=Category,y=Diff))+geom_boxplot()+theme_classic()
aggregate(data_empiric$Diff,by=list(data_empiric$Category),mean)
aggregate(data_empiric$Diff,by=list(data_empiric$Category),sd)

#Wilcox test for sum network size
wilcox.test(data_empiric$Sum~data_empiric$Category, data=data_empiric)
ggplot(data=data_empiric,aes(x=Category,y=Sum))+geom_boxplot()+theme_classic()
aggregate(data_empiric$Sum,by=list(data_empiric$Category),mean)
aggregate(data_empiric$Sum,by=list(data_empiric$Category),sd)


## >> mixed effect models ----


#for nestedness
mod_nestedness=lmer(data=data_empiric,
                     formula=Nestedness~as.character(data_empiric$Category)+(1|Type))
summary(mod_nestedness)
Anova(mod_nestedness)

#for modularity
mod_modularity=lmer(data=data_empiric,
                     formula=Modularity~as.character(data_empiric$Category)+(1|Type))
summary(mod_modularity)
Anova(mod_modularity)

#for connectance
mod_connectance=lmer(data=data_empiric,
                     formula=Connectance~as.character(data_empiric$Category)+(1|Type))
summary(mod_connectance)
Anova(mod_connectance)


#for difference in guild size
mod_diff=lmer(data=data_empiric,
              formula=Diff~as.character(data_empiric$Category)+(1|Type))
summary(mod_diff)
Anova(mod_diff)

#for sum in guild size
mod_sum=lmer(data=data_empiric,
              formula=Sum~as.character(data_empiric$Category)+(1|Type))
summary(mod_sum)
Anova(mod_sum)

## >> Anova + post hoc test for global metrics x type of interaction ----

#for connectance
anov_connectance=lm(Connectance~Type, data=data_empiric) 
emmeans(anov_connectance,pairwise~Type,adjust="bonferroni")

