for (name_csv in c("1","2A","2B","3A","3B")){
  
  
  test=read.table(paste0("../Figures/Supplementary_paper/Table/Table_S",name_csv,".csv"),sep=",",header = F)
  
  
  
  if (name_csv=="1"){
    N_m=216
    N_a=127
    
    
  }
  
  if (name_csv %in% c("3A","3A")){
    N_m=470

    N_a=705

  }
  
  if (name_csv %in% c("2A","2A")){
    N_m=23185
    
    N_a=27414
    
  }
  
  
  Percentage_classif=data.frame(Antagonistic=0,Mutualistic=0)
  sd_table=data.frame(Antagonistic=0,Mutualistic=0)
  for( i in 1:nrow(test)){
    for( j in seq(1,13,2)){
      percent=cbind(strsplit(test[i,j],split=" ")[[1]][1],
                    strsplit(test[i,j+1],split=" ")[[1]][1])
      sd=cbind(strsplit(test[i,j],split=" ")[[1]][3],
               strsplit(test[i,j+1],split=" ")[[1]][3])
      
      if (i==1 && j==1 && name_csv=="2B"){
        percent=cbind(strsplit(strsplit(test[i,j],split=" ")[[1]][1],split = "¿")[[1]][2],
                      strsplit(test[i,j+1],split=" ")[[1]][1])

        sd=cbind(strsplit(test[i,j],split=" ")[[1]][3],
                 strsplit(test[i,j+1],split=" ")[[1]][3])

      }
      colnames(percent)=c("Antagonistic",'Mutualistic')
      colnames(sd)=c("Antagonistic",'Mutualistic')
      
      Percentage_classif=rbind(Percentage_classif,percent)
      
      sd_table=rbind(sd_table,sd)
    }
    
  }
  Percentage_classif=Percentage_classif[-1,]
  sd_table=sd_table[-1,]
  
  F1_score=function(percent_AA,percent_MM,N_m=216,N_a=127){
    
    AA=percent_AA
    MM=percent_MM
    AM=100-percent_AA
    MA=100-percent_MM
    
    precision_a=(N_a*AA/100)/((N_a*AA/100)+N_m*MA/100)
    recall_a=(N_a*AA/100)/((N_a*AA/100)+N_a*AM/100)
    F1_a=2*precision_a*recall_a/(precision_a+recall_a)
    
    precision_m=(N_m*MM/100)/((N_m*MM/100)+N_a*AM/100)
    recall_m=(N_m*MM/100)/((N_m*MM/100)+N_m*MA/100)
    F1_m=2*precision_m*recall_m/(precision_m+recall_m)
    
    return(list(round(F1_a,2),round(F1_m,2)))
    
  }
  
  
  F1_table=data.frame(Antagonistic=0,Mutualistic=0)
  
  for (i in 1:nrow(Percentage_classif)){
    Fscore=cbind(F1_score(N_m=216,N_a=127,as.numeric(Percentage_classif[i,1]),
                    as.numeric(Percentage_classif[i,2]))[[1]],
             F1_score(N_m=216,N_a=127,as.numeric(Percentage_classif[i,1]),
                      as.numeric(Percentage_classif[i,2]))[[2]])
    colnames(Fscore)=c("Antagonistic",'Mutualistic')
    
    F1_table=rbind(F1_table,Fscore)
    
    
    
  }
  F1_table=F1_table[-1,]
  
  All_data_merged=c()
  
  for (i in 1:nrow(Percentage_classif)){
    for (k in 1:2){
      All_data_merged=c(All_data_merged,
                        paste0(Percentage_classif[i,k]," +/- ",
                               sd_table[i,k]," (",F1_table[i,k],")")
                               )
    }
    
    
    
  }  
    
    
  Final_table=matrix(All_data_merged,nrow = nrow(test),ncol=ncol(test),byrow = T)
  colnames(Final_table)=rep(c("global metrics",	"motifs","spectral density",
                              "global metrics + spectral density"	,
                              "global metrics + motifs",
                              "spectral density + motifs",
                              "spectral density + motifs + global metrics"),each=2)
  rownames(Final_table)=c("Probit","Lasso","Neural network")
  write.table(Final_table,paste0("../Figures/Supplementary_paper/test/Table_S",name_csv,".csv"),sep=";")
    
}
  
  
