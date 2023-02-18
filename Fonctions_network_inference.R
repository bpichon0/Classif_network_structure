
#packages, please make sure that there packages are well installed
library(RPANDA)
library(tidyverse)
library(FactoMineR) 
library(factoextra)
library(MLmetrics)
library(corrplot)
library(emmeans)
library(matrixStats)
library(maxnodf)
library(bipartite)
library(reshape2)
library(neuralnet)
library(lmerTest)
library(lme4)
library(car)
library(ggpubr)
library(igraph)
library(bmotif)



#Fonction####
#Spectral density calculation
spectR_network_abs <-
  function(network,
           weighted = NULL,
           bandwith = bw.ucv,#use unbiaised function to determine the bandwidth of the density
           n = 200) { #number of point used to describe the density
    # adjacency matrix
    matrix_network <- as.matrix(network)
    A = rbind(cbind(matrix(
      0, nrow = nrow(network), ncol = nrow(network)
    ), matrix_network),
    cbind(t(matrix_network), matrix(
      0, nrow = ncol(network), ncol = ncol(network)
    )))
    # degree matrix
    D = matrix(0, nrow = nrow(A), ncol = ncol(A))
    diag(D) <- c(rowSums(A))
    
    # Laplacian graph
    
    MGL <- D - A #laplacian matrix
    
    rowSums(MGL)
    colSums(MGL)
    
    
    sqrt_D <-
      eigen(D)$vectors %*% diag(1 / sqrt(eigen(D)$values)) %*% t(eigen(D)$vectors) #used to normalised the laplacian matrix
    nMGL = diag(nrow(D)) - sqrt_D %*% A %*% sqrt_D #normalized laplacian matrix
    
    skewness <- function(x, na.rm = FALSE) {#calculation of the skewness of the density
      if (is.matrix(x))
        apply(x, 2, skewness, na.rm = na.rm) 
      else if (is.vector(x)) {
        if (na.rm)
          x <- x[!is.na(x)]
        n <- length(x)
        (sum((x - mean(x)) ^ 3) / n) / (sum((x - mean(x)) ^ 2) / n) ^ (3 /
                                                                         2)
      }
      else if (is.data.frame(x))
        sapply(x, skewness, na.rm = na.rm)
      else
        skewness(as.vector(x), na.rm = na.rm)
    }
    
    e = eigen(nMGL , symmetric = T, only.values = F)#eigenvalue of the laplacian matrix
    x = abs(e$values) #absolute values because approximation of 0 can be negative
    d = dens_network(x, n = n, bw = bandwith) #calculation of the density
    bw = d$bw #getting the bandwith
    dsc = d$y / (integr(d$x, d$y))
    principal_eigenvalue <- max(x)
    skewness_v <- skewness(x)
    peak_height <- max(dsc)
    gaps <- abs(diff(x))
    gapMat <- as.matrix(gaps)
    modalities <- c(1:length(gapMat))
    gapMatCol <- cbind(modalities, gapMat)
    eigenGap <-
      subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 2]))
    res <-
      list(
        eigenvalues = x,
        principal_eigenvalue = principal_eigenvalue,
        asymmetry = skewness_v,
        peakedness = peak_height,
        eigengap = eigenGap[, 1],
        bw = bw,
        y = d$y,
        adj=A,
        dsc=dsc
      )
    
    class(res) <- "spectR"
    
    return(res)
  }

#Do the integral used for the spectral density
integr <- function(x, f) {
  if (!is.numeric(x)) {
    stop("The variable of integration \"x\" is not numeric.")
  }
  if (!is.numeric(f)) {
    stop("The integrand \"f\" is not numeric.")
  }
  if (length(x) != length(f)) {
    stop("The lengths of the variable of integration and the integrand do not match.")
  }
  n = length(x)
  integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] +
                                                    f[1:(n - 1)]))
  return(integral)
}

#dens_network
dens_network <- function(x,
                         bw = bw.ucv,
                         n = 120,
                         from = NULL,
                         to = NULL) {
  if (has.na <- any(is.na(x))) {
    x <- na.omit(x)
    if (length(x) == 0)
      stop("too infinite.")
  }
  kernelG <-
    function(x, mean = 0, sd = 0.1)
      dnorm(x, mean = mean, sd = sd)
  
  sd <- (if (is.numeric(bw))
    bw[1]
    else
      bw(x))
  
  if (is.null(from)) {
    from <- 0
  }
  if (is.null(to)) {
    to <- 2
  }
  
  X <- seq(from, to, len = n)
  M <- outer(X, x, kernelG, sd = sd)
  structure(
    list(
      x = X,
      y = rowMeans(M),
      bw = sd,
      call = match.call(),
      n = length(x),
      data.name = deparse(substitute(x)),
      has.na = has.na
    ),
    class = "density"
  )
}


 
#NESTEDNESS INDICE

#NODF index for empiric networks using NODF_c methods (https://doi.org/10.1101/2020.03.20.000612)

NODF_network <- function(network){
  m=1
  NODF_c_empiric=NULL
  for(i in network){
    i=as.matrix(i)
    NODF_c_empiric[m]=NODFc(i,quality=0)
    m=m+1
    
  }
  return(NODF_c_empiric)
  
}


modularity_network<- function (network){
  modules=NULL
  m=1
  for ( k in network){
    network = as.matrix(network)
    modules[m]=computeModules(network,method = "Beckett")
    m=m+1}
  return(modules)
}

#GENERATION OF NESTED MATRIX

nested_matrix <- function(rows, cols, coeff = 1) {
  row = max(rows, cols)
  col = min(rows, cols)
  nested <- matrix(0, row, col)
  m = nrow(nested)
  fun = function(x) { #function that reproduce nested pattern
    (1 / (row + 1)) ^ (coeff - 1) * (-x + (row + 1)) ^ coeff  + (coeff - 1)
  }
  m = floor(fun(1:col))
  for (i in 1:ncol(nested)) {
    if (m[i] == 0) {
      m[i] = 1
    }
    nested[1:m[i], i] <- 1
  }
  return(nested)
}

#MATRIX RANDOMINSATION WHILE CONSERVING CONNECTANCE

rando_matrix <- function(matrix, chgt) {
  matrix <- matrix
  #choose an interaction to remove
  for (j in 1:chgt) {
    row1 <- floor(runif(1, 1, nrow(matrix)))
    col1 <- floor(runif(1, 1, ncol(matrix)))
    #the interaction can't be in a column (row) with only one interaction 
    while (colSums(matrix)[col1] <= 1 |
           rowSums(matrix)[row1] <= 1  |
           matrix[row1, col1] == 0) {
      row1 <- floor(runif(1, 1, nrow(matrix)))
      col1 <- floor(runif(1, 1, ncol(matrix)))
    }
    matrix[row1, col1] <- 0
    #choose an interaction to add 
    row2 <- floor(runif(1, 1, nrow(matrix)))
    col2 <- floor(runif(1, 1, ncol(matrix)))
    while (matrix[row2, col2] == 1) {
      row2 <- floor(runif(1, 1, nrow(matrix)))
      col2 <- floor(runif(1, 1, ncol(matrix)))
    }
    matrix[row2, col2] <- 1
  }
  return(matrix)
}

#DELETE INTERACTION IN ORDER TO REDUCE CONNECTANCE

delete_matrix<- function( matrix, chgt) {
  matrix <- matrix
  #choose an interaction to remove
  for (j in 1:chgt) {
    row1 <- floor(runif(1, 1, nrow(matrix)))
    col1 <- floor(runif(1, 1, ncol(matrix)))
    #the interaction can't be in a column (row) with only one interaction 
    while (colSums(matrix)[col1] <= 1 |
           rowSums(matrix)[row1] <= 1  |
           matrix[row1, col1] == 0) {
      row1 <- floor(runif(1, 1, nrow(matrix)))
      col1 <- floor(runif(1, 1, ncol(matrix)))
    }
    matrix[row1, col1] <- 0
    
  }
  return(matrix)
}



#Generating a modular matrix
module_matrix <- function(row, col, nb_module, perfect = T) {
  if (row / nb_module < 2) {
    stop('too many modules')
  }
  module <- matrix(0, row, col)
  nb_module=nb_module+1
  step_row <- round(seq(0, row, length.out = nb_module))
  step_col = round(seq(0, col, length.out = nb_module))
  if (!perfect & nb_module!=2) {#asymmetric matrix, only if enough modules
    a <- round(runif(floor((nb_module - 2) / 2),-min(step_row[2:(nb_module-1)])+1, 0))
    b = round(runif(floor((nb_module - 2) / 2),-min(step_col[2:(nb_module-1)])+1, 0))
    #use a and -a to keep sum of change in the step equal to 0
    if (nb_module %% 2 == 1) { #if odd number of modules, one more step to have the same number of step
      a <- c(a,-a, 0)
      b = c(b,-b, 0)
    } else{
      a <- c(a,-a)
      b = c(b,-b)
    }
    #randomisation of the asymmetry
    a <- sample(a) 
    b = sample(b)
    step_row[2:(nb_module - 1)] <- step_row[2:(nb_module - 1)] + a
    step_col[2:(nb_module - 1)] = step_col[2:(nb_module - 1)] + b
  }
  for (i in 1:(nb_module - 1)) {
    module[(step_row[i]+1):(step_row[i + 1]), (step_col[i]+1):(step_col[i + 1])] <-
      1
  }
  return(module)
}


# 
# weighted_sample = function (network,taux=70){ 
#   'Fonction that takes a weighted network and return the network without taux*number of interaction
#   in the network interactions
#   network=interaction matrix
#   taux = percentage of interaction that is kept
#   '
#   
#   save=network
#   c=c()
#   dedouble=NULL
#   network=as.matrix(network)
#   colnames(network)=1:dim(network)[2]
#   rownames(network)=1:dim(network)[1]
#   row=dim(network)[1]
#   col=dim(network)[2]
#   if (dim(network)[1]>9 & dim(network)[2]>9){
#     f=taux/100   # % of interaction kept
#     nb_delete=floor((1-f)*sum(network))+1 #number of deleted interactions
#     
#     for (l in 1:floor(nb_delete)){
#       
#       melting_network=melt(network) #we melt the interaction
#       melting_network[,1]=1:dim(network)[1]
#       #melting_network$value=melting_network$value/min(melting_network$value[melting_network$value>0])
#       melting_network$value=floor(melting_network$value)
#       melting_network$Proba=(melting_network$value/sum(melting_network$value))# probability of being sample is proportional to the occurrence of interaction
#       colnames(melting_network)=c("Ligne","Colonne","Nb_inter","Proba")
#       c=c()
#       for (k in 1:dim(network)[2]){
#         c=c(c,rep(k,dim(network)[1]))
#       }
#       melting_network$Colonne=c # we create a matrix with the probability of being sampled for each interaction
#       reduce=melting_network[which(melting_network$Nb_inter>0),]      #We only consider the interaction, not the 0
#       reduce$Number=1:dim(reduce)[1]
#       interac_sample=sample(reduce$Number,prob=reduce$Proba,1) # we selection the interaction
#       row_sample=reduce$Ligne[which(reduce$Number==interac_sample)] #the corresponding row
#       col_sample=reduce$Colonne[which(reduce$Number==interac_sample)] #the corresponding col
#       network[row_sample,col_sample]=network[row_sample,col_sample]-1 #we delete one interaction
#       
#       network=network[rowSums(network)>0,]
#       network=network[,colSums(network)>0]
#     }
#     if (dim(network)[1]>9 & dim(network)[2]>9){
#       dedouble=network
#     }
#   }
#   return(dedouble)
# }
# 
# 
# 
# 
# 
# 

JSDtree_network <- function (networkEigen, names_network) {
  dist.JSD <- function(inMatrix, pseudocount = 1e-06, ...) {
    KLD <- function(x, y) sum(x * log(x/y))
    JSD <- function(x, y) sqrt(0.5 * KLD(x, (x + y)/2) + 
                                 0.5 * KLD(y, (x + y)/2))
    matrixColSize <- length(colnames(inMatrix))
    matrixRowSize <- length(rownames(inMatrix))
    colnames <- colnames(inMatrix)
    resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
    inMatrix = apply(inMatrix, 1:2, function(x) ifelse(x == 0, pseudocount, x))
    for (i in 1:matrixColSize) {
      for (j in 1:matrixColSize) {
        resultsMatrix[i, j] = JSD(as.vector(inMatrix[, 
                                                     i]), as.vector(inMatrix[, j]))
      }
    }
    rownames(resultsMatrix) <- colnames(resultsMatrix) <- colnames
    resultsMatrix <- as.dist(resultsMatrix)
    attr(resultsMatrix, "method") <- "dist"
    return(resultsMatrix)
  }
  
  m <- c()
  d <- c()
  for (i in 1:length(networkEigen)) {
    m[[i]] <- subset(networkEigen[[i]], networkEigen[[i]] >= 0) #=
    
  }
  
  # compute range x
  min_range <- c()
  max_range <- c()
  for (i in 1:length(networkEigen)) {
    br <- bw.bcv(m[[i]]) # use bw.nrd0 because bw.nrd0 not working
    min_range <- c(min_range, min(c(m[[i]]))-br*5 )
    max_range <- c(max_range, max(c(m[[i]]))+br*5 )
  }
  
  # compute density
  d <- c()
  for (i in 1:length(networkEigen)) {
    d[[i]] <- dens_network(m[[i]]  , from=min(min_range)  , to=max(max_range))#/length(m[[i]])) ### try normalization
    # compute integral
    d[[i]]$y = d[[i]]$y/(integr(d[[i]]$x, d[[i]]$y))
    
  }
  Ds <- c()
  for (i in 1:length(d)) {
    #Ds <- as.data.frame(cbind(Ds, d[[i]]$x))
    Ds <- as.data.frame(cbind(Ds, d[[i]]$y))
  }
  colnames(Ds) <- names_network
  #JSD <- as.matrix(JSDist(abs(Ds)))
  JSD <- as.matrix(dist.JSD(Ds))
  return(JSD)
}




Plot_motif_frequency=function(f_motif){
  
  print(ggplot(f_motif)+
    geom_point(aes(x=motif,y=normalise_sum,color=as.factor(nodes)),shape=8)+
    theme_classic()+
    theme(legend.position = "bottom")+
    labs(x="Motif id",y=expression(sqrt("Frequency motif")),color="# of species \n in the motif"))
  
}



plot_spectR_network <- function (spectR, bw = bw.ucv) {
  if (!inherits(spectR, "spectR"))
    stop("object \"spectR\" is not of class \"spectR\"")
  m = abs(spectR$eigenvalues)
  d <- dens_network(m, bw = bw)
  dint <- integr(d$x, d$y)
  dsc <- (d$y / dint)
  d2=tibble(eigen=d$x,dens_eigen=dsc)
  print(ggplot(d2)+
          geom_line(aes(x=eigen,y=dens_eigen),color="blue")+
          theme_classic()+
          labs(x="Eigenvalues",y="Density of eigenvalues"))
  
}
