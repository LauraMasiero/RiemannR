library(ks)
library(MASS)
library(expm)
#library(e1071)
library(SDMTools)
library(powerplus)
library(rlist)
library(CholWishart)
library(mvtnorm)
library(matrixcalc)
library(purrr)
#library(obliqueRF)
library(LncFinder)
#library(xgboost)
#####SCM#####

#Funzione per calcolare la matrice di covarianza campionaria di un trial 
#la matrice di input è NxT, N=numero di canali e T= numero di ossservazioni (lunghezza temporale della registrazione EEG)

SCM<- function(matrix){
  t <- ncol(matrix) 
  scm <- (t-1)^(-1)*matrix%*%t(matrix)
  return(scm)
}

#Supertrial--costruzione

#la lista deve avere una costruzione specifica:
#ogni elemento deve avere un $epochs che contiene le epochs provenienti da ogni canale
#e un $labels che conntiene l'etichetta dell'epoch corrispondente
#... indicano eventuali cose aggiuntive->messo per togliere il canale 17 che in realtà non è un canale nel mio caso
SuperTrial <- function(list, n=16, length_epochs=dim(list[[i]]$epochs)[3], tidy=F){  #n:numero canali. tidy=F significa che devo togliere il canale 17, altrimenti vuol dire che passo una matrice già a posto
  P <- length(list)
  mean_epoch <- matrix(NA, n, length_epochs)#129 dove c'è length epochs
  super_trial <- list()
  scm_st <- list()
  for(i in 1:P){
    new_trials <- list()
    indices<- which(list[[i]]$labels==1, arr.ind = T)
    if(tidy=F) epochs <- list[[i]]$epochs[indices[,2], -17 ,]
    else epochs <-  list[[i]]$epochs[indices[,2], ,]
    N <- sum(list[[i]]$labels) #dovrebbe essere sempre 80
    for(c in 1:n){
      for(j in 1:length_epochs){ #129 dove c'è length epochs
        mean_epoch[c,j] <- sum(epochs[,c,j])/N
      }
    }
    if(tidy=F){
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,-17,] ) 
        super_trial[[i]]<- new_trials
      }
    }
    else{
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,,] ) 
        super_trial[[i]]<- new_trials
      }
    }
  }
  return(list(mean_epoch=mean_epoch, super_trial=super_trial))
}


#####DISTANZA RIEMANNIANA#####

#Distanza in una varietà riemanniana con la metrica di Rao-Fisher 

RaoDist <- function(matrix1, matrix2){
  (sum(log(eigen(solve(matrix1)%*%matrix2)$values)^2))
}

#####VECT#####
#funzione per vettorizzare le matrici di covarianza con l'operatore vect
vect <- function(matrix){
  d <- diag(diag(matrix))
  m <- (matrix-diag(diag(matrix)))*sqrt(2)
  lav_matrix_vechr(m+d) #solo per matrici simmetriche: lower part stack rows equivalente a upper part per simmetria
}

#####MEDIA GEOMETRICA RIEMANNIANA#####
#MDRM
#C=centroide (baricentro)
log_map <- function(scm_list, C){
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  mapped_values <- list()
  for(i in 1:length(scm_list)) 
    mapped_values[[i]] <- C_star%*%logm(C_2star%*%scm_list[[i]]%*%C_2star)%*%C_star
  return(mapped_values)
}

exp_map <- function(mapped_values, C){
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  orig_space_values <- list()
  orig_space_values[[1]] <- C_star%*%expm(C_2star%*%mapped_values%*%C_2star)%*%C_star
  return(orig_space_values)
}

geometric_mean <- function(scm_list, eps=1/3){
  P <- length(scm_list)
  cat("Initializing...", "\n")
  geo_mean <- P^(-1)*Reduce("+", scm_list) #inizializzazione
  cat("Estimating geometric mean...", "\n")
  repeat{
    x <- geo_mean
    S_hat <-P^(-1)*Reduce("+", log_map(scm_list, x)) #media aritmetica nello spazio tangente
    geo_mean <- exp_map(S_hat, x)
    if(norm(S_hat, type="F")<eps) break
  }
  cat("Geometric mean estimated.", "\n")
  return(geo_mean)
}


#centroids_list: primo elemento distanza dal centroide degli eventi, secondo elemento: distanza dal centroide degli zeri
#SOLO PER CLASSIFICAZIONE BINARIA
MDRM <- function(scm_list, centroids_list){
  P <- length(scm_list)
  classes <- rep(0,P)
  distances <- array(NA, dim=c(P,length(centroids_list))) #centroids_list=2
  for(m in 1:P){
    for(c in 1:length(centroids_list)) 
      distances[m,c] <- RaoDist(centroids_list[[c]], scm_list[[m]])
  }
  classes<- apply(distances, 1, which.min)
  class <- rep(0, P)
  for(i in 1:P){
    if(classes[i]==1) class[i] <- 1
  }
  return(list(classes=class, distances=distances))
} 



#####TRASFORMAZIONE SCM PER SVM LINEARE######

scm_transform <- function(scm_list, centroid){
  C_2star <- Matpow(centroid, -0.5,1) 
  transformed_matrices <- list()
  cat("Transforming matrices...", "\n")
  for(i in 1:length(scm_list)) 
    transformed_matrices[[i]] <- logm(C_2star%*%scm_list[[i]]%*%C_2star)
  cat("Vectorizing matrices...", "\n")
  vect_matrices <- lapply(transformed_matrices, vect) 
  vect_matrices <- lapply(vect_matrices, as.vector) 
  final_matrix<- matrix(NA, length(vect_matrices), length(vect_matrices[[1]]))
  cat("Constructing final matrix...","\n")
  for(i in 1:length(vect_matrices)) 
    final_matrix[i,] <- unlist(vect_matrices[[i]])[1:length(vect_matrices[[i]])]
  cat("Done.","\n")
  return(final_matrix)
}


#####DISTRIBUZIONE RIEMANN-GAUSS#####
##Calcolo mu
#funzione di varianza empirica 

emp_var_function <- function(mu, scm_list) {
  P <- length(scm_list)
  P^(-1)*Reduce("+",lapply(scm_list, function(x) RaoDist(x, mu)))
}

eta <- function(sigma) -1/(2*sigma^2)

#stima monte carlo zeta (costante di normalizzazione)
#n_channels:numero canali
#R:numero replicazioni
zeta_est <- function(eta,n_channels,R){ 
  mu <- rep(0,n_channels)
  Sigma <- diag(-1/eta,n_channels)
  set.seed(123) 
  val <- rmvnorm(R,mu,Sigma)
  temp <- rep(1,R)
  for (i in 1:(n_channels-1))
    for (j in (i+1):n_channels)
      temp <- temp*(sinh(abs(val[,i]-val[,j])/2))
  omegap <- ((((pi^((n_channels^2)/2))*8^((n_channels*(n_channels-1))/4)))/(factorial(n_channels)*lmvgamma(n_channels/2,n_channels))) 
  res <- omegap*mean(temp)
  res <- res*(2*pi)^(n_channels/2)*(det(Sigma))^(1/2)
  return(res) 
  }

#funzione per calcolare il rapporto di verosimiglianza

#Funzione di densità di Riemann-Gauss
RG_distr <- function(mu,sigma, p, scm_matrix) zeta_est(eta(sigma),p,10000)^(-1)*exp(-RaoDist(scm_matrix,mu)/(2*sigma^2))


LR_classification <- function(mu, sigma, scm_list){
  cat("Constructing lik_events...", "\n")
  lik_events <- lapply(scm_list, function(x) RG_distr(mu[[1]],sigma[1],16, x))
  cat("Constructing lik_0...", "\n")
  lik_0 <- lapply(scm_list, function(x) RG_distr(mu[[2]],sigma[2],16, x))
  lik_ratios <- list()
  classes <- rep(0,length(scm_list))
  for(i in 1:length(scm_list)){ 
    lik_ratios[[i]] <- log(lik_events[[i]])-log(lik_0[[i]])
    if(lik_ratios[[i]]>=0) classes[i] <- 1
    #cat("Iteration: ", i, "\t")
  }
  return(list(lik_ratios=lik_ratios, classes=classes))
}


createAllSCM <- function(list, n=8, st=T) { #passo un oggetto creato con la funzione SuperTrial o non_adaptive_train/online e crea le SCM di tutti
  tmp <- list()
  final_list <- list() #lista di liste di scm calcolate coi supertrial o no
  if(st==T){
    for(j in 1:n){#solo i primi 8 elementi della lista del test set sono del soggetto 1 (per gli altri soggetti (da 8 a 24) devo mettere n=1 perchè hanno solo un elemento della lista)
      tmp <- list()
      for(i in 1:length(list$super_trial[[j]])){
        tmp[[i]] <- SCM(list$super_trial[[j]][[i]])
        final_list[[j]] <- tmp
      }
    }
  }
  else{
    for(j in 1:n){
      tmp <- list()
      for(i in 1:dim(list[[j]]$epochs)[1]){
        tmp[[i]] <- SCM(list[[j]]$epochs[i,-17,])
        final_list[[j]] <- tmp
      }
    }
  }
  return(final_list)
}

#PER DUE SOLI GRUPPI:1 E 0
extractGroups <- function(list, dataset, events=1){ 
  group_list <- list()
  if(events==1){
    for(t in 1:length(list)) {
      pos <- which(dataset[[t]]$labels==1, arr.ind=T)[,2]
      for(j in 1:length(list[[t]])) group_list[[t]] <- list[[t]][pos]
    }
  } 
  else {
    for(t in 1:length(list)) {
      pos <- which(dataset[[t]]$labels==0, arr.ind=T)[,2]
      for(j in 1:length(list[[t]])) group_list[[t]] <- list[[t]][pos]
    }
  }
  return(group_list)
}


#PER RISPOSTA BINARIA
extractResponse <- function(response_list, data_set, events=1){ #data_set: non_adaptive_train o non_adaptive_online (in questo caso specificare l'intervallo del soggetto considerato, anche su response list)
  if(events==1){
    labels <- as.vector(unlist(lapply(data_set, function(y) y$labels==1)))
    resp_list <- unlist(response_list)[labels==T]
  } 
  else {
    if(events==0){
      labels <- as.vector(unlist(lapply(data_set, function(y) y$labels==0)))
      resp_list <- unlist(response_list)[labels==T]
    }
  }
  resp_list <- resp_list
  return(resp_list)
}


VectorizeSCM <- function(scm_list){ 
  require(purrr)
  final_list <- lapply(flatten(scm_list), vect)
  final_list <- lapply(final_list, as.vector)
  return(final_list)
}


vectToMatrix <- function(scm_vect_list){
  mat <- matrix(NA, length(scm_vect_list), length(scm_vect_list[[1]]))
  for(i in 1:length(scm_vect_list)) mat[i,] <- unlist(scm_vect_list[[i]])
  return(mat)
}

#######################################################
extractBlock <- function(scm_st_list){
  dimension <- unlist(lapply(scm_st_list, function(x) dim(x)[1]))
  dim_block <- unlist(lapply(dimension, function(x) x/2))
  blocks <- list()
  for(d in 1:length(dim_block)){
    new_mat <- matrix(NA, dim_block[d], dim_block[d])
    for(i in 1:dim_block[d]){
      for(j in dim_block[d]:dimension[d]){
        new_mat[j-dim_block[d],i] <- scm_st_list[[d]][j,i]  
      }
    }
    blocks[[d]] <- new_mat
  }
  return(blocks)
}

