
#####SCM#####
SCM<- function(matrix){
  t <- ncol(matrix)
  scm <- (t-1)^(-1)*matrix%*%t(matrix)
  return(scm)
}

#####SUPER TRIAL#####
SuperTrial <- function(list, n=16, length_epochs=129, tidy=F, col=17){
  P <- length(list)
  mean_epoch <- matrix(NA, n, length_epochs)
  super_trial <- list()
  scm_st <- list()
  if(tidy==F){
    for(i in 1:P){
      new_trials <- list()
      indices<- which(list[[i]]$labels==1, arr.ind = T)
      epochs <- list[[i]]$epochs[indices[,2], -col ,]
      N <- sum(list[[i]]$labels)
      for(c in 1:n){
        for(j in 1:129){
          mean_epoch[c,j] <- sum(epochs[,c,j])/N
        }
      }
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,-col,] )
        super_trial[[i]]<- new_trials
      }
    }
  }
  else{
    for(i in 1:P){
      new_trials <- list()
      indices<- which(list[[i]]$labels==1, arr.ind = T)
      epochs <- list[[i]]$epochs[indices[,2], ,]
      N <- sum(list[[i]]$labels) 
      for(c in 1:n){
        for(j in 1:129){
          mean_epoch[c,j] <- sum(epochs[,c,j])/N
        }
      }
      for(t in 1:dim(list[[i]]$epochs)[1]){
        new_trials[[t]] <- rbind(mean_epoch,list[[i]]$epochs[t,,] )
        super_trial[[i]]<- new_trials
      }
    }
  }
  return(list(mean_epoch=mean_epoch, super_trial=super_trial))
}


##### RIEMANNIAN DISTANCE#####

RaoDist <- function(matrix1, matrix2){
  (sum(log(eigen(solve(matrix1)%*%matrix2)$values)^2))
}

#####VECT#####

vect <- function(matrix){
  require(lavaan)
  d <- diag(diag(matrix))
  m <- (matrix-diag(diag(matrix)))*sqrt(2)
  lav_matrix_vechr(m+d)
}

#####LOGARITHMIC MAPPING#####
log_map <- function(scm_list, C){
  require(powerplus)
  require(expm)
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  mapped_values <- list()
  for(i in 1:length(scm_list))
    mapped_values[[i]] <- C_star%*%logm(C_2star%*%scm_list[[i]]%*%C_2star)%*%C_star
  return(mapped_values)
}

#####EXPONENTIAL MAPPING#####
exp_map <- function(mapped_values, C){
  require(powerplus)
  require(expm)
  C_star <- Matpow(C, 0.5,1)
  C_2star <- Matpow(C, -0.5,1)
  orig_space_values <- list()
  for(i in 1:length(mapped_values))
    orig_space_values[[i]] <- C_star%*%expm(C_2star%*%mapped_values[[i]]%*%C_2star)%*%C_star
  return(orig_space_values)
}

#####GEOMETRIC MEAN#####
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


#####MINIMUM DISTANCE TO RIEMANNIAN MEAN#####

MDRM <- function(scm_list, centroids_list){
  P <- length(scm_list)
  classes <- rep(0,P)
  distances <- array(NA, dim=c(P,length(centroids_list))) #centroids_list= LENGTH 2
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



#####SCM TRANSFORMATION FOR LINEAR SVM APPLICATION######

scm_transform <- function(scm_list, centroid){
  require(powerplus)
  require(expm)
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


#####EMPIRICAL VARIANCE FUNCTION#####

emp_var_function <- function(mu, scm_list) {
  P <- length(scm_list)
  P^(-1)*Reduce("+",lapply(scm_list, function(x) RaoDist(x, mu)))
}

eta <- function(sigma) -1/(2*sigma^2)

######NORMALIZATION CONSTANT ESTIMATION####à
zeta_est <- function(eta,n_channels,R){
  require(mvtnorm)
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

#####RIEMANN-GAUSS DISTRIBUTION#####

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

#####CREATION OF SCMS#####

createAllSCM <- function(list, n=8, st=T, tidy=F, col=17) {
  tmp <- list()
  final_list <- list()
  if(tidy==F){
    if(st==T){
      for(j in 1:n){
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
          tmp[[i]] <- SCM(list[[j]]$epochs[i,-col,])
          final_list[[j]] <- tmp
        }
      }
    }
  }
  else{
    if(st==T){
      for(j in 1:n){
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
          tmp[[i]] <- SCM(list[[j]]$epochs[i,,])
          final_list[[j]] <- tmp
        }
      }
    }
  }
  return(final_list)
}

#####GROUPS EXTRACTION#####
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

#####SCM VECTORIZATION AND MISC#####

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
