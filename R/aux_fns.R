# Matrix ####
#### Generates a random uniform matrix with a given sparsity and max eigenvalue ####
rmatunif <- function(d, density, range=c(-1,1), seed=FALSE, rho_max=1) {
  stopifnot(density>0 && density<=1)
  R <- matrix(0,d,d)
  if(seed) set.seed(seed)
  pos <- sample(1:length(R),ceiling(density*length(R)))
  val <- runif(length(pos),min=range[1],max=range[2])
  R[pos] <- val
  lam_max <- max(abs(eigen(R)$values))
  if (lam_max>rho_max) R <- R*rho_max/lam_max
  return(R)
}

#### Generates a random symmetric uniform matrix with/out diagonal with given sparsity ####
rmatsym <- function(d, density, range=c(-1,1), diag=TRUE, seed=FALSE) {
  stopifnot(density>0 && density<=1)
  R <- matrix(0,d,d)
  if(seed) set.seed(seed)
  ind <- low.tri_vec(d,diag = diag)
  pos <- sample(ind,ceiling(density*length(ind)))
  val <- runif(length(pos),min=range[1],max=range[2])
  R[pos] <- val
  tR <- t(R); diag(tR) <- 0
  return(R+tR)
}

#### Generates a random normal matrix with given sparsity and max eigenvalue ####
rmatnorm <- function(d, density, sd=1, seed=FALSE, rho_max=1) {
  stopifnot(density>=0 && density<=1)
  R <- matrix(0,d,d)
  if(density>0){
  if(seed) set.seed(seed)
  pos <- sample(1:length(R),ceiling(density*length(R)))
  val <- rnorm(length(pos),sd=sd)
  R[pos] <- val
  lam_max <- max(abs(eigen(R)$values))
  if (lam_max>rho_max) R <- R*rho_max/lam_max
  }
  return(R)
}

#### random matix with max eigenvalue specified ####
A_rand <- function(d,ev_max,seed) {
  set.seed(seed)
  A <- array(runif(d^2,-1,1),c(d,d))
  A <- ev_max*A/max(abs(eigen(A)$values))
  return(A)
}

#### random sparse matrix with 0=empty, 1=full ####
A_S <- function(d,sp,seed) {
  A <- array(0,c(d,d))
  s <- ceiling(sp*d^2)
  A[sample(d^2,s)] <- runif(s,-1,1)
  return(A)
}

#### random matrix with with rank and max eigenvalue specified ####
A_LR <- function(d,rk,ev_max,seed) {
  set.seed(seed)
  A <- array(runif(d,-1,1),c(d,rk)) %*% t(array(runif(d,-1,1),c(d,rk)))
  A <- ev_max*A/max(abs(eigen(A)$values))
  return(A)
}

#### random matrix composition of specifed rank and sparsity where the infomation ratio and max eigenvalue are given ####
A_LRS <- function(d,rk,sp,IR,ev_max,seed) {
  set.seed(seed)
  A <- array(runif(d,-1,1),c(d,rk)) %*% t(array(runif(d,-1,1),c(d,rk)))
  S <- A_S(d=d,sp=sp,seed=seed+1)
  A <- A + IR*max(abs(A))*S/max(abs(S))
  A <- ev_max*A/max(abs(eigen(A)$values))
  return(A)
}

#### wrapper function for fixed plus changing matrix generation ####
fpc_mat <- function(m, FUN, FUN_args, seed, fixed=numeric()) {
  d <- FUN_args[["d"]]
  As <- list()
  if (length(fixed)==0) fixed <- array(0,c(d,d))
  for (i in 0:(m-1)) {
    A <- do.call(FUN,c(list(seed = seed + i), FUN_args))
    As <- c(As,list(fixed + A))
  }
  return(As)
}

#### rescales a set of matrices by the max eigenvalue across them ####
stable_VAR <- function(A, spectral_radius) {
  if ("list" %in% class(A)) {
    rhos <- c()
    for (Ai in A) {
      rho <- max(abs(eigen(Ai)$values))
      rhos <- c(rhos,rho)
    }
    for (i in 1:length(A)) {
      A[[i]] <- A[[i]]*spectral_radius/max(rhos) #rhos[i] to fix spectral radius
    }
  } else {
    A <- A*spectral_radius/max(abs(eigen(A)$values))
  }
  return(A)
}

#### sparsity of a given matrix of any dimension ####
A_sp <- function(A) {
  d <- dim(A)
  n_val <- length(which(abs(A)>0))
  n_val/prod(d)
}

#### rank of a given matrix of any dimension ####
A_rk <- function(A) {
  qr(A,LAPACK = FALSE)$rank
}

# #### symmetric matrix from low/upp tri list + diag list ####
# sym_mat <- function(Tri,D,diag=F) {
#   S <- matrix(0, nrows=length(D), ncols=length(D))
#   ####WRITE this!
#   S[lower.tri(S)] = t(S)[lower.tri(S)]
# }

# #### symmetric matrix with given rank ####
# smatLR <- function(D,rank,rho) {
#   D_eigen <- eigen(D)
# }

#### eigenvalues of matrix to 2dp ####
ev <- function(D) {
  print(round(eigen(D)$values,2))
}

#### max eigenvalue in an array of matrices ####
ev_max <- function(D) {
  if (length(dim(D))==3) {
    l_max <- 0
    for (i in 1:dim(D)[3]) {
      lam <- max(abs(eigen(D[,,i])$values))
      if (lam > l_max) l_max <- lam
    }
  } else {
    l_max <- max(abs(eigen(D)$values))
  }
  return(l_max)
}

#### sign of eigenvector dertermined by first entry >thr ####
ev_sign <- function(D,thr=0.1) {
  U <- eigen(D)$vectors
  sgn <- apply(U,2,
               FUN = function(v) {
                 i <- which(abs(v)>thr)[1]
                 if (is.na(i>0)) {i<-1}
                 if (v[i] < -thr) {return(-1)} else {return(1)}})
  return(U%*%diag(sgn))
}

#### sets the sign of the 1st entry in 1st eigenvector to be positive ####
eigen_graph <- function(G) {
  U <- eigen(G)
  U$vectors[,1] <- U$vectors[,1] * sign(U$vectors[1,1])
  return(U)
}

#### switches matrix sign so 1st entry positive ####
.esign <- function(A) {
      S <- diag(sign(A[1,]))
      X <- A %*% S
      return(X)
    }

#### computes the condition number of a matrix ####
cond_num <- function(D) {
  d <- svd(D)$d
  return(head(d,1)/tail(d,1))
}

#### Distance matrix pertubation tests ####
ev_dist <- function(D,kappa=c(1,0.5,0.1,0.01)) {
  D_eigenvalues <- c()
  D_eigenvectors <- c()
  for (k in kappa) {
    eigen_res <- eigen(exp(-k*D^2))
    D_eigenvalues <- c(D_eigenvalues,eigen_res$values)
    D_eigenvectors <- c(D_eigenvectors,eigen_res$vectors)
  }
  D_eigenvalues <- matrix(D_eigenvalues,nrow=dim(D)[1])
  D_eigenvectors <- array(D_eigenvectors,c(dim(D),length(kappa)))
  matplot(kappa,t(D_eigenvalues),type="l",pch=3,#log="y",
          xlab=expression(kappa),ylab=expression(paste(lambda,", eigenvalues")),
          main=bquote("Eigenvalues of"~exp*"(-"*kappa*D^2*")"~"for varying values of"~kappa))
  for (i in 1:dim(D)[1]) {
    D_eigenvectors[,i,] <- D_eigenvectors[,i,]%*%diag(sign(D_eigenvectors[which.max(abs(D_eigenvectors[,i,1])),i,]))
    matplot(kappa,t(D_eigenvectors[,i,]),type="l",
            xlab=expression(kappa),ylab="Eigenvector Component Value",
            main=bquote("Eigenvector"~.(i)~"of"~exp*"(-"*kappa*D^2*")"~"for varying values of"~kappa))
  }
  return(list(D_eval = D_eigenvalues,D_evec = D_eigenvectors))
}

#### Threshold abs matrix values strictly less than thr to val ####
th <- function(X,thr,val=0) {
  X[which(abs(X)<thr)] <- val
  return(X)
}

#### inverts abs matrix values greater than thr about thr #### 
.inv_th <- function(A,th) {
    a <- which(abs(A)>=th)
    A[a] <- th/A[a]
    return(A)
  }

#### Calculate sparsity of a given matrix ####
sparsity <- function(X, threshold=0) {
  if ("list" %in% class(X)) X <- unlist(X)
  s <- 0
  for (x in c(X)) {
    if (abs(x)>threshold) s <- s+1
  }
  return(s/length(X))
}

#### Computes a rank rk representation of a given matrix ####
lr_mat <- function(A, rk) {
  svd <- svd(A)
  A_lr <- svd$u[,1:rk] %*% diag(svd$d[1:rk],rk) %*% t(svd$v[,1:rk])
  return(A_lr)
}

#### Threshold singular values or entries in matrix for lr or sparse representation ####
pc_th <- function(X,pc_thr,lowrank=F) {
  if (lowrank){# Low rank aprroximation through svd truncation
    d <- min(dim(X))
    r <- ceiling(d*pc_thr)
    X_svd <- svd(X,nu=r,nv=r)
    A <- X_svd$u%*%diag(X_svd$d[1:r],nrow=r)%*%t(X_svd$v)
    return(list(A,list(r,X_svd$d[r])))
  } else {# strictly less than pc_thr percentile values
    rk <- sort(unique(as.vector(abs(X))))
    X[which(abs(X)<rk[ceiling(length(rk)*(1-pc_thr))])] <- 0
    #print(rk[ceiling(length(rk)*(1-pc_thr))])
    return(list(X,list(pc_thr,rk[ceiling(length(rk)*(1-pc_thr))])))
  }
}

#### Computes the rank approximation of a given matrix ####
lr_aprx <- function(X,rank,diff=F) {
  X_svd <- svd(X,nu=rank,nv=rank)
  A <- X_svd$u%*%diag(X_svd$d[1:rank],nrow=rank)%*%t(X_svd$v)
  if (diff) {
    print(paste("2-norm:",norm((A-X),"2")))
    print(paste("Frobemius-norm:",norm((A-X),"F")))
  }
  return(A)
}

#### matrix power function ####
mat_pow <- function(A,n) {
  E <- eigen(A)
  An <- E$vectors %*% diag(E$values)^n %*% solve(E$vectors)
  return(An)
}

#### low-rank resconstrunction of a matrix or its inverse ####
lr_eig <- function(X,rank,inv=F) {
  E <- eigen(X)
  d <- length(E$values)
  if (inv) {
    v <- tail(1:d,rank)
    D <- diag(1/E$values[v],rank)
  } else {
    v <- head(1:d,rank)
    D <- diag(E$values[v],rank)
  }
  L <- E$vectors[,v] %*% D %*% t(E$vectors[,v])
  return(L)
}

#### rescale matrix max abs value to 1 ####
.renorm <- function(A) {
    l <- max(abs(A), na.rm = T)
    return(A/l)
  }

# Data Structure ####
#### convert list of d-arrays into an d+1-array ####
lst_to_arr <- function(L) {
  d <- dim(L[[1]])
  arr_dim <- c(d, length(L))
  return(array(unlist(L), dim = arr_dim))
}

#### Convert list of matrices into one matrix ####
lst2mat <- function(list) {
  mat <- NULL
  for (l in list) {
    mat <- cbind(mat,l)
  }
  return(mat)
}

# Index ####
#### character vector of coordinates of the upper triangular matrix of a square matrix array inc/exc diagonal ####
upp.tri_ind <- function(d, diag=T) {
  if (diag){
    apply(cbind((matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = T)]-1)%%d+1,
                (matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = T)]-(matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = T)]-1)%%d-1)/d+1)
          ,1,function(x)paste(x[1],x[2],sep=", "))
  } else {
    apply(cbind((matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = F)])%%d,
                (matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = F)]-(matrix(c(1:d^2),ncol = d)[upper.tri(diag(d), diag = F)])%%d)/d+1)
          ,1,function(x)paste(x[1],x[2],sep=", "))
  }
}

#### character vector of coordinates for a square matrix dimension d by column or row ####
mat_ind <- function(d, byrow=F) {
  apply(cbind(c((matrix(c(1:d^2),ncol = d,byrow = byrow)-1)%%d+1),
              c((matrix(c(1:d^2),ncol = d,byrow = byrow)-(matrix(c(1:d^2),ncol = d,byrow = byrow)-1)%%d-1)/d+1))
        ,1,function(x)paste0(x[1],".",x[2]))
}

#### convert index of array to corresponding coordinate of square matrix of dimension d ####
ind_to_rc <- function(ind,d) {
  r <- (ind-1)%%d+1
  c <- (ind-1)%/%d+1
  return(c(r,c))
}

#### indices of d^2 vector for dxd lower triangular matrix ####
low.tri_vec <- function(d, diag=F, byrow=F) {
  if (byrow){
    s <- seq(1,d^2,d+1)
    f <- seq(d^2-d+1,d^2,1)
    by <- d
  } else {
    s <- seq(1,d^2,d+1)
    f <- seq(d,d^2,d)
    by <- 1
  }
  .fun_seq <- function(s,f,by) {
    sapply(1:min(length(s),length(f)),function(l) {seq(s[l],f[l],by)})
  }
  if (diag) {
    ind <- unlist(mapply(FUN = .fun_seq, s=s, f=f, by=by))
  } else {
    ind <- unlist(mapply(FUN = .fun_seq, s=s[-d]+by, f=f[-d], by=by))
  }
  return(ind)
}

#### index for upper-triangular array by column ####
upp.tri <- function(x) {
  return(x[2] * (x[2]-1) / 2 + x[1])
}

#### convert positive matrix indices lists to vector index ####
mat_vec_ind <- function(ii,jj,d,byrow=F) {
  .fun_lst <- function(ii,jj) {
    sapply(1:min(length(ii),length(jj)),function(l) {
      if(ii[l]>0 && ii[l]<=d && jj[l]>0 && jj[l] <=d) {
        ii[l] + (jj[l]-1)*d
      }
      })
  }
  if (byrow) {
    return(unlist(mapply(FUN = .fun_lst, ii=jj, jj=ii)))
  } else {
    return(unlist(mapply(FUN = .fun_lst, ii=ii, jj=jj)))
  }
}

#### convert n indices to matrix coordinates as a 2xn array ####
vec_mat_ind <- Vectorize(function(kk,d,byrow=F) {
  if (byrow) {c((kk-1)%/%d+1,(kk-1)%%d+1)}
  else {c((kk-1)%%d+1,(kk-1)%/%d+1)}
})

#### segment indicator for multiple change-points ####
CP_ind <- function(t1,t2,t0,t) {
  t <- sort(c(0,t0,t))
  t1 <- cut(t1,t,labels=1:(length(t0)+1))
  t2 <- cut(t2,t,labels=1:(length(t0)+1))
  return(upp.tri(sort(as.numeric(c(t1,t2)))))
}

#### order unique values in an array by size ####
val <- function(X) {
  return(sort(unique(as.vector(abs(X)))))
}

#### shift a vector right by k steps ####
shift <- function(x,k) {
  y <- rep(NA,length(x)+k)
  y[(k+1):(length(x)+k)] <- x
  return(y[1:length(x)])
}

shift_mat <- function(X,y) {
  for(i in c(1:dim(X)[1])) {
    X[i,] <- shift(X[i,],y[i])
  }
  return(X)
}

#### square lattice coordinates ####
gridPos_d <- function(d, lim=d/4) {
  matrix(
    c(rep(0:d,each=d+1),rep(0:d,times=d+1)),
    ncol = 2, nrow = (d+1)^2
  )/d*lim
}

# MISC ####
#### 2-norm ####
norm2 <- function(x){
  sqrt(sum(x^2))
}

#### l-infinity norm of matrix or list of matrices (max col sum) ####
Linf_norm <- function(X,print=TRUE) {
  if ("list" %in% class(X)) {
    x_old <- 0 
    for (Xi in X) {
      x <- max(colSums(abs(Xi)))
      x_new <- max(x,x_old)
      if(print) print(x_new)
      x_old <- x
    }
  } else {
    x_new <- max(colSums(abs(X)))
  }
  return(x_new)
}

#### list of points l1-norm distance p from points in dxd grid ####
l1_dist <- function(i,j,p,d,byrow=F) {
  if (p==0) {l <- mat_vec_ind(i,j,d,byrow)}
  else {
    il <- c(i+p,i-p,i,i)
    jl <- c(j,j,j+p,j-p)
    l <- mat_vec_ind(il,jl,d,byrow)
    if (p>1){
      for (q in 1:(p-1)) {
        il <- c(i+p-q,i+p-q,i-p+q,i-p+q)
        jl <- c(j+q,j-q,j+q,j-q)
        l <- append(l, mat_vec_ind(il,jl,d,byrow))
      }
    }
  }
  return(sort(l))
}

#### vectorizes l1_dist ####
vec_l1_dist <- function(l,p,d,byrow=F) {
  l1_dist(l[1],l[2],p,d,byrow)
}

#### d^2xd^2 Pairwise distance matrix from dxd points on lattice ####
spatial_dist <- function(x,y,d) {
  i1 <- x %% d + 1
  i2 <- y %% d + 1
  j1 <- x %/% d
  j2 <- y %/% d
  s <- sqrt((i1-i2)^2 + (j1-j2)^2)
  return(s)
}

#### returns mean of two vectors entrywise ####
mean_vec <- function(x,y) {rowMeans(cbind(x,y))}