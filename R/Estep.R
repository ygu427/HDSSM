## Compute the (expected) sufficient statistics for a single kalman filter sequence.
## ARMODE=1 specifies that C=I, R=0,i.e., a Gauss-Markov process. (Default 0).
##
## Latest Version on May.2016
## Written by Yu Gu

Estep <- function (y,A,C,Q,R,initx,initV,ARmode) {
  
  y <- as.matrix(y)
  os <- nrow(y)    ## observation size
  Tp <- ncol(y)
  
  ss <- nrow(A)    ## state size 
  
  if (is.array(A) & is.array(C) & is.array(Q) & is.array(R)){
    check <- is.matrix(A)+is.matrix(C)+is.matrix(Q)+is.matrix(R)
    if (check==4){
      A <- array(A,dim=c(dim(A),1))
      C <- array(C,dim=c(dim(C),1))
      Q <- array(Q,dim=c(dim(Q),1))
      R <- array(R,dim=c(dim(R),1))
    } else{stop("Inputs have unmatched dimension size")}
  } else {stop("Inputs must be either matrix or array")}
  
  if (ARmode) {
  	xsmooth <- y
  	Vsmooth <- array(0,dim=c(ss,ss,Tp))   ## no uncertainty about the hidden states
  	VVsmooth <- array(0,dim=c(ss,ss,Tp))
  	loglik <- 0
  } else {
  	kSmooth <- kalman_smoother(y,A,C,Q,R,initx,initV)
  	xsmooth <- kSmooth$xsmooth
  	Vsmooth <- kSmooth$Vsmooth
  	VVsmooth <- kSmooth$VVsmooth
  	loglik <- kSmooth$loglik
  } 
  
  delta <- matrix(0,os,ss)
  gamma <- matrix(0,ss,ss)
  beta <- matrix(0,ss,ss)
  
  for (t in 1:Tp) {
  	y_mat<-as.matrix(y[,t])
  	xsmooth_mat<-as.matrix(xsmooth[,t])
  	Vsmooth_mat<-as.matrix(Vsmooth[,,t])
  	xsmooth_former_mat<-as.matrix(xsmooth[,t-1])
  	VVsmooth_mat<-as.matrix(VVsmooth[,,t])
  	delta <- delta + y_mat %*% t(xsmooth_mat)
  	gamma <- gamma + xsmooth_mat %*% t(xsmooth_mat) + Vsmooth_mat
  	if (t>1){
  		beta <- beta + xsmooth_mat %*% t(xsmooth_former_mat) + VVsmooth_mat
  	}
  } 
  
  xsmooth_T_mat<-as.matrix(xsmooth[,Tp])
  Vsmooth_T_mat<-as.matrix(Vsmooth[,,Tp])
  xsmooth_1_mat<-as.matrix(xsmooth[,1])
  Vsmooth_1_mat<-as.matrix(Vsmooth[,,1])
  gamma1 <- gamma - xsmooth_T_mat %*% t(xsmooth_T_mat) - Vsmooth_T_mat
  gamma2 <- gamma - xsmooth_1_mat %*% t(xsmooth_1_mat) - Vsmooth_1_mat
  
  x1 <- xsmooth_1_mat
  V1 <- Vsmooth_1_mat
  xsmooth <- xsmooth
  
  # Output
  return(list(beta=beta,gamma=gamma,delta=delta,gamma1=gamma1,
              gamma2=gamma2,xsmooth=xsmooth,x1=x1,V1=V1,loglik=loglik))

}