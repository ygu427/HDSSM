## Kalman filter.
## INPUTS:
## y[,t]   - the observation at time t
## A - the system matrix
## C - the observation matrix
## Q - the system covariance
## R - the observation covariance
## init_x - the initial state (column) vector
## init_V - the initial state covariance
##
## OPTIONAL INPUTS (string/value pairs [default in brackets])
## 'model' - model[t] = m means using params from model m at time t
## In this case, all above matrices have an additional final dimension,
## i.e. A[,,m], C[,,m], Q[,,m], R[,,m]
## Note that init_x and init_V are independent on model[1].
## 'u'     - u[,t] the control signal at time t
## 'B'     - B[,,m] the input regression matrix for model m
##
##  Latest Version on May. 2016
##  Written by Yu Gu

kalman_filter<-function(y,A,C,Q,R,init_x,init_V, s.prop=.1^6, ...) {
  os <- nrow(y)    ## observation size
  Tp <- ncol(y)
  ss <- nrow(A)    ## state size

  ## Declare and set default value of input parameters
  model <- matrix(1,1,Tp)
  u <- array()
  B <- array()
  ndx <- list()

  ## Get the optional input arguments
  getArgs<-list(...)
  nargs<-length(getArgs)
  if (nargs!=0){
    for (i in seq(1,nargs,2)) {
      arg <- getArgs[[i]]
      if (arg == "model") {model <- getArgs[[i+1]]}
      else if (arg == "u") {u <- getArgs[[i+1]]}
      else if (arg == "B") {B <- getArgs[[i+1]]}
      else if (arg == "ndx") {ndx <- getArgs[[i+1]]}
      else {stop("Input parameter is unrecognizable!")}
    }
  }

  ## Declaration and Initialization
  x <- matrix(0,ss,Tp)
  V <- array(0,dim=c(ss,ss,Tp))
  VV <- array(0,dim=c(ss,ss,Tp))
  K <- array(0,dim=c(ss,os,Tp))
  Ain<-matrix()
  Cin<-matrix()
  Qin<-matrix()
  Rin<-matrix()
  loglik <- 0

  for (t in 1:Tp){
    m <- model[t]
    Ain<-A[,,m]
    Cin<-C[,,m]
    Qin<-Q[,,m]
    Rin<-R[,,m]

    if (t==1){
      prevx<-init_x
      prevV<-init_V
      initial <- TRUE
    }
    else {
      prevx <- x[,t-1]
      prevV <- V[,,t-1]
      initial <- FALSE
    }

    yin<-y[,t]
    if(is.na(u)){
      ## Ain, Cin, Qin, Rin, prevV are matrices
      ## yin, prevx are vectors
      onestep <- kalman_update(Ain,Cin,Qin,Rin,yin,prevx,prevV,"initial",initial)
      x[,t] <- onestep$xnew
      V[,,t] <- onestep$Vnew
      LL <- onestep$loglik
      VV[,,t] <- onestep$VVnew
      K[,,t] <- onestep$K
    } else {
      if(is.na(ndx)){
        ## Ain, Cin, Qin, Rin, prevV, B[,,m] are matrices
        ## yin, prevx, u[,t] are vectors
        onestep <- kalman_update(Ain,Cin,Qin,Rin,yin,prevx,prevV,
                                 "initial",initial,"u",u[,t],"B",B[,,m])
        x[,t] <- onestep$xnew
        V[,,t] <- onestep$Vnew
        LL <- onestep$loglik
        VV[,,t] <- onestep$VVnew
        K[,,t] <- onestep$K
      } else {
        i <- ndx[[t]]
        x[,t] <- prevx
        prevP <- Rinv(prevV, s.prop=s.prop)
        prevPsmall <- prevP[i,i]
        prevVsmall <- 1/prevPsmall
        ## Ain[i,i], Qin[i,i], prevx[i], prevVsmall is numeric
        ## yin, Cin[,i], u[,t], B[i,,m] are vectors
        ## Rin is matrix
        onestep <- kalman_update(Ain[i,i],Cin[,i],Qin[i,i],Rin,yin,prevx[i],
                                 prevVsmall,"initial",initial,"u",u[,t],"B",B[i,,m])
        x[i,t] <- onestep$xnew
        smallV <- onestep$Vnew
        LL <- onestep$loglik
        VV[i,i,t] <- onestep$VVnew
        K[i,,t] <- onestep$K

        smallP <- 1/smallV
        prevP[i,i] <- smallP
        V[,,t] <- Rinv(prevP, s.prop=s.prop)
      }
    }
    loglik <- loglik + LL
  }

  ## OUTPUT
  return(list(x=x,V=V,VV=VV,loglik=loglik,K=K))
}
