## Kalman/RTS smoother.
##
## The inputs are the same as for kalman_filter.
## The outputs are almost the same, except we condition on y[,1:Tp]
## (and u[,1:Tp] if specified), instead of on y[,1:t].
##
##  Latest Version on May.2016
##  Written by Yu Gu

kalman_smoother<-function(y,A,C,Q,R,init_x,init_V,...){
  os<-nrow(y)  ## observation size
  Tp<-ncol(y)
  ss <- nrow(A)  ## state size

  ## Set Default Params
  model <- matrix (1,1,Tp)
  u <- array()
  B <- array()
  ndx <- list()

  getArgs<-list(...)
  nargs<-length(getArgs)
  if (nargs!=0) {
    for (i in seq(1,nargs, by=2)) {
      arg <- getArgs[[i]]
      if (arg == "model") {model <- getArgs[[i+1]]}
      else if (arg == "u") {u <- getArgs[[i+1]]}
      else if (arg == "B") {B <- getArgs[[i+1]]}
      else if (arg == "ndx") {ndx <- getArgs[[i+1]]}
      else {stop("Input parameter is unrecognizable!")}
    }# end for
  }# end if

  xsmooth <- matrix(0,ss,Tp)
  Vsmooth <- array(0,dim = c(ss,ss,Tp))
  VVsmooth <- array(0,dim = c(ss,ss,Tp))

  ## Forward pass
  kf <- kalman_filter(y,A,C,Q,R,init_x,init_V,...)
  xfilt <- kf$x
  Vfilt <- kf$V
  VVfilt <- kf$VV
  loglik <- kf$loglik


  ## Backward pass
  xsmooth[,Tp] <- xfilt[,Tp]
  Vsmooth[,,Tp] <- Vfilt[,,Tp]

  for (t in (Tp-1):1) {
    m <- model [t+1]
    Ain<-A[,,m]
    Cin<-C[,,m]
    Qin<-Q[,,m]
    Rin<-R[,,m]

    xsmooth_mat <- xsmooth[,t+1]
    Vsmooth_mat <- Vsmooth[,,t+1]
    xfilt_mat<-xfilt[,t]
    Vfilt_mat<-Vfilt[,,t]
    Vfilt_future_mat<-Vfilt[,,t+1]
    VVfilt_mat<-VVfilt[,,t+1]

    if (is.na(B)){
      nextstep <- smooth_update(xsmooth_mat,Vsmooth_mat,xfilt_mat,Vfilt_mat,
                                Vfilt_future_mat,VVfilt_mat,Ain,Qin,Rin,B,u)
      xsmooth[,t] <- nextstep$xsmooth
      Vsmooth[,,t] <-nextstep$Vsmooth
      VVsmooth[,,t+1] <- nextstep$VVsmooth_future
    } else {
      nextstep <- smooth_update(xsmooth_mat,Vsmooth_mat,xfilt_mat,Vfilt_mat,
                                Vfilt_future_mat,VVfilt_mat,Ain,Qin,Rin,B[,,m],u[,t+1])
      xsmooth[,t] <- nextstep$xsmooth
      Vsmooth[,,t] <-nextstep$Vsmooth
      VVsmooth[,,t+1] <- nextstep$VVsmooth_future
    }
  }

  VVsmooth[,,1] <- matrix(0,ss,ss)

  ## Output
  return(list(xsmooth=xsmooth,Vsmooth=Vsmooth,VVsmooth=VVsmooth,
              loglik=loglik))
}# end function
