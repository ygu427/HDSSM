# KALMAN_UPDATE Do a one step update of the kalman filter
#
# Inputs:
# y[,t]   - the observation at time t
# A - the system matrix
# C - the observation matrix
# Q - the system covariance
# R - the observation covariance
# x - E[X|y[,1:t-1]] prior mean
# V - Cov[X|y[,1:t-1]] prior covariance

# OPTIONAL INPUTS (string/value pairs [default in brackets])
# "initial" - 1 means x and V are taken as initial conditions, then A and Q can be ignored [0]
# "u" -- the control signal at time t
# "B" -- the input regression matrix
#
# OUTPUTS (where X is the hidden state being estimated)
#  xnew =   E[ X | y[,1:t] ]
#  Vnew = Var[ X[,t] | y[,1:t] ]
#  VVnew = Cov[ X[,t], X[,t-1] | y[,1:t] ]
#  loglik = log P(y[,t] | y[,1:t-1]) log-likelihood of innovation
#
#   Latest Version on May. 2016
#   Written by Yu Gu

kalman_update<-function(A,C,Q,R,y,x,V,...) {
  # Set default params
  u<-array()
  B<-array()
  initial <- FALSE
  ## Get optional input parameters
  getArgs<-list(...)
  nargs<-length(getArgs)
  for (i in seq(1,nargs,2)) {
    arg <- getArgs[[i]]
    if (arg == "u") {u <- getArgs[[i+1]]}
    else if (arg == "B") {B <- getArgs[[i+1]]}
    else if (arg == "initial") {initial <- getArgs[[i+1]]}
    else {stop("Input parameter is unrecognizable!")}
  }
  #  xpred = E[X_t+1 | y[,1:t]
  #  Vpred = Cov[X_t+1 | y[,1:t]]
  if (initial){
    if (is.na(u)){
      xpred <- x;
    } else {xpred <- x + B %*% u}
    Vpred <- V
  } else {
    if(is.na(u)){
      xpred <- A %*% x
    } else {xpred <- A %*% x + B %*% u}
    Vpred <- A %*% V %*% t(A) + Q
  }

  e <- y - C %*% xpred    # error (innovation)
  ss <- nrow(A)
  S <- C %*% Vpred %*% t(C) + R  # covariance
  os <- nrow(R)
  invR <- diag(1/diag(R))
  m<- matrix(0,1,length(e))  # mean

  loglik <- gaussian_prob(e,m,S,TRUE)

  # Recursive formula for Vnew
  VR <- Vpred
  for (i in 1:os) {
    C_mat <- t(as.matrix(C[i,]))  ## Coerce vector as matrix
    temp <- as.numeric(1/(R[i,i]+ C_mat %*% VR %*% t(C_mat)))
    VR <- VR - temp * VR %*% t(C_mat) %*% C_mat %*% t(VR)
  }

  ## Note that due to numerical error, VR thus computed may not be
  ## *exactly* symmetric. Over many iterations, such small asymmetry
  ## can be accumulated into a large bias.  That's why I decide to
  ## manually symmetrize it.
  Vnew <- (VR + t(VR))/2
  K <- Vnew %*% t(C) %*% invR
  xnew <- xpred + K %*% e
  VVnew <- (diag(ss) - K %*% C) %*% A %*% V
  #Output
  return(list(xnew=xnew, Vnew=Vnew, loglik=loglik, VVnew=VVnew, K=K))
}# end function
