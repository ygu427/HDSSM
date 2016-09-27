## One step of the backwards RTS smoothing equations.
## 
## INPUTS
## xsmooth_future = E[X_t+1|T]
## Vsmooth_future = Cov[X_t+1|T]
## xfilt = E[X_t|t]
## Vfilt = Cov[X_t|t]
## Vfilt_future = Cov[X_t+1|t+1]
## VVfilt_future = Cov[X_t+1,X_t|t+1]
## A = system matrix for time t+1
## Q = system covariance for time t+1
## B = input matrix for time t+1 (or [] if none)
## u = input vector for time t+1 (or [] if none)
##
## OUTPUTS:
## xsmooth = E[X_t|T]
## Vsmooth = Cov[X_t|T]
## VVsmooth_future = Cov[X_t+1,X_t|T]
##
## Latest Version on May.2016
## Written by Yu Gu

smooth_update<-function(xsmooth_future,Vsmooth_future,xfilt,Vfilt,
                        Vfilt_future,VVfilt_future,A,Q,R,B,u) {

  if (is.na(B)){
  	xpred <- A %*% xfilt         ## xpred = E[X(t+1) | t]
  }
  else {
  	xpred <- A %*% xfilt + B %*% u
  }
  
  Vpred <- A %*% Vfilt %*% t(A) + Q   ## Vpred = Cov[X(t+1) | t]
  
  ## recursive formula for inv(A'invQA+invVfilt) for inv(pred)
  ss <- nrow(A)
  Vn <- Vfilt
  for (i in 1:ss){
  	Amax <- t(as.matrix(A[i,]))
  	temp <- as.numeric((1/(Q[i,i]+Amax %*% Vn %*% t(Amax))))
  	Vn <- Vn - temp * Vn %*% t(Amax) %*% Amax %*% t(Vn) 
  }
  
  invQ<-diag(1/diag(Q))
  invVpred <- invQ - invQ %*% A %*% Vn %*% t(A) %*% invQ  ## can further simplify when Q=I
  
  ## inv(Cov[X_t+1|t+1])= C'invR C + invVpred, C=I, R diag
  invVfilt_future <- diag(1/diag(R)) + invVpred
  
  ## smoother gain matrix
  J <- Vfilt %*% t(A) %*% invVpred ## smoother gain matrix
  xsmooth <- xfilt + J %*% (xsmooth_future - xpred)
  Vsmooth <- Vfilt + J %*% (Vsmooth_future - Vpred) %*% t(J)
  VVsmooth_future <- VVfilt_future + (Vsmooth_future - Vfilt_future) %*% invVfilt_future %*% VVfilt_future  
  
  ## Output
  return(list(xsmooth=xsmooth,Vsmooth=Vsmooth,VVsmooth_future=VVsmooth_future))
}# end function