## ssmal is row-based ERM algorithm. The difference between row-based and
## matrix based is in the R step.
##
##  Latest Version on May.2016
##  Written by Yu Gu

ssmal <- function(data,A,C,Q,R,initx,initV,max_iter=10,diagQ=FALSE,
                  diagR=FALSE,ARmode=FALSE, s.prop=.1^6, ...){
  verbose <- TRUE
  ## EM algorithm convergence likelihood func slope threshold
  thresh <- 5e-3
  if (length(dim(data))<3) {
    N <- 1
  } else {
    N <- dim(data)[3]
  }

  dataList <- list()

  if (N==1) {
    dataList[[N]]<-data
  } else {
    for (i in 1:N) {
      ## each elt of the 3rd dim gets its own cell
      dataList[[i]]<-data[,,i]
    }
  }

  ss <- nrow(A)     ## state size
  os <- nrow(C)     ## observation size

  alpha <- matrix(0,os,os)
  Tsum <- 0
  for (ex in 1:N) {
    y <- dataList[[ex]]
    Tp <- ncol(y)
    Tsum <- Tsum + Tp
    alpha_temp <- matrix(0,os,os)
    for (t in 1:Tp){
      y_mat<-as.matrix(y[,t])
      alpha_temp <- alpha_temp + y_mat %*% t(y_mat)
    }
    alpha <- alpha + alpha_temp
  }

  previous_loglik <- -Inf
  loglik <- 0
  converged <- FALSE
  num_iter <- 1
  LL <- vector()

  while (!converged & (num_iter <= max_iter)) {

    ## E step
    delta <- matrix(0,os,ss)
    gamma <- matrix(0,ss,ss)
    gamma1 <- matrix(0,ss,ss)
    gamma2 <- matrix(0,ss,ss)
    beta <- matrix(0,ss,ss)
    P1sum <- matrix(0,ss,ss)
    x1sum <- matrix(0,ss,1)
    xcurve <- array(0, dim=c(ss,Tp,N))
    loglik <- 0

    for (ex in 1:N) {
      y = dataList[[ex]]
      estep <- Estep(y,A,C,Q,R,initx,initV,ARmode)
      beta <- beta + estep$beta
      gamma <- gamma + estep$gamma
      delta <- delta + estep$delta
      gamma1 <- gamma1 + estep$gamma1
      gamma2 <- gamma2 + estep$gamma2
      P1sum <- P1sum + estep$V1 + estep$x1 %*% t(estep$x1)
      x1sum <- x1sum + estep$x1
      xcurve[,,ex] <- estep$xsmooth
      loglik <- loglik + estep$loglik
    }
    LL[num_iter] <- loglik

    if (verbose) {
      sprintf("interation %i, loglik = %f",num_iter,loglik)
    }

    num_iter <- num_iter + 1

    ## R step (row-based)
    Tsum1 <- Tsum -N

    ## 10/13/2016.  Xing added an option of LargeP
    svec <- svd(gamma1)$d
    if (min(svec)< max(svec)*.1^10 && length(gamma1)!=1) {
      ## This is the large p, small n case.
      Aols <- beta %*% Rinv(diag(gamma1), s.prop=s.prop)   ## Equation (18)
    } else {                            #small p case
      Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)   ## This is OLS/MLE
    }


    A <- array()

    XTX <- list()
    Tp <- ncol(estep$xsmooth)   # the N-th state sequence
    xlm <- t(matrix(xcurve[,1:(Tp-1),],ss,(Tp-1)*N))

    for (k in 1:ss){
      w <- abs(t(Aols[k,]))
      xs <- xlm * repmat(w,nrow(xlm),1)  ## x times weights
      ylm <- matrix(xcurve[k,2:Tp,],(Tp-1)*N,1)
      ww <- t(w) %*% w
      XTX$xtx <- gamma1 * ww
      XTX$xty <- t(beta[k,]*w)

      ## LARS-Lasso
      main <- emlars(ylm, xs, XTX, ...)
      sol <- main$history

      ## select min bic
      df <- vector()
      rss <- vector()
      for (d in 1:length(sol)) {
        df[d] <- sol[[d]]$df
      }

      for (r in 1:length(sol)){
        rss[r] <- sol[[r]]$rss
      }

      DF <- df[-is.na(df)]
      RSS <- rss[-is.na(rss)]

      ch <- rep(NA,length(sol)-1)

      for (i in 1:length(sol)-1){
        comb_num <- choose(ss,length(sol[[i+1]]$active_set))  # number of combinations
        ch[i] <- 2*log(comb_num)  # extra term in extended BIC
      }

      ga <- 1
      n <- nrow(xs)
      # RSS version; extended BIC by Jiajua Chen
      bic <- log(n)*DF + n*log(RSS/n) + ga*ch
      ks <- which(bic == min(bic))
      A <- rbind(A,sol[[ks+1]]$beta * w)
    }
    A <- A[-1,]

    ## M step
    Q <- as.matrix(Q)

    if (diagQ) {
      Q <- (gamma2 - A %*% t(beta)) / Tsum1
      Q <- diag(diag(Q))
    }

    if (!ARmode) {
      C <- diag(os)   # Let C=I
      R <- R
      if (diagR) {
        R <- (alpha - C %*% t(delta)) / Tsum
        R <- diag(diag(R))
      }
    }
    initx <- x1sum / N
    initV <- P1sum / N - initx %*% t(initx)

    emConverge <- em_converged(loglik,previous_loglik,thresh)
    converged <- emConverge$converged
    previous_loglik <- loglik
  }
  ## Output
  return(list(A=A,C=C,Q=Q,R=R,initx=initx,initV=initV,LL=LL,
              xcurve=xcurve,bic=bic,num_iter=num_iter))
}
