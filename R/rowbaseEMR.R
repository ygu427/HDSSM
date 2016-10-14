## rowbaseEMR is row-based ERM algorithm. The difference between row-based and
## matrix based is in the R step.
##
##  Latest Version on Oct.2016
##  Written by Yu Gu

rowbaseEMR <- function(y,initA,initC,initQ,initR,initx,initV,max_iter=100,diagQ=FALSE,
                  diagR=FALSE,ARmode=FALSE, s.prop=.1^6, ...){
  verbose <- TRUE
  ## EM algorithm convergence likelihood func slope threshold
  thresh <- 5e-3

  ss <- nrow(initA)     ## state size
  os <- nrow(initC)     ## observation size

  alpha <- matrix(0,os,os)
  Tp <- ncol(y)
  for (t in 1:Tp){
    y_mat<-as.matrix(y[,t])
    alpha <- alpha + y_mat %*% t(y_mat)
  }

  previous_loglik <- -Inf
  loglik <- 0
  converged <- FALSE
  num_iter <- 1
  LL <- vector()
  estA <- initA
  estC <- initC
  estQ <- initQ
  estR <- initR
  estx <- initx
  estV <- initV

  while (!converged & (num_iter <= max_iter)) {

    ## E step
    delta <- matrix(0,os,ss)
    gamma <- matrix(0,ss,ss)
    gamma1 <- matrix(0,ss,ss)
    gamma2 <- matrix(0,ss,ss)
    beta <- matrix(0,ss,ss)
    P1sum <- matrix(0,ss,ss)
    x1sum <- matrix(0,ss,1)
    xcurve <- matrix(0,ss,Tp)
    loglik <- 0


    estep <- Estep(y,estA,estC,estQ,estR,estx,estV,ARmode)
    beta <- estep$beta
    gamma <- estep$gamma
    delta <- estep$delta
    gamma1 <- estep$gamma1
    gamma2 <- estep$gamma2
    P1sum <- estep$V1 + estep$x1 %*% t(estep$x1)
    x1sum <- estep$x1
    xcurve <- estep$xsmooth
    loglik <- estep$loglik

    LL[num_iter] <- loglik

    if (verbose) {
      sprintf("interation %i, loglik = %f",num_iter,loglik)
    }

    num_iter <- num_iter + 1

    ## R step (row-based)
    Tp1 <- Tp - 1

    ## 10/13/2016.  Xing added an option of LargeP
    svec <- svd(gamma1)$d
    if (min(svec)< max(svec)*.1^10 && length(gamma1)!=1) {
      ## This is the large p, small n case.
      Aols <- beta %*% Rinv(diag(gamma1), s.prop=s.prop)   ## Equation (18)
    } else {                            #small p case
      Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)   ## This is OLS/MLE
    }


    A <- vector()

    XTX <- list()
    #Tp <- ncol(estep$xsmooth)   # the N-th state sequence
    xlm <- t(xcurve[,1:(Tp-1)])

    for (k in 1:ss){
      w <- abs(t(Aols[k,]))
      xs <- xlm * repmat(w,nrow(xlm),1)  ## x times weights
      ylm <- as.matrix(xcurve[k,2:Tp])
      ww <- t(w) %*% w
      XTX$xtx <- gamma1 * ww
      XTX$xty <- t(beta[k,]*w)

      ## LARS-Lasso
      main <- emlasso(ylm, xs, XTX,s.prop=s.prop, ...)
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
    estA <- A

    ## M step
    estQ <- as.matrix(estQ)

    if (diagQ) {
      estQ <- (gamma2 - estA %*% t(beta)) / Tp1
      estQ <- diag(diag(estQ))
    }

    if (!ARmode) {
      estC <- diag(os)   # Let C=I
      estR <- as.matrix(estR)
      if (diagR) {
        estR <- (alpha - estC %*% t(delta)) / Tp
        estR <- diag(diag(estR))
      }
    }
    estx <- x1sum
    estV <- P1sum - estx %*% t(estx)

    emConverge <- em_converged(loglik,previous_loglik,thresh)
    converged <- emConverge$converged
    previous_loglik <- loglik
  }
  ## Output
  return(list(estA=estA,estC=estC,estQ=estQ,estR=estR,estx=estx,estV=estV,LL=LL,
              xcurve=xcurve,bic=bic,num_iter=num_iter))
}
