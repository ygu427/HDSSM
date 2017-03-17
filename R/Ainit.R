### Ainit.R
###
### This R script is used to generate initial system matrix A
### using row-based ERM algorithm.
###
### Written by Yu GU
### Latest version: Feb, 2017


Ainit <- function(inputData, s.prop=.1^6, normalized=TRUE,...){
  
  ### If TRUE, Normalization is processed
  if (normalized){
    meanVector <- colMeans(inputData)
    meanMat <- matrix(rep(meanVector,nrow(inputData)),nrow = nrow(inputData),
                      byrow = TRUE)
    inputData <- inputData - meanMat
    std <- apply(inputData,2,sd)
    stdMat <- matrix(rep(std,nrow(inputData)),nrow = nrow(inputData),
                     byrow = TRUE)
    inputData <- inputData / stdMat
  }
  
  xsmooth <- t(inputData)      ## columns are time points now
  
  
  ## Check the variability of each row of xsmooth
  if(is.null(rownames(xsmooth))){
    vynames <- seq(1, nrow(xsmooth))
  } else {
    vynames <- rownames(xsmooth)
  }
  
  vys <- apply(xsmooth[, -1], 1, var)
  zero.vars <- which(vys==0)
  if (length(zero.vars)>0) {
    warning(paste("Variable", paste(vynames[zero.vars], collapse=", "), "has zero variance (excluding the first time point). We have to exclude this/these variables from network reconstruction."))
    xsmooth <- xsmooth[-zero.vars,]
  }
  
  Tp <- ncol(xsmooth)        ## number of time points
  os <- nrow(xsmooth)        ## object space size
  ss <- nrow(xsmooth)        ## state space size

  ### OLS for system matrix A
  gamma <- matrix(0,ss,ss)
  beta <- matrix(0,ss,ss)
  for (t in 1:Tp){
    gamma <- gamma + xsmooth[,t] %*% t(xsmooth[,t])
    if (t>1) {
      beta <- beta + xsmooth[,t] %*% t(xsmooth[,t-1])
    }
  }
  gamma1 <- gamma - xsmooth[,Tp] %*% t(xsmooth[,Tp])
  
  ### calculate Aols for two cases: large N small p or large p small N
  svec <- svd(gamma1)$d
  if (min(svec)< max(svec)*.1^10 && length(gamma1)!=1) {
    ## This is the large p, small n case.
    Aols <- beta %*% Rinv(diag(diag(gamma1)), s.prop=s.prop)   ## Equation (18)
  } else {                            #small p case
    Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)   ## This is OLS/MLE
  }

  ### Starting from M step for ssmal
  As <- array()
  XTX <- list()
  xlm <- t(xsmooth[,1:(Tp-1)])

  for (k in 1:ss){
    w <- abs(t(Aols[k,]))
    xs <- xlm * repmat(w,nrow(xlm),1)  ## x times weights
    ylm <- matrix(xsmooth[k,-1],ncol=1)
    ww <- t(w) %*% w
    XTX$xtx <- gamma1 * ww
    XTX$xty <- t(beta[k,]*w)

    ## LASSO
    main <- emlasso(ylm, xs, XTX, s.prop=s.prop, ...)
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
    As <- rbind(As,sol[[ks+1]]$beta * w)
  }
  As <- As[-1,]
  rownames(As) <- rownames(xsmooth)
  colnames(As) <- rownames(xsmooth)

  return(list(Ainit=As, norm.data=xsmooth))
}
