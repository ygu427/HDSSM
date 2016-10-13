### initialization.R
###
### This R script includes two functions: Arow.init and Amat.init
### The two functions are used to generate initial system matrix A
### for row-based and matrix-based ERM algorithm respectively.
###
### Written by Yu GU
### Latest version: Sep.30th, 2016


Arow.init <- function(inputData, s.prop=.1^6){
  ### Normalization
  meanVector <- colMeans(inputData)
  meanMat <- matrix(rep(meanVector,nrow(inputData)),nrow = nrow(inputData),
                    byrow = TRUE)
  inputData <- inputData - meanMat
  std <- apply(inputData,2,sd)
  stdMat <- matrix(rep(std,nrow(inputData)),nrow = nrow(inputData),
                   byrow = TRUE)
  inputData <- inputData / stdMat
  Tp <- nrow(inputData)
  os <- ncol(inputData)
  ss <- ncol(inputData)
  xsmooth <- t(inputData)

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
  Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)

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

    ## LARS-Lasso
    stopCriterion = list()
    stopCriterion[[1]] <- c("maxIterations",100)
    main <- emlars(yin = ylm,xin = xs,XTX = XTX,regressiontype = "lasso",
                   stopCriterion = stopCriterion)
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

  return(list(Arow=As, norm.data=xsmooth))
}


Amat.init <- function(inputData, s.prop=.1^6){
  ### Normalization
  meanVector <- colMeans(inputData)
  meanMat <- matrix(rep(meanVector,nrow(inputData)),nrow = nrow(inputData),
                    byrow = TRUE)
  inputData <- inputData - meanMat
  std <- apply(inputData,2,sd)
  stdMat <- matrix(rep(std,nrow(inputData)),nrow = nrow(inputData),
                   byrow = TRUE)
  inputData <- inputData / stdMat
  Tp <- nrow(inputData)
  os <- ncol(inputData)
  ss <- ncol(inputData)
  xsmooth <- t(inputData)

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
  Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)

  ### Starting from M step for learn_kalman
  w <- abs(matrix(t(Aols),1,ss^2))  ## raw adaptive lasso weights

  XTX <- list()
  xlm <- array()

  for (t in 1:(Tp-1)){
    xlm <- rbind(xlm,kronecker(diag(ss),t(xsmooth[,t])))
  }
  xlm <- xlm[-1,]

  ## x times weights
  xs <- xlm * repmat(w,nrow(xlm),1)
  ylm <- matrix(xsmooth[,2:Tp],ss*(Tp-1),1)

  ww <- t(w) %*% w
  XTX$xtx <- kronecker(diag(ss),gamma1) * ww
  XTX$xty <- matrix(t(beta),ss^2,1) * t(w)

  ## LARS-Lasso
  stopCriterion = list()
  stopCriterion[[1]] <- c("maxKernels",100)
  main <- emlars(yin = ylm,xin = xs,XTX = XTX,regressiontype = "lasso",
                 stopCriterion = stopCriterion)
  sol <- main$history

  ### select min bic
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
  n <- nrow(xs)
  ch <- matrix(NA,1,length(sol)-1)

  for (i in 1:length(sol)-1){
    comb_num <- choose(ss^2,length(sol[[i+1]]$active_set))  # number of combinations
    ch[i] <- 2*log(comb_num)
  }

  ga <- 1
  bic <- log(n)*DF + n*log(RSS/n) + ga*ch
  ks <- which(bic == min(bic))
  temp <- sol[[ks+1]]$beta * w
  Am <- matrix(temp,ss,ss,byrow = TRUE)

  return(list(Amat=Am,norm.data=xsmooth))
}
