## Used to generate system matrix and all other required input parameters
## inputdata is the reponse y
## funtype indicates which function will be called later
## It's fullfiled by running R step of ERM algorithm once.
##
## Latest version on May. 2016
## Written by Yu Gu

preprocessing<-function(inputdata,funtype){
  require("pracma")
  xsmooth<- as.matrix(inputdata)
  ss <- nrow(xsmooth)      ## state size
  os <- nrow(xsmooth)      ## observation size
  Tp <- ncol(xsmooth)       ## time point
     
  gamma <- matrix(0,ss,ss)
  beta <- matrix(0,ss,ss)
     
  for (t in 1:Tp){
    gamma <- gamma + xsmooth[,t] %*% t(xsmooth[,t])
    if (t>1) beta <- beta + xsmooth[,t] %*% t(xsmooth[,t-1])
  } 
     
  gamma1 <- gamma - xsmooth[,Tp] %*% t(xsmooth[,Tp])
  Aols <- beta %*% Rinv(gamma1)
  A <- array()

  XTX<-list()

     
  if (funtype == "ssmal"){
    xlm<-matrix()
    xlm<-t(xsmooth[,1:Tp-1])
    for (k in 1:ss){
      w <- abs(t(Aols[k,]))
    	xs <- xlm * repmat(w,nrow(xlm),1)  ## x times weights
    	ylm <- as.matrix(xsmooth[k,2:Tp])
    	ww <- t(w) %*% w
    	XTX$xtx <- gamma1 * ww
    	XTX$xty <- t(beta[k,]*w)
    	
    	## LARS-Lasso
    	stopCriterion = list()
    	stopCriterion[[1]] <- c("maxKernels",100)
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
    	A <- rbind(A,sol[[ks+1]]$beta * w)
    }
    A <- A[-1,]
  } else if (funtype == "learn_kalman"){
    w <- abs(matrix(t(Aols),1,ss^2)) ## raw adaptive lasso weights. size:1*ss^2
    
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
    #A <- t(matrix(temp,ss,ss))
    A <- matrix(temp,ss,ss,byrow = TRUE)
  }
  
  C<-diag(os)
  
  ## Default Q & R matrix
  Q <- diag(ss)
  dev <- 0.2
  Q <- Q + abs(Q) * matrix(runif(ss*ss,min=-dev,max=dev),ss,ss)
  R <- 0.1 * diag(os)
      
  ## Output
  return(list(A=A,C=C,Q=Q,R=R))
}
