## matbaseERM is the matrix-based ERM algorithm.
##
## matbaseERM(y, initA, initC, initQ, initR, initX, initV) fits
## the parameters which are defined as follows
##   x(t+1) = A*x(t) + w(t),  w ~ N(0, Q),  x(0) ~ N(init_x, init_V)
##   y(t)   = C*x(t) + v(t),  v ~ N(0, R)
## initA is the initial value, estA is the final value, etc.
## y[,t] is the observation vector at time t.
## LL is the "learning curve": a vector of the log lik. values at each iteration.
## LL might go positive, since prob. densities can exceed 1, although this probably
## indicates that something has gone wrong e.g., a variance has collapsed to 0.
##
## There are several optional arguments, that should be passed in the following order.
## matbaseERM(y, initA, initC, initQ, initR, initX, initV, max_iter, diagQ, diagR, ARmode)
## max_iter specifies the maximum number of EM iterations (default 10).
## diagQ=1 specifies that the Q matrix should be diagonal. (Default 0: fixed at true value).
## diagR=1 specifies that the R matrix should also be diagonal. (Default 0: same as above).
## ARmode=1 specifies that initC=I, initR=0. i.e., a Gauss-Markov process. (Default 0).
##
##  Latest Version on Oct.2016
##  Written by Yu Gu

matbaseERM <- function(y,initA,initC,initQ,initR,initx,initV,max_iter=100,
                         diagQ=0,diagR=0,ARmode=0, s.prop=.1^6, ...){

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

  	## R step (matrix-based)
  	Tp1 <- Tp - 1

  	svec <- svd(gamma1)$d
  	if (min(svec)< max(svec)*.1^10 && length(gamma1)!=1) {
  	  ## This is the large p, small n case.
  	  Aols <- beta %*% Rinv(diag(diag(gamma1)), s.prop=s.prop)   ## Equation (18)
  	} else {                            #small p case
  	  Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)   ## This is OLS/MLE
  	}

  	w <- abs(matrix(t(Aols),1,ss^2))  ## raw adaptive lasso weights

  	XTX <- list()
  	xlm <- array()

  	for (t in 1:(Tp-1)){
  		xlm <- rbind(xlm,kronecker(diag(ss),t(estep$xsmooth[,t])))
  	}
  	xlm <- xlm[-1,]

  	## x times weights
  	xs <- xlm * repmat(w,nrow(xlm),1)
  	ylm <- matrix(estep$xsmooth[,2:Tp],ss*(Tp-1),1)

  	ww <- t(w) %*% w
  	XTX$xtx <- kronecker(diag(ss),gamma1) * ww
  	XTX$xty <- matrix(t(beta),ss^2,1) * t(w)

  	## Lasso
  	main <- emlasso(ylm,xs,XTX,s.prop=s.prop,...)
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
  	estA <- matrix(temp,ss,ss,byrow = TRUE)

  	## M step
  	estQ <- as.matrix(estQ)		# Fix Q at true value

  	if (diagQ) {
  		estQ <- (gamma2 - estA %*% t(beta)) / Tp1
  		estQ <- diag(diag(estQ))
  	}

  	if (!ARmode) {
  		estC <- diag(os)		# Let C=I
  		estR <- estR			# Fix R at true value
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

