## learn_kalman is the matrix-based ERM algorithm. 
##
## [A, C, Q, R, INITX, INITV, LL] = LEARN_KALMAN(DATA, A0, C0, Q0, R0, INITX0, INITV0) fits
## the parameters which are defined as follows
##   x(t+1) = A*x(t) + w(t),  w ~ N(0, Q),  x(0) ~ N(init_x, init_V)
##   y(t)   = C*x(t) + v(t),  v ~ N(0, R)
## A0 is the initial value, A is the final value, etc.
## DATA[,t,l] is the observation vector at time t for sequence l. If the sequences are of
## different lengths, you can pass in a cell array.
## LL is the "learning curve": a vector of the log lik. values at each iteration.
## LL might go positive, since prob. densities can exceed 1, although this probably
## indicates that something has gone wrong e.g., a variance has collapsed to 0.
##
## There are several optional arguments, that should be passed in the following order.
## LEARN_KALMAN(DATA, A0, C0, Q0, R0, INITX0, INITV0, MAX_ITER, DIAGQ, DIAGR, ARmode)
## MAX_ITER specifies the maximum number of EM iterations (default 10).
## DIAGQ=1 specifies that the Q matrix should be diagonal. (Default 0: fixed at true value).
## DIAGR=1 specifies that the R matrix should also be diagonal. (Default 0: same as above).
## ARMODE=1 specifies that C=I, R=0. i.e., a Gauss-Markov process. (Default 0).
## 
##  Latest Version on May.2016
##  Written by Yu Gu

learn_kalman <- function(data,A,C,Q,R,initx,initV,max_iter=10,
                         diagQ=0,diagR=0,ARmode=0, s.prop=.1^6, ...){
  
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
  
  ss <- nrow(A)  ## state size
  os <- nrow(C)  ## observation size
  
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
  	
  	## R step (matrix-based)
  	Tsum1 <- Tsum -N
  	Aols <- beta %*% Rinv(gamma1, s.prop=s.prop)   ## This is OLS/MLE
  	w <- abs(matrix(t(Aols),1,ss^2))  ## raw adaptive lasso weights
  
  	XTX <- list()
  	Tp <- ncol(estep$xsmooth)
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
  
  	## LARS-Lasso
  	stopCriterion = list()
  	stopCriterion[[1]] <- c("manKernels",100)
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
  
  	## M step
  	Q <- as.matrix(Q)		# Fix Q at true value
  	
  	if (diagQ) {
  		Q <- (gamma2 - A %*% t(beta)) / Tsum1
  		Q <- diag(diag(Q))
  	}
  
  	if (!ARmode) {
  		C <- diag(os)		# Let C=I
  		R <- R			# Fix R at true value
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

