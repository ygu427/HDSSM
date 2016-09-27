## GAUSSIAN_PROB Evaluate a multivariate Gaussian density.
## p = gaussian_prob(X, m, C,use_log=FALSE)
## p[i] = N(X[,i], m, C) where C = covariance matrix
## and each COLUMN of x is a datavector
## p = gaussian_prob(X, m, C, use_log=TRUE) returns 
##        log N(X[,i], m, C) (to prevents underflow).
##
## If X has size dxN, then p has size Nx1, where N = number of examples 
##
## Latest Version on May.2016
## Written by Yu Gu

gaussian_prob <- function (x,m,C,use_log=FALSE) {


if (length(m) == 1) {
	x <- matrix (x,1,length(x))
}

d <- nrow(x)
N <- ncol(x)
m <- matrix(m,length(m),1)   ## length(m) == d
M <- m %*% matrix(1,1,N)   ## d * N
denom <- (2*pi)^(d/2)*sqrt(abs(det(C)))
mahal <- rowSums ((t(x-M) %*% solve(C)) * t(x-M)) ## Chris Bregler's trick
if (any(mahal<0)){
	warnings("mahal<0 => C is not psd")
}

if (use_log){
	p <- -0.5 * mahal - log(denom)
} else {
	p <- exp(-0.5 * mahal)/(denom + .Machine$double.eps)
}

return(p)
}