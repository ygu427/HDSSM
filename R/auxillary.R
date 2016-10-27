## All non-essential functions should be moved to this file

## Rinv stands for ridge-regression inverse. This is the SVD version.
Rinv <- function(Mat, s.prop=.1^6){
  ## Sometimes, Mat may be totally empty!!
  if (length(Mat)==0) {
    warning("Mat is an empty matrix!"); return(Mat)
  } else {
    ss <- svd(Mat); svec <- ss[["d"]]; U <- ss[["u"]]; V <- ss[["v"]]
    if( length(svec)==1){
      ## If Mat is just a number, svec will have length 1. We simply
      ## report the truncated inverse of this number. Note that in this
      ## case, we should compare the singular value (which is just
      ## abs(Mat)) to s.prop, not to sum(svec)*s.prop.
      truncated.svec <- max(svec, s.prop)
      inv.mat <- V * 1/truncated.svec * U
    } else {
      truncated.svec <- pmax(svec, sum(svec)*s.prop)
      inv.mat <- V %*% diag(1/truncated.svec) %*% t(U)
    }
    return(inv.mat)
  }
}

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
  mahal <- rowSums ((t(x-M) %*% Rinv(C)) * t(x-M)) ## Chris Bregler's trick
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


#-------------------------------------------------------------------------------------#
# Copyright (c) 2010-2011 Takeshi Arabiki                                             #
# Licensed under the terms of the MIT License (see LICENSE.txt)                       #
#-------------------------------------------------------------------------------------#



#-------------------------------------------------------------------------------------#
# Description:                                                                        #
#    Replicate and tile array                                                         #
#                                                                                     #
# Usage:                                                                              #
#    repmat(A, m, n = m)                                                              #
#                                                                                     #
# Arguments:                                                                          #
#    A:      a vector, factor, matrix, array, list, or data.frame.                    #
#    m:      a nonnegative integer or nonnegative integer vector                      #
#            to specify the number of replicatoin about each dimension.               #
#    n:      an integer ignored if the number of dimensions of m is larger than 1.    #
#-------------------------------------------------------------------------------------#
repmat <- function(A, m, n = m) {
  if (missing(m))
    stop("Requires at least 2 inputs.")

  if (is.vector(A))
    A <- matrix(A, 1)

  if (length(m) == 1)
    m[2] <- n

  if (any(as.integer(m) != m) || any(m <= 0)) {
    stop("'m' should be a nonnegative integer or a nonnegative integer vector!")
  }

  d <- dim(A)
  if (length(d) <= 2) {
    nr <- d[1L]
    nc <- d[2L]
    # kronecker(array(1, m), A) is slower
    tmpA <- matrix(t(A), nrow = nr * m[1], ncol = nc, byrow = TRUE)
    A <- matrix(tmpA, nrow = nr * m[1], ncol = nc * m[2])
    if (length(m) > 2) {
      A <- array(tmpA, c(nr * m[1], nc * m[2], m[-(1:2)]))
    }
  } else {
    A <- kronecker(array(1, m), A)
  }

  return(A)
}
