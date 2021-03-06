\name{matbaseERM}
\alias{matbaseERM}
\title{
matrix-based Expectation-Regularization-Maximization(ERM) algorithm
}
\description{
A high-dimensional linear State Space Model(SSM) with a new Expectation-Regularization-Maximization(ERM) algorithm to construct the dynamic Gene Regularization Network(GRN). The new ERM algorithm employs the idea of the adaptive LASSO-based variable selection method to preserve the sparsity property of GRN.
}
\usage{
matbaseERM(y, initA, initC, initQ, initR, initx, initV, max_iter = 100, diagQ = 0, diagR = 0, ARmode = 0, s.prop = 0.1^6, ...)
}

\arguments{
  \item{y}{
y[,t] the observation vector at time t
}
  \item{initA}{
the initial system matrix
}
  \item{initC}{
the initial observation matrix
}
  \item{initQ}{
the initial variance for normally distributed system noise
}
  \item{initR}{
the initial variance for normally distributed measurement noise
}
  \item{initx}{
mean value vector for initial state x0
}
  \item{initV}{
covariance matrix for initial state x0
}
  \item{max_iter}{
specifies the maximum number of EM iterations (default 100)
}
  \item{diagQ}{
boolean value. 1 specifies that the Q matrix should be diagonal. Default value is 0,indicating fixed at true value.
}
  \item{diagR}{
boolean value. 1 specifies that the R matrix should be diagonal. Default value is 0,indicating fixed at true value.
}
  \item{ARmode}{
boolean value. 1 specifies that C=I, R=0. i.e., a Gauss-Markov process. (Default 0).
}
  \item{s.prop}{
%%     ~~Describe \code{s.prop} here~~
}
  \item{\dots}{
more optional arguments
}
}

\value{
  \item{estA}{the estimated high-dimensional sparse system matrix A}
  \item{estC}{the estimated observation matrix}
  \item{estQ}{the estimated variance for normally distributed system noise}
  \item{estR}{the estimated variance for normally distributed measurement noise}
  \item{estX}{the mean value for state vector}
  \item{estV}{the covariance matrix for state vector}
  \item{LL}{the log likelihood vector}
  \item{xcurve}{smoothed real-valued hidden state variable vector}
  \item{bic}{BIC}
  \item{num_iter}{the number of iteration has been processed}
}
\references{
Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Green, P. J. (1990). On use of the EM for penalized likelihood estimation. Journal of the Royal Statistical Society. Series B (Methodological), 443-452.

Harrison, J., & West, M. (1999). Bayesian Forecasting & Dynamic Models. Springer.

Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American statistical association, 101(476), 1418-1429.
}
\author{
Yu Gu

Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}
