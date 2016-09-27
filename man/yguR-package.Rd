\name{HDSSM-package}
\alias{HDSSM-package}
\alias{HDSSM}
\alias{ssmal}
\alias{learn_kalman}
\docType{package}
\title{
High-dimensional linear SSM for dynamic GRN
}
\description{
A high-dimensional linear State Space Model(SSM) with a new Expectation-Regularization-Maximization(ERM) algorithm to construct the dynamic Gene Regularization Network(GRN). The new ERM algorithm employs the idea of the adaptive LASSO-based variable selection method to preserve the sparsity property of GRN.
}
\usage{
  ssmal(data,A,C,Q,R,initx,initV,max_iter=10,diagQ=FALSE,diagR=FALSE,ARmode=FALSE,...)
}
\usage{
  learn_kalman(data,A,C,Q,R,initx,initV,max_iter=10,diagQ=FALSE,diagR=FALSE,ARmode=FALSE,...)
}
\arguments{
  \item{data}{data[,t,l] the observation vector at time t for sequence l}
  \item{A}{the system matrix}
  \item{C}{the observation matrix}
  \item{Q}{the variance for normally distributed system noise}
  \item{R}{the variance for normally distributed measurement noise}
  \item{initx}{mean value vector for initial state x0}
  \item{initV}{covariance matrix for initial state x0}
  \item{max_iter}{specifies the maximum number of EM iterations (default 10)}
  \item{diagQ}{boolean value. 1 specifies that the Q matrix should be diagonal. Default value is 0,indicating fixed at true value.}
  \item{diagR}{boolean value. 1 specifies that the R matrix should be diagonal. Default value is 0,indicating fixed at true value.}
  \item{ARmode}{boolean value. 1 specifies that C=I, R=0. i.e., a Gauss-Markov process. (Default 0).}
  \item{...}{more optional arguments}
}
\details{
  \tabular{ll}{
    Package: \tab HDSSM\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2016-07-26\cr
    License: \tab GPL-3\cr
  }
To use this package, you need to update the R to at least version 3.1 and install the package "pracma". The most important functions are ssmal and learn_kalman, which are row-based ERM algorithm and matrix-based ERM algorithm respectively.

Keywords and phrases: adaptive LASSO, dynamic Gene Regulatory Network, L1 regularization, State Space Model, VAR model
}
\value{
  \item{A}{the estimated high-dimensional sparse system matrix A}
  \item{C}{the observation matrix}
  \item{Q}{the estimated variance for normally distributed system noise}
  \item{R}{the estimated variance for normally distributed measurement noise}
  \item{initX}{the mean value for state vector}
  \item{initV}{the covariance matrix for state vector}
  \item{LL}{the log likelihood vector}
  \item{xcurve}{smoothed real-valued hidden state variable vector}
  \item{bic}{BIC}
  \item{num_iter}{the number of iteration}
}
\author{
Yu Gu

Maintainer: Yu Gu <yu_gu@urmc.rochester.edu>
}
\references{
Chen, J., & Chen, Z. (2008). Extended Bayesian information criteria for model selection with large model spaces. Biometrika, 95(3), 759-771.

Green, P. J. (1990). On use of the EM for penalized likelihood estimation. Journal of the Royal Statistical Society. Series B (Methodological), 443-452.

Harrison, J., & West, M. (1999). Bayesian Forecasting & Dynamic Models. Springer.

Zou, H. (2006). The adaptive lasso and its oracle properties. Journal of the American statistical association, 101(476), 1418-1429.
}

\keyword{package}
