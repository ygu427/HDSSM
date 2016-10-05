# This program implements "LARS" algorithm and its lasso modification
# introduced by Efron et. al. 2003.  Read the paper to understand codes
# of this function.  Each line of this file has corresponding equation
# number in Efron et. al. 2003 for reader's convenience.
#
# *** CAUTION
# history[[1]]$mu_OLS contains original 'yin' to provide convenience in
# writing user defined stop criterion function. Actual mu_OLS of the first
# step should be just mean of yin which is a simple array of copy of
# history[[1]]$mu. So, if history[[1]]$mu_OLS contains that information, it is
# redundant. In fact, this is also contained in history[[1]]$b, which is a
# bias of the output. Therefore, to provide more information to user who
# want to write his/her own stop criterion function, history[[1]]$mu_OLS
# contains yin.
#
#
# Example: moderate size x.
#     stopCrioterion <- list()
#     stopCrioterion[[1]] = c('maxKernels',100)
#
# Note:
#     Users can add any kind of stop criterion by editing
#     the corresponding portion of this file.  See the code
#     for existing examples.
#
# Note 2:
#     This file does not implement routine for missing data.
#
# Latest Version on May. 2016
# Wriiten by Yu Gu


emlars <- function (yin,xin,XTX,regressiontype,stopCriterion,regularization=0,
                    trace=FALSE,quiet=FALSE, s.prop=.1^6){
  resolution_of_lars <- 1e-12
  ## Tikhonov regularization factor (or the ridge regression factor)
  regularization_factor <- 1e-11

  ## Notify the user what cause "stop".
  stopReason <- list()

  # Check Parameters
  if (length(yin) == 0 | length(xin) == 0){
    warning("\n Input or Output has zero length.\n")
    stopReason[[1]] <- "Parameter error"
    stopReason[[2]] <- 0
    stop(stopReason)
  }
  if (nrow(yin) != nrow(xin)){
    warning("\n Size of y does not match to that of x.\n")
    stopReason[[1]] <-"Parameter error"
    stopReason[[2]] <- 0
    stop(stopReason)
  }
  if ((regressiontype!="lasso")&(regressiontype!="lars")&(regressiontype!="forward_stepwise")){
    warning("\n Unknown type of regression.\n")
    stopReason[[1]] <- "Parameter error"
    stopReason[[2]] <- 0
    stop(stopReason)
  }
  if (regressiontype=="forward_stepwise"){
    warning("\n Forward_stepwise is not implemented.\n")
    stopReason[[1]] <- "Parameter error"
    stopReason[[2]] <- 0
    stop(stopReason)
  }


  if (!missing(regularization)){
    if (!is.na(regularization)) {
      regularization <-10
    }
  }

  if (!missing(quiet)) {
    if (quiet) {trace <- FALSE}
  }

  ##########################################################
  ###                                                    ###
  ###                Data preparation                    ###
  ###                                                    ###
  ##########################################################

  x <- xin   # nomalized xin
  n <- nrow(x)  # number of samples
  p <- ncol(x)  # number of predictors

  mx <- matrix(0,1,p)  #XTX.mx  # mean xin
  sx <- matrix(1,1,p)  #XTX.sx  # length of each column of xin

  all_candidate <- seq_len(p)    # indices for all possible columns
  xtx <- XTX$xtx  # x'x matrix
  xty <- XTX$xty  # x'y matrix

  my <- 0   ## mean of y, i.e., mean(y)
  y <- yin  ## y <- yin - mean(y)

  # Determine the maximum number of kernels
  existMaxKernels <- FALSE
  existMSE <- FALSE

  for(is in 1:length(stopCriterion)){
    if (stopCriterion[[is]][1]=="maxKernels"){
      existMaxKernels <-TRUE
      stopCriterion[[is]][2] <- min(qr(xin)$rank,length(all_candidate))
      if (as.numeric(stopCriterion[[is]][2])<1){
        warning("\n Max Kernel is less than 1. It must be larger than 0.\n")
        stopCriterion[[is]][2] <-1
      }
    }

    if (stopCriterion[[is]][1]=="MSE"){
      existMSE <-TRUE
      if (as.numeric(stopCriterion[[is]][2])<1e-10){
        warning("\n Maximum MSE is too small. Automatically set to 1.0e-10\n")
        stopCriterion[[is]][2]<-1e-10
      }
    }
  }

  if (!existMaxKernels) {
    is <- length(stopCriterion)
    stopCriterion[[is+1]]<-c("maxKernels",min(qr(xin)$rank,length(all_candidate)))
  }

  if (!existMSE){
    is <- length(stopCriterion)
    stopCriterion[[is+1]]<-c("MSE",1e-10)
  }

  ###########################################################
  ###                                                     ###
  ###                   Initialization                    ###
  ###                                                     ###
  ###########################################################

  active <- NA    # active set
  inactive <- all_candidate # inactive set

  mu_a <- matrix(0,n,1)   # current estimate
  mu_a_plus <- mu_a     # next estimate
  mu_a_OLS <- mu_a      # OLS estimate
  beta <- matrix(0,1,ncol(x))
  beta_new <- beta
  beta_OLS <- beta

  ### Create the list for storing history
  ### include: active_set,add,dropmatrix,beta_norm,beta,df,rss,b,
  ###          mu,beta_OLS_norm,beta_OLS,b_OLS,mu_OLS,MSE,R_square,
  ###          resolution_warning
  history<-list()
  history[[1]]<-list()

  history[[1]]$active_set <- NA
  history[[1]]$add <- NA
  history[[1]]$dropmatrix <- array(NA)
  history[[1]]$beta_norm <- array(NA)
  history[[1]]$beta <- array(NA)
  history[[1]]$df <- NA
  history[[1]]$rss <- NA
  history[[1]]$b <- my
  history[[1]]$mu <- my
  history[[1]]$beta_OLS_norm <- array(NA)
  history[[1]]$beta_OLS <- array(NA)
  history[[1]]$b_OLS <- my
  history[[1]]$mu_OLS <- my * array(1,dim=dim(y))
  history[[1]]$MSE <- colSums(y^2)/length(y)
  history[[1]]$R_square <-0
  history[[1]]$resolution_warning <- array(NA)

  if (diag(var(yin))==0){
    stopReason[[1]]<-"zeroVarY"
    stopReason[[2]]<-diag(var(yin))
    stop(stopReason)
  }

  ########################################################
  ###                                                  ###
  ###                  Main Loop                       ###
  ###                                                  ###
  ########################################################

  c<-NA
  C_max <- NA
  C_max_ind <- vector()
  C_max_ind_pl <- vector()
  dropmatrix <- array(NA)   #used for "lasso"
  k <- 1

  while (TRUE) {

    ############################################
    ###           Exit Criterions            ###
    ############################################

    if (!missing(stopCriterion)) {

      ## If there is any stop criterion
      ## Any of these is satisfied, algorithm stops.
      ## Other criterion: make your own cost function

      for(is in 1:length(stopCriterion)){
        if(stopCriterion[[is]][1]=="maxDrops"){
          # Default Criterions
          # Maximum number of consecutive drops or maximum number of drops within a window
          # drop_window: number of backward history to be considered
          # Zero means all the whole history
          drop_window <- as.numeric(stopCriterion[[is]][2])
          if (drop_window == 0){
            drop_window <-k
          }

          # drop_n: number of maximum drops within the window.
          drop_n<-min(drop_window,as.numeric(stopCriterion[[is]][3]))
          drop_vector<-array(NA)

          for (z in max(k-drop_window+1,1):k){
            drop_vector<-cbind(drop_vector,history[[k]]$dropmatrix)
          }

          if (length(drop_vector)>=drop_n){
            stopReason[[1]]<-"maxDrops"
            stopReason[[2]]<-drop_n
            break
          }
        }# end checking maximum number of drops

        ## # Maximum number of kernels.
        ## if (stopCriterion[[is]][1]=="maxKernels"){
        ##   if (length(active)>=min(as.numeric(stopCriterion[[is]][2]),nrow(xin)-1,length(all_candidate))){
        ##     stopReason[[1]]<-"maxKernels"
        ##     stopReason[[2]]<-length(active)
        ##     break
        ##   }
        ## } # end checking maximum number of kernels

        # Maximum number of iterations
        if (stopCriterion[[is]][1]=="maxIterations") {
          if (k>= as.numeric(stopCriterion[[is]][2])){
            stopReason[[1]]<-"maxIterations"
            stopReason[[2]]<-k
            break
          }
        }# end checking maximum number of iterations

        # MSE
        if (stopCriterion[[is]][1]=="MSE") {
          if (history[[k]]$MSE<=as.numeric(stopCriterion[[is]][2])){
            stopReason[[1]]<-"MSE"
            stopReason[[2]]<-history[[k]]$MSE
            break
          }
        }
      }
    }# end of stop criterion checking

    if (length(stopReason)>0){
      ## if there is any reason of stopping the algorithm, exit loop
      ## stop(stopReason)
      break
    }

    #################################################
    ###             LARS Algorithm                ###
    #################################################

    c <- xty - t(x) %*% mu_a  # x'(y-mu_a) the current correlation
    c1 <- abs(c[inactive])
    C_max <- max(c1)
    C_max_ind <- inactive[which(c1==C_max)]

    # Because of machine limit,there can be multiple new predictor.
    # This dramatically improves the overall precision of the result,
    # and speeds up the whole process.
    C_max_ind_pl <- c1 > (C_max - resolution_of_lars)
    C_max_ind_pl <- inactive[C_max_ind_pl]
    active <- sort(union(active,C_max_ind_pl), na.last = NA)
    inactive <- setdiff(all_candidate,active)

    if (regressiontype == "lasso"){
      ## If there is a drop, that must have the maximum correlation in inactive set.
      if (!is.na(dropmatrix)) {
        if (length(which(dropmatrix==C_max_ind))==0) {
          if ((!quiet) & (trace>=0)){
            warning("\n Dropped item and index of maximum correlation is not the same.
                    But it is being ignored here...\n")
          }
          active[which(active==C_max_ind)] <- NA
          sort(active, na.last = NA) ## eliminate NA value
        }
        C_max_ind <- NA
        C_max_ind_pl <- NA
      }
      active <- setdiff(active,dropmatrix)
      inactive <- sort(union(inactive,dropmatrix))
    }

    s <- sign(c[active])
    xa <- as.matrix(x[,active]) * repmat(t(s),n,1)
    ga <- xtx[active,active] * (s %*% t(s))

    if (regularization >2 ){
      # This routine will make the test below
      ga <- ga + diag(nrow(ga)) * regularization_factor
    }

    invga <- Rinv(ga, s.prop=s.prop)

    aa <- sum(invga)^(-0.5)
    wa <- aa * rowSums(invga)
    ua <- xa %*% wa   # equiangular vector

    ## test whether Xa'Ua=Aa1a and ||Ua||^2=1
    test_1 <- t(xa) %*% ua
    test_2 <- matrix(aa,length(active),1)
    test_1_2 <- sum(abs(test_1 - test_2))
    test_3 <- abs(norm(ua,type="2")-1)

    k <- k+1
    history[[k]] <- list()
    history[[k]]$resolution_warning <- 0
    thresh.tmp <- resolution_of_lars*100
    if ((test_1_2 > thresh.tmp)|(test_3 > thresh.tmp)) {
      if (regularization <=2 ) {
        if ((!quiet) & trace){
          warning("\n Eq 2.7 test failure. \n")
        }
        regularization <- regularization + 1
        if (regularization > 2) {
          if((!quiet) & trace){
            warning("\n Lots of Eq 2.7 test failure. Regularization will be applied from now on.\n")
          }
        }
      }
      history[[k]]$resolution_warning <- 1
    }

    a <- t(x) %*% ua
    tmp_1 <- (C_max - c[inactive])/(aa - a[inactive])
    tmp_2 <- (C_max + c[inactive])/(aa + a[inactive])
    tmp_3 <- cbind(tmp_1,tmp_2)
    tmp <- as.matrix(tmp_3[which(tmp_3>0)])

    if (length(tmp)==0) {
      # if this is the last step (i.e. length(active)==maxKernels)
      gamma = C_max/aa
    } else {
      gamma <- min(tmp)
    }

    d <- matrix(0,1,p)
    d[active] <- s * wa

    if (length(which(d[active]==0))) {
      warning("\n Something wrong with vector d: eq 3.4. \n")
    }

    tmp1 <- matrix(0,1,p)
    tmp1[active] <- -(beta[active] / d[active])   ##gamma_j
    tmp2 <- tmp1[which(tmp1>0)]     ##{gamma_j:gamma_j>0}

    dropmatrix <- array()
    gamma_tilde <- Inf
    if (length(tmp2)!=0) {
      if (gamma >= min(tmp2)){
        gamma_tilde <- min(tmp2)
        dropmatrix <- which(tmp1 == gamma_tilde)
      }
    }

    if (regressiontype == "lars") {
      mu_a_plus <- mu_a + gamma %*% ua
      beta_new <- beta + gamma %*% d
      dropmatrix <- array()
    } else if (regressiontype == "lasso") {
      mu_a_plus <- mu_a + min(gamma,gamma_tilde)*ua
      beta_new <- beta + min(gamma,gamma_tilde)*d
      active <- setdiff(active,dropmatrix)
      inactive <- setdiff(all_candidate,active)
      beta_new[dropmatrix] <- 0
    } else if (regressiontype == "forward_stepwise") {
      dropmatrix <- array()
      stop("Forward.stepwise has not been implemented yet.")
    }

    mu_a_OLS <- mu_a + C_max/aa * ua
    beta_OLS <- beta + C_max/aa * d
    MSE <- sum((y - mu_a_OLS)^2)/length(y)

    #########################################################
    ###               Update and Save                     ###
    #########################################################
    mu_a <- mu_a_plus
    beta <- beta_new

    res <- y - x %*% t(beta)
    rss <- sum(res^2)

    # history with scale correction
    history[[k]]$active_set <- active
    history[[k]]$df <- length(active)
    history[[k]]$dropmatrix <- dropmatrix
    history[[k]]$add <- C_max_ind_pl
    history[[k]]$beta <- beta
    history[[k]]$rss <- rss
    history[[k]]$b <- my - sum(mx/sx*beta)
    history[[k]]$mu <- xin %*% t(beta/sx) + history[[k]]$b
    history[[k]]$beta_OLS_norm <- beta_OLS[active]
    history[[k]]$beta_OLS <- beta_OLS[active]/sx[active]
    history[[k]]$b_OLS <- my - sum(mx/sx*beta_OLS)
    history[[k]]$mu_OLS <- xin %*% t(beta_OLS/sx) + history[[k]]$b_OLS
    history[[k]]$MSE <- MSE
    history[[k]]$R_square <- 1 - var(yin - history[[k]]$mu_OLS)/var(yin)

    ## exit if exact matching is achieved.
    if (abs(C_max/aa - min(gamma,gamma_tilde))<resolution_of_lars) {
      stopReason[[1]] <- "ExactMatching"
      stopReason[[2]] <- 0
      break
    }
  }

  ##Output
  return(list(history = history,stopReason = stopReason))
}
