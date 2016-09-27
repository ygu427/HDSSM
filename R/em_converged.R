##  Check Whether EM has been converged
##  EM is converged when the slope of the log-likelihood function falls below "threshold"
##  i.e.  |f(t)-f(t-1)| / avg < threshold
##  
##  If MAP estimation is processed (using priors), the likelihood can decrease, even though
##  the mode of the posterior is increasing
##  Input:
##    loglik -- current log-likelihood
##    pre_loglik -- previous log-likelihood
##    threshold -- determine whether EM is coverged. The default value is 1e-4.
##  Output:
##    converged -- boolean value to indicate whether converged or not
##    decrease -- boolean value to indicate whether log-likelihood decreases or not
##
##  Latest Verstion on May. 2015
##  Wriiten by Yu Gu

em_converged <- function(loglik,pre_loglik,threshold=1e-4,check_increased=TRUE){

# Vaviable declaration and Initialization
converged <- FALSE      ## Boolean value: whether converged
decrease <- FALSE       ## Boolean value: whether the likeloglihood is decreasing
delta_loglik<-NA        ## Difference between current and previous loglikelihood
avg_loglik<-NA          ## Average value of current and previous loglikelihood


## Determine whether trace the loglikelihood's variation, decrease or increase
if (check_increased)
{
	if((loglik - pre_loglik) < -1e-2) {
		sprintf("likelihood decreased from %6.4f to %6.4f !",pre_loglik,loglik)
		decrease <- TRUE
		converged <- FALSE
	} else {
		delta_loglik <- abs(loglik-pre_loglik)
		avg_loglik <- (abs(loglik)+abs(pre_loglik)+.Machine$double.eps)/2
  }
  
  temp <- delta_loglik/avg_loglik
  if(is.na(temp)) {temp <- Inf}
	if(temp < threshold) converged <- TRUE
} else {
  ## won't check the variation of likelihood
  delta_loglik <- abs(loglik-pre_loglik)
  avg_loglik <- (abs(loglik)+abs(pre_loglik)+.Machine$double.eps)/2

  temp<-delta_loglik/avg_loglik
  if(is.na(temp)) {temp <- Inf}
  if(temp < threshold) converged <- TRUE
}


## Output
if (check_increased) {
  return(list(converged = converged, decrease = decrease))
} else {
  return(list(converged = converged))
}

}