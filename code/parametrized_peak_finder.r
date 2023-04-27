# parametrixed_peak_finder.R
# Author: Ariella Fuzaylov
#     Secondary Author: Natalie Dupont
# Date Created: 2022-11-21
# Date last modified: 2023-04-20
# BIOL 490, Pedersen Lab, Concordia University 
# Description: This code creates useable function to perform the GAM derivative method.

library(dplyr)
library(mgcv)
library(tibble)

#Using functions as written from testing_maxima.R so that updating one function updates it here too.
  source("code/testing_maxima.R")
    #for Rmd
    #source("testing_maxima.R")


#simulating a model 
sim_model<- function(n=100, sigma=1, a=1, b=1, mid=2, low_lim=0.1, high_lim=4, step_size=0.1){
  
  # data frame with x, values, random y values, and the true curve
  training_data = tibble(x = seq(low_lim,high_lim, length=n),y= rnorm(n,fit_func(x,a,b,mid),sigma),
                             true_val = fit_func(x,a,b,mid) )
  
  # Extracts the actual minimum for the function (or the left boundary if the min
  # is outside the data range)
  true_min = optim(par = 1, fit_func, a=a,b=b,mid=mid,lower = low_lim,upper=high_lim,method = "L-BFGS-B")$par[1]
  
  
  "
  The fitted model, using a 20 basis function thin plate spline smoother with REML
  fitting criteria. m=3 specifies that the model should penalize squared third derivatives.
  This is important as if m=2 (the default) then prior simulations from the fit are too wiggly, and
  end up with too wide a range of 2nd derivatives
  "
  mod = gam(y~s(x, bs= "tp", k=20,m = 3), data=training_data,method="REML")
  
  params = list( "n" = n, "sigma" = sigma, "a" = a, "b"= b, "low_lim" = low_lim, "high_lim" = high_lim, "step_size" = 0.1)
  
  result <- list("training_data" = training_data, "true_min" = true_min, "mod"=mod, "params" = params)
  return(result)
}


"
note: can make more customization easily, keeping it simple for now
args:
      n_sims: number of simulations for posterior model
returns a list of:
        mod: gam model created from data set
"
#Taking a gam model object and formatting it for parametrized peak finder

#Alternate function to create an object passable to subsequent functions from external data
create_peakfinder_input <- function(mod, data, low_lim=NULL, high_lim=NULL, step_size=0.001){
  if(is.null(low_lim)) {
    low_lim <- min(mod$model$x, na.rm = TRUE)
  }
  if(is.null(high_lim)) {
    high_lim <- max(mod$model$x, na.rm = TRUE)
  }
  step_size <- step_size
  n <- length(mod$residuals)

  output <- list()
  
  #Training data is just the input data as a dataframe of x, y (and potentially true_vals)
  output$training_data <- data  
  output$mod <- mod
  output$params <- list("n"=n, 
                        "low_lim" = low_lim,
                        "high_lim" = high_lim,
                        "step_size" = step_size)
  return(output) #output follows the same structure as the output from sim_model()
  }

model_data <- function(true_model_summary, n_sims){
  low_lim = true_model_summary$params$low_lim
  high_lim = true_model_summary$params$high_lim
  step_size = true_model_summary$params$step_size
  mod = true_model_summary$mod
  #You need the multivariate normal RNG from the MASS package
  mvrnorm = MASS::mvrnorm
  
  # The test data, with one x per step unit across the range.
  test_data = tibble(x=seq(low_lim, high_lim,by= step_size))
  
  #Simulate new functions from the posterior distribution of functions, using the
  #test data and 500 simulations
  mod_coef = coef(mod) # mean values for all (intercept and k-1 basis functions) basis functions
  mod_vcov =vcov(mod)  # posterior variance-covariance matrix (20x20 matrix)s
  mod_sims = mvrnorm(n_sims, mod_coef,mod_vcov) #random parameter draws. Each column is n_sims number of  random draws of 1 of the basis function coefficients
  test_lp = predict.gam(mod,newdata = test_data,type = "lpmatrix") #value of each basis function at each step of test_data x
  test_sims = test_lp %*% t(mod_sims) #random parameters times basis functions. Actual n_sims number of simulated curves
  sim_summary <- list("gam"=mod, "coef" = mod_coef, "vcov" = mod_vcov, "gam_sims" = mod_sims)
  result <- list("test_data" = test_data, "test_sims"=test_sims, "test_lp" = test_lp, "sim_summary" = sim_summary)
  return(result)
}

find_peaks = function(deriv_1_bounds, 
                      deriv_2_bounds, 
                      test = c("two-deriv","crossing")){
  
  stopifnot(nrow(deriv_1_bounds)==nrow(deriv_2_bounds))
  stopifnot(ncol(deriv_1_bounds)==2&ncol(deriv_2_bounds)==2)
  
  n_vals = nrow(deriv_1_bounds)
  
  # calculate the sign of the derivatives. If both ends of the interval are pos/neg
  # at a point, this will return a value of +1/-1 for that point. It will be zero
  # iff the CI overlaps zero at that point.
  deriv1_sign = (sign(deriv_1_bounds[,1]) + sign(deriv_1_bounds[,2]))/2
  deriv2_sign  = (sign(deriv_2_bounds[,1]) + sign(deriv_2_bounds[,2]))/2
  
  #Derivative series start with na values
  #This removes that issue (to prevent a run of na values)
  deriv1_sign[c(1,n_vals)] = deriv1_sign[c(2, n_vals-1)] 
  deriv2_sign[c(1,n_vals)] = deriv2_sign[c(2, n_vals-1)] 
  
  if(test[1]=="two-deriv"){
    is_candidate = (deriv1_sign==0)&(deriv2_sign!=0)
  }else if(test[1] =="crossing"){
    deriv1_runs = rle(deriv1_sign)
    n_runs = length(deriv1_runs$lengths)
    deriv1_run_start = c(1, cumsum(deriv1_runs$lengths)[-n_runs]+1)
    deriv1_run_end = cumsum(deriv1_runs$lengths)
    deriv1_run_sign = deriv1_runs$values
    is_candidate = rep(FALSE, times=nrow(deriv_1_bounds))
    for(i in 1:n_runs){
      if(i>1&i<n_runs & deriv1_run_sign[i]==0){
        if(deriv1_sign[deriv1_run_start[i-1]]*deriv1_sign[deriv1_run_end[i+1]]==-1){
          is_candidate[deriv1_run_start[i]:deriv1_run_end[i]] = TRUE
        }
      }
    }
  }else {
    stop("this is not an implemented test")
  }
  return(is_candidate)
}

peak_finder <- function (model_coefs, true_model_summary, test){
  test_sims = model_coefs$test_sims
  step_size = true_model_summary$params$step_size
  
  x_vals <- seq(from=true_model_summary$params$low_lim, to=true_model_summary$params$high_lim, by=true_model_summary$params$step_size)
  #Calculates estimated first and second derivatives 
  test_1st_deriv = apply(test_sims,MARGIN = 2,calc_1st_deriv, delta= step_size)
  test_2nd_deriv = apply(test_sims,MARGIN = 2,calc_2nd_deriv, delta= step_size)
  
  # 95% confidence intervals for the function, 1st, and 2nd derivatives
  test_CI = t(apply(test_sims,
                    MARGIN = 1,
                    FUN = quantile,
                    probs=c(0.025,0.5,0.975),
                    na.rm=T))
  test_1st_deriv_CI = t(apply(test_1st_deriv ,
                              MARGIN = 1,
                              FUN = quantile,
                              probs=c(0.025,0.5,0.975),
                              na.rm=T))
  test_2nd_deriv_CI = t(apply(test_2nd_deriv ,
                              MARGIN = 1,
                              FUN = quantile,
                              probs=c(0.025,0.5, 0.975),
                              na.rm=T))
  
  candidate_peaks = as.vector(find_peaks(test_1st_deriv_CI[,c(1,3)], 
                                                   test_2nd_deriv_CI[,c(1,3)], test = test))
  
  candidate_peaks = ifelse(is.na(candidate_peaks), F, candidate_peaks)
  CI <- list("test_CI" =  test_CI, "test_1st_deriv_CI" = test_1st_deriv_CI, "test_2nd_deriv_CI" = test_2nd_deriv_CI)
  result <- list("x_val"=x_vals, "test_1st_deriv" = test_1st_deriv, "test_2nd_deriv" = test_2nd_deriv, "CI" = CI, "candidate_peaks" = candidate_peaks)
}

plot_peaks_against_true_vals <- function(possible_peaks, model_coefs, true_model_summary, plot_2nd_deriv=TRUE, poisson=FALSE){
  training_data=true_model_summary$training_data
  test_data=model_coefs$test_data
  test_CI= possible_peaks$CI$test_CI
  candidate_peaks = possible_peaks$candidate_peaks
  test_CI_only_CI <- test_CI
  test_CI_only_CI[!candidate_peaks,] <- NA
  #true_min=true_model_summary$true_min
  if(plot_2nd_deriv){
    par(mfrow=c(3,1))}
  else{par(mfrow=c(2,1))}
  
  # Plot of raw data and model fit, with true function in blue and 
  # estimated minima in red. Vertical blue line is the true minimum
  # par(mar=c(1,1,1,1))
  plot(y~x, data= training_data)
  points(true_val~x, data=training_data, col="blue",type="l",lwd=1)
  if(!poisson){
  matplot(test_data$x, test_CI,type="l",col="black",lty=c(2,1,2),add = T)
  matplot(test_data$x, test_CI_only_CI,
          type="l",col="red",lty=c(2,1,2),lwd=2,add = T)
}
  if(poisson=="TRUE"){
    matplot(test_data$x, exp(test_CI),type="l",col="black",lty=c(2,1,2),add = T)
    matplot(test_data$x, exp(test_CI_only_CI),
            type="l",col="red",lty=c(2,1,2),lwd=2,add = T)
  }
  
  #abline(v= true_min, col="blue",lty=2)
  
}

plot_derivs <- function(model_coefs, possible_peaks, plot_2nd_deriv=TRUE, poisson=FALSE){
  test_data = model_coefs$test_data
  candidate_peaks = possible_peaks$candidate_peaks
  test_1st_deriv_CI = possible_peaks$CI$test_1st_deriv_CI
  test_1st_deriv_CI_only_CI <- test_1st_deriv_CI
  test_1st_deriv_CI_only_CI[!candidate_peaks,] <- NA
    
  test_2nd_deriv_CI = possible_peaks$CI$test_2nd_deriv_CI
  test_2nd_deriv_CI_only_CI <- test_2nd_deriv_CI
  test_2nd_deriv_CI_only_CI[!candidate_peaks,] <- NA
  
  
  #plot of first derivatives plus CI
  if(!poisson){
  matplot(test_data$x,test_1st_deriv_CI,type="l",col="black",lty=c(2,1,2))
  matplot(test_data$x, test_1st_deriv_CI_only_CI,
          type="l",col="red",lty=c(2,1,2),lwd=2,add = T)
  }
  if(poisson){
    matplot(test_data$x,exp(test_1st_deriv_CI),type="l",col="black",lty=c(2,1,2))
    matplot(test_data$x, exp(test_1st_deriv_CI_only_CI),
            type="l",col="red",lty=c(2,1,2),lwd=2,add = T)
  }
  abline(h=0,lty=3,col="red")
  
  if(plot_2nd_deriv){
    
    #Plot of estimated 2nd derivative plus CI
    matplot(test_data$x, test_2nd_deriv_CI,type="l",col="black",lty=c(2,1,2))
    matplot(test_data$x, test_2nd_deriv_CI_only_CI,
            type="l",col="red",lty=c(2,1,2),add = T)
    if(poisson){
      matplot(test_data$x, exp(test_2nd_deriv_CI),type="l",col="black",lty=c(2,1,2))
      matplot(test_data$x, exp(test_2nd_deriv_CI_only_CI),
              type="l",col="red",lty=c(2,1,2),add = T)
    }
    abline(h=0, lty=3, col="red")
  }
}


