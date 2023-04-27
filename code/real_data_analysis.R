#real_data_analysis.R
# Author: Natalie Dupont
# Code for functions peaktest_quad, reg2 and peaktest_twolines, and all derivatives thereof are written by Uri Simonsohn (2018)
# Source: https://osf.io/wdbmr
# Date created: 2023-03-23
# Date last modified: 2023-04-25
# BIOL 490, Pedersen Lab, Concordia University, Montreal
# Description: This script serves to analyse real data from published literature on diversity-disturbance. 
# data is processed by data_cleaning.R beforehand. We apply the quadratic test, the two-lines test, the scam test and the gam derivative test

#outline
#1: Loading packages
#2 Collins 1995
#3 Hiura 1995
#4 Halpern and Spies 1995
#5 Marra 2014
#6 Wilson and Keddy 1988
#7 Nilsson

# Loading Packages ####
  library(scam)
  library(ggplot2)
  library(dplyr) #for data handling
  library(gratia)
  source("code/parametrized_peak_finder.r")
  dev.off()

#Reloading functions used in simulation tests ####
  peaktest_scam=function(x, y, trend="increasing", timeLimit="yes", scamType="m21"){
    #Setting trend to default increasing and returning a warning if not specified, ensuring it is either increasing or decreasing
    if(missing(trend)) {
      warning("no trend selected; default to using increasing - bs=mpi in constrained smooth")
      trend=match.arg(trend)
    }
    else if(!(trend=="increasing"|| trend=="inc"||trend=="increase"||trend=="decreasing"||trend=="dec"||trend=="decrease")){
      stop("invalid trend selected")
    }
    
    #Setting time limit to default increasing and returning a warning if not specified
    if(missing(timeLimit)){
      warning("time limit not specified; default to using a time limit (20 minutes)")
      timeLimit=match.arg(timeLimit)
    }
    else if (!(timeLimit=="Yes" ||timeLimit=="yes" ||timeLimit=="Y" ||timeLimit=="YES"||timeLimit=="No" ||timeLimit=="no" ||timeLimit=="N" ||timeLimit=="NO")){
      stop("invalid time limit request selected")
    }
    
    if(missing(scamType)){
      warning("scamType not specified, default to using m=c(2, 1")
      scamType=match.arg(scamType)
    }
    else if(!(scamType=="m22"||scamType=="m21")){
      stop("invalid scam penalization type selected")
    }
    
    # Determining n/2 to assign k (number of basis functions per predictor)
    half_length <- length(x)/2
    
    # Setting time limit 
    #timeout limit for functions to skip non-convergence
    timeout <- 60*20
    
    
    #fitting a shape-constained model with mpi or mpd based on user input. 
    #Null hypothesis is mpi/mpd constrained smooth predictor and alternate hypothesis is the unconstrained ps predictor
    #each predictor (constrained and unconstrained) have n/2 basis functions
    
    #if-else statements have been checked and work appropriately
    #withTimout has been checked (separately) and is confirmed to return an object but populated by NULL
    if(scamType=="m21"){
      if (timeLimit=="Yes" ||timeLimit=="yes" ||timeLimit=="Y" ||timeLimit=="YES"){
        if (trend=="increasing"|| trend=="inc"||trend=="increase"){
          scam = withTimeout(scam(y~s(x,k=half_length, bs="mpi", m=2)+s(x, k=half_length, bs="ps",  m =c(2,1))),
                             timeout=timeout, 
                             onTimeout="warning")
        }
        else{
          scam = withTimeout(scam(y~s(x,k=half_length, bs="mpd", m=2)+s(x, k=half_length, bs="ps",  m =c(2,1))),
                             timeout=timeout, 
                             onTimeout="warning")
        }
      }
      else {
        if (trend=="increasing"|| trend=="inc"||trend=="increase"){
          scam = scam(y~s(x,k=half_length, bs="mpi", m=2)+s(x, k=half_length, bs="ps",  m =c(2,1)))
        }
        else{
          scam = scam(y~s(x,k=half_length, bs="mpd", m=2)+s(x, k=half_length, bs="ps",  m =c(2,1)))
        }
      }
    }
    #for m=c(2, 2)
    else {
      if (timeLimit=="Yes" ||timeLimit=="yes" ||timeLimit=="Y" ||timeLimit=="YES"){
        if (trend=="increasing"|| trend=="inc"||trend=="increase"){
          scam = withTimeout(scam(y~s(x,k=half_length, bs="mpi", m=2)+s(x, k=half_length, bs="ps",  m =c(2,2))),
                             timeout=timeout, 
                             onTimeout="warning")
        }
        else {
          scam = withTimeout(scam(y~s(x,k=half_length, bs="mpd", m=2)+s(x, k=half_length, bs="ps",  m =c(2,2))),
                             timeout=timeout, 
                             onTimeout="warning")
        }
      }
      else {
        if (trend=="increasing"|| trend=="inc"||trend=="increase"){
          scam = scam(y~s(x,k=half_length, bs="mpi", m=2)+s(x, k=half_length, bs="ps",  m =c(2,2)))
        }
        else {
          scam = scam(y~s(x,k=half_length, bs="mpd", m=2)+s(x, k=half_length, bs="ps",  m =c(2,2)))
        }
      }
    }
    
    
    
    #determining if there is a significant peak and returning it as a boolean 0 or 1
    #if function timed out and scam is null, return a vector or 2 NAs. 
    #withTimeout reassigns NULL even if the object is already populated (no concern of reassigning the previous run's results)
    if(is.null(scam)){
      return(c(NA, NA))
    }
    
    #if function ran correctly and scam object is populated, extract the p-value of nonmono smooth and the boolean significance
    summary <- summary(scam)
    
    p_nonmono <- summary$s.table[,"p-value"][[2]]
    u <- 0 #initializing to 0 should be fine since a timed-out analysis will already return a value and exit the function
    if(p_nonmono<=0.05){
      u <- 1
    }
    
    #return  whether the test result is significant or not and nonmono p-value
    return(c(u, p_nonmono))
  }
  
  
  ### 2-1 b: loading quadratic regression function peaktest_quad (flood.u in source code) source: https://osf.io/wdbmr ####
  # **** RETURNS NaN sometimes, and I don't know why... (happened with skewed p=8 beta curve)
  #happens in the j1 and j2 functions (sqrt(-)
  peaktest_quad=function(x,y) {
    #Sample size 
    n=length(x)
    #Square x to include in the regression as predictor
    x2=x^2
    #Run quadratic regression
    qlm=lm(y~x+x2)   
    #Run floodlight on results
    coef <- coefficients(qlm) #point estimates
    vc <- vcov(qlm)           #VAR matrix
    t <- qt(.975,df = n-3)    #critical t-stat
    #Critical points where slope is significantly of each sign
    jn1 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) - 
             sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
    jn2 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) + 
             sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
    
    #if both are within the range of values, we have a u-shape
    u=0
    if (!is.na(jn1) & !is.na(jn2)) if (jn1>min(x) & jn1<max(x) & jn2>min(x) & jn2<max(x) & jn1>jn2) u=1
    #return u=1 if significant u (both ts have opposite sign and are p<.05), =0 otherwise
    return(u)
  }
  
  ### 2-1c:  Loading two-lines function peaktest_twolines (reg2hood in source code) source: https://osf.io/wdbmr ####
  #reg2 fits the two-line regression based on a given break point xc
  reg2=function(x,y,xc)
  { 
    #xc is included in the first line
    xlow1=ifelse(x<=xc,x-xc,0)     #xlow=x-xc when x<xc, 0 otherwise
    xhigh1=ifelse(x>xc,x-xc,0)     #xhigh=x when x<xmax, 0 otherwise
    high1=ifelse(x>xc,1,0)         #high dummy, allows interruption
    
    #Now include xc in second line
    xlow2=ifelse(x<xc,x-xc,0)     
    xhigh2=ifelse(x>=xc,x-xc,0)     
    high2=ifelse(x>=xc,1,0)         #high dummy, allows interruption
    #Run the regressions (they differe only on whetehr xc is 'high' or 'low')
    lm1=lm(y~xlow1+xhigh1+high1)     #estimate regression
    lm2=lm(y~xlow2+xhigh2+high2)     #estimate regression
    #Fitted values
    yhat1=fitted(lm1)              
    yhat2=fitted(lm2)              
    #Regression results
    lmc1=summary(lm1)$coefficients 
    lmc2=summary(lm2)$coefficients 
    #Turn into individual results (scalars)
    #Results for 1st line come from 1st regression, including xc1
    b1=lmc1[2,1]
    t1=lmc1[2,3]
    p1=lmc1[2,4]
    #Results for 2nd line come from 2nd regression, including xc1 as well
    b2=lmc2[3,1]
    t2=lmc2[3,3]
    p2=lmc2[3,4]
    #Is the u-shape significnat?
    u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)                     
    #All results
    res=list(b1=b1,p1=p1,b2=b2,p2=p2,u.sig=u.sig,xc=xc,t1=t1,t2=t2,yhat1=yhat1,yhat2=yhat2)  #Output list with all those parameters, betas, t-values, p-values and significance for u
    return(res)
  }
  
  #peaktest_twolines (reg2hood in source code) determines break points based on different criteria
  #rmax: xc=max(yhat)
  #rmed: xc=median(yflat)- yflat defined as the range where yhat+1SE>=max(yhat)
  #rprop: uses robin hood algorithm on the rmed break point to determine a new breakpoint
  #u.3lines: uses 3-lines method (not sure beyond that...)
  peaktest_twolines=function(x,y)
  {
    #Syntax:
    #1 Run gam()
    g=gam(y~s(x,bs="cr"))  
    #2 Get fitted values
    g.fit=predict.gam(g,se.fit=T)
    y.hat=g.fit$fit
    y.se =g.fit$se.fit
    #3 Focus on the middle 80% of the x-values (smoothed values are not reliable near the end)
    x10=quantile(x,.1)
    x90=quantile(x,.9)
    middle=(x>x10 & x<x90)       #Don't consider extreme values for cutoff
    y.ub=y.hat+1*y.se            #+1 SE is the flat region
    x.middle=x[middle]           #X values associated with flat region
    xc.max=x.middle[match(max(y.hat[middle]),y.hat[middle])]       #find value of x associated with highest predicted value
    #4 Find flat maximum  
    flat=(y.ub>max(y.hat) & middle)
    xflat=x[flat] 
    #5 If empty, use median x
    if (length(xflat)==0) xflat=median(x) 
    
    #6 Regression split based on predicted maximum
    rmax=reg2(x,y,xc=xc.max) 
    
    #7  Regression split based median of xflat
    rmed=reg2(x,y,xc=median(xflat))  #At the median of  xflat
    #8 Adjust split point based on ratio of t1,t2, move split away from less significant slope, to give more of a chance
    t1=abs(rmed$t1)             
    t2=abs(rmed$t2)             
    xc.prop=quantile(xflat,t2/(t1+t2))  #Choose the percentile value in flat region proportional to t1,t2
    #For example, if t2=t1 it stays at median flat, if t2>t1, it moves lower
    #9 Regression split based on adjusted based on t1,t2    
    rprop=reg2(x,y,xc=xc.prop) 
    
    #10 Run 3-semgent regression
    xc1=min(x.middle)
    xc2=max(x.middle)
    
    #10.1 Gen variables for twice interrupted regression
    x1=ifelse(x<xc1,x-xc1,0)
    x2=ifelse(x>xc1 & x<xc2,x-xc1,0)
    x3=ifelse(x>xc2, x-xc2,0)
    high1=ifelse(x<xc1,0,1)
    high2=ifelse(x<xc2,0,1)
    
    #10.2 Run the regressions (they differe only on whetehr xc is 'high' or 'low')
    lm1=lm(y~x1+x2+x3+high1+high2)      #estimate regression
    lmc=summary(lm1)$coefficients  
    
    #10.3 Turn into individual results (scalars)
    #Results for 1st line come from 1st regression, including xc1
    b1=lmc[2,1]
    t1=lmc[2,3]
    p1=lmc[2,4]
    
    #10.4 Results for 3rd line come 
    b3=lmc[4,1]
    t3=lmc[4,3]
    p3=lmc[4,4]
    
    #10.5 Is the u-shape significant?
    u.3lines =ifelse(b1*b3<0 & p1<.05 & p3<.05,1,0)      
    
    #11 Output is four 1/0 dummies for significant/not of each of the four regressions
    res=c(rmax$u.sig,rmed$u.sig, rprop$u.sig,u.3lines)  
    return(res)
    #to find significance values, rmax$p1 and p2 wwill be the p-values for each line
    #u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)      
    #considered significant if slopes are of opposite signs and both slope p-values are significant
  }
  
  
  
# Altered functions ####
  ##Quadratic test returns model object ####
  peaktest_quad_return=function(x,y) {
    #Sample size 
    n=length(x)
    #Square x to include in the regression as predictor
    x2=x^2
    #Run quadratic regression
    qlm=lm(y~x+x2)   
    #Run floodlight on results
    coef <- coefficients(qlm) #point estimates
    vc <- vcov(qlm)           #VAR matrix
    t <- qt(.975,df = n-3)    #critical t-stat
    #Critical points where slope is significantly of each sign
    jn1 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) - 
             sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
    jn2 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) + 
             sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
    
    #if both are within the range of values, we have a u-shape
    u=0
    if (!is.na(jn1) & !is.na(jn2)) if (jn1>min(x) & jn1<max(x) & jn2>min(x) & jn2<max(x) & jn1>jn2) u=1
    #return u=1 if significant u (both ts have opposite sign and are p<.05), =0 otherwise
    return(list(u=u, qlm=qlm))
  }
  
  ##Quadratic test but poisson ####
    peaktest_quad_pois_return=function(x,y) {
      #Sample size 
      n=length(x)
      #Square x to include in the regression as predictor
      x2=x^2
      #Run quadratic regression
      qlm=glm(y~x+x2, family=poisson)   
      #Run floodlight on results
      coef <- coefficients(qlm) #point estimates
      vc <- vcov(qlm)           #VAR matrix
      t <- qt(.975,df = n-3)    #critical t-stat
      #Critical points where slope is significantly of each sign
      jn1 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) - 
               sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
      jn2 = (-(4*t^2*vc[2,3]-4*coef[2]*coef[3]) + 
               sqrt((4*t^2*vc[2,3] - 4*coef[2]*coef[3])^2 - 4*(4*t^2*vc[3,3]-4*coef[3]^2)*(t^2*vc[2,2]-coef[2]^2)))/(2*(4*t^2*vc[3,3]-4*coef[3]^2))
      
      #if both are within the range of values, we have a u-shape
      u=0
      if (!is.na(jn1) & !is.na(jn2)) if (jn1>min(x) & jn1<max(x) & jn2>min(x) & jn2<max(x) & jn1>jn2) u=1
      #return u=1 if significant u (both ts have opposite sign and are p<.05), =0 otherwise
      return(list(u=u, qlm=qlm))
    }
    
    ##Two_lines-Return model ####
    reg2_returnModel=function(x,y,xc)
    { 
      #xc is included in the first line
      xlow1=ifelse(x<=xc,x-xc,0)     #xlow=x-xc when x<xc, 0 otherwise
      xhigh1=ifelse(x>xc,x-xc,0)     #xhigh=x when x<xmax, 0 otherwise
      high1=ifelse(x>xc,1,0)         #high dummy, allows interruption
      
      #Now include xc in second line
      xlow2=ifelse(x<xc,x-xc,0)     
      xhigh2=ifelse(x>=xc,x-xc,0)     
      high2=ifelse(x>=xc,1,0)         #high dummy, allows interruption
      #Run the regressions (they differe only on whetehr xc is 'high' or 'low')
      lm1=lm(y~xlow1+xhigh1+high1)     #estimate regression
      lm2=lm(y~xlow2+xhigh2+high2)     #estimate regression
                  
      #All results
      res=list(lm1=lm1,lm2=lm2, xc) 
      return(res)
    }
    
    #peaktest_twolines (reg2hood in source code) determines break points based on different criteria
    #rmax: xc=max(yhat)
    #rmed: xc=median(yflat)- yflat defined as the range where yhat+1SE>=max(yhat)
    #rprop: uses robin hood algorithm on the rmed break point to determine a new breakpoint
    #u.3lines: uses 3-lines method (not sure beyond that...)
    peaktest_twolines_returnModel=function(x,y)
    {
      #Syntax:
      #1 Run gam()
      g=gam(y~s(x,bs="cr"))  
      #2 Get fitted values
      g.fit=predict.gam(g,se.fit=T)
      y.hat=g.fit$fit
      y.se =g.fit$se.fit
      #3 Focus on the middle 80% of the x-values (smoothed values are not reliable near the end)
      x10=quantile(x,.1)
      x90=quantile(x,.9)
      middle=(x>x10 & x<x90)       #Don't consider extreme values for cutoff
      y.ub=y.hat+1*y.se            #+1 SE is the flat region
      x.middle=x[middle]           #X values associated with flat region
      xc.max=x.middle[match(max(y.hat[middle]),y.hat[middle])]       #find value of x associated with highest predicted value
      #4 Find flat maximum  
      flat=(y.ub>max(y.hat) & middle)
      xflat=x[flat] 
      #5 If empty, use median x
      if (length(xflat)==0) xflat=median(x) 
      
      #6 Regression split based on predicted maximum
      rmax=reg2(x,y,xc=xc.max) 
      
      #7  Regression split based median of xflat
      rmed=reg2(x,y,xc=median(xflat))  #At the median of  xflat
      #8 Adjust split point based on ratio of t1,t2, move split away from less significant slope, to give more of a chance
      t1=abs(rmed$t1)             
      t2=abs(rmed$t2)             
      xc.prop=quantile(xflat,t2/(t1+t2))  #Choose the percentile value in flat region proportional to t1,t2
      #For example, if t2=t1 it stays at median flat, if t2>t1, it moves lower
      #9 Regression split based on adjusted based on t1,t2    
      rprop=reg2_returnModel(x,y,xc=xc.prop) 
      
      #10 Run 3-semgent regression
      res=rprop
      return(res)
      #to find significance values, rmax$p1 and p2 wwill be the p-values for each line
      #u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)      
      #considered significant if slopes are of opposite signs and both slope p-values are significant
    }
    
    ## Two-lines-Pois####
    reg2_pois=function(x,y,xc)
    { 
      #xc is included in the first line
      xlow1=ifelse(x<=xc,x-xc,0)     #xlow=x-xc when x<xc, 0 otherwise
      xhigh1=ifelse(x>xc,x-xc,0)     #xhigh=x when x<xmax, 0 otherwise
      high1=ifelse(x>xc,1,0)         #high dummy, allows interruption
      
      #Now include xc in second line
      xlow2=ifelse(x<xc,x-xc,0)     
      xhigh2=ifelse(x>=xc,x-xc,0)     
      high2=ifelse(x>=xc,1,0)         #high dummy, allows interruption
      #Run the regressions (they differe only on whetehr xc is 'high' or 'low')
      lm1=glm(y~xlow1+xhigh1+high1, family=poisson)     #estimate regression
      lm2=glm(y~xlow2+xhigh2+high2, family=poisson)     #estimate regression
      #Fitted values
      yhat1=fitted(lm1)              
      yhat2=fitted(lm2)              
      #Regression results
      lmc1=summary(lm1)$coefficients 
      lmc2=summary(lm2)$coefficients 
      #Turn into individual results (scalars)
      #Results for 1st line come from 1st regression, including xc1
      b1=lmc1[2,1]
      t1=lmc1[2,3]
      p1=lmc1[2,4]
      #Results for 2nd line come from 2nd regression, including xc1 as well
      b2=lmc2[3,1]
      t2=lmc2[3,3]
      p2=lmc2[3,4]
      #Is the u-shape significnat?
      u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)                     
      #All results
      res=list(b1=b1,p1=p1,b2=b2,p2=p2,u.sig=u.sig,xc=xc,t1=t1,t2=t2,yhat1=yhat1,yhat2=yhat2, lm1=lm1, lm2=lm2)  #Output list with all those parameters, betas, t-values, p-values and significance for u
      return(res)
    }
    
    #peaktest_twolines (reg2hood in source code) determines break points based on different criteria
    #rmax: xc=max(yhat)
    #rmed: xc=median(yflat)- yflat defined as the range where yhat+1SE>=max(yhat)
    #rprop: uses robin hood algorithm on the rmed break point to determine a new breakpoint
    #u.3lines: uses 3-lines method (not sure beyond that...)
    peaktest_twolines_pois=function(x,y)
    {
      #Syntax:
      #1 Run gam()
      g=gam(y~s(x,bs="cr"), family=poisson)  
      #2 Get fitted values
      g.fit=predict.gam(g,se.fit=T)
      y.hat=g.fit$fit
      y.se =g.fit$se.fit
      #3 Focus on the middle 80% of the x-values (smoothed values are not reliable near the end)
      x10=quantile(x,.1)
      x90=quantile(x,.9)
      middle=(x>x10 & x<x90)       #Don't consider extreme values for cutoff
      y.ub=y.hat+1*y.se            #+1 SE is the flat region
      x.middle=x[middle]           #X values associated with flat region
      xc.max=x.middle[match(max(y.hat[middle]),y.hat[middle])]       #find value of x associated with highest predicted value
      #4 Find flat maximum  
      flat=(y.ub>max(y.hat) & middle)
      xflat=x[flat] 
      #5 If empty, use median x
      if (length(xflat)==0) xflat=median(x) 
      
      #6 Regression split based on predicted maximum
      rmax=reg2_pois(x,y,xc=xc.max) 
      
      #7  Regression split based median of xflat
      rmed=reg2_pois(x,y,xc=median(xflat))  #At the median of  xflat
      #8 Adjust split point based on ratio of t1,t2, move split away from less significant slope, to give more of a chance
      t1=abs(rmed$t1)             
      t2=abs(rmed$t2)             
      xc.prop=quantile(xflat,t2/(t1+t2))  #Choose the percentile value in flat region proportional to t1,t2
      #For example, if t2=t1 it stays at median flat, if t2>t1, it moves lower
      #9 Regression split based on adjusted based on t1,t2    
      rprop=reg2_pois(x,y,xc=xc.prop) 
      
      #10 Run 3-semgent regression
      xc1=min(x.middle)
      xc2=max(x.middle)
      
      #10.1 Gen variables for twice interrupted regression
      x1=ifelse(x<xc1,x-xc1,0)
      x2=ifelse(x>xc1 & x<xc2,x-xc1,0)
      x3=ifelse(x>xc2, x-xc2,0)
      high1=ifelse(x<xc1,0,1)
      high2=ifelse(x<xc2,0,1)
      
      #10.2 Run the regressions (they differe only on whetehr xc is 'high' or 'low')
      lm1=glm(y~x1+x2+x3+high1+high2, family=poisson)      #estimate regression
      lmc=summary(lm1)$coefficients  
      
      #10.3 Turn into individual results (scalars)
      #Results for 1st line come from 1st regression, including xc1
      b1=lmc[2,1]
      t1=lmc[2,3]
      p1=lmc[2,4]
      
      #10.4 Results for 3rd line come 
      b3=lmc[4,1]
      t3=lmc[4,3]
      p3=lmc[4,4]
      
      #10.5 Is the u-shape significant?
      u.3lines =ifelse(b1*b3<0 & p1<.05 & p3<.05,1,0)      
      
      #11 Output is four 1/0 dummies for significant/not of each of the four regressions
      res=c(rmax$u.sig,rmed$u.sig, rprop$u.sig,u.3lines)  
      return(list(res=res, lm1=rprop$lm1, lm2=rprop$lm2, xc=xc.prop))
      #to find significance values, rmax$p1 and p2 wwill be the p-values for each line
      #u.sig =ifelse(b1*b2<0 & p1<.05 & p2<.05,1,0)      
      #considered significant if slopes are of opposite signs and both slope p-values are significant
    }

# 1: Collins 1995####
  Collins1995 <- read.csv("data/processed/CLEAN_Collins1995.csv")
  head(Collins1995)
  plot(n_sp~y_pb, data=Collins1995)
  
  ## a: Peak detection ####
    ### Scam ####
      #decisded on increasing
      Collins1995_forscam <- Collins1995
      colnames(Collins1995_forscam) <- c("Years_since_burn", "Species_Richness")
      
      Collins1995_scam <- scam(Species_Richness~s(Years_since_burn,k=10, bs="mpi", m=2)+
                                 s(Years_since_burn, k=10, bs="ps",  m =c(2,1)),
                               data=Collins1995_forscam, family=poisson(link="log"))
      scam.check(Collins1995_scam)
      plot.scam(Collins1995_scam, pages=1, 
                residuals=TRUE, 
                pch=16, cex=.7, rug=FALSE,
                shift = coef(Collins1995_scam)[1], seWithMean=FALSE)
      
      Collins1995_scam_p <- summary(Collins1995_scam)$s.table[,"p-value"][[2]]
      
      lines(fitted(Collins1995_scam)~Collins1995$y_pb, col="red", lwd=2, lty=2)
      
      #peaked
      
    ###Quadtratic ####
      quad_model_Collins1995 <- glm(n_sp~y_pb+I(y_pb^2), data=Collins1995, family=poisson)
      plot(quad_model_Collins1995)
      summary(quad_model_Collins1995)
      Collins1995_quad<- peaktest_quad_pois_return(x = Collins1995$y_pb, y=Collins1995$n_sp)
      summary(Collins1995_quad$qlm)
        #Peak detected
      
      Collins_quad_conf <- predict(Collins1995_quad$qlm, se.fit = TRUE, type="response")
      Collins_quad_plot <- ggplot(data=Collins1995)+
        geom_point(aes(x=y_pb, y=n_sp), size=1)+
        geom_line(aes(y=fitted(Collins1995_quad$qlm), x=y_pb))+
        geom_line(aes(y=I(Collins_quad_conf$fit+1.96*Collins_quad_conf$se.fit), x=y_pb), lty=2, colour="gray")+
        geom_line(aes(y=I(Collins_quad_conf$fit-1.96*Collins_quad_conf$se.fit), x=y_pb), lty=2, colour="gray")+
        theme_bw()+theme(axis.line = element_line(color='black'),
                         plot.background = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.grid.major = element_blank() )+
        labs(title="a) Quadratic Regression",
             x="Years since burn", y="Species Richness")
      Collins_quad_plot
      
    ###Two-Lines ####
      peaktest_twolines(x = Collins1995$y_pb, y=Collins1995$n_sp)
      Collins_twoLines_pois <- peaktest_twolines_pois(x = Collins1995$y_pb, y=Collins1995$n_sp)
      summary(Collins_twoLines_pois$lm1)
      summary(Collins_twoLines_pois$lm1)

      Collins1995_TL1 <- fitted(Collins_twoLines_pois$lm1)
      Collins1995_TL1 <- Collins1995_TL1[which(Collins1995$y_pb<=6)]
      Collins1995_TL2 <- fitted(Collins_twoLines_pois$lm2)
      Collins1995_TL2 <- Collins1995_TL2[which(Collins1995$y_pb>=6)]
      
      Collins_TL1_conf <- predict(Collins_twoLines_pois$lm1, se.fit = TRUE, type="response")
      Collins_TL1_conf$fit <- Collins_TL1_conf$fit[which(Collins1995$y_pb<=6)]
      Collins_TL1_conf$se.fit <- Collins_TL1_conf$se.fit[which(Collins1995$y_pb<=6)]

      Collins_TL2_conf <- predict(Collins_twoLines_pois$lm2, se.fit = TRUE, type="response")
      Collins_TL2_conf$fit <- Collins_TL2_conf$fit[which(Collins1995$y_pb>=6)]
      Collins_TL2_conf$se.fit <- Collins_TL2_conf$se.fit[which(Collins1995$y_pb>=6)]
      
      Collins_twoLines_plot <- ggplot()+
        geom_point(aes(x=Collins1995$y_pb, y=Collins1995$n_sp), size=1)+
        geom_line(aes(y=Collins1995_TL1, x=Collins1995$y_pb[which(Collins1995$y_pb<=6)]))+
        geom_line(aes(y=Collins1995_TL2, x=Collins1995$y_pb[which(Collins1995$y_pb>=6)]))+
        geom_line(aes(y=I(Collins_TL1_conf$fit+1.96*Collins_TL1_conf$se.fit), x=Collins1995$y_pb[which(Collins1995$y_pb<=6)]), lty=2, colour="gray")+
        geom_line(aes(y=I(Collins_TL1_conf$fit-1.96*Collins_TL1_conf$se.fit), x=Collins1995$y_pb[which(Collins1995$y_pb<=6)]), lty=2, colour="gray")+
        geom_line(aes(y=I(Collins_TL2_conf$fit+1.96*Collins_TL2_conf$se.fit), x=Collins1995$y_pb[which(Collins1995$y_pb>=6)]), lty=2, colour="gray")+
        geom_line(aes(y=I(Collins_TL2_conf$fit-1.96*Collins_TL2_conf$se.fit), x=Collins1995$y_pb[which(Collins1995$y_pb>=6)]), lty=2, colour="gray")+ 
        theme_bw()+
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank() )+
        labs(title="b) Two-Lines Regression",
             x="Years since burn", y="Species Richness")
      Collins_twoLines_plot
    
    ###Plotting predicted curve for each approach ####
      Collins_quad_plot+Collins_twoLines_plot
      
      predicted_Collins1995 <- ggplot(data=Collins1995, aes(x=y_pb))+
        geom_point(aes(y=n_sp))+
        geom_line(aes(y=fitted(Collins1995_scam)))+
        geom_line(aes(y=fitted(glm(n_sp~y_pb+I(y_pb^2), data=Collins1995, family=poisson))))+
        geom_line(aes(y=fitted(Collins_twoLines_pois$lm1)), colour="red")
      predicted_Collins1995
      
      
  ##Peak quantification ####
    data_Collins_mod <- tibble(x = Collins1995$y_pb,y= Collins1995$n_sp,
                                 true_val = NA)
    Collins_mod <- gam(y~s(x, bs="ad", m=3, k=10),method="REML", data=data_Collins_mod, family=poisson)
    appraise(Collins_mod)
    summary(Collins_mod)
    plot(Collins_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Collins1995_scam)[1])
    plot(Collins1995$n_sp~Collins1995$y_pb)
    lines(fitted(Collins_mod)~Collins1995$y_pb)
    # lines(fitted(Collins_mod)+exp(predict.gam(Collins_mod)$se.fit)~Collins1995$y_pb, lty=2)
    # lines(fitted(Collins_mod)-exp(predict.gam(Collins_mod)$se.fit)~Collins1995$y_pb, lty=2)

   # lines(predict(Collins_mod,newdata =Collins1995$y_pb,  interval="confidence", level=0.95), type="l")
    
    model_summary <- create_peakfinder_input(mod = Collins_mod, data = data_Collins_mod,step_size = diff(range(Collins1995$y_pb))/1000)
    set.seed(1)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    #No peak found
  
# 2: Hiura 1995 #####
    Hiura1995 <- read.csv("data/processed/CLEAN_Hiura1995.csv")
    head(Hiura1995)
    plot(H_prime~ln_MWI, data=Hiura1995)
    Hiura1995$sqrt_MWI <- sqrt(Hiura1995$ln_MWI)
    
    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    Hiura1995_scam <- scam(H_prime~s(sqrt_MWI,k=8, bs="mpi", m=2)+
                               s(sqrt_MWI, k=8, bs="ps",  m =c(2,1)),
                             data=Hiura1995)
    scam.check(Hiura1995_scam)
      #Not normally distributed. right-skewed
      #Sqrt instead of ln-transformed looks better...
    plot.scam(Hiura1995_scam, pages=1, 
              residuals=TRUE, 
              pch=1, cex=1,
              shift = coef(Hiura1995_scam)[1])
    
    Hiura1995_scam_p <- summary(Hiura1995_scam)$s.table[,"p-value"][[2]]
      #not peaked
    
    ###Quadtratic ###
    quad_model_Hiura1995 <- glm(H_prime~ln_MWI+I(ln_MWI^2), data=Hiura1995)
      plot(quad_model_Hiura1995)
      #not too bad
      peaktest_quad_return(x = Hiura1995$ln_MWI,y=Hiura1995$H_prime )
      #not peaked
    
    ###Two-Lines ####
    peaktest_twolines(x = Hiura1995$ln_MWI,y=Hiura1995$H_prime )
      #peaked
    Hiura1995_twoLinesModel <- peaktest_twolines_returnModel(x = Hiura1995$ln_MWI,y=Hiura1995$H_prime )

    Hiura1995_TL1 <- fitted(Hiura1995_twoLinesModel$lm1)
    Hiura1995_TL1 <- Hiura1995_TL1[which(Hiura1995$ln_MWI<=1.876325)]
    Hiura1995_TL2 <- fitted(Hiura1995_twoLinesModel$lm2)
    Hiura1995_TL2 <- Hiura1995_TL2[which(Hiura1995$ln_MWI>=1.876325)]
    
    plot(H_prime~ln_MWI, data=Hiura1995)
    lines(Hiura1995_TL1~Hiura1995$ln_MWI[which(Hiura1995$ln_MWI<=1.876325)],col="red", lwd=2, lty=3)
    lines(Hiura1995_TL2~Hiura1995$ln_MWI[which(Hiura1995$ln_MWI>=1.876325)],col="red", lwd=2, lty=3)
    
    ###Plotting predicted curve for each approach ####
    predicted_Hiura1995 <- ggplot(data=Hiura1995, aes(x=ln_MWI))+
      geom_point(aes(y=H_prime))+
      geom_line(aes(y=fitted(Hiura1995_scam)))+
      geom_line(aes(y=fitted(glm(H_prime~ln_MWI+I(ln_MWI^2), data=Hiura1995))))+
      geom_line(aes(y=fitted(Hiura1995_twoLinesModel$lm1)), colour="red")
    predicted_Hiura1995
    
    ##Peak quantification ####
    data_Hiura_mod <- tibble(x = Hiura1995$ln_MWI,y= Hiura1995$H_prime,
                               true_val = NA)
    Hiura_mod <- gam(y~s(x, bs="ad", k=10, m=3),method="REML", data=data_Hiura_mod)
    gam.check(Hiura_mod)
    plot(Hiura_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Collins1995_scam)[1])
    
    model_summary <- create_peakfinder_input(mod = Hiura_mod, data = data_Hiura_mod,step_size = diff(range(Hiura1995$ln_MWI))/1000)
    set.seed(2)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    #No peak found

    
# 3: Halpern and Spies 1995, W1 ####
  HalpernSpies1995_w1 <- read.csv("data/processed/CLEAN_HalpernSpies1995_Fig1d_w1.csv")
  plot(HalpernSpies1995_w1$heterogeneity~HalpernSpies1995_w1$t_disturbance)
  
  ## a: Peak detection ####
  ### Scam ####
  #decisded on increasing
  Halpernw1_scam <- scam(heterogeneity~s(t_disturbance,k=6, bs="mpd", m=2)+
                           s(t_disturbance, k=8, bs="ps",  m =c(2,1)),
                         data=HalpernSpies1995_w1)
  scam.check(Halpernw1_scam)
  #weird residuals, but just 1-2 on the ends
  plot.scam(Halpernw1_scam, pages=1, 
            residuals=TRUE, 
            pch=16, cex=0.7,
            shift = coef(Halpernw1_scam)[1])
  
  Halpernw1_scam_p <- summary(Halpernw1_scam)$s.table[,"p-value"][[2]]
  Halpernw1_scam_p
  #PEAKED
  
  ###Quadtratic ####
  quad_model_Halpernw1 <- lm(heterogeneity~t_disturbance+I(t_disturbance^2), data=HalpernSpies1995_w1)
  plot(quad_model_Halpernw1)
  
  Halpern1995_quad <- peaktest_quad_return(x = HalpernSpies1995_w1$t_disturbance,y=HalpernSpies1995_w1$heterogeneity )
  summary(Halpern1995_quad$qlm)
  # peaked
  
  Halpern1_quad_conf <- predict(Halpern1995_quad$qlm, se.fit = TRUE, type="response")
  Halpern1_quad_plot <- ggplot(data=HalpernSpies1995_w1)+
    geom_point(aes(x=t_disturbance, y=heterogeneity), size=1)+
    geom_line(aes(y=fitted(Halpern1995_quad$qlm), x=t_disturbance))+
    geom_line(aes(y=I(Halpern1_quad_conf$fit+1.96*Halpern1_quad_conf$se.fit), x=t_disturbance), lty=2, colour="gray")+
    geom_line(aes(y=I(Halpern1_quad_conf$fit-1.96*Halpern1_quad_conf$se.fit), x=t_disturbance), lty=2, colour="gray")+
    theme_bw()+theme(axis.line = element_line(color='black'),
                     plot.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major = element_blank() )+
    labs(title="a) Quadratic Regression",
         x="Years since burn", y="Species Richness")
  Halpern1_quad_plot
  
  
  ###Two-Lines ####
  peaktest_twolines(x = HalpernSpies1995_w1$t_disturbance,y=HalpernSpies1995_w1$heterogeneity )
  #peaked
  HalperntwoLinesModel <- peaktest_twolines_returnModel(x = HalpernSpies1995_w1$t_disturbance,y=HalpernSpies1995_w1$heterogeneity )
  summary(HalperntwoLinesModel$lm1)
  summary(HalperntwoLinesModel$lm2)
  HalperntwoLinesModel
  
  Halpern11995_TL1 <- fitted(HalperntwoLinesModel$lm1)
  Halpern11995_TL1 <- Halpern11995_TL1[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]
  Halpern11995_TL2 <- fitted(HalperntwoLinesModel$lm2)
  Halpern11995_TL2 <- Halpern11995_TL2[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]
  
  Halpern1_TL1_conf <- predict(HalperntwoLinesModel$lm1, se.fit = TRUE, type="response")
  Halpern1_TL1_conf$fit <- Halpern1_TL1_conf$fit[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]
  Halpern1_TL1_conf$se.fit <- Halpern1_TL1_conf$se.fit[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]
  
  Halpern1_TL2_conf <- predict(HalperntwoLinesModel$lm2, se.fit = TRUE, type="response")
  Halpern1_TL2_conf$fit <- Halpern1_TL2_conf$fit[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]
  Halpern1_TL2_conf$se.fit <- Halpern1_TL2_conf$se.fit[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]
  
  Halper1_twoLines_plot <- ggplot()+
    geom_point(aes(x=HalpernSpies1995_w1$t_disturbance, y=HalpernSpies1995_w1$heterogeneity), size=1)+
    geom_line(aes(y=Halpern11995_TL1, x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]))+
    geom_line(aes(y=Halpern11995_TL2, x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]))+
    geom_line(aes(y=I(Halpern1_TL1_conf$fit+1.96*Halpern1_TL1_conf$se.fit), x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]), lty=2, colour="gray")+
    geom_line(aes(y=I(Halpern1_TL1_conf$fit-1.96*Halpern1_TL1_conf$se.fit), x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance<=4.631718)]), lty=2, colour="gray")+
    geom_line(aes(y=I(Halpern1_TL2_conf$fit+1.96*Halpern1_TL2_conf$se.fit), x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]), lty=2, colour="gray")+
    geom_line(aes(y=I(Halpern1_TL2_conf$fit-1.96*Halpern1_TL2_conf$se.fit), x=HalpernSpies1995_w1$t_disturbance[which(HalpernSpies1995_w1$t_disturbance>=4.631718)]), lty=2, colour="gray")+ 
    theme_bw()+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank() )+
    labs(title="b) Two-Lines Regression",
         x="Years since burn", y="Species Richness")
  Halper1_twoLines_plot
  
  
  ##Peak quantification ####
  data_halperm_mod <- tibble(x = HalpernSpies1995_w1$t_disturbance,y=HalpernSpies1995_w1$heterogeneity ,
                           true_val = NA)
  Halpern_mod <- gam(y~s(x, bs="ad", k=10, m=3),method="REML", data=data_halperm_mod)
  appraise(Halpern_mod)
  summary(Halpern_mod)
    #Doesn't look great
  plot(Halpern_mod,
       residuals=TRUE, 
       pch=1, cex=1,
       shift = coef(Halpern_mod)[1])
  
  model_summary <- create_peakfinder_input(mod = Halpern_mod, data = data_halperm_mod,step_size = diff(range(data_halperm_mod$x))/1000)
  set.seed(3)
  model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
  possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
  
  #Confidence interval:
  CI_Halpern <- possible_peaks$x_val[which(possible_peaks$candidate_peaks)]
  
  plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                               plot_2nd_deriv = FALSE ,
                               poisson=FALSE)
  plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
  #Peak found
  
    peak_halpernw1 <- data.frame(x=seq(from=min(data_halperm_mod$x), to=max(data_halperm_mod$x), by= diff(range(data_halperm_mod$x))/1000),
                                 peaks=possible_peaks$candidate_peaks)
    peak_CI_Halpernw2 <- c("2.5%"=min(peak_halpernw1$x[which(peak_halpernw1$peaks)]),
                           "97.5"=max(peak_halpernw1$x[which(peak_halpernw1$peaks)]))
    
# 3: Halpern and Spies 1995, W2 ####
    HalpernSpies1995_w2 <- read.csv("data/processed/CLEAN_HalpernSpies1995_Fig1d_w2.csv")
    plot(HalpernSpies1995_w2$heterogeneity~HalpernSpies1995_w2$t_disturbance)
    
    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    Halpernw2_scam <- scam(heterogeneity~s(t_disturbance,k=8, bs="mpd", m=2)+
                             s(t_disturbance, k=8, bs="ps",  m =c(2,1)),
                           data=HalpernSpies1995_w2)
    scam.check(Halpernw2_scam)
    #weird residuals
    plot.scam(Halpernw2_scam, pages=1, 
              residuals=TRUE, 
              pch=1, cex=1,
              shift = coef(Halpernw2_scam)[1])
    
    Halpernw2_scam_p <- summary(Halpernw2_scam)$s.table[,"p-value"][[2]]
    Halpernw2_scam_p
    #PEAKED
    
    ###Quadtratic ###
    quad_model_Halpernw2 <- lm(heterogeneity~t_disturbance+I(t_disturbance^2), data=HalpernSpies1995_w2)
    plot(quad_model_Halpernw1)
    peaktest_quad(x = HalpernSpies1995_w2$t_disturbance,y=HalpernSpies1995_w2$heterogeneity )
    # peaked
    
    ###Two-Lines ####
    peaktest_twolines(x = HalpernSpies1995_w2$t_disturbance,y=HalpernSpies1995_w2$heterogeneity )
    #NOT PEAKED
    Halpernw2_twoLinesModel <- peaktest_twolines_returnModel(x = HalpernSpies1995_w2$t_disturbance,y=HalpernSpies1995_w2$heterogeneity )
    
    ##Peak quantification ####
    data_halperm_w2_mod <- tibble(x = HalpernSpies1995_w2$t_disturbance,y=HalpernSpies1995_w2$heterogeneity ,
                               true_val = NA)
    Halpern_w2_mod <- gam(y~s(x, bs="ad", k=15, m=3),method="REML", data=data_halperm_w2_mod)
    gam.check(Halpern_w2_mod)
    summary(Halpern_w2_mod)
    #Doesn't look great
    plot(Halpern_w2_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Halpern_w2_mod)[1])
    
    model_summary <- create_peakfinder_input(mod = Halpern_w2_mod, data = data_halperm_w2_mod,step_size = diff(range(data_halperm_w2_mod$x))/1000)
    set.seed(4)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    # No Peak found

# 4: Marra 2014 ####
    Marra2014 <- read.csv("data/processed/CLEAN_Marra2014.csv")
    plot(Marra2014$sp_rarefaction~Marra2014$mort_mean)
    
    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    Marra2014_scam <- scam(sp_rarefaction~s(mort_mean,k=10, bs="mpd", m=2)+
                             s(mort_mean, k=10, bs="ps",  m =c(2,1)),
                           data=Marra2014)
    scam.check(Marra2014_scam)
    #looks really good
    plot.scam(Marra2014_scam, pages=1, 
              residuals=TRUE, 
              pch=16, cex=1,
              shift = coef(Marra2014_scam)[1])
    
    Marra2014_scam_p <- summary(Marra2014_scam)$s.table[,"p-value"][[2]]
    Marra2014_scam_p
    #not peaked (but really close tho)
    
    ###Quadtratic ###
    quad_model_Marra2014 <- lm(sp_rarefaction~mort_mean+I(mort_mean^2), data=Marra2014)
    plot(quad_model_Marra2014)
    summary(quad_model_Marra2014)
    peaktest_quad(x = Marra2014$mort_mean,y=Marra2014$sp_rarefaction )
    # peaked
    
    
    ###Two-Lines ####
    peaktest_twolines(x = Marra2014$mort_mean,y=Marra2014$sp_rarefaction )
    #NOT PEAKED
    Marra2014_twoLinesModel <- peaktest_twolines_returnModel(x = Marra2014$mort_mean,y=Marra2014$sp_rarefaction )
    Marra2014_lm1_predict <- predict(Marra2014_twoLinesModel$lm1, Marra2014, interval="confidence", level=0.95)
    Marra2014_lm2_predict <- predict(Marra2014_twoLinesModel$lm2, Marra2014, interval="confidence", level=0.95)
    
    
    Marra2014_TL1 <- fitted(Marra2014_twoLinesModel$lm1)
    Marra2014_TL1 <- Marra2014_TL1[which(Marra2014$mort_mean<=24.77648 )]
    Marra2014_TL2 <- fitted(Marra2014_twoLinesModel$lm2)
    Marra2014_TL2 <- Marra2014_TL2[which(Marra2014$mort_mean>=24.77648 )]
    
    plot(sp_rarefaction~mort_mean, data=Marra2014)
    lines(Marra2014_TL1~Marra2014$mort_mean[which(Marra2014$mort_mean<=24.77648 )],col="red", lwd=2, lty=3)
    lines(Marra2014_lm1_predict[which(Marra2014$mort_mean<=24.77648 ),2]~Marra2014$mort_mean[which(Marra2014$mort_mean<=24.77648 )],col="grey", lwd=2, lty=3)
    lines(Marra2014_lm1_predict[which(Marra2014$mort_mean<=24.77648 ),3]~Marra2014$mort_mean[which(Marra2014$mort_mean<=24.77648 )],col="grey", lwd=2, lty=3)
    lines(Marra2014_TL2~Marra2014$mort_mean[which(Marra2014$mort_mean>=24.77648 )],col="red", lwd=2, lty=3)
    lines(Marra2014_lm2_predict[which(Marra2014$mort_mean>=24.77648 ),2]~Marra2014$mort_mean[which(Marra2014$mort_mean>=24.77648 )],col="grey", lwd=2, lty=3)
    lines(Marra2014_lm2_predict[which(Marra2014$mort_mean>=24.77648 ),3]~Marra2014$mort_mean[which(Marra2014$mort_mean>=24.77648 )],col="grey", lwd=2, lty=3)
    
    ##Peak quantification ####
    data_Marra2014_mod <- tibble(x = Marra2014$mort_mean , y=Marra2014$sp_rarefaction ,
                                  true_val = NA)
    Marra_mod <- gam(y~s(x, bs="ad", k=20, m=3),method="REML", data=data_Marra2014_mod)
    appraise(Marra_mod)
    summary(Marra_mod)
    #Not the best...
    plot(Marra_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Marra_mod)[1])
    
    model_summary <- create_peakfinder_input(mod = Marra_mod, data = data_Marra2014_mod,step_size = diff(range(data_Marra2014_mod$x))/1000)
    set.seed(5)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=FALSE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    # No Peak found

# 5: Wilson and Keddy 1988 ####
    WilsonKeddy1988 <- read.csv("data/processed/CLEAN_WilsonKeddy1988.csv")
    WilsonKeddy1988$ln_somc_perc <- log(WilsonKeddy1988$somc_perc)
    plot(WilsonKeddy1988$n_sp~WilsonKeddy1988$ln_somc_perc)

    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    WilsonKeddy1988_scam <- scam(n_sp~s(ln_somc_perc,k=10, bs="mpi", m=2)+
                             s(ln_somc_perc, k=10, bs="ps",  m =c(2,1)),
                           data=WilsonKeddy1988)
    scam.check(WilsonKeddy1988_scam)
      #looks normally distributed, even if should be poisson, but its a pain in the butt to do poisson
    #looks worse if i do poisson
    plot.scam(WilsonKeddy1988_scam, pages=1, 
              residuals=TRUE, 
              pch=1, cex=1,
              shift = coef(WilsonKeddy1988_scam)[1])
    
    WilsonKeddy1988_scam_p <- summary(WilsonKeddy1988_scam)$s.table[,"p-value"][[2]]
    WilsonKeddy1988_scam_p
    #not peaked (but really close tho)
    
    ###Quadtratic ###
    quad_model_WilsonKeddy1988 <- glm(n_sp~ln_somc_perc+I(ln_somc_perc^2), data=WilsonKeddy1988)
    summary(quad_model_WilsonKeddy1988)
    
    plot(quad_model_WilsonKeddy1988)
    plot((WilsonKeddy1988$n_sp~WilsonKeddy1988$ln_somc_perc))
    lines((quad_model_WilsonKeddy1988$fitted.values)~WilsonKeddy1988$ln_somc_perc, col="red", lwd=2)
    peaktest_quad(x = WilsonKeddy1988$ln_somc_perc,y=WilsonKeddy1988$n_sp )
    # peaked
    
    plot(WilsonKeddy1988$n_sp~WilsonKeddy1988$somc_perc)
    lines(predict(quad_model_WilsonKeddy1988)~(WilsonKeddy1988$somc_perc))
    
    ###Two-Lines ####
    peaktest_twolines(x = WilsonKeddy1988$ln_somc_perc,y=WilsonKeddy1988$n_sp )
    
    #NOT PEAKED
    WilsonKeddy1988_twoLinesModel <- peaktest_twolines_returnModel(x = WilsonKeddy1988$ln_somc_perc,y=WilsonKeddy1988$n_sp )
    summary(WilsonKeddy1988_twoLinesModel$lm1)
    summary(WilsonKeddy1988_twoLinesModel$lm2)
    
    
    plot(WilsonKeddy1988$n_sp~WilsonKeddy1988$somc_perc)
    lines(predict(WilsonKeddy1988_twoLinesModel$lm1)[1:78]~(WilsonKeddy1988$somc_perc[1:78]), col="red")
    lines(predict(WilsonKeddy1988_twoLinesModel$lm2)[79:197]~(WilsonKeddy1988$somc_perc[79:197]), col="red")
    
    
    
    
    ##Peak quantification ####
    data_WilsonKeddy1988_mod <- tibble(x = WilsonKeddy1988$ln_somc_perc,y=WilsonKeddy1988$n_sp ,
                                 true_val = NA)
    WilsonKeddy_mod <- gam(y~s(x, bs="ad", k=20, m=3),method="REML", data=data_WilsonKeddy1988_mod[-1, ], family=nb)
    appraise(WilsonKeddy_mod)
    #Not bad
      #normal distribution shows better qqplot but has very heteroscedastic residuals
      #NB is better for residuals, but worse for qqplot, and first value is not appropriate 
    plot(WilsonKeddy_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(WilsonKeddy_mod)[1])
    
    model_summary <- create_peakfinder_input(mod = WilsonKeddy_mod, data = data_WilsonKeddy1988_mod,step_size = diff(range(data_WilsonKeddy1988_mod$x))/1000)
    set.seed(6)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    # No Peak found
 
# 6: Nilsson 1987 ####
    
    Nilsson1987 <- read.csv("data/processed/CLEAN_Nilsson1987_Fig2a.csv")
    plot(Nilsson1987$num_spp~Nilsson1987$vel_m_s)
    
    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    Nilsson1987_scam <- scam(num_spp~s(vel_m_s,k=10, bs="mpd", m=2)+
                                   s(vel_m_s, k=10, bs="ps",  m =c(2,1)),family=poisson,
                                 data=Nilsson1987)
    scam.check(Nilsson1987_scam)
    #looks really good, even if should be poisson, but its a pain in the butt to do poisson
    #looks worse if i do poisson
    plot.scam(Nilsson1987_scam, pages=1, 
              residuals=TRUE, 
              pch=1, cex=1,
              shift = coef(Nilsson1987_scam)[1])
    
    Nilsson1987_scam_p <- summary(Nilsson1987_scam)$s.table[,"p-value"][[2]]
    Nilsson1987_scam_p

    ###Quadtratic ###
    quad_model_Nilsson1987 <- glm(num_spp~vel_m_s+I(vel_m_s^2), data=Nilsson1987)
    plot(quad_model_Nilsson1987)
    Nilsson1987_quadpeaktest <- peaktest_quad_pois_return(x = Nilsson1987$vel_m_s,y=Nilsson1987$num_spp )
    # peaked
    
    ###Two-Lines ####
    Nilsson1987_twolines <- peaktest_twolines_pois(x = Nilsson1987$vel_m_s,y=Nilsson1987$num_spp )
    #NOT PEAKED
    peaktest_twolines_pois(x = Nilsson1987$vel_m_s,y=Nilsson1987$num_spp )

    ##Peak quantification ####
    data_Nilsson1987_mod <- tibble(x = Nilsson1987$vel_m_s,y=Nilsson1987$num_spp ,
                                       true_val = NA)
    Nilsson_mod <- gam(y~s(x, bs="ad",k=20, m=3),method="REML", data=data_Nilsson1987_mod, family=poisson)
    gam.check(Nilsson_mod)
    #Not bad
    plot(Nilsson_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Nilsson_mod)[1])
    
    model_summary <- create_peakfinder_input(mod = Nilsson_mod, data = data_Nilsson1987_mod,step_size = diff(range(data_Nilsson1987_mod$x))/1000)
    set.seed(7)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    # No Peak found
    
# 7: Keddy 1983 ####   
    
    Keddy1983 <- read.csv("data/processed/CLEAN_Keddy1983_Fig4.csv")
    plot(Keddy1983$n_spp~Keddy1983$exp)
    
    ## a: Peak detection ####
    ### Scam ####
    #decisded on increasing
    Keddy1983_scam <- scam(n_spp~s(exp,k=10, bs="mpd", m=2)+
                               s(exp, k=10, bs="ps",  m =c(2,1)),family=poisson,
                             data=Keddy1983)
    scam.check(Keddy1983_scam)
    #looks really good, even if should be poisson, but its a pain in the butt to do poisson
    #looks worse if i do poisson
    plot.scam(Keddy1983_scam, pages=1, 
              residuals=TRUE, 
              pch=1, cex=1,
              shift = coef(Keddy1983_scam)[1])
    
    Keddy1983_scam_p <- summary(Keddy1983_scam)$s.table[,"p-value"][[2]]
    Keddy1983_scam_p
    #not peaked (but really close tho)
    
    ###Quadtratic ###
    quad_model_Keddy1983 <- glm(n_spp~exp+I(exp^2), data=Keddy1983, family=poisson)
    plot(quad_model_Nilsson1987)
    Keddy1983_quadpeaktest <- peaktest_quad_pois_return(x = Keddy1983$exp,y=Keddy1983$n_spp )
    # not peaked
    
    ###Two-Lines ####
    Keddy1983_twolines <- peaktest_twolines_pois(x = Keddy1983$exp,y=Keddy1983$n_spp )
    #NOT PEAKED
    
    ##Peak quantification ####
    data_Keddy1983_mod <- tibble(x = Keddy1983$exp,y=Keddy1983$n_spp ,
                                   true_val = NA)
    Keddy_mod <- gam(y~s(x, bs="ad", k=20, m=3),method="REML", data=data_Keddy1983_mod, family=poisson)
    gam.check(Keddy_mod)
    summary(Keddy_mod)
    #Not bad
    plot(Keddy_mod,
         residuals=TRUE, 
         pch=1, cex=1,
         shift = coef(Keddy_mod)[1])
    
    model_summary <- create_peakfinder_input(mod = Keddy_mod, data = data_Keddy1983_mod,step_size = diff(range(data_Keddy1983_mod$x))/1000)
    set.seed(8)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=TRUE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    # No Peak found 
    
################################################################################################
    # END OF CODE ####  
################################################################################################
    