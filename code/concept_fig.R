# concept_fig.R
# Author: Natalie Dupont
# Code for functions peaktest_quad, reg2 and peaktest_twolines written by Uri Simonsohn (2018)
  # Source: https://osf.io/wdbmr
# Date created: 2023-02-02
# Date last modified: 2023-04-27
# BIOL 490, Pedersen Lab, Concordia University
# Description: This code is used to generate conceptual figures to highlight the 
# importance of this study

########################################################################################
#Outline

# 0: Loading necessary function
# 1: Quadratic method applied to a saturating relationship
# 2: Quadratic method with skewed data
# 3: Two-Lines regression for both scenario 1 and 2


########################################################################################
# 0: Loading necessary functions for two-lines models
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


set.seed(12)

x <- 1:100
par(mfrow=c(1, 2))

# 1 : Quadratic method applied to a saturating relationship ####
  #True relationship
  y_log <- log(x)

  #normally-distributed random variation
  y_log_noise <- y_log+rnorm(n=length(y_log), mean=0, sd=0.2)
  plot(y_log_noise~x, cex=0.1,
       ylab="", xlab="", 
       main="a) False Positive",
       xaxt='n', yaxt='n')
  title(ylab="Outcome", mgp=c(0.5,0,0),cex.lab=1.2)
  title(xlab="Predictor", mgp=c(0.5,0,0),cex.lab=1.2)
  
  lines(y_log~x)
  
  lmquad_log <- lm(y_log_noise~x+I(x^2))
  summary(lmquad_log)
  lmquad_log
  #y=-0.004472x^2+0.0738282x+1.244808
  # dy/dt=-0.004472*2*x+0.0738282
  # dy/dt=0=-0.004472*2*x+0.0738282
  vertex_x <- (-lmquad_log$coefficients[2])/(lmquad_log$coefficients[3]*2)
  vertex_x <- data.frame(x=vertex_x)
  vertex_y <- predict(lmquad_log,newdata =vertex_x )
  
  lines(lmquad_log$fitted.values~x, col="firebrick")
  points(vertex_y~vertex_x$x[1], col="firebrick", pch=17, cex=1.4)
  # legend("bottomright",
  #        legend=c("True Relationship", "Quadratic Regression Model"),
  #        col=c("black", "firebrick"),
  #        lty=1, cex=1)

# 2: Quadratic method with skewed data #### 
  #True relationship
    y_beta <- dbeta(x=x/100, shape1=8, shape2=2)
  #Adding residuals
    y_beta_noise <- y_beta+rnorm(n=length(y_log), mean=0, sd=0.2)
    
    plot(y_beta_noise~x, cex=0.1,
         ylab="", xlab="", 
         main="b) False Negative",
         xaxt='n', yaxt='n')
    title(ylab="Outcome", mgp=c(0.5,0,0),cex.lab=1.2)
    title(xlab="Predictor", mgp=c(0.5,0,0),cex.lab=1.2)
    lines(y_beta~x)
    
    #Fitting quadratic model
    lmquad_beta <- lm(y_beta_noise~x+I(x^2))
    lmquad_beta
    summary(lmquad_beta)
    lines(lmquad_beta$fitted.values~x, col="firebrick")
    points(y=max(y_beta), x=x[which(y_beta==max(y_beta))], pch=17, cex=1.4)
    # legend("topleft",
    #        legend=c("True Relationship", "Quadratic Regression Model"),
    #        col=c("black", "firebrick"),
    #        lty=1,cex=1)

#3: Two-Lines regression for both scenario 1 and 2 ####
  # Saturating
    plot(y_log_noise~x, cex=0.1,
         ylab="", xlab="", 
         main="a) No Peak Detected",
         xaxt='n', yaxt='n')
    title(ylab="Outcome", mgp=c(0.5,0,0),cex.lab=1.2)
    title(xlab="Predictor", mgp=c(0.5,0,0),cex.lab=1.2)
    lines(y_log~x)
    
    two_lines_log <- peaktest_twolines_returnModel(x=x, y=y_log_noise)
    y_log_fitted1 <- fitted(two_lines_log$lm1)
    y_log_fitted1 <- y_log_fitted1[1:88]
    y_log_fitted2 <- fitted(two_lines_log$lm2)
    y_log_fitted2 <- y_log_fitted2[89:100]
    lines(y_log_fitted1~x[1:88],col="blue", lwd=2)
    lines(y_log_fitted2~x[89:100],col="red", lwd=2)
    
  # Skewed
    plot(y_beta_noise~x, cex=0.1,
         ylab="", xlab="", 
         main="b) Peak Detected",
         xaxt='n', yaxt='n')
    title(ylab="Outcome", mgp=c(0.5,0,0),cex.lab=1.2)
    title(xlab="Predictor", mgp=c(0.5,0,0),cex.lab=1.2)
    lines(y_beta~x)
    
    two_lines_beta <- peaktest_twolines_returnModel(x=x, y=y_beta_noise)
    y_beta_fitted1 <- fitted(two_lines_beta$lm1)
    y_beta_fitted1 <- y_beta_fitted1[1:86]
    y_beta_fitted2 <- fitted(two_lines_beta$lm2)
    y_beta_fitted2 <- y_beta_fitted2[87:100]
    lines(y_beta_fitted1~x[1:86],col="blue", lwd=2)
    lines(y_beta_fitted2~x[87:100],col="red", lwd=2)

########################################################################################################
#End of Code
########################################################################################################
