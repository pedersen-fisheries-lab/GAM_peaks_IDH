# Title: test_sims 
# Author: Natalie Dupont
# Code for functions peaktest_quad, reg2 and peaktest_twolines written by Uri Simonsohn (2018)
  # Source: https://osf.io/wdbmr
# Date created: 2022-09-30
# Date last modified: 2023-04-20
# BIOL 490, Pedersen Lab, Concordia.

# Outline ####
  # 1: Loading packages
  # 2: loading functions
    # 2-1: Peak detection test performance functions
      # 2-1a: Loading scam ANOVA function
      # 2-1b: loading quadratic regression function peaktest_quad (flood.u in source code) source: https://osf.io/wdbmr
      # 2-1c: Loading two-lines function peaktest_twolines (reg2hood in source code) source: https://osf.io/wdbmr
    # 2-2: Peak quantification functions
      # 2-2 a: Pedersen gam function
  # 3: importing sim data
  # 4: Running simulations for peak detection
    # 4a: Testing peak detection methods on each sample using parallel processing
    # 4b: Exporting sim analysis results to a csv
  # 5: Data analysis for peak quantification
    # Measure all peaked relationships, test non-peaked, but cannot assess accuracy and precision
    # 5a: running peak quantifying function on each sample using loops
    # 5b: Analysing peak quantification performance
      # 5 b 1: Analyzing precision
      # 5b 2: analyzing accuracy
 
# 1: Loading Packages ####
  #library(mgcv)#builds off of nlme not necessary since scam is beaing loaded
  library(scam)# builds off of mgcv-nlme, for shape-constrained gams
  library(ggplot2) #for visualizing
  library(dplyr) #for data handling
  library(R.utils) #for timeout function
  library(pracma)#for tic toc
  library(foreach) #parallelized for loop
  library(parallel) #to set up cores
  library(doParallel) #to register cores
  library(purrr)

# 2:  loading functions ####
## 2-1: Peak detection functions ####
  ### 2-1a: Scam ANOVA function ####
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
  ###2-1 d: parallelized analysis of all peak detection methods ####
    #parallel analysis of all three peak detection methods on multiple simulated datasets
    #pass it a dataframe object of samples organized by sample ID
    #outputs a dataframe of sampleID, param, error, scam_p, scam_u, quad_u, rmax_u, rmed_u, rprop_u, lines3_u, and func
    #
    #parallelized analysis for each sample. Results output as a list of numerical vectors
    #pass it a dataframe with multiple samples, each with a unique ID number, and a log name in quotation marks
    peaktest_all_prll <- function(sim_data, logname, scamType="m21"){ 
      sim_data <- sim_data
      scamType <- scamType
      
      #using the user-defined logname to set up a new parallel processing logging file
      reportFileName <- paste0("data/simulation/PP_log-", logname, ".txt")
      file.create(reportFileName)
      
      startTime <- as.numeric(Sys.time())

      #Parallel processing using the sim_data passed in the function
      sim_data_detect_results_list <- foreach (i = unique(sim_data$sampleID),    #iterating for each ID number
                                               .packages=c('scam', 'dplyr', 'R.utils','pracma'),        #Packages necessary for loop processing
                                               .combine="rbind")%dopar% {                      #merges list output into a dataframe, where each iteration output is a row
         
          cat(paste0(Sys.time(), ",Iteration Start,", "Sample ID,", i,  ",Time since foreach loop start,", as.numeric(Sys.time())-startTime, "\n" ),
              file=reportFileName, append=TRUE)
                                                                                       
          ##### selecting only the test data for iteration i
           test_data <- sim_data%>%
             filter(sampleID==i)
           x <- test_data$x
           y <- test_data$y
         
         #### scam anova test (2-1 a)
           #scam_results <- NA #need to initialize. In regular loop, the previous iteration's values would stay in scam_results, and if processing times out, would be assigned to a new section
           #no need to initialize with scam_results being external
           tic()
            scam_results <-  peaktest_scam(x=x, y=y, trend="increasing", scamType=scamType)
           #return value vector: u (boolean significance), and p-value of unconstrained ps spline
           time_Scam <- toc()
           
           cat(paste0(Sys.time(), ",Scam test,", "Sample ID,", i,  ",Duration of scam analysis,", time_Scam, "\n" ),
               file=reportFileName, append=TRUE)
         
         #### quadratic regression (2-1 b)
           tic()
           quad_results <- peaktest_quad(x=x, y=y)
           #returns a numerical value u (boolean significance) under the criteria 
           #"if (!is.na(jn1) & !is.na(jn2)) if (jn1>min(x) & jn1<max(x) & jn2>min(x) & jn2<max(x) & jn1>jn2) u=1
           #FIGURE OUT HOW TO UNDERSTAND THIS
         
         #### two_lines (2-1 c) 
           twolines_results <- peaktest_twolines(x=x, y=y)
           # returns boolean u significance for 3 different rHood methods and 3lines test
           #res=c(rmax$u.sig,rmed$u.sig, rprop$u.sig,u.3lines)  
           #rmax: xc=max(yhat)
           #rmed: xc=median(yflat)- yflat defined as the range where yhat+1SE>=max(yhat)
           #rprop: uses robin hood algorithm on the rmed break point to determine a new breakpoint
           #u.3lines: uses 3-lines method (not sure beyond that...)
           
         
         #returing a vector with all the values in the correct order
         #func will be added later based on sample ID
           results <- c(test_data$sampleID[1], 
                        test_data$param[1],
                        test_data$error[1],
                        scam_results[2] ,  #scam_p
                        scam_results[1],  #scam_u
                        quad_results,     #quad_u
                        twolines_results[1],   #rmax_u
                        twolines_results[2],    #rmed_u
                        twolines_results[3], #rprop_u
                        twolines_results[4]) #lines3_u
         time_simon <- toc()
         gc()
         cat(paste0(Sys.time(), ",All analyses complete,", "Sample ID,", i,  ",Duration of Simon analyses,", time_simon, "\n" ),
             file=reportFileName, append=TRUE)
         
         #returing concatenated analysis results for sample ID=i
         return(results)
                                               }
      
      #reassigning column names
      colnames(sim_data_detect_results_list) <- c("sampleID", "param", "error", 
                                                  "scam_p", "scam_u", 
                                                  "quad_u",
                                                  "rmax_u",
                                                  "rmed_u",
                                                  "rprop_u",
                                                  "lines3_u")
      #adding func
      func_id_mergedf <- sim_data%>%
        group_by(func, param, rep, error, sampleID)%>%
        summarize()%>%
        ungroup()
      
      sim_data_detect_results_df <- merge(sim_data_detect_results_list, func_id_mergedf[,c(1, 5)], by="sampleID")
      
      #returing results for all samples in sim_data passed through the function
      return(sim_data_detect_results_df)
    }

##2-2 Peak quantification functions ####
    source("code/parametrized_peak_finder.r")
    dev.off()
    
# 3: importing sim data  ####
  sim_data <- read.csv("data/simulation/CLEAN_simulationData.csv")
  head(sim_data)

# 4: Running simulations for peak detection####
  ## 4a: Testing peak detection methods on each sample with parallel processing####
    #set number of cores
      my.cluster <- makeCluster(detectCores()-2)

    #register cores
      registerDoParallel(my.cluster)
      
    #check if it is registered 
      foreach::getDoParRegistered()
      
    #how many workers are available?
      foreach::getDoParWorkers()
      
      # Exporting functions necessary for parallel processing within peaktest_all_prll
      clusterExport(my.cluster, c("peaktest_scam", "peaktest_quad", "peaktest_twolines", "reg2"))
      
      #Analysing each function type separately 
      #analysing using m=c(2, 1)
      for(s in unique(sim_data$scenario)){
          sim_data_subset_subset_param <- sim_data%>%
                filter(scenario==s)
          #setting variables for object naming convention
          f <- sim_data_subset_subset_param$func[1]
          e <- sim_data_subset_subset_param$error[1]
          p <-sim_data_subset_subset_param$param[1]
          
          tic()
          #Perfomring paralell analysis of all unique samples within the f/e subset. 
          sim_data_subset_detect_results <- peaktest_all_prll(sim_data_subset_subset_param, logname=paste0(f, "_e", e, "_p", p), scamType = "m21")
          toc()
       
          #saving the file under the same name just in case
          saveRDS(sim_data_subset_detect_results, file=paste0("data/simulation/detect_results_F", f, "_E", e, "_P", p, ".RDS"))
      }
      
      #importing and merging RDS for m=c(2, 1)
      detect_results_file_names <- list.files(path="./data/simulation", pattern="^detect_results_F", full.names=TRUE)
      detect_results_list <- lapply(detect_results_file_names, function(x) readRDS(x))
      detect_results_list
      
      detect_results_list <- detect_results_file_names%>%
        map_dfr(readRDS)
      
      detect_results_list <-lapply(detect_results_file_names, function(x) readRDS(paste0("data/simulation/",x)))
      detect_results_all <- do.call(rbind,detect_results_list)
      
      #merging all dataframes into a single dataframe
      detect_results_all <- do.call(rbind, mget(ls(pattern="^detect_results_F")))
      head(detect_results_all)
      summary(detect_results_all)
      
      #exporting detect_results_all to csv
      write.csv(detect_results_all,"data/simulation/detect_results.csv", row.names = FALSE)
      
# 5: Data analysis for peak quantification ####
  ## 5a: selecting necessary datasets ####
    sim_data_quant <- sim_data
    head(sim_data_quant)
    length(sim_data_quant$func)
    
  ## 5b: running peak quantifying function on each sample using loops####
    #set number of cores
    my.cluster <- makeCluster(detectCores()-2)
    
    #register cores
    registerDoParallel(my.cluster)
    
    #check if it is registered 
    foreach::getDoParRegistered()
    
    #how many workers are available?
    foreach::getDoParWorkers()
    
    #parallel loop for true peaks
    for (s in unique(sim_data_quant$scenario[1])){
      sim_data_quant_subset <- sim_data_quant%>%
        filter(scenario==s)
      reportFileName <- paste0("data/simulation/PP_quant_full_",s,".txt")
      
      file.create(reportFileName)
      
      iterationStartTime <- as.numeric(Sys.time())
      
      tic()
      sim_data_quant_results <- foreach(i=unique(sim_data_quant_subset$sampleID),
                                        .packages=c('mgcv','dplyr','R.utils', 'tibble','pracma'))%dopar%{
          #setting seed as the ID number of the sample
          set.seed(i)
                                          
          #Subsetting just the dataset of interest
          test_data <- sim_data_quant_subset%>%
            filter(sampleID==i)
          
          cat(paste0(Sys.time(), ",Iteration Start,", "Sample ID,", i,  ",Time since foreach loop start,", as.numeric(Sys.time())-iterationStartTime, "\n" ),
              file=reportFileName, append=TRUE)
          
          #reformatting the dataset for it to work with parametrized peakfinder code
          data <- tibble(x = test_data$x,y= test_data$y,
                         true_val = test_data$y_mean)
          
          #running gam
          tic()
          mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
          time_mod <- toc()
          
          cat(paste0(Sys.time(), ",Mod constructed,", "Sample ID,", i,  ",Duration of mod construction,", time_mod, "\n" ),
              file=reportFileName, append=TRUE)
          
          #formatting data to be compatible with parametrized peakfinder code
          tic()
          model_summary <- create_peakfinder_input(mod = mod, data = data)
          
          model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
          
          possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
          
          #Peaks true-false series
          peak <- possible_peaks$candidate_peaks
          
          time_peakfinder <- toc()
          
          #Characterizing extrema as peaks or troughs
            rle <- rle(peak)
            CIs_index <- which(rle$values)
            CI_first_deriv <- possible_peaks$CI$test_1st_deriv_CI
            
            # counting the Number of extrema CI
            n_extr <- length(CIs_index)
            
            if(n_extr!=0){
              CIs_lengths <-rle$lengths[c(CIs_index)]
              CI_width <- CIs_lengths*0.001
              
              CIs_start <- cumsum(rle$lengths)[CIs_index-1]+1
              CIs_end <- cumsum(rle$lengths)[CIs_index]
              CI_start_deriv <- CI_first_deriv[c(CIs_start-1),c(1, 3)]
            } else {
              CI_width <- 0
              CIs_start <- 0
              CIs_end <- 0
            }
            
            if(n_extr==0){
              peak_or_trough <- 0
            } else if (n_extr==1){
              peak_or_trough <- sign(CI_start_deriv[1])
            } else {
              peak_or_trough <- sign(CI_start_deriv[,1])
            }
          
          #creating a list of final output objects
          gam_output <- list()
          gam_output$sampleID <- i
          mod_summary <- summary(mod)
          gam_output$mod_summary <- c(mod_summary$s.table[[4]], #Smooth p-value
                                      mod_summary$r.sq,   #R-squared
                                      mod_summary$edf,    #edf
                                      mod_summary$dev.expl) #deviance explained

          gam_output$peaks <- peak
          
          gam_output$first_deriv <- CI_first_deriv
          
          gam_output$n_CI<- n_extr
          gam_output$CI_summary <- data.frame(peak_or_trough=peak_or_trough, 
                                              CI_start_ID=CIs_start, 
                                              CI_end_ID=CIs_end, 
                                              CI_width=CI_width)
            
          cat(paste0(Sys.time(), ",peak finder done,", "Sample ID,", i,  ",Duration of peak finding, ", time_peakfinder, "\n" ),
              file=reportFileName, append=TRUE)
          return(gam_output)
        }
      toc()
      
      #organizing results from all replicates of a scenario together
      sampleID<- unlist(lapply(sim_data_quant_results, function(l) l[[1]]), recursive=TRUE, use.names=TRUE)
      
      #Model summary results
      mod_summary <- as.data.frame(do.call(rbind,lapply(sim_data_quant_results, function(l) l[[2]])))
      colnames(mod_summary) <-  c("p_val","Rsq_adj", "edf","dev_expl")
      mod_summary <- cbind(data.frame(sampleID=sampleID),mod_summary)
      
      #Peaks results
      peaks_logi <- as.data.frame(do.call(cbind, lapply(sim_data_quant_results, function(l) l[[3]])))
      colnames(peaks_logi) <- as.character(sampleID)
      
      #First derivative CI results
      first_deriv_list <- lapply(sim_data_quant_results, function(l) l[[4]])
      names(first_deriv_list) <- sampleID
      
      #CI number results
      n_CI <- data.frame(sampleID=sampleID,
                         n_CI=unlist(lapply(sim_data_quant_results, function(l) l[[5]]), recursive=TRUE, use.names=TRUE))
      
      #CI summary result
      CI_summary <- lapply(sim_data_quant_results, function(l) l[[6]])
      names(CI_summary) <- sampleID
      CI_summary <- map_df(CI_summary, data.frame, .id = 'name')
      rownames(CI_summary) <- NULL
      colnames(CI_summary)[1] <- "sampleID"
      
      #saving all these results as RDS files for later use
      saveRDS(mod_summary, paste0("data/simulation/mod_sum_", s,".RDS"))
      saveRDS(peaks_logi, paste0("data/simulation/peak_pres_", s, ".RDS"))
      saveRDS(first_deriv_list, paste0("data/simulation/first_deriv_list_", s,".RDS"))
      saveRDS(n_CI, paste0("data/simulation/n_CI_", s,".RDS"))
      saveRDS(CI_summary, paste0("data/simulation/CI_summary_", s,".RDS"))
    }
    parallel::stopCluster(cl = my.cluster)

    ###importing and merging RDS ####
      #For mod_summary table
      mod_sum_filenames <- list.files(path="./data/simulation", pattern="^mod_sum_", full.names=TRUE)
      mod_sum_all <- mod_sum_filenames%>%
        map_dfr(readRDS)
      View(mod_sum_all)
  
        #exporting to csv
        write.csv(mod_sum_all,"data/simulation/quant_mod_sum_results.csv", row.names = FALSE)
        
      #For peak detection table
      quant_peaks_filenames <- list.files(path="./data/simulation", pattern="^peak_pres_", full.names=TRUE)
      
      quant_peaks <- quant_peaks_filenames%>%
        map_dfc(readRDS)
      
      #adding x values for reference. Ideally, should be added earlier on...
      quant_peaks_all <- as.data.frame(cbind(data.frame(x=c(seq(0, 1, by=0.001))), 
                                             quant_peaks))
      View(quant_peaks_all)
      
        #exporting  to csv
        write.csv(quant_peaks_all,"data/simulation/quant_peaks_results.csv", row.names = FALSE)
    
      #n_CI exporting to csv
      n_CI_filenames <- list.files(path="./data/simulation", pattern="^n_CI_", full.names=TRUE)
      n_CI_all <- n_CI_filenames%>%
        map_dfr(readRDS)
      View(n_CI_all)

        #exporting to csv
        write.csv(n_CI_all,"data/simulation/quant_n_CI_results.csv", row.names = FALSE)
        
      #CIsummary exporting to csv
      CI_summary_filenames <- list.files(path="./data/simulation", pattern="^CI_summary_", full.names=TRUE)
      CI_summary_all <- CI_summary_filenames%>%
        map_dfr(readRDS)
      View(CI_summary_all)
      
      #exporting to csv
      write.csv(CI_summary_all,"data/simulation/quant_CI_summary_results.csv", row.names = FALSE)
 
###################################################################################################      
  # END OF CODE ####         
###################################################################################################