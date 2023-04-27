# Title: imulations 
# Author: Natalie Dupont
# Date created:2022-09-30
# Date last modified: 2023-04-06
# BIOL 490, Pedersen Lab, Concordia University
# Description: This code serves to generate random samples for simulation analysis of 
# the scam peak detection, and gam derivative peak quantification methods. It also makes
# conceptual figures of this suite of random samples
########################################################################################
#Outline

# 0 : Loading packages
# 1 : Preparing simulation data
#   1a: Setting up dataframe parameters
#   1b: Generating a dataframe 
#   1c: scaling response variable names
# 2: exporting data
# 3: visualizing data


########################################################################################
#0: Loading packages ####
  library(dplyr) #for data management
  library(tidyr) #for data management
  library(ggplot2) #For visualization
  library(patchwork) #For visualization

########################################################################################
# 1: Preparing a simulation dataframe using expand_grid ####
  # Sample size of 20, x evenly distributed between 0 and 1
  # 6 curve types: beta, normal-varying mean, normal varying sd, saturating, linear, polynomial
  # 4 error levels: 0.2, 0.3, 0.4, 0.5
    #All error structures will be normally distributed
  # Parameter values: 2 for beta, 3 for normal-varying sd, 2 for normal varying mean, 3 for saturating, 1 for linear, 3 for polynomial
    #total situations to simulate: 56 scenarios to simulate
  #500 replicates per scenario

  ## 1 a: Setting up dataframe skeleton with parameters ####
    set.seed(20221110)

    n <- 20 
    x = seq(0,1,length = n)
      #Domain distribution is restricted to 0<=x<=1, and uniformely distributed
    error <- c(.2,.3,.4, 0.5)
    reps <- 500
    func=c("norm-var-sd", "norm-var-mean", "beta", "saturating", "linear", "polynomial") 
      #norm-var-sd: Normal distribution with constant mean but varying sd
        #tests the effect of width
      #norm-var-mean : Normal distribution with constant sd but varying mean
        #tests the effect of a non-centered peak
      #beta: varying tradeoff between shape 1 and shape 2
        #tests the effect of skew on test power
      #Saturating: asymptotic curve with varying slope- decelerating rate of change
      #Linear: 1 linear line 
      #Polynomial: Polynomial curves with varying exponents- accelerating rate of change
    
    param <- seq(2, 8, by=2)
      #Setting up a range of parameters to apply to all curve types
    
  ## 1 b: generating simulation dataset and computing response #### 
    # Generating all possible combinations of these parameters
    sim_data <-  expand_grid(func = func,
                           param=param,
                           rep = 1:reps,
                           error = error,
                           x=x)%>%
     #filtering out redundant curves
       filter(!(func=="norm-var-sd"&param==6))%>% #curve for param 6 and 8 are very close
       filter(!(func=="norm-var-mean"&(param==2|param==4)))%>% # symmetric, therefore redundant
       filter(!(func=="beta"&(param==2|param==4)))%>% # symmetric, therefore redundant
       filter(!(func=="saturating"&param==6))%>% #curve for param 6 and 8 are very close
       filter(!(func=="linear"&(param==4 | param==6 |param==8)))%>%  #All curves overlap once standardized, only need 1
       filter(!(func=="polynomial"&(param==6 )))%>% #curve for param 6 and 8 are very close
     #Calculating true relationship response variables using parameter values. 
    mutate(y_mean_unscaled = case_when(func=="norm-var-sd"~10*dnorm(x,mean=.5,sd =param/20),
                                       func=="norm-var-mean"~10*dnorm(x,mean=(param/10),sd =4/20),
                                       func=="beta"~dbeta(x, shape1=param, shape2=10-param),
                                       func=="saturating"~param*x/(1+(param-1)*2*x), 
                                       func=="linear"~(param)/2*x, 
                                       func=="polynomial"~x^(param)))
                                
    
  ## 1 c: Scaling all functions to 0<=y<=1 ####
    sim_data <- sim_data%>%
      group_by(func, param, error)%>%
      mutate(y_mean=(y_mean_unscaled-min(y_mean_unscaled)))%>% #shifting the whole range down vertically to start at 0
      mutate(y_mean=y_mean/max(y_mean), #scaling the whole range to have a maximum of 1
           y = y_mean + rnorm(n = n(),mean = 0, sd = error))%>%  #Generating response data with incorporated error levels
      ungroup()%>%
      as.data.frame(sim_data)#not necessary
    
    #removing unscaled values
    sim_data <- select(sim_data,-y_mean_unscaled)

    #adding sample ID, Curve ID and scenario (curve-error combination) ID
    sim_data <- sim_data%>%
      group_by( func, param, rep, error)%>%
      mutate(sampleID=cur_group_id(),
             CurveID=paste0(func, param), 
             scenario=paste0(CurveID, error))%>%
      ungroup()
    
    head(sim_data)
    
# 2:  exporting simulation data to data folder ####
    write.csv(sim_data,"data/simulation/CLEAN_simulationData.csv", row.names = FALSE)
    
# 3:  Visualizing sim data ####
    #filtering a subset for visualization
    sim_data_names <- filter(sim_data, rep==1, error==0.5)
    
    #Adding more appropriate curve ID for publication
    sim_data_names$CurveID2 <- NA
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="saturating8")] <- "D"    
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="saturating4")] <- "D"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="saturating2")] <- "D"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="linear2")] <- "C"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="polynomial2")] <- "A"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="polynomial4")] <- "A"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="polynomial8")] <- "A"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="norm-var-sd8")] <- "W"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="norm-var-sd4")] <- ""
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="norm-var-sd2")] <- "W"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="norm-var-mean6")] <- "L"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="norm-var-mean8")] <- "L"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="beta6")] <- "S"
    sim_data_names$CurveID2[which(sim_data_names$CurveID=="beta8")] <- "S"
    
    sim_data_names$CurveID3 <- NA
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="saturating8")] <- "3"    
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="saturating4")] <- "2"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="saturating2")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="linear2")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="polynomial2")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="polynomial4")] <- "2"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="polynomial8")] <- "3"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="norm-var-sd8")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="norm-var-sd4")] <- ""
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="norm-var-sd2")] <- "2"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="norm-var-mean6")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="norm-var-mean8")] <- "2"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="beta6")] <- "1"
    sim_data_names$CurveID3[which(sim_data_names$CurveID=="beta8")] <- "2"
    
    #Adding a marker for which curve shapes will be reported on
    sim_data_names$ofInterest <- NA
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="saturating8")] <- "Yes"    
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="saturating4")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="saturating2")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="linear2")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="polynomial2")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="polynomial4")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="polynomial8")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="norm-var-sd8")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="norm-var-sd4")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="norm-var-sd2")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="norm-var-mean6")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="norm-var-mean8")] <- "Yes"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="beta6")] <- "No"
    sim_data_names$ofInterest[which(sim_data_names$CurveID=="beta8")] <- "Yes"

    sim_data_names$refCurve <- NA
    sim_data_names$refCurve[which(sim_data_names$CurveID=="norm-var-sd4")] <- "Reference"
    refx <- filter(sim_data_names, refCurve=="Reference")$x
    refy <- filter(sim_data_names, refCurve=="Reference")$y_mean
    #true relationship curves
    
    
    #Plotting curves peaked
    plot_true_peaked <- ggplot(data=filter(sim_data_names,func=="beta"|func=="norm-var-mean"|func=="norm-var-sd"), aes(x=x, y=y_mean))+
      geom_line(aes(group=paste(param), colour=as.factor(CurveID3), alpha=as.factor(ofInterest)))+
      annotate(geom="line", x=refx, y=refy, colour="black")+
      annotate("text", x=0.22, y=0.5, label= "Ref")+
      facet_wrap(~factor(func, levels=c("norm-var-sd", "norm-var-mean", "beta")),
                 labeller=as_labeller(c("norm-var-sd"="W: Varying peak width", "norm-var-mean"="L: Varying peak location", "beta"="S: Varying skew")))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position="none")+
      geom_text(data=filter(sim_data_names, func=="beta"&y_mean<=0.6&y_mean>=0.4&x<0.7|((func=="norm-var-mean"|func=="norm-var-sd")&y_mean<=0.55&y_mean>=0.4&x<0.6)),
                aes(label=paste0(CurveID2, CurveID3), colour=as.factor(CurveID3), alpha=as.factor(ofInterest)), nudge_x=-0.04, y=0.5)+
      scale_x_continuous(breaks=c(0,0.5,1))+
      scale_y_continuous(breaks=c(0,0.5,1))+
      scale_alpha_discrete(range=c(.35,1))+
      labs(y="True Response Curve",
           x="Predictor", 
           title="b) Peaked Curves")
    plot_true_peaked
    
    #plotting curves not peaked
    plot_true_notpeaked <- ggplot(data=filter(sim_data_names,func=="polynomial"|func=="linear"|func=="saturating"), 
                                  aes(x=x, y=y_mean))+
      geom_line(aes(group=paste(param), colour=CurveID3, alpha=as.factor(ofInterest)))+
      labs(y="True Response Curve",
           x="Predictor",
           title="a) Non-Peaked Curves",
           colour="Curve ID")+
      facet_wrap(~factor(func, levels=c("saturating", "linear", "polynomial")),
                 labeller=as_labeller(c("saturating"="D: Decelerating rate of change", "linear"="C: Constant Rate of Change", "polynomial"="A: Accelerating Rate of Change")))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position="none")+
      geom_text(data=filter(sim_data_names, rep==1, error==0.5, func=="linear"&y_mean<=0.51&y_mean>=0.45&x<0.7|(func=="saturating"&y_mean<=0.75&y_mean>=0.7|func=="polynomial"&y_mean<=0.33&y_mean>=0.25)),
                aes(label=paste0(CurveID2, CurveID3), colour=CurveID3, alpha=as.factor(ofInterest)), nudge_x=-0.08)+
      scale_x_continuous(breaks=c(0,0.5,1))+
      scale_y_continuous(breaks=c(0,0.5,1))+
      scale_alpha_discrete(range=c(.35,1))
    plot_true_notpeaked
    
    #combing peaked and non-peaked curves to make 1 figure
    plot_true_notpeaked+plot_true_peaked+plot_layout(nrow=2)
    
    #plotting error levels
    plot_errors <- ggplot(data=filter(sim_data, rep==2, func=="norm-var-sd", param==4, error==0.2|error==0.5), aes(x=x, y=y_mean))+
      geom_line(aes(group=paste(param), colour="True Response Curve"))+
      geom_point(aes(group=paste(param), y=y))+
      # geom_ribbon(aes(ymin=y_mean-error, ymax=y_mean+error, group=paste(param), colour="Residual Error"), alpha=0.2, 
      #             fill=NA, lty=2)+
      labs(y="Response", 
           x="Predictor")+
      facet_wrap(~factor(error),
                         labeller=as_labeller(c("0.2"="Residual Error 0.2","0.3"="Residual Error 0.3", "0.4"="Residual Error 0.4", "0.5"="Residual Error 0.5")))+      
      theme_bw()+
      scale_colour_manual(name="",labels=c("True Response Curve"),
                          values=c("black"))+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            legend.position="bottom")+
      scale_x_continuous(breaks=c(0,0.5,1))+
      scale_y_continuous(breaks=c(0,0.5,1, 1.5, 2))
    plot_errors
    
    
  plot_simulations <- ggplot(data=filter(sim_data, scenario=="saturating40.3"))+
    geom_point(aes(x=x, y=y, colour=as.factor(rep)))+
    theme_bw()+
    theme(legend.position = "none")
  plot_simulations
  
  #Plotting curves peaked
  plot_true_peaked_mini <- ggplot(data=filter(sim_data_names,CurveID=="norm-var-sd4"), aes(x=x, y=y_mean))+
    geom_line(aes(group=paste(param), colour=as.factor(CurveID3), alpha=as.factor(ofInterest)))+
    annotate(geom="line", x=refx, y=refy, colour="black", linewidth=1.8)+
    facet_wrap(~factor(func, levels=c("norm-var-sd", "norm-var-mean", "beta")),
               labeller=as_labeller(c("norm-var-sd"="W: Varying peak width", "norm-var-mean"="L: Varying peak location", "beta"="S: Varying skew")))+
    theme_bw()+
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position="none")+
    scale_x_continuous(breaks=c(0,0.5,1))+
    scale_y_continuous(breaks=c(0,0.5,1))+
    scale_alpha_discrete(range=c(.35,1))+
    labs(y="True Response Curve",
         x="Predictor", 
         title="b) Peaked Curves")
  plot_true_peaked_mini
  
########################################################################################################
#End of Code
########################################################################################################