# sim_result_analysis.R
# Author: Natalie Dupont
# Date created: 2023-02-02
# Date last modified: 2023-04-27
# BIOL 490, Pedersen Lab, Concordia.
# Description: Analysis of results of methods applied to simulated data in 
# tests_sims_Servercopy.R. This test the type I error rate and power of the scam 
# detection method compared to the two-lines robin hood and quadratic regression methods.
# It also tests the performance of the gam peak quantification method as both a peak
# detection method and a quantification method, by assessing coverage.


#####################################################################################################################################
# Outline
  # 0: Loading packages and data
#Scam peak detection results
  #1 : Type 1 error analysis
    #1a: Rate of type I error 
    #1b: plotting type I error
    #1c Analysing example cases  
  #2: : Power analysis
    #2a: Rate of true positive
    #2b: plotting Power
    #2c: Analysing example cases
#Gam Peak quantification results 
  #3: Peak accuracy

#####################################################################################################################################
# 0: Loading packages and data####
  ##loading packages ####
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(scam)
    library(gratia)

  ##loading data ####
    #original simulation data
    sim_data <- read.csv("data/simulation/CLEAN_simulationData.csv")
    
    #scam peak detection results
    sim_data_detect_results<- read.csv("data/simulation/detect_results.csv")
    
    #adding IDs for easier handling
    sim_data_detect_results$scenario <- paste0(sim_data_detect_results$func, sim_data_detect_results$param, sim_data_detect_results$error)
    sim_data_detect_results$CurveID <- paste0(sim_data_detect_results$func, sim_data_detect_results$param)
    View(sim_data_detect_results)
    
      #Search for NA and NaN values
        sum(is.na(sim_data_detect_results$scam_p))
        #104 unconverged Scam results
        sum(is.nan(sim_data_detect_results$lines3_u))
    
      #Extraction only values that converged
        detect_results_noNA <- sim_data_detect_results%>%
          filter(!(is.na(scam_p)))
      #extraction only values that did not converge
        detect_results_NA <- sim_data_detect_results%>%
          filter((is.na(scam_p)))
        summary(detect_results_NA)

# Scam Peak detection results ####
# 1 Type I error ####
  ## 1a Type I error rate ####
    #Filtering only peakless curves  calculating the rate of type I error for each curve
    detect_results_noNA_curvetype_noPeak<- detect_results_noNA %>% 
      filter(func=="linear" | func=="polynomial"|func=="saturating") 
    
    detect_results_noNA_curvetype_noPeak$SN <- 1/detect_results_noNA_curvetype_noPeak$error
    detect_results_noNA_curvetype_noPeak$scenario <- paste0(detect_results_noNA_curvetype_noPeak$CurveID,detect_results_noNA_curvetype_noPeak$error)
    
    #calculating the rate of type I error for each curve
    typeIerrorRate <-detect_results_noNA_curvetype_noPeak%>%
      group_by(scenario, CurveID, func, param, error, SN)%>%
      summarise(sample_size=n(),
                scam_T1R=sum(scam_u)/length(scam_u),      
                quad_T1R=sum(quad_u)/length(quad_u),
                rmax_T1R=sum(rmax_u)/length(rmax_u),
                rmed_T1R=sum(rmed_u)/length(rmed_u),
                rprop_T1R= sum(rprop_u)/length(rprop_u),
                lines3_T1R= sum(lines3_u)/length(lines3_u))%>%                              
      ungroup()
      
  ## 1b: Plotting Type I error ####
    ###histogram of p-values for scam method ####
      histo_noPeak <- ggplot(data=detect_results_noNA_curvetype_noPeak)+
        geom_histogram(aes(x=scam_p, y=after_stat(density)), bins=20)+
        facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~SN)+
        labs(x="Scam test p-value",y="Density")+
      theme_bw()
      histo_noPeak
        #larger signal gives larger TypeI error inflation
      
    ### ECDF of p-values for scam method ####
      ecdf_noPeak <- ggplot(data=detect_results_noNA_curvetype_noPeak, aes(x=scam_p))+
        stat_ecdf()+
        facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~SN)+
        geom_abline(aes(intercept=0, slope=1), colour="firebrick", linetype=2)+
        theme_bw()+
        labs(x="Scam test p-value", y="CDF")
      ecdf_noPeak
      
      ecdf_noPeak_2 <- ggplot(data=detect_results_noNA_curvetype_noPeak, aes(x=scam_p))+
        stat_ecdf(aes(group=SN, colour=as.factor(round(SN, 2))))+
        facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~.)+
        geom_abline(aes(intercept=0, slope=1), colour="black", linetype=2)+
        theme_bw()+
        labs(x="Scam test p-value", y="CDF", colour="Signal-to-Noise")
        #theme(legend.position="bottom")
      ecdf_noPeak_2
      
      ecdf_noPeak_onlyextremes <- ggplot(data=filter(detect_results_noNA_curvetype_noPeak,
                                                     CurveID=="saturating8"|CurveID=="linear2"|CurveID=="polynomial8"), aes(x=scam_p))+
        stat_ecdf(aes(group=error, colour=as.factor(error)))+
        facet_grid(factor(CurveID, levels=c("saturating8","linear2","polynomial8"))~.,
                   labeller=as_labeller(c("saturating8"="NP1","linear2"="NP4","polynomial8"="NP7")))+
        geom_abline(aes(intercept=0, slope=1), colour="black", linetype=2)+
        theme_bw()+
        labs(x="Scam test p-value", y="ECDF", colour="Residual Error",
             title="a) ECDF of unconstrained smooth p-values")+
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(), 
              legend.position="bottom")+
        scale_colour_brewer(palette="Reds", direction=-1)+
        scale_x_continuous(breaks=c(0,0.5,1))+
        scale_y_continuous(breaks=c(0,0.5,1))
      #theme(legend.position="bottom")
      ecdf_noPeak_onlyextremes
    
    ###TypeI error rate by curve type ####
      plot_T1error <- ggplot(data=typeIerrorRate, aes(x=error))+
        geom_line(aes(y=scam_T1R, colour="Scam test"))+
        geom_line(aes(y=quad_T1R, colour="Quadratic Test"))+
        geom_line(aes(y=rprop_T1R, colour="Two-Lines test"))+    
        #geom_line(aes(y=falsePosRate, x=error, colour="Quant Peak Pres"), data=falsePosRate_quant)+
        geom_hline(yintercept=0.05, lty=2)+
        facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~.)+
        theme_bw()+
        scale_color_manual(values=c('red', "blue", "darkgreen", "purple"))+
        labs(x="Signal-to-Noise Ratio", y="Type I Error Rate")+
        scale_y_continuous(trans='log10') +
        scale_x_continuous(trans='log10')
      plot_T1error
      
      plot_T1error_3tests <- ggplot(data=filter(typeIerrorRate,
                                         CurveID=="saturating8"|CurveID=="linear2"|CurveID=="polynomial8"), aes(x=error))+
        geom_line(aes(y=scam_T1R, colour="Scam test"))+
        geom_line(aes(y=quad_T1R, colour="Quadratic Test"))+
        geom_line(aes(y=rprop_T1R, colour="Two-Lines test"))+
        geom_hline(yintercept=0.05, lty=5, colour="black", linewidth=0.75)+
        facet_grid(factor(CurveID, levels=c("saturating8", "linear2","polynomial8"))~.,
                   labeller=as_labeller(c("saturating8"="Decelerating Rate of Change","linear2"="Constant Rate of Change","polynomial8"="Accelerating Rate of Change")))+
        theme_bw()+
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor=element_line(colour="darkgrey"),
              panel.grid.major=element_line(colour="darkgrey"),
              legend.position="bottom")+
        scale_color_manual(values=c('purple', "red", "blue", "forestgreen"))+
        labs(x="Residual error", y="Rate of False-Positive Errors", colour="")+
        scale_y_continuous(trans='log10', breaks=c(0.005,  0.05, 0.5),
                           minor_breaks = c(0.005, 0.01,0.02,0.03,0.04,0.05, 0.1,0.2,0.3,0.4,0.5))+
        guides(colour=guide_legend(nrow=2,byrow=TRUE))
      plot_T1error_3tests
        #may be limited by the splines available (may work better if adaptive splines were present)
          #but its struggling with the polynomial, where the derivative remains constnat???
          #struggles a lot with saturating, which is in line with issues with adaptive splines

    #Reference plot of curves
      true_curves_noPeak <- ggplot(data=filter(sim_data, rep==1, error==0.5, func=="linear" | func=="polynomial"|func=="saturating"), aes(x=x, y=y_mean))+
        geom_line(aes(group=paste(param)))+
        labs(y="True curve")+
        facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~.)+
        theme_bw()
      true_curves_noPeak
      
      #don't really need to be put together
      ecdf_noPeak_onlyextremes+plot_T1error_publish

      
  ## 1c: Examples cases #### 
   ### False postives ####
    #filtering out cases of inflated type I error
      #selecting any datasets where scam_u=1, and only selecting saturating8 and 2, and polynomial8 and 2 options
      sim_data_falsePos <- sim_data%>%
        filter(sampleID%in%sim_data_detect_results[which(sim_data_detect_results$scam_u==1),]$sampleID)%>%
        filter(func=="saturating" | func=="polynomial",
               param==8 | param==2,
               sampleID%in%c(22815, 16509,23613,17485,24122,26875,18698))
      
      scam22815 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data_falsePos, sampleID==22815))
      summary(scam22815)
      plot(y~x, data=filter(sim_data_falsePos, sampleID==22815))
      plot(scam22815, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam22815)[1])
      plot(scam22815, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam22815)[1])
      scam.check(scam22815)
      appraise(scam22815) #histogram kinda sucks
      
      scam23613 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data_falsePos, sampleID==23613))
      summary(scam23613)
      plot(y~x, data=filter(sim_data_falsePos, sampleID==23613))
      plot(scam23613, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam23613)[1])
      plot(scam23613, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam23613)[1])
      scam.check(scam23613)
      appraise(scam23613) #histogram kinda sucks

      scam26875 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data_falsePos, sampleID==26875))
      summary(scam26875)
      plot(y~x, data=filter(sim_data_falsePos, sampleID==26875))
      plot(scam26875, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam26875)[1], all.terms = TRUE)
      plot(scam26875, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam26875)[1])
      scam.check(scam26875)
      appraise(scam26875) 
      
      scam16509 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data_falsePos, sampleID==16509))
      summary(scam16509)
      plot(y~x, data=filter(sim_data_falsePos, sampleID==16509))
      lines(y_mean~x, data=filter(sim_data_falsePos, sampleID==16509), lty=2)
      plot(scam16509, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam16509)[1])
      plot(scam16509, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam16509)[1])
      scam.check(scam16509)
      appraise(scam16509) 
      
      scam17485 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data_falsePos, sampleID==17485))
      summary(scam17485)
      plot(y~x, data=filter(sim_data_falsePos, sampleID==17485))
      plot(scam17485, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam17485)[1])
      plot(scam17485, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam17485)[1])
      scam.check(scam17485)
      appraise(scam17485) 
      
   ### True Negatives ####
      scam3306 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==3306))
      summary(scam3306)
      plot(y~x, data=filter(sim_data, sampleID==3306))
      plot(scam3306, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam3306)[1])
      plot(scam3306, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam3306)[1])
      scam.check(scam3306)
      appraise(scam3306)
      
      scam4205 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==4205))
      summary(scam4205)
      plot(y~x, data=filter(sim_data, sampleID==4205))
      plot(scam4205, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam4205)[1])
      plot(scam4205, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam4205)[1])
      scam.check(scam4205)
      appraise(scam4205)
      
      scam27756 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==27756))
      summary(scam27756)
      plot(y~x, data=filter(sim_data, sampleID==27756))
      plot(scam27756, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam27756)[1])
      plot(scam27756, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam27756)[1])
      scam.check(scam27756)
      appraise(scam27756)
    
      scam19619 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==19619))
      summary(scam19619)
      plot(y~x, data=filter(sim_data, sampleID==19619))
      plot(scam19619, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam19619)[1])
      plot(scam19619, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam19619)[1])
      scam.check(scam19619)
      appraise(scam19619)
      
      scam16667 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==16667))
      summary(scam16667)
      plot(y~x, data=filter(sim_data, sampleID==16667))
      plot(scam16667, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam16667)[1])
      plot(scam16667, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam16667)[1])
      scam.check(scam16667)
      appraise(scam16667)
      
      #finally a sensible fit?
      scam26692 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==26692))
      summary(scam26692)
      plot(y~x, data=filter(sim_data, sampleID==26692))
      plot(scam26692, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam26692)[1])
      plot(scam26692, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam26692)[1])
      scam.check(scam26692)
      appraise(scam26692)
      
      #Strange that the curve is flat...
      scam27855 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==27855))
      summary(scam27855)
      plot(y~x, data=filter(sim_data, sampleID==27855))
      plot(scam27855, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam27855)[1])
      plot(scam27855, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam27855)[1])
      scam.check(scam27855)
      appraise(scam27855)
      
      scam23509 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==23509))
      summary(scam23509)
      plot(y~x, data=filter(sim_data, sampleID==23509))
      plot(scam23509, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam23509)[1])
      appraise(scam23509)
      
# 2: power analysis ####
  ##2a: Rate of true positive ####
    #filtering results from peaked curves
    detect_results_noNA_curvetype_Peak <- detect_results_noNA %>% 
      filter(func=="beta" | func=="norm-var-mean"|func=="norm-var-sd")%>% #selecting only peaked curves
      group_by( func, param)%>%
      mutate(CurveID=paste0(func, param))%>%    #Curve ID
      mutate(SN=1/error)%>%                      #Signal-to noise ratio
      ungroup()
    
    power <- detect_results_noNA_curvetype_Peak%>%
      group_by(func, param, error)%>%
      summarise(sample_size=n(),
                scam_Pr=sum(scam_u)/length(scam_u),
                quad_Pr=sum(quad_u)/length(quad_u),
                rmax_Pr=sum(rmax_u)/length(rmax_u),
                rmed_Pr=sum(rmed_u)/length(rmed_u),
                rprop_Pr= sum(rprop_u)/length(rprop_u),
                lines3_Pr= sum(lines3_u)/length(lines3_u))%>%
      mutate(CurveID=paste0(func, param),
             SN=1/error,
             scenario=paste0(func, param, error))%>%
      ungroup()
    View(power)

  ## 2b:Plotting Power ####
    ###histogram of p-values for scam method ####
    histo_Peak <- ggplot(data=detect_results_noNA_curvetype_Peak)+
      geom_histogram(aes(x=scam_p, y=after_stat(density)), bins=20)+
      facet_grid(CurveID~SN)+
      labs(x="Scam test p-value",y="Density")
    histo_Peak
    #larger signal gives larger TypeI error inflation

    ### ECDF of p-values for scam method ####
    ecdf_Peak <- ggplot(data=detect_results_noNA_curvetype_Peak, aes(x=scam_p))+
      stat_ecdf()+
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~SN)+
      geom_abline(aes(intercept=0, slope=1), colour="firebrick", linetype=2)+
      theme_bw()+
      labs(x="Scam test p-value", y="CDF")
    ecdf_Peak

    ###power rate by curve type ####
    plot_power_3tests <- ggplot(data=filter(power, CurveID=="norm-var-sd4"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+ 
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.)+
      #facet_grid(func~as.factor(param))+
      theme_bw()+
      scale_color_manual(values=c('red', "blue", "darkgreen"))+
      labs(x="error", y="Power")
    plot_power_3tests
  
    #may be limited by the splines available (may work better if adaptive splines were present)
    #but its struggling with the polynomial, where the derivative remains constnat???
    #struggles a lot with saturating, which is in line with issues with adaptive spline
    plot_power_var_sd <- ggplot(data=filter(power, CurveID=="norm-var-sd2"|CurveID=="norm-var-sd8"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+   
      # annotate(geom="line", y=filter(power, CurveID=="norm-var-sd4")$scam_Pr,x=error,colour="blue",alpha=0.4)+
      # annotate(geom="line",y=filter(power, CurveID=="norm-var-sd4")$quad_Pr, x=error,colour="red",alpha=0.4)+
      # annotate(geom="line",y=filter(power, CurveID=="norm-var-sd4")$rprop_Pr,x=error, colour="darkgreen",alpha=0.4)+   
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.,
                 labeller=as_labeller(c("norm-var-sd2"="P3: Narrower Peak","norm-var-sd8"="P1: Wider peak")))+
      theme_bw()+
      theme(legend.position="none")+
      scale_color_manual(values=c('red', "blue", "darkgreen"))+
      labs(x="error", y="Power", title="b) Varying peak width")
    plot_power_var_sd
    
    plot_power_var_mean <- ggplot(data=filter(power, CurveID=="norm-var-mean6"|CurveID=="norm-var-mean8"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+    
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.,
                 labeller=as_labeller(c("norm-var-mean6"="P4: Peak shifted to 0.6","norm-var-mean8"="P5: Peak shifted to 0.8")))+
      theme_bw()+
      theme(legend.position="none")+
      scale_color_manual(values=c('red', "blue", "darkgreen"))+
      labs(x="error", y="Power", title="Varying peak location")
    plot_power_var_mean
    
    plot_power_var_skew <- ggplot(data=filter(power, CurveID=="beta6"|CurveID=="beta8"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+    
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.,
                 labeller=as_labeller(c("beta6"="P6: Skewed curve","beta8"="P7: More Skewed Curve")))+
      theme_bw()+
      theme(legend.position="none")+
      scale_color_manual(values=c('red', "blue", "darkgreen"))+
      labs(x="error", y="Power", title="Varying Skew")
    plot_power_var_skew
    
    plot_power_var_sd+plot_power_var_mean+plot_power_var_skew
    
    plot_power_two_lines <- ggplot(data=power, aes(x=error))+
      geom_line(aes(y=rprop_Pr, colour=CurveID))+ 
      #annotate(geom="line", x=error, y=filter(power, CurveID=="norm-var-sd4")$rprop_Pr)+
      facet_grid(factor(func)~.)+
      theme_bw()+
      theme(legend.position="none")+
      labs(x="error", y="Power", title="a) Two-Lines Test")
    plot_power_two_lines
    
    plot_power_quad <- ggplot(data=power, aes(x=error))+
      geom_line(aes(y=quad_Pr, colour=CurveID))+ 
      #annotate(geom="line", x=error, y=filter(power, CurveID=="norm-var-sd4")$quad_Pr)+
      facet_grid(factor(func)~.)+
      theme_bw()+
      theme(legend.position="bottom")+
      labs(x="error", y="Power", title="b) Quadatic Test")
    plot_power_quad
    
    plot_power_scam <- ggplot(data=power, aes(x=error))+
      geom_line(aes(y=scam_Pr, colour=CurveID))+ 
      #annotate(geom="line", x=error, y=filter(power, CurveID=="norm-var-sd4")$scam_Pr)+
      facet_grid(factor(func)~.)+
      theme_bw()+
      theme(legend.position="none")+
      labs(x="error", y="Power", title="c) Scam Test")
    plot_power_scam
    plot_power_two_lines+plot_power_scam
    
    #Reference plot of curves
    true_curves_Peak <- ggplot(data=filter(sim_data, rep==1, error==0.5, func=="beta" | func=="norm-var-mean"|func=="norm-var-sd"), aes(x=x, y=y_mean))+
      geom_line(aes(group=paste(param)))+
      labs(y="True curve")+
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.)
    true_curves_Peak
    
    true_curves_Peak+ecdf_Peak
  
    
  ## 2c: example cases ####
    ### True Positives ####
      #Case where quad did not detect peak
      scam61 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==61))
      summary(scam61)
      plot(y~x, data=filter(sim_data, sampleID==61))
      plot(scam61, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam61)[1])
      plot(scam61, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam61)[1])
      scam.check(scam61)
      appraise(scam61) 
  
      #Case where only scam detected peak
      scam2145 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==2145))
      summary(scam2145)
      plot(y~x, data=filter(sim_data, sampleID==2145))
      plot(scam2145, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam2145)[1])
      plot(scam2145, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam2145)[1])
      scam.check(scam2145)
      appraise(scam2145) 
      
      #lowest p-value case
      scam13049 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==13049))
      summary(scam13049)
      plot(y~x, data=filter(sim_data, sampleID==13049))
      plot(scam13049, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam13049)[1])
      plot(scam13049, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam13049)[1])
      scam.check(scam13049)
      appraise(scam13049) 
      
      #lowest p-value where another test did not identify peak
      scam2557 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==2557))
      summary(scam2557)
      plot(y~x, data=filter(sim_data, sampleID==2557))
      plot(scam2557, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam2557)[1])
      plot(scam2557, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam2557)[1])
      scam.check(scam2557)
      appraise(scam2557) 

    ### False Negatives ####
      #All methods but Scam found peak
        #Just above cutoff, looks well-fit tho
      scam13691 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==13691))
      summary(scam13691)
      plot(y~x, data=filter(sim_data, sampleID==13691))
      plot(scam13691, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam13691)[1])
      plot(scam13691, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam13691)[1])
      scam.check(scam13691)
      appraise(scam13691) 
      
      #no Method found Peak, p=1
      scam3306 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==3306))
      summary(scam3306)
      plot(y~x, data=filter(sim_data, sampleID==3306))
      plot(scam3306, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam3306)[1])
      plot(scam3306, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam3306)[1])
      scam.check(scam3306)
      appraise(scam3306) 
      
      #only 2-lines method found peak
      scam2476 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==2476))
      summary(scam2476)
      plot(y~x, data=filter(sim_data, sampleID==2476))
      plot(scam2476, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam2476)[1])
      plot(scam2476, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam2476)[1])
      scam.check(scam2476)
      appraise(scam2476) 
      
      #All methods but Scam found peak
      scam14598 = scam(y~s(x,k=10, bs="mpi", m=2)+s(x, k=10, bs="ps",  m =c(2,1)), data=filter(sim_data, sampleID==14598))
      summary(scam14598)
      plot(y~x, data=filter(sim_data, sampleID==14598))
      plot(scam14598, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           #seWithMean=TRUE, # WTF??
           shift = coef(scam14598)[1])
      plot(scam14598, pages=1, residuals=TRUE,
           pch = 1, cex = 1,
           shade=TRUE,
           seWithMean=TRUE, # WTF??
           shift = coef(scam14598)[1])
      scam.check(scam14598)
      appraise(scam14598) 
    
# 3 Gam Peak Quantification analysis ####
    #Importing peak quantification results 
      quant_CI_summary <- read.csv("data/simulation/quant_CI_summary_results.csv")
      quant_n_CI <- read.csv("data/simulation/quant_n_CI_results.csv")
      
#### peaks_results File too large to upload to GIT, please skip to quant_Sample_results ####
    peaks_results <- read.csv("data/simulation/quant_peaks_results.csv")
      
    #Importing raw simulation data in case not done before
    sim_data <- read.csv("data/simulation/CLEAN_simulationData.csv")
    
    ##3a Computing true peak location for each scenario####
      # true Reference peaks
        true_peaks <- sim_data%>%
          group_by(scenario, CurveID,func, param, error, sampleID)%>%
          summarise(x_true_peak=NA)%>%
          ungroup()
        
        true_peaks$x_true_peak[true_peaks$CurveID=="beta6"] <- 0.625
        true_peaks$x_true_peak[true_peaks$CurveID=="beta8"] <- 0.875
        true_peaks$x_true_peak[true_peaks$CurveID=="linear2"] <-NA
        true_peaks$x_true_peak[true_peaks$CurveID=="norm-var-mean6"] <-0.6
        true_peaks$x_true_peak[true_peaks$CurveID=="norm-var-mean8"] <-0.8
        true_peaks$x_true_peak[true_peaks$func=="norm-var-sd"] <-0.5
        true_peaks$x_true_peak[true_peaks$func=="polynomial"] <-NA
        true_peaks$x_true_peak[true_peaks$func=="saturating"] <-NA
        
        #Adding and relating the column names
        true_peaks$colnames <- paste0("X",true_peaks$sampleID)
        
        #Adding booleand for true-peaked or not peaked
        true_peaks$true_peaked <- 1
        true_peaks$true_peaked[which(is.na(true_peaks$x_true_peak))] <- 0
    
    ##3b Summarizing peak quantification results by sample ####
      #Identifying whether each critical maxima is identified as within a CI
        crit_withinCI_test <- function(x){
          return(c(x[626], x[875], x[601], x[801], x[501]))
        }
      true_peaks_within_CI <- apply(peaks_results, MARGIN=2,FUN=crit_withinCI_test )
      true_peaks_within_CI <- as.data.frame(t(true_peaks_within_CI[,-1]))
      colnames(true_peaks_within_CI) <- c("C_0.625", "C_0.875","C_0.600", "C_0.800", "C_0.500")
      true_peaks_within_CI$colnames <- rownames(true_peaks_within_CI)

      true_peaks_within_CI <- merge(true_peaks, true_peaks_within_CI, by="colnames")

      #Function to identify which of the critical maxima is the appropriate one
      peak_withinCI_test <- function(x){
        if(is.na(x[8]))
          return(NA)
        else if(x[8]=="0.625")
          return(as.numeric(x[10]))
        else if(x[8]=="0.875")
          return(as.numeric(x[11]))
       else if(x[8]=="0.600")
          return(as.numeric(x[12]))
        else if(x[8]=="0.800")
          return(as.numeric(x[13]))
        else if(x[8]=="0.500")
          return(as.numeric(x[14]))
      }
      
      #Determining on a per-sample basis if the CI includes the true peak, from the list of potential maximia
      true_peaks_within_CI$C_incl_peak <- apply(true_peaks_within_CI, MARGIN=1,FUN=peak_withinCI_test )

      #adding per-sample number of CI
      quant_Sample_results <- merge(true_peaks_within_CI, quant_n_CI, by="sampleID")
      write.csv(quant_Sample_results, "data/simulation/quant_Sample_results.csv",row.names = FALSE)
      

      #adding per-sample number fo peak sand troughs
      quant_n_PorT <- quant_CI_summary%>%
        group_by(sampleID)%>%
        summarise(n_peak=sum(peak_or_trough==1),
                  n_trough=sum(peak_or_trough==-1))%>%
        ungroup()

      quant_Sample_results <- merge(quant_Sample_results, quant_n_PorT, by="sampleID")
      
      #Calculating total CI width for each sample to add to quant_Sample_results
      quant_CI_summary_peaks <- filter(quant_CI_summary, peak_or_trough==1)
      
      total_CI_width_peaks <- quant_CI_summary_peaks%>%
        group_by(sampleID)%>%
        summarize(total_width_peak=sum(CI_width),
                  n_peaks=n())%>%
        ungroup()
      
      quant_Sample_results <- merge(quant_Sample_results, total_CI_width_peaks[,c(1, 2)], by="sampleID", all=TRUE)
      
      total_CI_width <- quant_CI_summary%>%
        filter(peak_or_trough==1| peak_or_trough==-1)%>%
        group_by(sampleID)%>%
        summarize(total_width_CI=sum(CI_width),
                  n_CI=n())%>%
        ungroup()
      
      quant_Sample_results <- merge(quant_Sample_results, total_CI_width[,c(1, 2)], by="sampleID", all=TRUE)

      #### Skip to here to avoid code blocks using peaks_results####
      quant_Sample_results <-  read.csv("data/simulation/quant_Sample_results.csv")
      
      
    #Summary of number of CIs
      quant_Sample_results%>%
        group_by(n_CI)%>%
        summarize(tally_CI=n())%>%
        ungroup()
      
      quant_Sample_results%>%
        group_by(n_peak)%>%
        summarize(tally_Peak=n())%>%
        ungroup()
      
      quant_Sample_results%>%
        group_by(n_trough)%>%
        summarize(tally_trough=n())%>%
        ungroup()
      
      quant_Sample_results%>%
        filter(n_CI==1)%>%
        summarize(n_peak=sum(n_peak),
                  n_trough=sum(n_trough))


      n_CI_by_scenario <- quant_Sample_results%>%
        group_by(n_CI, scenario)%>%
        summarize(tally_CI=n())%>%
        ungroup()
      n_CI_by_scenario
      #Computing a scenario-specific 'accuracy rate' (defined as the rate of the CI including the true peak)
        #only performed for cases with a true peak

##3c: Summarizing results by scenario ####
        #Power rate , accuracy rate, trough rate and peak rate 
        #(ie, # peaked curves where a CI was detected)
        # curves where at least 1 CI includes mean, 
        # # curves where troughs were detected
        # # curves where peaks were detected
      quant_scenario_results <- quant_Sample_results%>%
        filter(!(is.na(x_true_peak)))%>% 
        group_by(scenario, error, param, func, x_true_peak, CurveID)%>%
        summarize(accRate=sum(C_incl_peak)/length(C_incl_peak),
                  powerRate=sum(n_CI!=0)/n(),
                  peakRate=sum(n_peak!=0)/n(),
                  troughRate=sum(n_trough!=0)/n(),
                  coverage=sum(C_incl_peak)/sum(n_peak!=0),
                  multiCIErrorRate=sum(n_CI>1)/n(),
                  mean_peakCI_width=mean(total_width_peak, na.rm=TRUE),
                  sd_peakCI_width=sd(total_width_peak, na.rm=TRUE),
                  n_sampleWpeak=sum(!is.na(total_width_peak)),
                  mean_CI_width=mean(total_width_CI, na.rm=TRUE),
                  sd_CI_width=sd(total_width_CI, na.rm=TRUE),
                  n_sampleWCI=sum(!is.na(total_width_CI)),
                  n_samples=n())%>%
        ungroup()
      quant_scenario_results
      View(quant_scenario_results)
      
      #Plotting the accuracy rate
      accuracyRatePlot <- ggplot(data=quant_scenario_results, aes(x=error))+
        geom_line(aes(y=accRate, colour=CurveID, lty="Accuracy Rate"))+
        geom_line(aes(y=powerRate, colour=CurveID, lty="Power Rate"))+
        geom_line(aes(y=coverage, colour=CurveID, lty="coverage"))+
        #geom_line(aes(y=troughRate, colour=CurveID, lty="Trough Rate"))+
        annotate(geom="line", x=filter(quant_scenario_results, CurveID=="norm-var-sd4")$error, 
                 y=filter(quant_scenario_results, CurveID=="norm-var-sd4")$accRate, colour="grey", lty=1)+
       # geom_point(aes(y=accRate))+
        facet_grid(factor(func)~.)+
        scale_linetype_manual(values=c(1, 2, 3),
                              name="Rate", labels=c("Accuracy Rate", "coverage", "Power Rate"))+
        theme_bw()+
        labs(x="Residual Error", y="Accuracy")
      accuracyRatePlot
      
      accuracyRate_refplot <- ggplot(data=filter(quant_scenario_results, CurveID=="norm-var-sd4"), aes(x=error))+
        geom_line(aes(y=accRate, colour=CurveID))+
        geom_line(aes(y=powerRate, colour=CurveID), lty=2)+
        # geom_point(aes(y=accRate))+
        facet_grid(factor(func)~., labeller=as_labeller(c("norm-var-sd"="P3: Reference Plot")))+
        theme_bw()+
        labs(x="Residual Error", y="Accuracy")
      accuracyRate_refplot
      
      accuracyRate_beta <- ggplot(data=filter(quant_scenario_results, func=="beta"), aes(x=error))+
        geom_line(aes(y=accRate, colour=CurveID))+
        geom_line(aes(y=powerRate, colour=CurveID), lty=2)+
        # geom_point(aes(y=accRate))+
        facet_grid(factor(func)~.)+
        theme_bw()+
        labs(x="Residual Error", y="Accuracy")
      accuracyRate_beta
      
      accuracyRate_varsd <- ggplot(data=filter(quant_scenario_results, func=="beta"), aes(x=error))+
        geom_line(aes(y=accRate, colour=CurveID))+
        geom_line(aes(y=powerRate, colour=CurveID), lty=2)+
        # geom_point(aes(y=accRate))+
        facet_grid(factor(func)~.)+
        theme_bw()+
        labs(x="Residual Error", y="Accuracy")
      accuracyRate_varsd
      
      #plotting coverage 
      quant_scenario_results

      quant_scenario_results$CurveID2 <- NA
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="saturating8")] <- "D"    
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="saturating4")] <- "D"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="saturating2")] <- "D"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="linear2")] <- "C"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="polynomial2")] <- "A"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="polynomial4")] <- "A"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="polynomial8")] <- "A"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="norm-var-sd8")] <- "W"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="norm-var-sd4")] <- ""
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="norm-var-sd2")] <- "W"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="norm-var-mean6")] <- "L"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="norm-var-mean8")] <- "L"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="beta6")] <- "S"
      quant_scenario_results$CurveID2[which(quant_scenario_results$CurveID=="beta8")] <- "S"
      
      quant_scenario_results$CurveID3 <- NA
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="saturating8")] <- "3"    
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="saturating4")] <- "2"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="saturating2")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="linear2")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="polynomial2")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="polynomial4")] <- "2"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="polynomial8")] <- "3"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="norm-var-sd8")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="norm-var-sd4")] <- "Ref"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="norm-var-sd2")] <- "2"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="norm-var-mean6")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="norm-var-mean8")] <- "2"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="beta6")] <- "1"
      quant_scenario_results$CurveID3[which(quant_scenario_results$CurveID=="beta8")] <- "2"
      
      coverage_plot <- ggplot(data=filter(quant_scenario_results,
                                          CurveID=="beta8"|CurveID=="norm-var-mean8"|func=="norm-var-sd"), aes(x=error))+
        geom_line(aes(y=coverage, colour=CurveID3))+
        geom_hline(yintercept = 0.95, lty=2, linewidth=0.2)+
        annotate(geom="line", x=filter(quant_scenario_results, CurveID=="norm-var-sd4")$error, 
                 y=filter(quant_scenario_results, CurveID=="norm-var-sd4")$coverage, colour="black")+
        facet_grid(factor(func)~.,
                   labeller=as_labeller(c("beta"="S: Varying Skew", "norm-var-mean"="L: Varying Peak Location", "norm-var-sd"="W: Varying Peak Width")))+
        theme_bw()+
        theme(axis.line = element_line(color='black'),
              plot.background = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position="bottom")+
        labs(x="Residual Error", y="Capture Rate", colour="")+
        scale_color_manual(values=c('#00BA38', "#619CFF", "black"),
                           labels=c("W1", "W2, L2, S2", "Reference"))
      coverage_plot
     
    ## 3d: redoing Type I error and power with gam derivative test ####
    ### Fale Positives ####
    falsePosRate_quant <- quant_Sample_results%>%
      filter(is.na(x_true_peak))%>%
      group_by(scenario, CurveID, error, func, param)%>%
      summarize(falsePosRate=sum(n_CI!=0)/n(),
                troughErrorRate=sum(n_trough!=0)/n(),
                multiCIErrorRate=sum(n_CI!=0)/n())%>%
      ungroup()

      typeIerrorRate <- merge(typeIerrorRate,falsePosRate_quant[, c(1, 6)], by="scenario")
      
    plot_T1error_quant <- ggplot(data=falsePosRate_quant, aes(x=error))+
      geom_line(aes(y=falsePosRate))+
      geom_hline(yintercept=0.05, lty=2)+
      facet_grid(factor(CurveID, levels=c("saturating8", "saturating4","saturating2","linear2","polynomial2","polynomial4","polynomial8"))~.)+
      theme_bw()+
      labs(x="error", y="Type I Error Rate")
    #theme(legend.position = "bottom")
    plot_T1error_quant
    
    plot_T1error_publish <- ggplot(data=filter(typeIerrorRate,
                                               CurveID=="saturating8"|CurveID=="linear2"|CurveID=="polynomial8"), aes(x=error))+
      geom_line(aes(y=scam_T1R, colour="Scam test"))+
      geom_line(aes(y=quad_T1R, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_T1R, colour="Two-Lines test"))+
      geom_line(aes(y=falsePosRate, x=error, colour="GAM Derivative \nTest"), data=filter(falsePosRate_quant,
                                                                                          CurveID=="saturating8"|CurveID=="linear2"|CurveID=="polynomial8"))+
      # geom_line(aes(y=scam_T1R, linetype="Scam test"), linewidth=1)+
      # geom_line(aes(y=quad_T1R,linetype="Quadratic Test"))+
      # geom_line(aes(y=rprop_T1R, linetype="Two-Lines test"))+  
      # geom_line(aes(y=falsePosRate, x=error,linetype="GAM"),linewidth=1, data=filter(falsePosRate_quant,
      #                                                                                               CurveID=="saturating8"|CurveID=="linear2"|CurveID=="polynomial8"))+
      geom_hline(yintercept=0.05, lty=5, colour="black", linewidth=0.75)+
      facet_grid(factor(CurveID, levels=c("saturating8", "linear2","polynomial8"))~.,
                 labeller=as_labeller(c("saturating8"="D3: Decelerating \nRate of Change","linear2"="C1: Constant \nRate of Change","polynomial8"="A3: Accelerating \nRate of Change")))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor=element_line(colour="grey87"),
            legend.position="bottom")+
      scale_color_manual(values=c('purple', "red", "blue", "forestgreen"))+
      #scale_linetype_manual(values=c(2, 3, 4, 5))+
      labs(x="Residual error", y="Rate of False-Positive Errors", colour="")+
      scale_y_continuous(trans='log10', breaks=c(0.005,  0.05, 0.5),
                         minor_breaks = c(0.005, 0.01,0.02,0.03,0.04,0.05, 0.1,0.2,0.3,0.4,0.5))+
      guides(colour=guide_legend(nrow=2,byrow=TRUE))
    plot_T1error_publish
    
    plot_T1error_publish_pres <- ggplot(data=filter(typeIerrorRate,
                                               CurveID=="saturating8"), aes(x=error))+
      geom_line(aes(y=scam_T1R, colour="Scam test"))+
      geom_line(aes(y=quad_T1R, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_T1R, colour="Two-Lines test"))+
      geom_line(aes(y=falsePosRate, x=error, colour="GAM Derivative \nTest"), data=filter(falsePosRate_quant,
                                                                                          CurveID=="saturating8"))+
      geom_hline(yintercept=0.05, lty=5, colour="black", linewidth=0.75)+
      facet_grid(factor(CurveID, levels=c("saturating8", "linear2","polynomial8"))~.,
                 labeller=as_labeller(c("saturating8"="Decelerating Rate of Change","linear2"="Constant Rate of Change","polynomial8"="Accelerating Rate of Change")))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor=element_line(colour="grey87"),
            panel.grid.major=element_line(colour="grey80"),
            legend.position="bottom")+
      scale_color_manual(values=c('purple', "red", "blue", "forestgreen"))+
      #scale_linetype_manual(values=c(2, 3, 4, 5))+
      labs(x="Residual error", y="Rate of False-Positive Errors", colour="")+
      scale_y_continuous(trans='log10', breaks=c(seq(0.01, 0.05, by=0.01), seq(0.1, 0.3, by=0.1)),
                         minor_breaks = c(seq(0.01, 0.5, by=0.01)))+
      guides(colour=guide_legend(nrow=2,byrow=TRUE))
    plot_T1error_publish_pres
    
  ###Power rate ####
    
    power_master <- merge(power, quant_scenario_results[, c(1, 7, 8, 9, 10,11, 19, 20)])
    
    power_master_plot <- ggplot(data=filter(power_master,
                                        CurveID=="beta8"|CurveID=="norm-var-mean8"|CurveID=="norm-var-sd2"|CurveID=="norm-var-sd8"|CurveID=="norm-var-sd4"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+
      geom_line(aes(y=peakRate, colour="Gam Derivative test"))+
      facet_grid(factor(CurveID, levels=c("norm-var-sd4", "beta8", "norm-var-mean8", "norm-var-sd8", "norm-var-sd2"))~.,
                 labeller=as_labeller(c("norm-var-sd4"="Ref: Reference Curve", "beta8"="S2: Skew", "norm-var-mean8"="L2: Peak Location", "norm-var-sd8"="W1: Wide Peak", "norm-var-sd2"="W2: Narrow Peak")))+
      ylim(c(0,1))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="bottom")+
      scale_linetype_manual(values=c(1, 2, 3, 4))+
      labs(x="Residual Error", y="Power", colour="", lty="")+
     scale_color_manual(values=c('purple', "red", "blue", "darkgreen"))+
      guides(colour=guide_legend(nrow=2,byrow=TRUE),
             lty=guide_legend(nrow=2,byrow=TRUE))
    power_master_plot
    
    power_master_plot_pres <- ggplot(data=filter(power_master,
                                            CurveID=="beta8"|CurveID=="norm-var-sd4"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+
      geom_line(aes(y=peakRate, colour="Gam Derivative test"))+
      facet_grid(factor(CurveID, levels=c("norm-var-sd4", "beta8", "norm-var-mean8", "norm-var-sd8", "norm-var-sd2"))~.,
                 labeller=as_labeller(c("norm-var-sd4"="Ref: Reference Curve", "beta8"="S2: Skew", "norm-var-mean8"="L2: Peak Location", "norm-var-sd8"="W1: Wide Peak", "norm-var-sd2"="W2: Narrow Peak")))+
      ylim(c(0,1))+
      theme_bw()+
      theme(axis.line = element_line(color='black'),
            plot.background = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position="bottom")+
      scale_linetype_manual(values=c(1, 2, 3, 4))+
      labs(x="Residual Error", y="Power", colour="", lty="")+
      scale_color_manual(values=c('purple', "red", "blue", "darkgreen"))+
      guides(colour=guide_legend(nrow=2,byrow=TRUE),
             lty=guide_legend(nrow=2,byrow=TRUE))
    power_master_plot_pres

    
    plot_power <- ggplot(data=filter(power_master, CurveID=="norm-var-sd4"), aes(x=error))+
      geom_line(aes(y=scam_Pr, colour="Scam test"))+
      geom_line(aes(y=quad_Pr, colour="Quadratic Test"))+
      geom_line(aes(y=rprop_Pr, colour="Two-Lines test"))+ 
      geom_line(aes(y=peakRate, x=error, colour="GAM"))+
      facet_grid(factor(CurveID, levels=unique(detect_results_noNA_curvetype_Peak$CurveID))~.)+
      #facet_grid(func~as.factor(param))+
      theme_bw()+
      scale_color_manual(values=c('red', "blue", "darkgreen", "orange"))+
      labs(x="error", y="Power")
    plot_power
    
   ###3e example curves ####
    source("code/parametrized_peak_finder.r")
    dev.off()
    
    test_data <- sim_data%>%
      filter(sampleID==1368)
    
    #reformatting the dataset for it to work with parametrized peakfinder code
    data <- tibble(x = test_data$x,y= test_data$y,
                   true_val = test_data$y_mean)
    
    #running gam
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
    
    model_summary <- create_peakfinder_input(mod = mod, data = data)
    
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 1000)
    
    set.seed(5)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=FALSE)
    points(x=0.625, y=1, pch=16, col="blue", cex=1.2)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    
    
    #nonCaptured peak
    test_data <- sim_data%>%
      filter(sampleID==21)
    
    #reformatting the dataset for it to work with parametrized peakfinder code
    data <- tibble(x = test_data$x,y= test_data$y,
                   true_val = test_data$y_mean)
    
    #running gam
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
    
    model_summary <- create_peakfinder_input(mod = mod, data = data)
    
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 1000)
    
    set.seed(5)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE ,
                                 poisson=FALSE)
    points(x=0.625, y=1, pch=16, col="blue", cex=1.2)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    
    
    #Too many peaks
    test_data <- sim_data%>%
      filter(sampleID==27295)
    
    #reformatting the dataset for it to work with parametrized peakfinder code
    data <- tibble(x = test_data$x,y= test_data$y,
                   true_val = test_data$y_mean)
    
    #running gam
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
    
    model_summary <- create_peakfinder_input(mod = mod, data = data)
    
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 1000)
    
    set.seed(27295)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    
    
    
    #2 CIs
    test_data <- sim_data%>%
      filter(sampleID==1936)
    
    #reformatting the dataset for it to work with parametrized peakfinder code
    data <- tibble(x = test_data$x,y= test_data$y,
                   true_val = test_data$y_mean)
    
    #running gam
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
    
    model_summary <- create_peakfinder_input(mod = mod, data = data)
    
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 1000)
    
    set.seed(1936)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                 plot_2nd_deriv = FALSE)
    plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    
    #Investigating low coverage of normvarmean 8 error level 0.3 
    IDnormvarsd80.3 <- filter(quant_Sample_results, scenario=="norm-var-mean80.3", n_CI>0, n_peak>0, C_incl_peak==0)$sampleID
    for (i in IDnormvarsd80.3){
      test_data <- sim_data%>%
        filter(sampleID==i)
      
      #reformatting the dataset for it to work with parametrized peakfinder code
      data <- tibble(x = test_data$x,y= test_data$y,
                     true_val = test_data$y_mean)
      
      #running gam
      mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data)
      
      model_summary <- create_peakfinder_input(mod = mod, data = data)
      
      model_coefs <- model_data(true_model_summary = model_summary, n_sims = 1000)
      
      set.seed(i)
      possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
      plot_peaks_against_true_vals(possible_peaks = possible_peaks, model_coefs =model_coefs,true_model_summary = model_summary,
                                   plot_2nd_deriv = FALSE)
      points(x=0.8, y=1, pch=16, col="blue", cex=1.2)
      
      plot_derivs(model_coefs =model_coefs,possible_peaks = possible_peaks,plot_2nd_deriv = FALSE)
    }
    
    #bad curve (most CIs identified)
    set.seed(4024)
    test_data <- sim_data%>%filter(sampleID==4024)
    data4024 <- tibble(x = test_data$x,y= test_data$y,
                               true_val = test_data$y_mean)
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data4024)
    model_summary <- create_peakfinder_input(mod = mod, data = data)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks=possible_peaks, model_coefs = model_coefs,true_model_summary = model_summary,plot_2nd_deriv = FALSE)
    plot_derivs(possible_peaks=possible_peaks, model_coefs = model_coefs, plot_2nd_deriv = FALSE)
    
    
    set.seed(15023)
    test_data <- sim_data%>%filter(sampleID==15023)
    data15023 <- tibble(x = test_data$x,y= test_data$y,
                       true_val = NA)
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data15023)
    model_summary <- create_peakfinder_input(mod = mod, data = data15023)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks=possible_peaks, model_coefs = model_coefs,true_model_summary = model_summary,plot_2nd_deriv = FALSE)
    plot_derivs(possible_peaks=possible_peaks, model_coefs = model_coefs, plot_2nd_deriv = FALSE)
      
    
    
    set.seed(7014)
    test_data <- sim_data%>%filter(sampleID==7014)
    data7014 <- tibble(x = test_data$x,y= test_data$y,
                        true_val = test_data$y_mean)
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data7014)
    model_summary <- create_peakfinder_input(mod = mod, data = data7014)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks=possible_peaks, model_coefs = model_coefs,true_model_summary = model_summary,plot_2nd_deriv = FALSE)
    plot_derivs(possible_peaks=possible_peaks, model_coefs = model_coefs, plot_2nd_deriv = FALSE)
    
    
    set.seed(7102)
    test_data <- sim_data%>%filter(sampleID==7102)
    data7102 <- tibble(x = test_data$x,y= test_data$y,
                       true_val = test_data$y_mean)
    mod <- gam(y~s(x, k=20, bs="ad", m=3),method="REML", data=data7102)
    model_summary <- create_peakfinder_input(mod = mod, data = data7102)
    model_coefs <- model_data(true_model_summary = model_summary, n_sims = 5000)
    possible_peaks= peak_finder(model_coefs = model_coefs, true_model_summary = model_summary, test = "crossing")
    plot_peaks_against_true_vals(possible_peaks=possible_peaks, model_coefs = model_coefs,true_model_summary = model_summary,plot_2nd_deriv = FALSE)
    plot_derivs(possible_peaks=possible_peaks, model_coefs = model_coefs, plot_2nd_deriv = FALSE)
    

################################################################################################
  # END OF CODE ####
################################################################################################