# data_cleaning.R
# Author: Natalie Dupont
# Date Created: 2022-09-30
# Date last modified: 2023-04-15 
# BIOL 490, Pedersen Lab, Concordia University 
# Description: This code serves to clean literature data from webplot digitizer: 
# This mainly involves ensuring the plots generated in R closely match the original
# figure, and rounding to the nearest integer for data known to be integer 
# (species richness, years since disturbance, etc.). Data are then rewritten as CSV for 
# later use in other scripts (real_data_analysis.R)

########################################################################################
#Outline

#0: Loading packages 
# 1: Collins 1995
# 2: Halpern and Spies 1995
# 3: Hiura 1995
# 4: Marra 2014
# 5: Lazarina 2019
# 6: Wilsson and Keddy 1988
# 7: Keddy 1983


########################################################################################
#0: Loading packages ####
library(mgcv)

# 1: Collins 1995 ####
  Collins1995 <- read.csv("data/Collins1995_Fig3b.csv")
  head(Collins1995)
  plot(n_sp~y_pb, data=Collins1995)
    #both variables are integer data, so should be rounded to the nearest digit
  #sorting in x-ascending
  Collins1995 <- Collins1995[order(Collins1995$y_pb), ]
  
  Collins1995$y_pb <- round(Collins1995$y_pb, 0)#integer value
  Collins1995$n_sp <- round(Collins1995$n_sp, 0)#integer value
  
  write.csv(Collins1995,"data/processed/CLEAN_Collins1995.csv",row.names = FALSE)

#2: Halpern and Spies 1995####
  HalpernSpies1995_w1 <- read.csv("data/HalpernSpies1995_Fig1d_w1.csv")
    plot(heterogeneity~t_disturbance, data=HalpernSpies1995_w1)
    
    #sorting in x-ascending
    HalpernSpies1995_w1 <- HalpernSpies1995_w1[order(HalpernSpies1995_w1$t_disturbance), ]
    
  HalpernSpies1995_w1$t_disturbance <- round(HalpernSpies1995_w1$t_disturbance, 0) #integer value
    #heterogeneity remains unchanged since it is continuous
  head(HalpernSpies1995_w1)
  
  write.csv(HalpernSpies1995_w1,"data/processed/CLEAN_HalpernSpies1995_Fig1d_w1.csv",row.names = FALSE)
  
    
  HalpernSpies1995_w2 <- read.csv("data/HalpernSpies1995_Fig1d_w2.csv")
    plot(heterogeneity~t_disturbance, data=HalpernSpies1995_w2)
    
    #sorting in x-ascending
    HalpernSpies1995_w2 <- HalpernSpies1995_w2[order(HalpernSpies1995_w2$t_disturbance), ]
  HalpernSpies1995_w2$t_disturbance <- round(HalpernSpies1995_w2$t_disturbance, 0) #integer value
  head(HalpernSpies1995_w2)
  
  write.csv(HalpernSpies1995_w2,"data/processed/CLEAN_HalpernSpies1995_Fig1d_w2.csv",row.names = FALSE)
  
#3: Hiura 1995 ####
  Hiura1995 <- read.csv("data/Hiura1995_Fig1.csv")
  head(Hiura1995)
  View(Hiura1995)
  plot(H_prime~ln_MWI, data=Hiura1995, col=StandID) 
    #both variables are continuous, may lose too much precision by rounding in introduce more error
  # Hiura1995$H_prime <- round(Hiura1995$H_prime, 1)
  # Hiura1995$ln_MWI <- round(Hiura1995$ln_MWI, 1)
  
  #sorting in x-ascending
  Hiura1995 <- Hiura1995[order(Hiura1995$ln_MWI), ]
  
  write.csv(Hiura1995,"data/processed/CLEAN_Hiura1995.csv",row.names = FALSE)
  
#4: Marra 2014 ####
  Marra2014 <- read.csv("data/Marra2014_Fig4c.csv")
  #sorting in x-ascending
  Marra2014 <- Marra2014[order(Marra2014$mort_mean), ]
  
  head(Marra2014)
  View(Marra2014)
  plot(Marra2014$sp_rarefaction~Marra2014$mort_mean)
  
  write.csv(Marra2014,"data/processed/CLEAN_Marra2014.csv",row.names = FALSE)
  
   
#5: Lazarina 2019 ####
  Lazarina2019_fig2a <- read.csv("data/Lazarina2019_Fig2a.csv")
  head(Lazarina2019_fig2a)
  plot(richness~fire_severity, data=Lazarina2019_fig2a)
  #sorting in x-ascending
  Lazarina2019_fig2a <- Lazarina2019_fig2a[order(Lazarina2019_fig2a$fire_severity), ]
  
  Lazarina2019_fig2b <- read.csv("data/Lazarina2019_Fig2b.csv")
  head(Lazarina2019_fig2b)
  plot(expShannon~fire_severity, data=Lazarina2019_fig2b)
  
  #plotting to compare to article
  par(mfrow=c(1, 2))
  plot.gam(gam(richness~s(fire_severity, bs="tp", k=17),data=Lazarina2019_fig2a, method="REML"),
           residuals = TRUE,
           pch = 1, cex = 1)
  plot.gam(gam(expShannon~s(fire_severity, bs="tp", k=17),data=Lazarina2019_fig2b, method="REML"),
           residuals = TRUE,
           pch = 1, cex = 1)

#6 Wilsson and Keddy 1988 ####
  WilsonKeddy1988<- read.csv("data/WilsonKeddy1988_Fig2.csv")
  plot(WilsonKeddy1988$n_sp~WilsonKeddy1988$somc_perc)
  
  #sorting in x-ascending
  WilsonKeddy1988 <- WilsonKeddy1988[order(WilsonKeddy1988$somc_perc), ]
  
  #n_sp is integer data
  WilsonKeddy1988$n_sp <- round(WilsonKeddy1988$n_sp, 0)#integer value
  write.csv(WilsonKeddy1988,"data/processed/CLEAN_WilsonKeddy1988.csv",row.names = FALSE)
  
#7 Nilsson 1987 ####
  Nilsson1987<- read.csv("data/Nilsson1987_Fig2a.csv")
  plot(Nilsson1987$num_spp~Nilsson1987$vel_m_s)
  
  Nilsson1987$num_spp <- round(Nilsson1987$num_spp, 0)#integer value
  
  #sorting in x-ascending
  Nilsson1987 <- Nilsson1987[order(Nilsson1987$vel_m_s), ]
  
  write.csv(Nilsson1987,"data/processed/CLEAN_Nilsson1987_Fig2a.csv",row.names = FALSE)
  
#8: Keddy 1983####
  Keddy1983<- read.csv("data/Keddy1983_Fig4.csv")
  plot(Keddy1983$n_spp~Keddy1983$exp)
  
  Keddy1983$n_spp <- round(Keddy1983$n_spp, 0)#integer value
  Keddy1983$exp <- round(Keddy1983$exp, 0)#integer value
  
  #sorting in x-ascending
  Keddy1983 <- Keddy1983[order(Keddy1983$exp), ]
  
  write.csv(Keddy1983,"data/processed/CLEAN_Keddy1983_Fig4.csv",row.names = FALSE)
########################################################################################################
#End of Code
########################################################################################################