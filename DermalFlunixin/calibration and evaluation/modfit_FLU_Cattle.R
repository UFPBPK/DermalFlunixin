## Load libraries
library(mrgsolve)    # Needed for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    # The pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(ggplot2)     # Needed for plot
library(FME)         # Package for MCMC simulation and model fitting
library(minpack.lm)  # Package for model fitting
library(reshape)     # Package for melt function to reshape the table
library(truncnorm)   # Package for the truncated normal distribution function   
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    # Package for inverse gamma distribution function
library(foreach)     # Package for parallel computing
library(doParallel)  # Package for parallel computing
library(bayesplot)   # Package for MCMC traceplot
library(tidyr)       # R-package for tidy messy data
library(tidyverse)   # R-package for tidy messy data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(truncnorm)   # R package for Truncated normal distribution
library(EnvStats)    # Package for Environmental Statistics, Including US EPA Guidance
library(ggpubr)      # R package for plotting the data
library(ggprism)
library(ggplot2)
library(tibble)
library(rlang)
library(patchwork)

# Input the PBPK model
source (file = "C:/Users/wuxue/OneDrive - University of Florida/cattle mrg pbpk general.R") # Loading the generic PBPK model code

## Build mrgsolve-based PBPK Model
mod <- mcode ("MrgSolve", PBPK)

## Model fitting
data_flu <- read.csv("C:/Users/wuxue/OneDrive - University of Florida/Data_FLU.csv")
head(data_flu)

## Read the data set and later used in model calibration
## Study 1: 2016 Kleinhenz,M.D.et al
Obs_A1     <- data_flu %>% filter(Study ==1)
Obs_A1_P   <- Obs_A1 %>% filter(Matrix == "P")%>%select(Time = Time, CV = Conc.)

## Study 3: 2016 Kleinhenz,M.D.et al Effect of age
Obs_A3     <- data_flu %>% filter(Study ==3)
Obs_A3_P   <- Obs_A1 %>% filter(Matrix == "P")%>%select(Time = Time, CV = Conc.)

## Study 4: 2018 Kleinhenz,M.D.et al The impact of pain
Obs_A4     <- data_flu %>% filter(Study ==4)
Obs_A4_P   <- Obs_A4 %>% filter(Matrix == "P")%>%select(Time = Time, CV = Conc.)

## Study 4.1: 2018 Kleinhenz,M.D.et al The impact of pain
Obs_A4.1     <- data_flu %>% filter(Study ==4.1)
Obs_A4.1_P   <- Obs_A4.1 %>% filter(Matrix == "P")%>%select(Time = Time, CV = Conc.)

## Study 7: 2017 FDA
Obs_A7     <- data_flu %>% filter(Study ==7)
Obs_A7_L   <- Obs_A7 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc.)
Obs_A7_M   <- Obs_A7 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc.)
Obs_A7_K   <- Obs_A7 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc.)
Obs_A7_F   <- Obs_A7 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc.)

# Define the prediction function
## Define the prediction function (for least squres fit using levenberg-marquart algorithm)
pred <- function(pars, BW, area, Dtimes) {
  
  ## Get out of log domain
  pars <- pars ## 
  
  ## Define the three exposure scenario
  ## Exposure scenario 
  tinterval     = 24
  End_time      = 9
  Dtimes        = 1
  TDOSE         = Dtimes
  BW            = 120 # KG
  PDOSEde       = 3.3 # mg/kg/day
  MW            = 296.24
  DOSE          = PDOSEde*BW/MW
  area          = 1500
  
  
  
  ev_1.B <- ev (ID   = 1, amt  = DOSE, ii = tinterval,
                
                addl = TDOSE-1, cmt  = "AS", replicate = FALSE)
  
  
  
  
  ## set up the exposure time
  tsamp  = tgrid(0, tinterval*(TDOSE-1) + tinterval*End_time, 0.1)
  
  
  ## Simulation
  out <- 
    mod %>% 
    param(pars) %>%
    update(atol = 1E-4, rtol=1E-2) %>%
    mrgsim_d(data = ev_1.B, tgrid=tsamp)
  
  CHECK<-as.data.frame(out)
  
  output <-cbind.data.frame(Time = out$time,
                            CV = out$conc,
                            CL = out$Liver,
                            CK = out$Kidney,
                            CM = out$Muscle,
                            CF = out$Fat,
                            CV1= out$Venous1,
                            CL1= out$Liver1,
                            CK1= out$Kidney1,
                            AUC_V = out$AUCCV,
                            AUC_L = out$AUCCL,
                            AUC_K = out$AUCCK,
                            AUC_M = out$AUCCM,
                            AUC_F = out$AUCCF,
                            Jf1= out$jf1,
                            Jf2=out$jf2,
                            Jf3=out$jf3,
                            Jf4=out$jf4,
                            Bal=out$bal,
                            solventvol=out$solventvol2,
                            AS=out$as,
                            ASC=out$asc,
                            AVS=out$avs)
  
  return(output)
  
}





#-----------------------------------------------------
##balance check

pars <-  c(logP=1.09,logPsw =2.952,MTf=6,D=12)

OUT<- pred(pars)

plot(OUT$Time,OUT$CV,ylim=c(0,0.1),type="l")
plot(OUT$Time,OUT$Bal,ylim=c(-0.00000000001,0.000000001),type="l")

#################################################################################################################################
#################################################################################################################################


### Make cost function
Cost <- function (pars,BW,area,Dtimes,w) {
 
  # Prediction
  #out_A1      <- pred (pars=pars, BW = 60.2, area=1200, Dtimes = 1)
  out_A3      <- pred (pars=pars, BW = 222,  area=3600,Dtimes = 1)
  #out_A4      <- pred (pars=pars, BW = 61.3, area=1200,Dtimes = 1)
  #out_A4.1    <- pred (pars=pars, BW = 61.3, area=1200,Dtimes = 1)
  out_A7      <- pred (pars=pars, BW = 309 , area=4800,Dtimes = 1)
  
  cost <- NULL
  
  ## Cost function
  #cost<- modCost  (model = out_A1,   obs = Obs_A1_P,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A3,   obs = Obs_A3_P,   x ="Time", cost = cost, weight = w)
  #cost<- modCost  (model = out_A4,   obs = Obs_A4_P,   x ="Time", cost = cost, weight = w)
  #cost<- modCost  (model = out_A4.1, obs = Obs_A4.1_P, x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A7,   obs = Obs_A7_L,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A7,   obs = Obs_A7_M,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A7,   obs = Obs_A7_K,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A7,   obs = Obs_A7_F,   x ="Time", cost = cost, weight = w)
  
  return(cost)
}

Cost(pars, w="mean")

## Initial parameters
pars1 = c(logP=1.09,logPsw =2.952,MTf=6,
          D=12,Prest=10.5,PK=4,solventevap=0.024,
          PL=10.52,PF=0.6,PM=0.5)

## Sensitivity function (FME)
## Check the Sensitivity parameters in the model
Sns <- sensFun(func = Cost,w = "mean",
               parms = pars1, varscale = 1)

Sen_1 <- summary(Sns)
plot(summary(Sns))


# ## Selected senstive parameters;
newpars <- pars[abs(Sen_1$Mean) > mean(abs(Sen_1$Mean))]
newpars

## Customized parameter value
pars2 = c(logP=1.09,logPsw =2.952,PK=4,PF=0.6)

# ## Selected sensitive parameters; 
Fit <- modFit(f = Cost, 
              p = pars2,
              lower =pars2*0.2, upper = pars2*2,
              method ="Marq",w = "mean", control = nls.lm.control(nprint = 1))

summary(Fit) 
Fit$par

preout<- pred(Fit$par)
## Make the data for the plot
Cost_new <-Cost(Fit$par, w="mean")

## Make a plot
PDat   <- cbind.data.frame (OBS = Cost_new$residuals$obs,
                            PRE = Cost_new$residuals$mod,
                            RES = Cost_new$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10))


## Estimate the r-squared values
fit <- lm(Log.OBS ~ Log.PRE, data = PDat)
summary(fit)

##MAPE
MAPE=100*mean(abs((PDat$OBS-PDat$PRE)/PDat$OBS))
summary(MAPE)

#segments(data$Time, data$low,data$Time, data$high)
plot(PDat$PRE,PDat$OBS)
plot(PDat$RES)



plot1<- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)

p1 <-
  ggplot(plot1, aes(Log.PRE, Log.OBS)) +
  geom_point  (colour = "steelblue4", size = 4)  +
  geom_abline (intercept = 0,
               slope     = 1,
               color     ="steelblue4", size = 1, alpha = 0.8) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-4,3),labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-4,3),labels = scales::math_format(10^.x))
p1

Data_Cal <-rbind.data.frame(
  
  Obs_A3_1  <-Obs_A3 %>% filter(Matrix == 'P')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A3_1'),
  Obs_A7_1  <-Obs_A7 %>% filter(Matrix == 'L')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A7_1'),
  Obs_A7_2  <-Obs_A7 %>% filter(Matrix == 'M')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A7_2'),
  Obs_A7_3  <-Obs_A7 %>% filter(Matrix == 'K')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A7_3'),
  Obs_A7_4  <-Obs_A7 %>% filter(Matrix == 'F')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A7_4')
)

Out_Cal <- rbind.data.frame(
  #out_A1_1  <- pred (pars=Fit$par, BW = 60.2, area=1200,  Dtimes = 1)%>%
    #select(Time = Time, Conc = CV)%>%mutate(STUDY = 'A1_1'),
  out_A3_1  <- pred (pars=Fit$par, BW = 222, area=3600,Dtimes = 1)%>%
    select(Time = Time, Conc = CV)%>%mutate(STUDY = 'A3_1'),
  out_A7_1  <- pred (pars=Fit$par, BW = 309, area=4800,Dtimes = 1)%>%
    select(Time = Time, Conc = CL)%>%mutate(STUDY = 'A7_1'),
  out_A7_2  <- pred (pars=Fit$par, BW = 309, area=4800,Dtimes = 1)%>%
    select(Time = Time, Conc = CM)%>%mutate(STUDY = 'A7_2'),
  out_A7_3  <- pred (pars=Fit$par, BW = 309, area=4800,Dtimes = 1)%>%
    select(Time = Time, Conc = CK)%>%mutate(STUDY = 'A7_3'),
  out_A7_4  <- pred (pars=Fit$par, BW = 309, area=4800,Dtimes = 1)%>%
    select(Time = Time, Conc = CF)%>%mutate(STUDY = 'A7_4')
)

## Plot
Levels <- c( 'A3_1','A7_1','A7_2','A7_3','A7_4')

Data_Cal$Study <- factor(Data_Cal$STUDY, levels=Levels)
Out_Cal$Study <- factor(Out_Cal$STUDY, levels=Levels)

p2<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
  geom_point(shape = 21, colour = "red", fill = "white", size = 1.5, stroke = 1.2) + 
  geom_line(data = data.frame(Out_Cal), aes(x = Time, y = Conc), 
            size = 0.6, colour = "black") +
  
  scale_y_log10 (labels = function(x) format(x, scientific = TRUE))  + 
  facet_wrap(~Study, scales = "free",ncol=6) +
  theme_bw() +
  theme (
    panel.background = element_rect(fill = "white"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    strip.background        = element_blank(),
    strip.text              = element_blank(),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 10, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 10, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "",  y = "")

p2


ggsave("p2.tiff",scale = 1.5,
       plot = p2,
       width = 25, height = 15, units = "cm", dpi=320)

##//////////////////////////////////////////////////////////////////////////////
## Model evaluation
## Study 11: Kissell et al. (2016)

### Make cost function
Cost_Cattle_eva <- function (pars,BW,area,Dtimes,w) {
  
  # Prediction
  #out_A1      <- pred (pars=pars, BW = 60.2, area=600,Dtimes = 1)
  out_A4      <- pred (pars=pars, BW = 61.3, area=600,Dtimes = 1)
  out_A4.1    <- pred (pars=pars, BW = 61.3, area=600,Dtimes = 1)
  
  cost <- NULL
  
  #cost<- modCost  (model = out_A1,   obs = Obs_A1_P,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A4,   obs = Obs_A4_P,   x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_A4.1, obs = Obs_A4.1_P, x ="Time", cost = cost, weight = w)
  
  return(cost)
}


Cost_Cattle_eva (Fit$par, w="mean")

## Make the data for the plot
Cattle_FLU_eva <-Cost_Cattle_eva (Fit$par, w="mean")


## Make a plot
PDat_eva1   <- cbind.data.frame (OBS1 = Cattle_FLU_eva$residuals$obs,
                                PRE1 = Cattle_FLU_eva$residuals$mod,
                                RES1 = Cattle_FLU_eva$residuals$res)

PDat_eva1 <- PDat_eva1 %>% mutate (Log.OBS1 = log(OBS1, 10), Log.PRE1 = log(PRE1, 10), Species = "Cattle")


## Estimate the r-squared values
fit_eva1 <- lm(Log.OBS1 ~ Log.PRE1, data = PDat_eva1)
summary(fit_eva1)

MAPE=100*mean(abs((PDat_eva1$OBS1-PDat_eva1$PRE1)/PDat_eva1$OBS1))
summary(MAPE)

## Plot the results of model evaluation
PlotDat_eva1 <- PDat_eva1 %>% mutate(prediction = predict(fit_eva1), OPR1 = PRE1/OBS1)


p3 <-
  ggplot(PlotDat_eva1, aes(Log.PRE1, Log.OBS1)) +
  geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
  geom_abline (intercept = 0,
               slope     = 1,
               color     ="steelblue4", size = 1, alpha = 0.8) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-4,2), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-4,2),labels = scales::math_format(10^.x))

p3

Data_Cal1 <-rbind.data.frame(
  
  Obs_A4_1  <-Obs_A4 %>% filter(Matrix == 'P')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A4_1'),
  Obs_A4_2  <-Obs_A4.1 %>% filter(Matrix == 'P')%>%
    select(Time = Time, Conc = Conc.)%>%mutate(STUDY = 'A4_2')
)

Out_Cal1 <- rbind.data.frame(
 
  out_A4_1  <- pred (pars=Fit$par, BW = 61.3,area=1200, Dtimes = 1)%>%
    select(Time = Time, Conc = CV)%>%mutate(STUDY = 'A4_1'),
  out_A4_2  <- pred (pars=Fit$par, BW = 61.3,area=1200,Dtimes = 1)%>%
    select(Time = Time, Conc = CV)%>%mutate(STUDY = 'A4_2')
)

## Plot
Levels <- c('A4_1', 'A4_2')

Data_Cal1$Study <- factor(Data_Cal1$STUDY, levels=Levels)
Out_Cal1$Study <- factor(Out_Cal1$STUDY, levels=Levels)

p4<-ggplot(data = Data_Cal1, aes(x=Time, y=Conc)) + 
  geom_point(shape = 21, colour = "red", fill = "white", size = 1.5, stroke = 1.2) + 
  geom_line(data = data.frame(Out_Cal1), aes(x = Time, y = Conc), 
            size = 0.6, colour = "black") +
  
  scale_y_log10 (labels = function(x) format(x, scientific = TRUE))  + 
  facet_wrap(~Study, scales = "free",ncol=6) +
  theme_bw() +
  theme (
    panel.background = element_rect(fill = "white"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    strip.background        = element_blank(),
    strip.text              = element_blank(),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    axis.text               = element_text (size   = 10, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 10, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
  labs (x = "",  y = "")

p4


ggsave("p4.tiff",scale = 1.5,
       plot = p,
       width = 25, height = 15, units = "cm", dpi=320)

## Save the fitting results
saveRDS(Cattle_FLU_eva, file = 'Cattle_FLU_eva.rds')
