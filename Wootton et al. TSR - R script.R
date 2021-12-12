# This code is associated with the publication "Reproductive trade-offs, not metabolism, make adult fish smaller at warmer temperatures".
# It analyses metabolic data recorded using intermittent-closed respirometry over multiple generations. 
# Further data of other life-history reposes that were analyzed in the manuscript can be found in the github repository.


rm(list = ls()) # clear R environment
#setwd("C:Users/...") # Set the directory to your location of choice
# Or, set the working directory to the location of this Rmd file
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#list.dirs()                    # look at directory structure
library(lme4)
library(ggplot2)
library(ggforce)
library(MASS)
library(FSA)
library(plyr)
library(dplyr)
library(effects)
library(MuMIn)
library(bbmle)
library(emmeans)

# Weight and length analyses ----

# Weight analysis ####
# first, read out data - data should be stored in the same folder as this script file.

weight_data <- read.csv(file = "Wootton et al. TSR - Size data.csv")

glimpse(weight_data)

# check for missing data
apply(weight_data, 2, function(x) any(is.na(x)))
weight_data <- weight_data[!is.na(weight_data$Weight), ]
apply(weight_data, 2, function(x) any(is.na(x)))

# format data
weight_data$Temperature <- as.factor(weight_data$Temperature)
weight_data$Temperature <- relevel(weight_data$Temperature , ref="L")
weight_data$Stage <- as.factor(weight_data$Stage)
weight_data$Stage <- relevel(weight_data$Stage , ref="Juv")

# determine the best model to predict total body weight
WGTM <- lmer(Weight ~ Temperature*Generation*Stage  + (1|Population)+ (1|Generation), data = weight_data, REML = T)

summary(WGTM)

plot(WGTM) # check resid spread 

qqnorm(residuals(WGTM))

# transform the response
WGTM1 <- lmer(sqrt(Weight) ~ Temperature*Generation*Stage + (1|Population) + (1|Generation), data = weight_data, REML = T)

summary(WGTM1)

plot(WGTM1) # check resid spread 

qqnorm(residuals(WGTM1))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=weight_data$Generation,pearson=residuals(WGTM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Stage=weight_data$Stage,pearson=residuals(WGTM1,type="pearson")),aes(x=Stage,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=weight_data$Temperature,pearson=residuals(WGTM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
WGTM1a <- lmer(sqrt(Weight) ~ Temperature*Generation*Stage  + (1|Population) + (1|Generation), data = weight_data, REML = F) 

# specify the model with 2-way interactions
WGTM1b <- lmer(sqrt(Weight) ~ Temperature*Generation + Temperature*Stage + Generation*Stage + (1|Population) + (1|Generation), data = weight_data, REML = F) 

# test the models using anova
anova(WGTM1a, WGTM1b)
#Data: weight_data
#Models:
#WGTM1b: sqrt(Weight) ~ Temperature * Generation + Temperature * Stage + 
#WGTM1b:     Generation * Stage + (1 | Population) + (1 | Generation)
#WGTM1a: sqrt(Weight) ~ Temperature * Generation * Stage + (1 | Population) + 
#WGTM1a:     (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)    
#WGTM1b   13 -3963.3 -3880.3 1994.7  -3989.3                         
#WGTM1a   15 -3973.3 -3877.6 2001.7  -4003.3 14.009  2  0.0009078 ***
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
WGTMbest <- lmer(sqrt(Weight) ~ Temperature*Generation*Stage  + (1|Population) + (1|Generation), data = weight_data, REML = T) 

summary(WGTMbest)
#Fixed effects:
#                                     Estimate Std. Error t value
#(Intercept)                          0.342201   0.030163  11.345
#TemperatureH                        -0.021879   0.017758  -1.232
#Generation                           0.015076   0.007525   2.003
#StageAdult                           0.713114   0.018158  39.272
#StageMature                          0.380662   0.018158  20.963
#TemperatureH:Generation              0.022304   0.003754   5.941
#TemperatureH:StageAdult             -0.077473   0.025665  -3.019
#TemperatureH:StageMature            -0.005708   0.025665  -0.222
#Generation:StageAdult               -0.026710   0.004679  -5.708
#Generation:StageMature              -0.003390   0.004679  -0.725
#TemperatureH:Generation:StageAdult  -0.020255   0.006614  -3.062
#TemperatureH:Generation:StageMature -0.019989   0.006614  -3.022

confint(WGTMbest)
#                                       2.5 %       97.5 %
#.sig01                               0.0046904120  0.024656610
#.sig02                               0.0143131032  0.049490001
#.sigma                               0.1495198961  0.155928149
#(Intercept)                          0.2855136418  0.398929451
#TemperatureH                        -0.0559240968  0.012161042
#Generation                           0.0007641275  0.029377570
#StageAdult                           0.6775957990  0.748718781
#StageMature                          0.3451438103  0.416266793
#TemperatureH:Generation              0.0149534732  0.029657002
#TemperatureH:StageAdult             -0.1277273911 -0.027204782
#TemperatureH:StageMature            -0.0559626747  0.044559934
#Generation:StageAdult               -0.0358825336 -0.017555875
#Generation:StageMature              -0.0125627645  0.005763894
#TemperatureH:Generation:StageAdult  -0.0332095719 -0.007302539
#TemperatureH:Generation:StageMature -0.0329438540 -0.007036821

plot(allEffects(WGTMbest)) # gives good visual of parameter estimates


# Length analysis ####

Length_data <- read.csv(file = "Wootton et al. TSR - Size data.csv")

glimpse(Length_data)

# check and delete any missing data
apply(Length_data, 2, function(x) any(is.na(x)))
Length_data <- Length_data[!is.na(Length_data$Length), ]
apply(Length_data, 2, function(x) any(is.na(x)))

# input data from random size selection episode in each generation. Details of this selection can be found in Wootton et al. (2021).
# note that we did not record weight during selection so do not replicate this approach in the weight analysis
selection_data <- read.csv(file = "Wootton et al. TSR - Selection data.csv")

# now write function to randomise order of cases, then select first 20 cases
Seldataselfunct <- function(data){
  set.seed(35)
  if(length(data$Length)>= 20){
    data1 <- data[sample(nrow(data), 20), ] 
  } else{
    
    data1 <- data
  }
  return(data1)
}

Seldat     <-  ddply(selection_data, .( Compa.No., Generation), Seldataselfunct)

# rbind (stack) growth data and the sampled selection data
Length_data <- rbind(Length_data, Seldat)

# format data
Length_data$Stage <- as.factor(Length_data$Stage)
Length_data$Stage <- relevel(Length_data$Stage , ref="Juv")
Length_data$Temperature <- as.factor(Length_data$Temperature)
Length_data$Temperature <- relevel(Length_data$Temperature , ref="L")

# determine the best model to predict total length
LengthM <- lmer(Length ~ Temperature*Generation*Stage  + (1|Population)+ (1|Generation), data = Length_data, REML = T)

summary(LengthM)

plot(LengthM) # check resid spread 
# outlier is a small fish, not an erroneous measurement
qqnorm(residuals(LengthM))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Length_data$Generation,pearson=residuals(LengthM,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Stage=Length_data$Stage,pearson=residuals(LengthM,type="pearson")),aes(x=Stage,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=Length_data$Temperature,pearson=residuals(LengthM,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
LengthMa <- lmer(Length ~ Temperature*Generation*Stage  + (1|Population) + (1|Generation), data = Length_data, REML = F) 

# specify the model with 2-way interactions
LengthMb <- lmer(Length ~ Temperature*Generation + Temperature*Stage + Generation*Stage + (1|Population) + (1|Generation), data = Length_data, REML = F) 

# test the models using anova
anova(LengthMa, LengthMb)
#Data: Length_data
#Models:
#LengthMb: Length ~ Temperature * Generation + Temperature * Stage + Generation * 
#LengthMb:     Stage + (1 | Population) + (1 | Generation)
#LengthMa: Length ~ Temperature * Generation * Stage + (1 | Population) + 
#LengthMa:     (1 | Generation)
#npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
#LengthMb   13 30159 30244 -15066    30133                         
#LengthMa   15 30130 30228 -15050    30100 32.459  2  8.948e-08 ***
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
LengthMbest <- lmer(Length ~ Temperature*Generation*Stage  + (1|Population) + (1|Generation), data = Length_data, REML = T) 

summary(LengthMbest)
#Fixed effects:
#                                    Estimate Std. Error t value
#(Intercept)                          21.4246     1.1324  18.920
#TemperatureH                         -0.3171     0.4700  -0.675
#Generation                            0.6340     0.2893   2.192
#StageAdult                           21.8465     0.5484  39.838
#StageMature                          12.8945     0.4676  27.579
#TemperatureH:Generation               0.8541     0.1133   7.541
#TemperatureH:StageAdult              -1.5740     0.7754  -2.030
#TemperatureH:StageMature              1.0504     0.6609   1.589
#Generation:StageAdult                -0.8972     0.1414  -6.346
#Generation:StageMature               -0.2141     0.1207  -1.774
#TemperatureH:Generation:StageAdult   -0.7058     0.1999  -3.531
#TemperatureH:Generation:StageMature  -0.9343     0.1706  -5.476

confint(LengthMbest)
#                                         2.5 %      97.5 %
#.sig01                               0.0000000  0.48551212
#.sig02                               0.5752961  1.92032047
#.sigma                               4.5304843  4.71000743
#(Intercept)                         19.2926824 23.55744226
#TemperatureH                        -1.2230026  0.58885819
#Generation                           0.0878373  1.17998749
#StageAdult                          20.7735471 22.92168600
#StageMature                         11.9799026 13.81142236
#TemperatureH:Generation              0.6322392  1.07587114
#TemperatureH:StageAdult             -3.0926910 -0.05538689
#TemperatureH:StageMature            -0.2440658  2.34478511
#Generation:StageAdult               -1.1743279 -0.62052364
#Generation:StageMature              -0.4507579  0.02202099
#TemperatureH:Generation:StageAdult  -1.0973444 -0.31425793
#TemperatureH:Generation:StageMature -1.2685216 -0.60011611

plot(allEffects(LengthMbest)) # gives good visual of parameter estimates

# Metabolism analyses ----
# first, read out and visualise the data - data should be stored in the same folder as this script file.

met_data <- read.csv(file = "Wootton et al. TSR - Metabolic data.csv")

# visualise and determine the oulier selection criteria for SMR
plot(met_data$Weight, met_data$SMR, pch = 19)
abline(h = 0.15) 
abline (v = 1) 
which(met_data$Weight > 1 & met_data$SMR < 0.15) # 1, 2 and 180 - We went back and checked the data sheets and the weights seem to have been recorded incorrectly here (e.g. large weight placed into a small chamber) - so we remove the whole row of data
# The missing values for SMR in the dataframe are where individuals died during metabolic trails

# visualize and determine the outlier selection criteria for MMR
plot(met_data$Weight, met_data$MMR, pch = 19)
abline(h=6)
abline(v=1)
abline(v=.4)
abline(h = 0.4)
which(met_data$MMR > 6) # 120

# remove outliers
met_data <- met_data[-which(met_data$Weight > 1 & met_data$SMR < 0.15),]
met_data <- met_data[-which(met_data$MMR > 6),]

# re-order the dataframe index
row.names(met_data) <- NULL

# convert the data into the correct format and order the data by generation
met_data$Temperature <- as.factor(met_data$Temperature)
met_data$Size <- as.factor(met_data$Size)
met_data$Stage <- as.factor(met_data$Stage)

# check to see if the weight relationship to SMR is linear
plot(met_data$Weight, met_data$SMR, pch = 19)
# check to see if the weight relationship to SMR is log-log
plot(log(met_data$Weight), log(met_data$SMR), pch = 19) 

# convert SMR to mass specific SMR in mg O2 kg-1 h-1
met_data$SMRms <- (met_data$SMR/met_data$Weight)*1000 

# check that the data seems reasonable
mean(met_data$SMRms, na.rm = T) #379 mg O2 kg-1 h-1; this compares to ca 700 in damselfish and 100 in perch
median(met_data$SMRms, na.rm = T) #338
max(met_data$SMRms, na.rm = T) #1444
min(met_data$SMRms, na.rm = T) #2.7 - this is unexpectedly low 

# find and remove the outlier
which.min(met_data$SMRms) # tiny SMR for a large fish - exclude this SMR value
met_data$SMR[met_data$SMRms < 3] <- NA # this removes the SMRvalue 
met_data$SMRms[met_data$SMRms < 3] <- NA # this removes the mass specific value

# check for normal distribution
hist(met_data$SMRms, breaks = 50) 

# convert MMR to mass specific MMR in mg O2 kg-1 h-1
met_data$MMRms <- (met_data$MMR/met_data$Weight)*1000

# check that the data seems reasonable
mean(met_data$MMRms, na.rm = T) #2348 mg O2 kg-1 h-1; this compares to ca 2000 in damselfish and 320 in perch 
median(met_data$MMRms, na.rm = T) #2185
max(met_data$MMRms, na.rm = T) #7864
min(met_data$MMRms, na.rm = T) #250

# check for normal distribution
hist(met_data$MMRms, breaks = 50)

# calculate aerobic scope and factorial aerobic scope
met_data$absolScope <- met_data$MMR - met_data$SMR
met_data$factScope <- met_data$MMR/met_data$SMR

# visualise and determine the oulier selection criteria for aerobic scope
plot(met_data$Weight, met_data$absolScope, pch = 19)
abline(h=0.3)
abline(v=0.4)
abline(h = 0.4)

which(met_data$Weight > 0.4 & met_data$absolScope < 0.3) # 3 23 108 115 340 
# low scope is plausible, however the MMR values look weirdly small so I will remove them
which(met_data$absolScope <= 0) # 23 108 
# These are obviously errors, MMR is lower than SMR so we will remove the MMR values
which(met_data$Weight < 1 & met_data$absolScope > 2.5) # 114 171 173 175 223
# these have high scope due to high (relative to weight) MMR readings. This is plausible so we will not exclude

# remove the chosen outliers
met_data$MMR[met_data$Weight > 0.4 & met_data$absolScope < 0.3] <- NA 
met_data$absolScope[met_data$Weight > 0.4 & met_data$absolScope < 0.3] <- NA 
met_data$MMR[met_data$absolScope <= 0] <- NA # this removes the MMR values for rows that have a negative slope
met_data$absolScope[met_data$absolScope <= 0] <- NA # this removes the negative absol scope values (preserves SMR which we assume to be ok)

# lets look at the spread of weight data analysed here (i.e. lets check the data is comparable across temperatures)
WGTViolin<-ggplot(met_data, aes(x= Temperature, y = Weight, colour = Temperature)) + geom_violin(aes(color = Temperature, fill = Temperature), alpha = 0.03) +
  geom_sina(aes( color = Temperature), alpha = 0.3)+
  labs(x = "Temperature", y = "Somatic weight (grams)")+ 
  theme_bw()  + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), legend.position = "right") + 
  scale_colour_manual(values = c( "deepskyblue", "red" )) + scale_fill_manual(values = c( "deepskyblue", "red"  ))+  facet_grid(Stage~Generation, labeller = label_both)

WGTViolin

# SMR analysis #####

# determine the best model to predict SMR
# create a copy of the dataframe 
met_data_Full_SMR <- met_data

# remove missing data from the dataframe
met_data_Full_SMR <- met_data_Full_SMR %>%dplyr:: filter (SMR != "NA")

# format data
met_data_Full_SMR$Temperature <- as.factor(met_data_Full_SMR$Temperature)
met_data_Full_SMR$Temperature <- relevel(met_data_Full_SMR$Temperature , ref="L")

# check to see if SMR should be transformed
SMRM <- lmer(SMR ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_SMR, REML = T) 
# is the model (almost/near) singular
isSingular(SMRM, tol = 1e-07) # this is driven by the population random effect, which we retain due to its structural importance
# check the summary table
summary(SMRM)
# view the residual spread
plot(SMRM)
# view the normal Q-Q plot
qqnorm(residuals(SMRM))

# repeat the above model with SMR transformed (square rooted)
SMRM1 <- lmer(sqrt(SMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_SMR, REML = T) 
# is the model (almost/near) singular
isSingular(SMRM1, tol = 1e-07) 
# check the summary table
summary(SMRM1)
# view the residual spread
plot(SMRM1) 
# view the normal Q-Q plot
qqnorm(residuals(SMRM1))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=met_data_Full_SMR$Generation,pearson=residuals(SMRM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=met_data_Full_SMR$Weight,pearson=residuals(SMRM1,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=met_data_Full_SMR$Temperature,pearson=residuals(SMRM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

# now check if a random effect of chamber position is important for the model.
# (i.e. the position of each respirometry chamber within the waterbath may have been closer to a bubbler etc. and so we test the inclusion of this random effect)
SMRM1Chambertest <- lmer(sqrt(SMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_SMR, REML = T) 
summary(SMRM1Chambertest) #Chamber explains no variation so we proceed without fitting it

# now run the model selection with SMR square rooted. Across analyses we use backwards stepwise regression to test for the significance of interactions of interest. 
# specify the full model
SMRM1a <- lmer(sqrt(SMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_SMR, REML = F) 
# specify the model with 2-way interactions
SMRM1b <- lmer(sqrt(SMR) ~ Temperature*Generation + Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_SMR, REML = F) 

# test the models using anova
anova(SMRM1a, SMRM1b)
#Data: met_data_Full_SMR
#Models:
#SMRMF1b: sqrt(SMR) ~ Temperature * Generation + Temperature * Weight + 
#SMRMF1b:     Generation * Weight + (1 | Population) + (1 | Generation)
#SMRMF1a: sqrt(SMR) ~ Temperature * Generation * Weight + (1 | Population) + 
#SMRMF1a:     (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)   
#SMRMF1b   10 -840.74 -802.11 430.37  -860.74                        
#SMRMF1a   11 -846.38 -803.88 434.19  -868.38 7.6324  1   0.005733 **
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
SMRM1best <- lmer(sqrt(SMR) ~Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_SMR, REML = T) 

# check the model summary
summary(SMRM1best)
#Fixed effects:
#                                 Estimate Std. Error         df t value   
#(Intercept)                      0.250019   0.033373  13.764617   7.492 
#Generation                      -0.003870   0.008575  13.814322  -0.451    
#TemperatureH                    -0.059418   0.034381 340.349231  -1.728   
#Weight                           0.263172   0.035571 341.579724   7.399 
#Generation:TemperatureH          0.019036   0.008975 340.317579   2.121   
#Generation:Weight                0.009368   0.009593 340.972263   0.977  
#TemperatureH:Weight              0.165541   0.051886 340.590097   3.190  
#Generation:TemperatureH:Weight  -0.038539   0.013988 340.474437  -2.755  

# check the confidence intervals
confint(SMRM1best)
#                                        2.5 %     97.5 %
#.sig01                          0.009768404  0.042003495
#.sig02                          0.000000000  0.012916094
#.sigma                          0.064583243  0.074967410
#(Intercept)                     0.188601346  0.311423642
#Generation                     -0.019648676  0.011912491
#TemperatureH                   -0.126347926  0.007620258
#Weight                          0.194036802  0.332556183
#Generation:TemperatureH         0.001540405  0.036510403
#Generation:Weight              -0.009348041  0.028019745
#TemperatureH:Weight             0.064414425  0.266551545
#Generation:TemperatureH:Weight -0.065760914 -0.011265765

plot(allEffects(SMRM1best)) # gives good visual of parameter estimates

# MMR analysis #####

# determine the best model to predict MMR
# create a copy of the dataframe 
met_data_Full_MMR <- met_data

# remove missing data from the dataframe
met_data_Full_MMR <- met_data_Full_MMR %>%dplyr:: filter (MMR != "NA")

# format data
met_data_Full_MMR$Temperature <- as.factor(met_data_Full_MMR$Temperature)
met_data_Full_MMR$Temperature <- relevel(met_data_Full_MMR$Temperature , ref="L")

# check to see if MMR should be transformed
MMRM <- lmer(MMR ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = T) 
# is the model (almost/near) singular
isSingular(MMRM, tol = 1e-07) 
# check the model summary
summary(MMRM)
# view the residual spread
plot(MMRM)
# view the normal Q-Q plot
qqnorm(residuals(MMRM))

#repeat the above model with MMR square rooted
MMRM1 <- lmer(sqrt(MMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = T) 
# is the model (almost/near) singular
isSingular(MMRM1, tol = 1e-07) 
# check the model summary
summary(MMRM1)
# view the residual spread
plot(MMRM1)
# view the normal Q-Q plot
qqnorm(residuals(MMRM1))
# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=met_data_Full_MMR$Generation,pearson=residuals(MMRM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=met_data_Full_MMR$Weight,pearson=residuals(MMRM1,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=met_data_Full_MMR$Temperature,pearson=residuals(MMRM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

# now check if a random effect of chamber position is important for the model.
# (i.e. the position of each respirometry chamber within the waterbath may have been closer to a bubbler etc. and so we test the inclusion of this random effect)
MMRM1Chambertest <- lmer(sqrt(MMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_MMR, REML = T) 
summary(MMRM1Chambertest) #Chamber explains no variation so we proceed without fitting it

# specify the full model
MMRM1a <- lmer(sqrt(MMR) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = F) 
# specify the model with 2-way interactions
MMRM1b <- lmer(sqrt(MMR) ~ Temperature*Generation + Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = F) 

# test the models using anova
anova(MMRM1a, MMRM1b)
#Data: met_data_Full_MMR
#Models:
#MMRM1b: sqrt(MMR) ~ Temperature * Generation + Temperature * Weight + 
#MMRM1b:     Generation * Weight + (1 | Population) + (1 | Generation)
#MMRM1a: sqrt(MMR) ~ Temperature * Generation * Weight + (1 | Population) + 
#MMRM1a:     (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#MMRM1b   10 -233.79 -195.19 126.90  -253.79                     
#MMRM1a   11 -232.05 -189.59 127.03  -254.05 0.2604  1     0.6098

# now test model fits with 2-way interactions of interest dropped (see methods in the manuscript for descriptions of interactions of interest)
# drop Temperature*Generation from the model
MMRM2a <- lmer(sqrt(MMR) ~  Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = F) 
anova(MMRM1b, MMRM2a)
#Data: met_data_Full_MMR
#Models:
#MMRM2a: sqrt(MMR) ~ Temperature * Weight + Generation * Weight + (1 | 
#MMRM2a:     Population) + (1 | Generation)
#MMRM1b: sqrt(MMR) ~ Temperature * Generation + Temperature * Weight + 
#MMRM1b:     Generation * Weight + (1 | Population) + (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#MMRM2a    9 -235.00 -200.25  126.5  -253.00                     
#MMRM1b   10 -233.79 -195.19  126.9  -253.79 0.7965  1     0.3721

# drop Temperature*Weight from the model
MMRM2b <- lmer(sqrt(MMR) ~ Temperature*Generation  + Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = F) 
anova(MMRM1b, MMRM2b)
#Data: met_data_Full_MMR
#Models:
#MMRM2b: sqrt(MMR) ~ Temperature * Generation + Generation * Weight + 
#MMRM2b:     (1 | Population) + (1 | Generation)
#MMRM1b: sqrt(MMR) ~ Temperature * Generation + Temperature * Weight + 
#MMRM1b:     Generation * Weight + (1 | Population) + (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#MMRM2b    9 -235.74 -200.99 126.87  -253.74                     
#MMRM1b   10 -233.79 -195.19 126.90  -253.79 0.0563  1     0.8124

# dropping the interactions of interest did not degrate fit so we fit and interpret the null model
MMRMnull <- lmer(sqrt(MMR) ~ Temperature + Generation + Weight  + (1|Population) + (1|Generation), data = met_data_Full_MMR, REML = T) 

# check the model summary
summary(MMRMnull)
#Fixed effects:
#             Estimate Std. Error t value
#(Intercept)  0.529728   0.035538  14.906
#TemperatureH 0.102694   0.026119   3.932
#Generation   0.003101   0.006806   0.456
#Weight       0.822006   0.027961  29.398  

# check the confidence intervals
confint(MMRMnull)
#             2.5 %     97.5 %
#.sig01        0.000000000 0.04153088
#.sig02        0.000000000 0.04257526
#.sigma        0.156793531 0.18220426
#(Intercept)   0.467040825 0.59180230
#TemperatureH  0.055066646 0.15063232
#Generation   -0.009923092 0.01610100
#Weight        0.766975848 0.87678174

plot(allEffects(MMRMnull)) # gives good visual of parameter estimates

# AAS analysis #####

# determine the best model to predict AAS
# create a copy of the dataframe 
met_data_Full_AAS <- met_data

# remove missing data from the dataframe
met_data_Full_AAS <- met_data_Full_AAS %>%dplyr:: filter (absolScope != "NA")

# format data
met_data_Full_AAS$Temperature <- as.factor(met_data_Full_AAS$Temperature)
met_data_Full_AAS$Temperature <- relevel(met_data_Full_AAS$Temperature , ref="L")

# check to see if AAS should be transformed
AASM <- lmer(absolScope ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_AAS, REML = T) 
# is the model (almost/near) singular
isSingular(AASM, tol = 1e-07) 
# check the model summary
summary(AASM)
# view the residual spread
plot(AASM)
# view the normal Q-Q plot
qqnorm(residuals(AASM))

#repeat the above model with MMR square rooted
AASM1 <- lmer(sqrt(absolScope) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation), data = met_data_Full_AAS, REML = T) 
# is the model (almost/near) singular
isSingular(AASM1, tol = 1e-07) 
# check the model summary
summary(AASM1)
# view the residual spread
plot(AASM1) 
# view the normal Q-Q plot
qqnorm(residuals(AASM1))
# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=met_data_Full_AAS$Generation,pearson=residuals(AASM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=met_data_Full_AAS$Weight,pearson=residuals(AASM1,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=met_data_Full_AAS$Temperature,pearson=residuals(AASM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

# now check if a random effect of chamber position is important for the model.
# (i.e. the position of each respirometry chamber within the waterbath may have been closer to a bubbler etc. and so we test the inclusion of this random effect)
AASM1Chambertest <- lmer(sqrt(absolScope) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = T) 
summary(AASM1Chambertest) #Chamber explains some variation so we will fit it

#repeat the above model with MMR square rooted and a random effect of chamber
AASM2 <- lmer(sqrt(absolScope) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = T) 
# is the model (almost/near) singular
isSingular(AASM2, tol = 1e-07) 
# check the model summary
summary(AASM2)
# view the residual spread
plot(AASM2) 
# view the normal Q-Q plot
qqnorm(residuals(AASM2))
# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=met_data_Full_AAS$Generation,pearson=residuals(AASM2,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=met_data_Full_AAS$Weight,pearson=residuals(AASM2,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=met_data_Full_AAS$Temperature,pearson=residuals(AASM2,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

# specify the full model
AASM2a <- lmer(sqrt(absolScope) ~ Temperature*Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = F) 
# specify the model with 2-way interactions
AASM2b <- lmer(sqrt(absolScope) ~ Temperature*Generation + Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = F) 

# test the models using anova
anova(AASM2a, AASM2b)
#Data: met_data_Full_AAS
#Models:
#AASM2b: sqrt(absolScope) ~ Temperature * Generation + Temperature * Weight + 
#AASM2b:     Generation * Weight + (1 | Population) + (1 | Generation) + 
#AASM2b:     (1 | Chamber)
#AASM2a: sqrt(absolScope) ~ Temperature * Generation * Weight + (1 | Population) + 
#AASM2a:     (1 | Generation) + (1 | Chamber)
#       npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AASM2b   11 -194.58 -152.23 108.29  -216.58                     
#AASM2a   12 -193.03 -146.83 108.51  -217.03 0.4484  1     0.5031

# now test model fits with 2-way interactions of interest dropped (see methods in the manuscript for descriptions of interactions of interest)
# drop Temperature*Generation from the model
AASM3a <- lmer(sqrt(absolScope) ~  Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = F) 
anova(AASM2b, AASM3a)
#Data: met_data_Full_AAS
#Models:
#AASM2a: sqrt(absolScope) ~ Temperature * Weight + Generation * Weight + 
#AASM2a:     (1 | Population) + (1 | Generation) + (1 | Chamber)
#AASM1b: sqrt(absolScope) ~ Temperature * Generation + Temperature * Weight + 
#AASM1b:     Generation * Weight + (1 | Population) + (1 | Generation) + 
#AASM1b:     (1 | Chamber)
#       npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AASM2a   10 -195.27 -156.78 107.64  -215.27                     
#AASM1b   11 -194.58 -152.23 108.29  -216.58 1.3034  1     0.2536

# drop Temperature*Weight from the model
AASM3b <- lmer(sqrt(absolScope) ~ Temperature*Generation  + Generation*Weight + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = F) 
anova(AASM2b, AASM3b)
#Data: met_data_Full_AAS
#Models:
#AASM2b: sqrt(absolScope) ~ Temperature * Generation + Generation * Weight + 
#AASM2b:     (1 | Population) + (1 | Generation) + (1 | Chamber)
#AASM1b: sqrt(absolScope) ~ Temperature * Generation + Temperature * Weight + 
#AASM1b:     Generation * Weight + (1 | Population) + (1 | Generation) + 
#AASM1b:     (1 | Chamber)
#       npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AASM2b   10 -196.57 -158.08 108.29  -216.57                     
#AASM1b   11 -194.58 -152.23 108.29  -216.58 0.0021  1     0.9633

# dropping the interactions of interest did not degrade fit so we fit and interpret the null model
AASMnull <- lmer(sqrt(absolScope) ~ Temperature + Generation + Weight  + (1|Population) + (1|Generation) + (1|Chamber), data = met_data_Full_AAS, REML = T)  

# check the model summary
summary(AASMnull)
#Fixed effects:
#             Estimate Std. Error t value
#(Intercept)  0.488161   0.033770  14.455
#TemperatureH 0.093701   0.028352   3.305
#Generation   0.001632   0.005625   0.290
#Weight       0.759372   0.029706  25.563 

# check the confidence intervals
confint(AASMnull)
#                  2.5 %     97.5 %
#.sig01        0.000000000 0.03774942
#.sig02        0.000000000 0.03306982
#.sig03        0.000000000 0.04522784
#.sigma        0.165443294 0.19215812
#(Intercept)   0.426414058 0.54863759
#TemperatureH  0.042833095 0.14455193
#Generation   -0.009847959 0.01310397
#Weight        0.702001868 0.82116083

plot(allEffects(AASMnull)) # gives good visual of parameter estimates

# Maturity models ----

# Weight at maturity ####

W50_data <- read.csv(file = "Wootton et al. TSR - Weight at maturity data.csv")

glimpse(W50_data)

# format data
W50_data$Temperature <- as.factor(W50_data$Temperature)
W50_data$Temperature <- relevel(W50_data$Temperature , ref="L")

# determine the best model to predict weight at maturity
W50M <- lmer(W50 ~ Temperature*Generation + (1|Population)+ (1|Generation), data = W50_data, REML = T)

summary(W50M)

plot(W50M) # check resid spread 

qqnorm(residuals(W50M))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=W50_data$Generation,pearson=residuals(W50M,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=W50_data$Temperature,pearson=residuals(W50M,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
W50Ma <- lmer(W50 ~ Temperature*Generation + (1|Population) + (1|Generation), data = W50_data, REML = F) 

# specify the null model
W50Mb <- lmer(W50 ~ Temperature + Generation + (1|Population) + (1|Generation), data = W50_data, REML = F) 

# test the models using anova
anova(W50Ma, W50Mb)
#Data: W50_data
#Models:
#W50Mb: W50 ~ Temperature + Generation + (1 | Population) + (1 | Generation)
#W50Ma: W50 ~ Temperature * Generation + (1 | Population) + (1 | Generation)
#           npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#W50Mb    6 -92.029 -82.528 52.015  -104.03                     
#W50Ma    7 -90.600 -79.516 52.300  -104.60 0.5713  1     0.4497

# as the interaction of interest is not supported, we fit and interpret the null model
W50Mnull <- lmer(W50 ~  Temperature + Generation + (1|Population) + (1|Generation), data = W50_data, REML = T) 

summary(W50Mnull)
#Fixed effects:
#             Estimate Std. Error t value
#(Intercept)   0.26956    0.03000   8.984
#TemperatureH -0.05424    0.01912  -2.836
#Generation    0.01298    0.00730   1.778

confint(W50Mnull)
#                    2.5 %      97.5 %
#.sig01        0.000000000  0.03336605
#.sig02        0.000000000  0.04441368
#.sigma        0.044444300  0.07316350
#(Intercept)   0.213384710  0.32572999
#TemperatureH -0.093356773 -0.01511715
#Generation   -0.000846669  0.02680876

plot(allEffects(W50Mnull)) # gives good visual of parameter estimates

# Age at maturity ####

A50_data <- read.csv(file = "Wootton et al. TSR - Age at maturity data.csv")

glimpse(A50_data)

# format data
A50_data$Temperature <- as.factor(A50_data$Temperature)
A50_data$Temperature <- relevel(A50_data$Temperature , ref="L")

# determine the best model to predict age at maturity
A50M <- lmer(A50 ~ Temperature*Generation + (1|Population)+ (1|Generation), data = A50_data, REML = T)

summary(A50M)

plot(A50M) # check resid spread 

qqnorm(residuals(A50M))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=A50_data$Generation,pearson=residuals(A50M,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=A50_data$Temperature,pearson=residuals(A50M,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
A50Ma <- lmer(A50 ~ Temperature*Generation + (1|Population) + (1|Generation), data = A50_data, REML = F) 

# specify the null model
A50Mb <- lmer(A50 ~ Temperature + Generation + (1|Population) + (1|Generation), data = A50_data, REML = F) 

# test the models using anova
anova(A50Ma, A50Mb)
#Data: A50_data
#Models:
#A50Mb: A50 ~ Temperature + Generation + (1 | Population) + (1 | Generation)
#A50Ma: A50 ~ Temperature * Generation + (1 | Population) + (1 | Generation)
#           npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)   
#A50Mb    6 233.58 243.08 -110.79   221.58                        
#A50Ma    7 228.31 239.39 -107.15   214.31 7.2722  1   0.007003 **
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# fit the supported model
A50Mbest <- lmer(A50 ~  Temperature*Generation + (1|Population) + (1|Generation), data = A50_data, REML = T) 

summary(A50Mbest)
#Fixed effects:
#                        Estimate Std. Error t value
#(Intercept)              69.8512     4.9219  14.192
#TemperatureH              0.4155     3.2596   0.127
#Generation                0.8506     1.2526   0.679
#TemperatureH:Generation  -2.2861     0.8024  -2.849

confint(A50Mbest)
#                            2.5 %     97.5 %
#.sig01                   0.000000  3.3381475
#.sig02                   1.828841  7.9798384
#.sigma                   3.130990  5.3931484
#(Intercept)             60.668326 79.0341229
#TemperatureH            -5.902381  6.7333264
#Generation              -1.500245  3.2013726
#TemperatureH:Generation -3.894298 -0.6779594

plot(allEffects(A50Mbest)) # gives good visual of parameter estimates

# Gonad weight analyses ----

# Recently mature female gonad weight ####

Mat_Gonad_data <- read.csv(file = "Wootton et al. TSR - Recently matured female reproductive data.csv")

# format data
Mat_Gonad_data$Temperature <- as.factor(Mat_Gonad_data$Temperature)
Mat_Gonad_data$Temperature <- relevel(Mat_Gonad_data$Temperature , ref="L")

# determine the best model to predict recently matured gonad weight
MatGM <- lmer(Gonad.weight ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation), data = Mat_Gonad_data, REML = T)

summary(MatGM)

plot(MatGM) # check resid spread 

qqnorm(residuals(MatGM))

# transform response
MatGM1 <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation), data = Mat_Gonad_data, REML = T)

summary(MatGM1)

plot(MatGM1) # check resid spread 

qqnorm(residuals(MatGM1))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Mat_Gonad_data$Generation,pearson=residuals(MatGM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=Mat_Gonad_data$Weight,pearson=residuals(MatGM1,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=Mat_Gonad_data$Temperature,pearson=residuals(MatGM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
MatGM1a <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation) , data = Mat_Gonad_data, REML = F)

# specify the model with 2-way interactions
MatGM1b <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation + Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = Mat_Gonad_data, REML = F) 

# test the models using anova
anova(MatGM1a, MatGM1b)
#Data: Mat_Gonad_data
#Models:
#MatGM1b: sqrt(Gonad.weight) ~ Temperature * Generation + Temperature * 
#MatGM1b:     Weight + Generation * Weight + (1 | Population) + (1 | Generation)
#MatGM1a: sqrt(Gonad.weight) ~ Temperature * Generation * Weight + (1 | 
#                                                                        MatGonam1a:     Population) + (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)  
#MatGM1b   10 -537.71 -506.35 278.85  -557.71                       
#MatGM1a   11 -539.79 -505.29 280.89  -561.79 4.0802  1    0.04339 *
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now fit the supported model
MatGMbest <- lmer(sqrt(Gonad.weight) ~   Temperature*Generation*Weight + (1|Population) + (1|Generation) , data = Mat_Gonad_data, REML = T) 

summary(MatGMbest)
#Fixed effects:
#                                 Estimate Std. Error t value
#(Intercept)                      -0.02410    0.04489  -0.537
#TemperatureH                     -0.02478    0.06081  -0.408
#Generation                       -0.01119    0.01330  -0.842
#Weight                            0.37064    0.07824   4.737
#TemperatureH:Generation           0.02687    0.01710   1.572
#TemperatureH:Weight               0.11639    0.10432   1.116
#Generation:Weight                 0.02434    0.02268   1.073
#TemperatureH:Generation:Weight   -0.05363    0.02857  -1.877 

confint(MatGMbest)
#                                        2.5 %       97.5 %
#.sig01                            0.000000000  0.015484501
#.sig02                            0.000000000  0.019831916
#.sigma                            0.041445007  0.051447252
#(Intercept)                      -0.109666663  0.062036427
#TemperatureH                     -0.146794849  0.087603951
#Generation                       -0.037488128  0.013449006
#Weight                            0.216531577  0.519580323
#TemperatureH:Generation          -0.004280244  0.062003333
#TemperatureH:Weight              -0.076156077  0.327138756
#Generation:Weight                -0.017592465  0.069976246
#TemperatureH:Generation:Weight   -0.112515007 -0.001761474

plot(allEffects(MatGMbest)) # gives good visual of parameter estimates

# Older adult (post spawn) female gonad weight ####
# Note that we combine fresh and stored samples here. Comparison of these stored samples can be found in the manuscript

Adult_Gonad_dataFull <- read.csv(file = "Wootton et al. TSR - Older adult female reproductive data.csv")

Adult_Gonad_data <- Adult_Gonad_dataFull %>% filter (Sample.type == "Wet")

Adult_Gonad_dataEtoH <- Adult_Gonad_dataFull %>% filter (Sample.type == "ETOH")

# calculate the correction for stored samples
Adult_Gonad_dataEtoH$Gonad.weight <- (Adult_Gonad_dataEtoH$Gonad.weight*1.39)+0.012935
Adult_Gonad_dataEtoH$Weight <- (Adult_Gonad_dataEtoH$Weight*1.23481)+0.07660

# combine datasets 
Adult_Gonad_data <- rbind(Adult_Gonad_data,  Adult_Gonad_dataEtoH)

# check and delete any missing data
apply(Adult_Gonad_data, 2, function(x) any(is.na(x)))
Adult_Gonad_data <- Adult_Gonad_data[!is.na(Adult_Gonad_data$Gonad.weight), ]
apply(Adult_Gonad_data, 2, function(x) any(is.na(x)))

# format data
Adult_Gonad_data$Temperature <- as.factor(Adult_Gonad_data$Temperature)
Adult_Gonad_data$Temperature <- relevel(Adult_Gonad_data$Temperature , ref="L")

# determine the best model to predict older adult gonad weight
AdultGM <- lmer(Gonad.weight ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = T)

summary(AdultGM)

plot(AdultGM) # check resid spread 

qqnorm(residuals(AdultGM))

# transform the response
AdultGM1 <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = T)

summary(AdultGM1)

plot(AdultGM1) # check resid spread 

qqnorm(residuals(AdultGM1))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=Adult_Gonad_data$Generation,pearson=residuals(AdultGM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Weight=Adult_Gonad_data$Weight,pearson=residuals(AdultGM1,type="pearson")),aes(x=Weight,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=Adult_Gonad_data$Temperature,pearson=residuals(AdultGM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()

#now go ahead with model selection
AdultGM1a <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation*Weight  + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = F)

# specify the model with 2-way interactions
AdultGM1b <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation + Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = F) 

# test the models using anova
anova(AdultGM1a, AdultGM1b)
#Data: Adult_Gonad_data
#Models:
#AdultGM1b: sqrt(Gonad.weight) ~ Temperature * Generation + Temperature * 
#AdultGM1b:     Weight + Generation * Weight + (1 | Population) + (1 | Generation)
#AdultGM1a: sqrt(Gonad.weight) ~ Temperature * Generation * Weight + (1 | 
#AdultGM1a:     Population) + (1 | Generation)
#                 npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AdultGM1b   10 -279.67 -250.54 149.83  -299.67                     
#AdultGM1a   11 -277.76 -245.72 149.88  -299.76 0.0949  1     0.7581

# now test model fits with 2-way interactions of interest dropped 
# drop Temperature*Generation from the model
AdultGM2a <- lmer(sqrt(Gonad.weight) ~  Temperature*Weight + Generation*Weight + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = F) 
anova(AdultGM1b, AdultGM2a)
#Data: Adult_Gonad_data
#Models:
#AdultGM2a: sqrt(Gonad.weight) ~ Temperature * Weight + Generation * Weight + 
#AdultGM2a:     (1 | Population) + (1 | Generation)
#AdultGM1b: sqrt(Gonad.weight) ~ Temperature * Generation + Temperature * 
#AdultGM1b:     Weight + Generation * Weight + (1 | Population) + (1 | Generation)
#          npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AdultGM2a    9 -281.00 -254.79 149.50  -299.00                     
#AdultGM1b   10 -279.67 -250.54 149.83  -299.67 0.6637  1     0.4153

# drop Temperature*Weight from the model
AdultGM2b <- lmer(sqrt(Gonad.weight) ~ Temperature*Generation  + Generation*Weight + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = F) 
anova(AdultGM1b, AdultGM2b)
#Data: Adult_Gonad_data
#Models:
#AdultGM2b: sqrt(Gonad.weight) ~ Temperature * Generation + Generation * 
#AdultGM2b:     Weight + (1 | Population) + (1 | Generation)
#AdultGM1b: sqrt(Gonad.weight) ~ Temperature * Generation + Temperature * 
#AdultGM1b:     Weight + Generation * Weight + (1 | Population) + (1 | Generation)
#          npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#AdultGM2b    9 -281.13 -254.91 149.56  -299.13                     
#AdultGM1b   10 -279.67 -250.54 149.83  -299.67 0.5401  1     0.4624

# as no interactions of interest are supported, we fit and interpret the null model
AdultGMnull <- lmer(sqrt(Gonad.weight) ~ Temperature + Generation + Weight + (1|Population) + (1|Generation), data = Adult_Gonad_data, REML = T) 

summary(AdultGMnull)
#Fixed effects:
#               Estimate Std. Error t value
#(Intercept)  -0.0474147  0.0512788  -0.925
#TemperatureH -0.0005284  0.0248762  -0.021
#Generation    0.0180463  0.0084871   2.126
#Weight        0.3468753  0.0206955  16.761

confint(AdultGMnull)
#                   2.5 %     97.5 %
#.sig01        0.000000000 0.04662683
#.sig02        0.000000000 0.04033334
#.sigma        0.070329662 0.09010582
#(Intercept)  -0.143663921 0.04980424
#TemperatureH -0.048918511 0.04625656
#Generation    0.002066069 0.03490466
#Weight        0.304958349 0.38615262

plot(allEffects(AdultGMnull)) # gives good visual of parameter estimates

# Feed rate analysis ----

PG2F5Feedratedat <- read.csv(file = "Wootton et al. TSR - Feed rate data.csv")

# format data
PG2F5Feedratedat$Temperature <- as.factor(PG2F5Feedratedat$Temperature)
PG2F5Feedratedat$Temperature <- relevel(PG2F5Feedratedat$Temperature , ref="L")

PG2F5Feedratedat$Generation <- as.factor(PG2F5Feedratedat$Generation)
PG2F5Feedratedat$Generation <- relevel(PG2F5Feedratedat$Generation , ref="Gen1")

#  determine the best model to predict feed rate per fish per week
FeedRM <- lmer(Feed.rate.per.fish ~ Temperature*Generation*ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = T)

# check the summary table
summary(FeedRM)
# view the residual spread
plot(FeedRM) 
# view the normal Q-Q plot
qqnorm(residuals(FeedRM))

# repeat the above model with feed rate transformed (log)
FeedRM1 <- lmer(log(Feed.rate.per.fish) ~ Temperature*Generation*ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = T)

summary(FeedRM1)

plot(FeedRM1) # check resid spread 

qqnorm(residuals(FeedRM1))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=PG2F5Feedratedat$Generation,pearson=residuals(FeedRM1,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=PG2F5Feedratedat$Temperature,pearson=residuals(FeedRM1,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(ave.weight.of.population=PG2F5Feedratedat$ave.weight.of.population,pearson=residuals(FeedRM1,type="pearson")),aes(x=ave.weight.of.population,y=pearson)) + geom_point() + theme_bw()

# log transform the average population weight variable

PG2F5Feedratedat$Log.ave.weight.of.population <- log(PG2F5Feedratedat$ave.weight.of.population)

FeedRM2 <- lmer(log(Feed.rate.per.fish) ~ Temperature*Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = T)

summary(FeedRM2)

plot(FeedRM2) # check resid spread 

qqnorm(residuals(FeedRM2))

# check linearity of fixed terms term (randomly distributed data indicates linearity)
ggplot(data.frame(Generation=PG2F5Feedratedat$Generation,pearson=residuals(FeedRM2,type="pearson")),aes(x=Generation,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Temperature=PG2F5Feedratedat$Temperature,pearson=residuals(FeedRM2,type="pearson")),aes(x=Temperature,y=pearson)) + geom_point() + theme_bw()
ggplot(data.frame(Log.ave.weight.of.population=PG2F5Feedratedat$Log.ave.weight.of.population,pearson=residuals(FeedRM2,type="pearson")),aes(x=Log.ave.weight.of.population,y=pearson)) + geom_point() + theme_bw()
# this looks ok now

#now go ahead with model selection
FeedRM2a <- lmer(log(Feed.rate.per.fish) ~ Temperature*Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = F)

# specify the model with 2-way interactions
FeedRM2b <- lmer(log(Feed.rate.per.fish)  ~ Temperature*Generation + Temperature*Log.ave.weight.of.population + Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = F) 

# test the models using anova
anova(FeedRM2a, FeedRM2b)
#Data: PG2F5Feedratedat
#Models:
#FeedRM2b: log(Feed.rate.per.fish) ~ Temperature * Generation + Temperature * 
#FeedRM2b:     Log.ave.weight.of.population + Generation * Log.ave.weight.of.population + 
#FeedRM2b:     (1 | Population) + (1 | Generation)
#FeedRM2a: log(Feed.rate.per.fish) ~ Temperature * Generation * Log.ave.weight.of.population + 
#FeedRM2a:     (1 | Population) + (1 | Generation)
#         npar     AIC    BIC logLik deviance  Chisq Df Pr(>Chisq)  
#FeedRM2b   10 -17.563 10.312 18.781  -37.563                       
#FeedRM2a   11 -19.379 11.284 20.689  -41.379 3.8161  1    0.05076 .
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# now test model fits with 2-way interactions of interest dropped 
# drop Temperature*Generation from the model
FeedRM3a <- lmer(log(Feed.rate.per.fish) ~  Temperature*Log.ave.weight.of.population + Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = F) 
anova(FeedRM2b, FeedRM3a)
#Data: PG2F5Feedratedat
#Models:
#FeedRM3a: log(Feed.rate.per.fish) ~ Temperature * Log.ave.weight.of.population + 
#FeedRM3a:     Generation * Log.ave.weight.of.population + (1 | Population) + 
#FeedRM3a:     (1 | Generation)
#FeedRM2b: log(Feed.rate.per.fish) ~ Temperature * Generation + Temperature * 
#FeedRM2b:     Log.ave.weight.of.population + Generation * Log.ave.weight.of.population + 
#FeedRM2b:     (1 | Population) + (1 | Generation)
#         npar      AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#FeedRM3a    9   8.9444 34.032  4.5278   -9.056                         
#FeedRM2b   10 -17.5626 10.312 18.7813  -37.563 28.507  1  9.336e-08 ***
#  ---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# drop Temperature*Weight from the model
FeedRM3b <- lmer(log(Feed.rate.per.fish) ~ Temperature*Generation  + Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = F) 
anova(FeedRM2b, FeedRM3b)
#Data: PG2F5Feedratedat
#Models:
#FeedRM3b: log(Feed.rate.per.fish) ~ Temperature * Generation + Generation * 
#FeedRM3b:     Log.ave.weight.of.population + (1 | Population) + (1 | Generation)
#FeedRM2b: log(Feed.rate.per.fish) ~ Temperature * Generation + Temperature * 
#FeedRM2b:     Log.ave.weight.of.population + Generation * Log.ave.weight.of.population + 
#FeedRM2b:     (1 | Population) + (1 | Generation)
#        npar     AIC     BIC logLik deviance  Chisq Df Pr(>Chisq)
#FeedRM3b    9 -17.412  7.6751 17.706  -35.412                     
#FeedRM2b   10 -17.563 10.3123 18.781  -37.563 2.1503  1     0.1425

# fit the supported model
FeedRMbest <- lmer(log(Feed.rate.per.fish) ~   Temperature*Generation + Generation*Log.ave.weight.of.population + (1|Population) + (1|Generation), data = PG2F5Feedratedat, REML = T) 

summary(FeedRMbest)

#Fixed effects:
#                                            Estimate Std. Error t value
#(Intercept)                                 -1.75618    0.05239 -33.519
#TemperatureH                                -0.10306    0.05274  -1.954
#GenerationGen6                              -0.23499    0.07991  -2.941
#Log.ave.weight.of.population                 0.47129    0.03013  15.642
#TemperatureH:GenerationGen6                  0.42017    0.07874   5.336
#GenerationGen6:Log.ave.weight.of.population  0.26461    0.05764   4.591

confint(FeedRMbest)  
#                                                 2.5 %       97.5 %
#.sig01                                       0.0000000  0.049961137
#.sig02                                       0.0000000  0.065249743
#.sigma                                       0.1849098  0.238256665
#(Intercept)                                 -1.8433229 -1.669049051
#TemperatureH                                -0.2046568 -0.001453932
#GenerationGen6                              -0.3710180 -0.098968194
#Log.ave.weight.of.population                 0.4132699  0.529313945
#TemperatureH:GenerationGen6                  0.2685404  0.571794918
#GenerationGen6:Log.ave.weight.of.population  0.1536148  0.375604336

plot(allEffects(FeedRMbest)) # gives good visual of parameter estimates
