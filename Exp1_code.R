knitr::opts_chunk$set(echo = TRUE)

library(tidyverse) #for all data wrangling
library(cowplot) #for manuscript ready figures
library(lme4) #for lmer & glmer models
library(sjPlot) #for plotting lmer and glmer mods
library(sjmisc) 
library(effects)
library(sjstats) #use for r2 functions
library(lmerTest)
library("scales")
library(ggplot2)

data <- read_csv("S1.csv")

data$Instr<- factor(data$Instr)
is.factor(data$Instr)

data$ID_n<- factor(data$ID_n)
is.factor(data$ID_n)

data$Task1Seek<- factor(data$Task1Seek)
is.factor(data$Task1Seek)


data$r_RealMC <- rescale(data$RealMC) 
data$r_uncertainty <- rescale(data$uncertainty) 

#control parameters to use in all models
control_params = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))


#ONLY share
mod <- lmer(seek/share~ Instr +r_RealMC +r_uncertainty+
              (Instr+r_RealMC+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==2), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ r_uncertainty+
              (r_uncertainty|ID_n),
            data=subset(data,Task1Seek==2), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ Instr +r_RealMC+
              (Instr+r_RealMC|ID_n),
            data=subset(data,Task1Seek==2), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ Instr +r_uncertainty+
              (Instr+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==2), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ r_RealMC +r_uncertainty+
              (r_RealMC+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==2), REML=FALSE, control = control_params)
summary(mod)

#ONLY  SEEK
mod <- lmer(seek/share~ Instr +r_RealMC +r_uncertainty+
              (Instr+r_RealMC+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)


mod <- lmer(seek/share~ r_RealMC+
              (r_RealMC|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)


mod <- lmer(seek/share~ Instr+
              (Instr|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)



mod <- lmer(seek/share~ r_uncertainty+
              (r_uncertainty|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ Instr +r_RealMC+
              (Instr+r_RealMC|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ Instr +r_uncertainty+
              (Instr+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)

mod <- lmer(seek/share~ r_RealMC +r_uncertainty+
              (r_RealMC+r_uncertainty|ID_n),
            data=subset(data,Task1Seek==1), REML=FALSE, control = control_params)
summary(mod)



mod <- lmer(seek~ (Instr +r_Market_change +r_uncertainty)*Task1Seek+
              (Instr +r_Market_change +r_uncertainty|ID_n),
            data=data, REML=FALSE, control = control_params)

summary(mod)


#ONLY SHARE OR SEEK
mod <- lmer(seek~ Instr +r_Market_change +r_uncertainty+
              (Instr +r_Market_change +r_uncertainty|ID_n),
            data=subset(data,Task1Seek==0), REML=FALSE, control = control_params)

summary(mod)
