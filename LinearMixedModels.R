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

data <- read_csv("yourdata.csv")

data$Instr<- factor(data$Instr)
is.factor(data$Instr)


data$r_mc_model <- rescale(data$mc_model) 
data$r_uncertainty <- rescale(data$uncertainty) 

#control parameters to use in all models
control_params = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000))


#ONLY SHARE OR SEEK
mod <- lmer(seek~ Instr +r_mc_model +r_uncertainty+
              (Instr +r_mc_model +r_uncertainty|ID),
            data=data,REML=FALSE, control = control_params)

