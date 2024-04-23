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

data1 <- read_csv("Share2.csv")

data1$Cl_Seek<- factor(data1$Cl_Seek)
is.factor(data1$Cl_Seek)

library(rgl)
library(car)
library(mgcv)
library(nlme)
library(plotly)
library(plot3D)
library(scatterD3)
library(car)
library(dplyr)

data1$Cl_Seek <- as.character(data1$Cl_Seek)

# Replace numbers with desired label
data1$Cl_Seek[data1$Cl_Seek == "1"] <- "."
data1$Cl_Seek[data1$Cl_Seek == "2"] <- "`"
data1$Cl_Seek[data1$Cl_Seek == "3"] <- ","

# 3D plot
scatter3d(x = data1$SeekUnc, y = data1$SeekMC, z = data1$SeekInstr, groups = factor(data1$Cl_Seek),
          ellipsoid = TRUE, grid = FALSE, surface = FALSE, axis.ticks = TRUE, axis.scales = FALSE, level = 0.5,
          text.col = c("white", "White", "White"), axis.col = c("Black", "Black", "Black"),
          surface.col = c('RED1','#0066ff', 'springgreen4'),
          xlab = "", ylab = "", zlab = "")

# Adjust the view angle (you can change these values)
par3d(windowRect = c(10, 10, 1200, 1200))
par3d(userMatrix = rotate3d(par3d("userMatrix"), pi/3.3, 1, 0, 0))
par3d(userMatrix = rotate3d(par3d("userMatrix"), -pi/4, 0, 1, 0))


# Save a snapshot with higher resolution
rgl.snapshot("Share2.png", fmt = "png")


