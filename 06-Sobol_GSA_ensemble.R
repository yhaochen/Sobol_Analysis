##==============================================================================
##
## Script creates the ensemble and the parameter sets required for GSA
##
## Authors: Iman Hosseini-Shakib (ishakib@gmail.com)
##          Klaus Keller (kzk10@psu.edu)
##
##==============================================================================
## Copyright 2022 Iman Hosseini-Shakib
## This file is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This file is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this file.  If not, see <http://www.gnu.org/licenses/>.
##==============================================================================
## Instructions to run:
## 1. If you have not already done so, change the working directory to the main 
##    folder 
##    To do so:
##      1. If on RStudio, open the README.md file. Then on the menu bar, go to 
##         Session-->Set Working Directory-->To Source File Location
##      2. If on RStudio, on the lower right box, open the main folder
##         Then, click on More --> Set as Working Directory
##      3. On console, type the following command: 
##         setwd("~/../MAIN_FOLDER") 
## 2. To Run:
##      1. Click on Source button or, on console, type: Source("../../....R")
##==============================================================================
#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")

library("sensobol")
library("TruncatedDistributions")

# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)

params <- c("Q", "Z", "W", "n_ch", "n_fp", "DEM","V","X") # name of parameters
N<-2000 # number of samples
k <- length(params) # number of parameters
R <- 1e4 # number of bootstrap samples
order<-"second"
type <- "norm" # method to compute the confidence intervals of Sobol indices
conf <- 0.95 # confidence interval
matrices=c("A", "B", "AB")
save(params,N,k,R,order,type,conf,matrices,file="Sobol_setup.RData")
#ensemble of parameters between 0 and 1
set.seed(1)
mat <- sobol_matrices(N = N, params = params,order=order,type="LHS"
                      ,matrices = matrices)
save(mat,file = './ensemble_uniform.RData')

# loading samples of discharge data from MCMC posterior
load('Q_sample_A.RData')
load('./Q_sample_B.RData')

# finding the pattern of matrices A and B for the discharge values
a=mat[1,1] # first element of matrix A
b=mat[(N+1),1] # first element of matrix B
index_a<-(which(mat[,1]==a)-1)/N+1 # where matrix A is used
index_b<-(which(mat[,1]==b)-1)/N+1 # where matrix B is used
  #pattern
pattern<-rep(NA,nrow(mat)/N)
pattern[index_a]<-'A'
pattern[index_b]<-'B'
  # replacing matrices A and B values based on the pattern
Q<-rep(NA,nrow(mat))
for(i in 1:length(pattern)){
  if(pattern[i]=="A") {Q[((i-1)*N+1):(i*N)]=Q_samp_A}
  else {Q[((i-1)*N+1):(i*N)]=Q_samp_B}
} 

# function map between [0,1] and a bounded parameter range
map_range <- function(x, bdin, bdout) {
  bdout[1] + (bdout[2] - bdout[1]) * ((x - bdin[1]) / (bdin[2] - bdin[1]))
}

# ensemble
mat[, "Q"]<-map_range(Q,range(Q),c(0,1))

set.seed(2)
mat[, "Z"] <- qtbeta(mat[, "Z"], alpha=5, beta=5, a=0, b=1)

set.seed(3)
mat[, "W"] <- qtbeta(mat[, "W"], alpha = 5, beta = 5, a=0, b=1)

set.seed(4)
mat[, "n_ch"] <- qtnorm(mat[, "n_ch"], mean=(0.03-0.02)/(0.1-0.02), sd=0.5,a=0, b=1)

set.seed(5)
mat[, "n_fp"] <- qtnorm(mat[, "n_fp"], mean=(0.12-0.02)/(0.2-0.02), sd=0.5,a=0, b=1)

set.seed(7)
mat[, "V"] <- qtbeta(mat[, "V"], alpha = 5, beta = 5, a=0, b=1)

set.seed(8)
mat[, "X"] <- qtbeta(mat[, "X"], alpha = 5, beta = 5, a=0, b=1)

save(mat,file = './ensemble.RData')

#parameter sets for model run
para<-mat
para[,'Q']<-Q
para[,'Z']<-map_range(para[,'Z'],c(0,1),c(-5,+5))
para[,'W']<-map_range(para[,'W'],c(0,1),c(-0.1,+0.1))
para[,'n_ch']<-map_range(para[,'n_ch'],c(0,1),c(0.02,0.1))
para[,'n_fp']<-map_range(para[,'n_fp'],c(0,1),c(0.02,0.2))
10->para[, "DEM"][para[, "DEM"]<1/3]
30->para[, "DEM"][para[, "DEM"]<2/3]
50->para[, "DEM"][para[, "DEM"]<3/3]
para[,'V']<-map_range(para[,'V'],c(0,1),c(-0.4,+0.4))
para[,'X']<-map_range(para[,'X'],c(0,1),c(-0.4,+0.4))

save(para,file = './parameter_set.RData')
rm(list = ls())
