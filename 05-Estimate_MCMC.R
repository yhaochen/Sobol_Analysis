##==============================================================================
##
## Script estimates GEV distribution parameters using MCMC 
##
## Authors: Iman Hosseini-Shakib (ishakib@gmail.com)
##          Klaus Keller (kzk10@psu.edu)
##
##  Modified from a code by Mahkameh Zarekarizi available at:
## https://github.com/scrim-network/Zarekarizi-flood-home-elavate/blob/master/Source_Code/S10_Estimate_MCMC.R
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
## 3. Outputs:
##      1. output includes a data file MCMC_Discharge_CMS.RData
##==============================================================================

wd<-getwd()
setwd(wd)

# Load the libraries required in this script  
library(evdbayes)
library(ismev)
library(extRemes)

# Load Annual Maximum Data  
load("./annual_maxima_cms.RData")
load('./GEV_Parameters.RData')
# Estimate parameters 
set.seed(1)

cov.mat <- diag(c(1000,100,1))
pn <- prior.norm(mean = c(0,0,0), cov = cov.mat)
pos<-posterior(50000, init = GEV_params, prior = pn,
               lh = "gev",data = annu_max_Q[,2],
               psd = c(0.2,.1,.1))
mu_chain<-(pos[,1])
sigma_chain<-(pos[,2])
xi_chain<-(pos[,3])

# Save MCMC chains 
save(mu_chain,xi_chain,sigma_chain,file="GEV_Parameters_MCMC.RData")

# produce discharge values from MCMC results
#load('./GEV_Parameters_MCMC.RData')
load('./PMF.RData') # probable maximum flood (10,000-yr flood)
n=1 # number of random discharge samples from each MCMC parameter set
discharge_df<-as.data.frame(NULL)
set.seed(1)
for (i in 1:length(mu_chain)) {
  discharge_df[i,1]<-revd(n,loc=mu_chain[i],scale=sigma_chain[i],shape=xi_chain[i])
  print(i)
}
colnames(discharge_df)<-c('discharge_cms')
discharge_df<-discharge_df[discharge_df$discharge_cms<=PMF,]
discharge_df<-data.frame(discharge_df)
save(discharge_df,file="MCMC_Discharge_CMS.RData")
rm(list = ls())
