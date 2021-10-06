wd<-getwd()
setwd(wd)

# Load the libraries required in this script  
library(evdbayes)
library(ismev)

# Load Annual Maximum Data  
#load("./floods.RData")
load("./annual_maxima_cfs.RData")
load('./GEV_Parameters.RData')
# Estimate parameters 
set.seed(1)

cov.mat <- diag(c(1000,100,1))
pn <- prior.norm(mean = c(0,0,0), cov = cov.mat)
pos<-posterior(50000, init = GEV_params, prior = pn, lh = "gev",data = log10(annu_max_Q[,2]*0.3048^3),psd = c(0.2,.1,.1))
mu=(pos[,1])
logsigma=(pos[,2])
xi=(pos[,3])

# Choose the last 10,000 iterations
#mu_chain <-mu[(length(xi)-10000+1):length(xi)]
#xi_chain <- xi[(length(xi)-10000+1):length(xi)]
#sigma_chain <- logsigma[(length(logsigma)-10000+1):length(logsigma)]
mu_chain<-mu
xi_chain<-xi
sigma_chain<-logsigma

# Save MCMC chains 
save(mu_chain,xi_chain,sigma_chain,file="GEV_Parameters_MCMC.RData")
