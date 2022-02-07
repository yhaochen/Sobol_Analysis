rm(list=ls())
wd<-getwd()
setwd(wd)
load("./annual_maxima_cms.RData")
library(extRemes)

####################################################################################################
####################################################################################################
dat<-annu_max_Q[,2]
####################################################################################################
####################################################################################################
# Maximum Likelihood Approach
####################################################################################################
####################################################################################################
# The following code fits the GEV model using the fevd function from the extRemes package
fit<-fevd(dat,type = 'GEV', method = 'MLE', initial = list(c(0,1,0.5)))
# GEV location, shape, and scale parameters
location<-fit$results$par[1]
scale<-fit$results$par[2]
shape<-fit$results$par[3]
GEV_params=c(location,scale,shape)
GEV_params
GEV_mle<-ci(fit, alpha = 0.05, type = c("parameter")) # From the fExtremes Packages
# location        scale        shape 
# 5056.1595463 1619.6022021    0.1927981 

# #estimate PMF from 10,000 flood
PMF<-qevd(p=1-1e-4, loc = GEV_params[1], scale = GEV_params[2],
          shape = GEV_params[3], type = c("GEV"))

####################################################################################################
# Alternative Method: Maximum Likelihood without the extRemes Package
####################################################################################################
# Define the likelihood function (use the likelihood function from extRemes Package)
llhd<-function(par,dat){
  if(par[2]<=0){return(-1e10)}else{
    loglik<-sum(devd(x=dat, loc = par[1], scale = par[2], shape = par[3], log = TRUE,
             type = c("GEV")))  
    if(loglik==-Inf){return(-1e10)}else{return(loglik)}
  }
}

par.init<-c(0,1,2)
llhd(par=par.init , dat=dat)

optimOutput<-optim(par = par.init ,fn = llhd, dat=dat,  control=list(fnscale=-1),hessian=TRUE,
                   method = "L-BFGS-B" , lower = c(-Inf,0,-Inf), upper = c(Inf,Inf,Inf))
GEV_est<-optimOutput$par

# Use the Delta method to construct the 95% confidence intervals
fisher_info<-solve(-optimOutput$hessian)
prop_sigma<-sqrt(diag(fisher_info))
upper<-GEV_est+1.959964*prop_sigma
lower<-GEV_est-1.959964*prop_sigma
interval<-data.frame(value=GEV_est, upper=upper, lower=lower)
round(interval,3)


#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Bayesian Approach
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
# Log Likelihood Function
llhd<-function(par,dat){
  sum(devd(x=dat, loc = par[1], scale = exp(par[2]), shape = par[3], log = TRUE,
       type = c("GEV")))
}

# Log Prior Density Function
lprior<-function(par){
  dnorm(x=par[1],mean=0,sd=sqrt(1e+10), log = TRUE) + dnorm(x=par[2],mean=0,sd=sqrt(100), log = TRUE) + dnorm(x=par[3],mean=0,sd=sqrt(1), log = TRUE)
}
# Log Posterior Density Function
lposterior<-function(par,dat){
  llhd(par,dat)+lprior(par)
}

# Bayesian Approach
library(adaptMCMC)
accept.mcmc = 0.234										# Optimal acceptance rate as # parameters->infinity
#	(Gelman et al, 1996; Roberts et al, 1997)
niter.mcmc = 200000										# number of iterations for MCMC
gamma.mcmc = 0.7										# rate of adaptation (between 0.5 and 1, lower is faster adaptation)
burnin = round(niter.mcmc*0.5)				# how much to remove for burn-in
dat=annu_max_Q[,2]
# Run MCMC algorithm
par.init<-c(GEV_params[1],log(GEV_params[2]),GEV_params[3])

amcmc.out = MCMC(p=lposterior, n=niter.mcmc, init=par.init, adapt=TRUE, 
                 acc.rate=accept.mcmc,
                 gamma=gamma.mcmc, list=TRUE, n.start=round(0.01*niter.mcmc), dat=dat)

amcmc.out$acceptance.rate
mcmcSamples<-amcmc.out$samples[-(1:burnin),]
mcmcSamples[,2]<-exp(mcmcSamples[,2])

#####################################
# Highest Posterior Density Function
## Using Ming-Hui Chen's paper in Journal of Computational and Graphical Stats.
hpd <- function(samp,p=0.05){
  ## to find an approximate (1-p)*100% HPD interval from a
  ## given posterior sample vector samp
  
  r <- length(samp)
  samp <- sort(samp)
  rang <- matrix(0,nrow=trunc(p*r),ncol=3)
  dimnames(rang) <- list(NULL,c("low","high","range"))
  for (i in 1:trunc(p*r)) {
    rang[i,1] <- samp[i]
    rang[i,2] <- samp[i+(1-p)*r]
    rang[i,3] <- rang[i,2]-rang[i,1]
  }
  hpd <- rang[order(rang[,3])[1],1:2]
  return(hpd)
}

######################################
bayesEstimator<-apply(mcmcSamples,2,function(x){c(mean(x), hpd(x))})
rownames(bayesEstimator)<-c("Posterior Mean" , "95%CI-Low", "95%CI-High")
# Results from Bayesian Approach (MCMC)
bayesEstimator
# Results from MLE Approach
GEV_mle

# Comparative plots for the results
library(evd)
xSeq<-seq(from=0,to=20000, length.out=100000)
yDensity<-dgev(x = xSeq , loc = GEV_params[1] , scale = GEV_params[2] , shape = GEV_params[3])
yDensityBayes<-dgev(x = xSeq , loc = bayesEstimator[1,1] , scale = bayesEstimator[1,2] , shape = bayesEstimator[1,3])
par(mfrow=c(2,2))
plot(x=annu_max_Q[,1], y=annu_max_Q[,2], pch=16, main = "Observations")
plot(density(annu_max_Q[,2]) , main= "Observation Density")
plot(x=xSeq , y=yDensity, typ="l",col="red" , main="Comparison of Densities: MLE (red) vs. Bayes (black)")
lines(x=xSeq , y=yDensityBayes)

# Trace plots for MCMC
par(mfrow=c(2,3), mar=c(2,2,2,2))
plot.ts(mcmcSamples[,1], main="location", ylim=range(mcmcSamples[,1], GEV_params[1])) ; abline(h=GEV_params[1], col="red")
plot.ts(mcmcSamples[,2], main="scale", ylim=range(mcmcSamples[,2], GEV_params[2])); abline(h=GEV_params[2], col="red")
plot.ts(mcmcSamples[,3], main="shape", ylim=range(mcmcSamples[,3], GEV_params[3])); abline(h=GEV_params[3], col="red")
plot(density(mcmcSamples[,1]), main="location"); abline(v=GEV_params[1], col="red")
plot(density(mcmcSamples[,2]), main="scale"); abline(v=GEV_params[2], col="red")
plot(density(mcmcSamples[,3]), main="shape"); abline(v=GEV_params[3], col="red")
