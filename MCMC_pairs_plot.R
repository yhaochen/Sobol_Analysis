library(GGally)
wd<-getwd()
setwd(wd)

load('./GEV_Parameters_MCMC.RData')
myplot<-ggpairs(data.frame(mu=mu_chain,sigma=sigma_chain,xi=xi_chain))
ggsave('./MCMC_pairs.pdf')
rm(list=ls())
