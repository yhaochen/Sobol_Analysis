library(DEoptim)
wd<-getwd()
setwd(wd)

load('./annual_maxima_cms.RData')
load('./GEV_Parameters.RData')
load('./GEV_Parameters_MCMC.RData')

# define empirical probability searching function
median.auxiliary.func <- function(p, e, n){
        out <- abs((1-pbinom(e-1, n, p))-0.5)
        return(out)
}

# Numerical median probability return period formula
median.prob <- function(obs){
        l <- length(obs)
        # define variables
        e <- 1:l # the ranks of the events
        n <- l # sample size for the events
        pb <- txtProgressBar(min = 0, max = l, initial = 0, char = '=', style = 1) # loading bar
        prob <- vector(mode = 'numeric', length = l)
        for (i in 1:l) {
                setTxtProgressBar(pb, i) # loading bar
                fit <- DEoptim(median.auxiliary.func, lower = 0, upper = 1, e = e[i], n = n, control = DEoptim.control(trace = FALSE))
                prob[i] <- fit$optim$bestmem
        }
        close(pb)
        out <- sort(prob, decreasing = FALSE)
        return(out)
}

cdf_fun<- function(x,mu,sigma,xi) {
        if(xi!=0) t<-(1+xi*((x-mu)/sigma))^(-1/xi) else t<-exp(-(x-mu)/sigma)
        cdf_value<-exp(-t)
        if (is.nan(cdf_value)==T) cdf_value=0
        return(cdf_value)
}


Qs<-as.array(sort(annu_max_Q$peak_va))
# Observation points CDF
CDF_Qs<-median.prob(sort(Qs))

CDF_MLE<-cdf_fun(Qs,mu=GEV_params[1],sigma=GEV_params[2],xi=GEV_params[3])
CDF_Bayes<-cdf_fun(Qs,mu=mean(mu_chain),sigma=mean(sigma_chain),
                   xi=mean(xi_chain))

MLE_residuals<-CDF_MLE-CDF_Qs
Bayes_residuals<-CDF_Bayes-CDF_Qs

par(mfrow=c(2,2))
plot(Qs,MLE_residuals,main="MLE Residuals")
plot(Qs,Bayes_residuals,main="Bayes Residuals")
plot(Qs,CDF_Qs,pch=16,main="CDFs of MLE (red) and Bayes (blue)")
lines(Qs,CDF_MLE,col='red')
lines(Qs,CDF_Bayes,col='blue')
