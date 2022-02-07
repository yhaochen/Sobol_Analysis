library(DEoptim)
wd<-getwd()
setwd(wd)

C.I<-c(0.75,0.9,0.95,0.99)

load('./annual_maxima_cms.RData')
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

surp_fun<-function(C.I){
        up=0.5+C.I/2
        low<-0.5-C.I/2
        upper_limit<-sapply(1:length(Qs), function (x){max(1e-20,quantile(as.matrix(CDF[x,]),up))})
        lower_limit<-sapply(1:length(Qs), function (x){max(1e-20,quantile(as.matrix(CDF[x,]),low))})
        count=0
        for(i in 1:length(Qs)) if(CDF_Qs[i]>upper_limit[i]) count=count+1
        for(i in 1:length(Qs)) if(CDF_Qs[i]<lower_limit[i]) count=count+1
        surprise<-round(count/length(Qs),3)
        return(c(C.I,surprise))
}


Qs<-as.array(sort(annu_max_Q$peak_va))
# Observation points CDF
CDF_Qs<-median.prob(sort(Qs))

CDF<-data.frame()
for (k in 1:50000) {
        print(k)
        for (l in 1:length(Qs)) {
                CDF[l,k]<-cdf_fun(x=Qs[l],mu=mu_chain[k+50000],
                                  sigma=sigma_chain[k+50000],xi=xi_chain[k+50000])
        }
}
save(CDF,file='./annual_maxima_uncertainty_cdf.RData')
# load('./annual_maxima_uncertainty_cdf.RData')

result<-sapply(1:length(C.I), function (x){surp_fun(C.I[x])})
row.names(result)<-c("C.I.","Surprise")
save(result,file='./surprise_index.RData')

