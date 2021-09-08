wd<-getwd()
setwd(wd)
N=2000 #number of samples
load('./Q100_MCMC.RData')
load('./ensemble.RData')
Q100<-10^Q100_MCMC
Q100_sample<- 10^Q100_MCMC[seq(2,length(Q100_MCMC),by=floor(length(Q100_MCMC)/N))]

#plotting the pdf
pdf("Q100_PDF.pdf",width =6, height =4)
plot(density(Q100),col="blue",lwd=3,main="100-yr Discharge PDF for Susquehanna River at Sunbury",
     xlab="Discharge (cms)", xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(density(Q100_sample),col="red",lwd=3)
legend(45000,0.00018, legend=c("MCMC Results", "Sample"),
       col=c("blue", "red"), lty=1, cex=0.8)
dev.off()

#plotting the cdf
pdf("Q100_CDF.pdf",width =6, height =4)
plot(ecdf(Q100),col="blue",lwd=3,main="100-yr Discharge CDF for Susquehanna River at Sunbury",
     xlab="Discharge (cms)", xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(ecdf(Q100_sample),col="red",lwd=3)
legend(45000,0.8, legend=c("MCMC Results", "Sample"),
       col=c("blue", "red"), lty=1, cex=0.8)
dev.off()

#plotting the survival function
Q100_sort<-Q100[order(Q100)]
Q100_sample_sort<-Q100_sample[order(Q100_sample)]
cdf<-ecdf(Q100_sort)
cdf_samp=ecdf(Q100_sample_sort)
pdf("Q100_SURVIVAL.pdf",width =6, height =4)
plot(Q100_sort,1-cdf(Q100_sort),type="l",col="blue",lwd=3,cex.main=0.95,main="100-yr Discharge Survival Function for Susquehanna River at Sunbury",
     xlab="Discharge (cms)",ylab="Survival", xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(Q100_sample_sort, 1-cdf(Q100_sample_sort),col="red",lwd=3)
legend(45000,1, legend=c("MCMC Results", "Sample"),
       col=c("blue", "red"), lty=1, cex=0.8)
dev.off()
