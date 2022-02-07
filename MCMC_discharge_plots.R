library(extRemes)
load('./GEV_Parameters.RData')
load('./GEV_Parameters_MCMC.RData')
load('./annual_maxima_cfs.RData')
load('./MCMC_Discharge_CMS.RData') #MCMC discharge(50k)
load('./MCMC_Q_samp_random.RData')
load('./PMF.RData')

n<-1e5
MLE<-10^revd(n, loc = GEV_params[1], scale = GEV_params[2],
             shape = GEV_params[3], type = c("GEV"))
mle_order<-MLE[order(MLE)]
cdf_mle<-ecdf(mle_order)
annu_Q_order<-annu_max_Q[order(annu_max_Q[,2]),2]*0.3048^3
cdf_Q<-ecdf(annu_Q_order)

mcmc_50k_order<-discharge_df[order(discharge_df[,1]),]
q_samp_order<-Q_samp[order(Q_samp)]

cdf_50k<-ecdf(mcmc_50k_order)
cdf_q_samp<-ecdf(q_samp_order)

#pdf('.pdf',width =8.5, height =11)
#par(cex.lab=1.5,cex.axis=1.5,mai=c(0.6,0.6,0.5,0.5), mfrow=c(3,1))

# plot histogram and best fit GEV (max likelihood) pdf of annual maxima of discharge


pdf('ann_max_hist.pdf',width =11, height =8.5)
hist(annu_max_Q$peak_va*0.3048^3,freq = F,
     cex.lab=1.5,cex.axis=1,
     main="",xlab="Discharge (cms)",ylab='', xaxt="n",yaxt='n',col='white')
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(density(MLE),lwd=3,col="red")
legend(12000,0.00020, legend=c("Best Fit GEV", "Normalized Histogram"),
       col=c("red", "black"), pch=c(NA_integer_, 22),
       lty=c(1,0),lwd=c(3,0),  cex=1.5)
dev.off()

# add MCMC samples on the previous plot
n_samp<-10
df=matrix(0,nrow=n,ncol=n_samp)
set.seed(666)
rows<-sample(seq(1:length(mu_chain)),n_samp,replace=F)
for (i in 1:n_samp) {
  df[(1:n),i]<-10^revd(n, loc = mu_chain[rows[i]], scale = sigma_chain[rows[i]],
                       shape = xi_chain[rows[i]], type = c("GEV"))
}

pdf('ann_max_hist_MCMC.pdf',width =11, height =8.5)
hist(annu_max_Q$peak_va*0.3048^3,freq = F,
     cex.lab=1.5,cex.axis=1,ylim=c(0,0.00025),
     main="",xlab="Discharge (cms)",ylab='', xaxt="n",yaxt="n",col='white')
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
for (i in 1:n_samp) {
 lines(density(df[,i]),lwd=0,col="lightgrey")
}
lines(density(MLE),lwd=3,col="red")

legend(12000,0.00020, legend=c("MCMC Samples", "Best Fit GEV", "Normalized Histogram"),
       col=c("lightgrey", "red", "black"), pch=c(NA_integer_,NA_integer_, 22),
       lty=c(1,1,0),lwd=c(0,3,0),  cex=1.5)
dev.off()

# showing the survival function of the previous plot


pdf('ann_max_survival.pdf',width =11, height =8.5)
plot(mle_order,1-cdf_mle(mle_order),lwd=3,col="red",type='l',
     cex.lab=1.5,cex.axis=1,xlim=c(0,20000),
     main="",xlab="Discharge (cms)",ylab='Survival', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
for (i in 1:n_samp) {
  data<-df[,i]
  data_order<-data[order(data)]
  cdf<-ecdf(data_order)
  lines(data_order,1-cdf(data_order),type='l',lwd=0,col="lightgrey")
}
lines(mle_order,1-cdf_mle(mle_order),lwd=3,col="red",type='l')
points(annu_Q_order,1-cdf_Q(annu_Q_order),col="black",pch=16)

legend(13000,1, legend=c("MCMC Samples", "Best Fit GEV","Observations"),
       col=c("lightgrey", "red","black"), pch=c(NA_integer_,NA_integer_, 16),
       lty=c(1,1,0),lwd=c(0,3,0),  cex=1.5)
dev.off()

# showing the 95% CI boundaries on the survival function

pdf('ann_max_survival_95CI.pdf',width =11, height =8.5)
plot(mle_order,1-cdf_mle(mle_order),lwd=3,col="red",type='l',
     cex.lab=1.5,cex.axis=1,xlim=c(0,20000),
     main="",xlab="Discharge (cms)",ylab='Survival', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
polygon(x = c(mle_order,rev(mle_order)), 
        y = c((1.1-cdf_mle(mle_order)),rev(0.9-cdf_mle(mle_order))), border = NA , col = "lightgreen")
lines(mle_order,1-cdf_mle(mle_order),lwd=3,col="red",type='l')
points(annu_Q_order,1-cdf_Q(annu_Q_order),col="black",pch=16)
legend(13000,1, legend=c("95% CI", "Best Fit GEV","Observations"),
       col=c("lightgreen", "red","black"), pch=c(15,NA_integer_, 16),
       lty=c(0,1,0),lwd=c(3,3,0),  cex=1.5)

dev.off()

# showing the survival functions of MCMC discharge (50k), and 2k samples


pdf('mcmc_comparison_survival.pdf',width =11, height =8.5)
plot(mcmc_50k_order,1-cdf_50k(mcmc_50k_order),lwd=3,col="red",type='l',
     cex.lab=1.5,cex.axis=1,
     main="",xlab=expression("Discharge (m"^3*"/s)"),ylab='Survival (1-CDF)', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(q_samp_order,1-cdf_q_samp(q_samp_order),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)

legend(10000,1, legend=c("Full MCMC Results", "MCMC Sample","Probable Maximum Flood"),
       col=c("red", "blue","black"), lty=c(1,2,1), lwd=3, cex=1.5)
dev.off()

# same plot with log survival
pdf('mcmc_comparison_survival_log.pdf',width =11, height =8.5)
plot(mcmc_50k_order,1-cdf_50k(mcmc_50k_order),lwd=3,col="red",type='l',
     cex.lab=1.5,cex.axis=1,log="y",
     main="",xlab=expression("Discharge (m"^3*"/s)"),ylab='Survival (1-CDF)', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(q_samp_order,1-cdf_q_samp(q_samp_order),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)
points(annu_Q_order,1-cdf_Q(annu_Q_order),col="black",pch=16)

legend(20000,1, legend=c("Full MCMC Results", "MCMC Sample","Probable Maximum Flood","Observed"),
       col=c("red", "blue","black","black"), lty=c(1,2,1,NA_integer_),
       pch=c(NA_integer_,NA_integer_,NA_integer_,16), lwd=3, cex=1.5)
dev.off()

# showing the survival functions of MCMC discharge (50k), and 2k samples
pdf('mcmc_comparison_PDFs.pdf',width =11, height =8.5)
plot(density(discharge_df[,1]),lwd=3,col="red",type='l',
     cex.lab=1.5,cex.axis=1,
     main="",xlab=expression("Discharge (m"^3*"/s)"), xaxt="n")
axis(side=1, at=axTicks(1),
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(density(Q_samp),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)

legend(12000,0.0002, legend=c("Full MCMC Results","MCMC Sample","Probable Maximum Flood"),
       col=c("red", "blue","black"), lty=c(1,2,1), lwd=3, cex=1.5)
dev.off()


################################
### Combined plots

# combined pdf and survival

pdf('mcmc_comparison_two_plots.pdf',width =11, height =11)
par(cex.axis=1.5,cex.lab=1.5, mai=c(1,1,1,1), mfrow=c(2,1))
plot(density(discharge_df[,1]),lwd=3,col="red",type='l',
    main="",xlab=expression("Discharge (m"^3*"/s)"), xaxt="n")
axis(side=1, at=axTicks(1),
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(density(Q_samp),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)

legend(20000,0.0002, legend=c("Full MCMC Results","MCMC Sample","Probable Maximum Flood"),
       col=c("red", "blue","black"), lty=c(1,2,1), lwd=3, cex=1.5)

plot(mcmc_50k_order,1-cdf_50k(mcmc_50k_order),lwd=3,col="red",type='l',
     main="",xlab=expression("Discharge (m"^3*"/s)"),ylab='Survival (1-CDF)', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(q_samp_order,1-cdf_q_samp(q_samp_order),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)
points(annu_Q_order,1-cdf_Q(annu_Q_order),col="black",pch=16)

legend(20000,1, legend=c("Full MCMC Results", "MCMC Sample","Probable Maximum Flood","Observed"),
       col=c("red", "blue","black","black"), lty=c(1,2,1,NA_integer_),
       pch=c(NA_integer_,NA_integer_,NA_integer_,16), lwd=3, cex=1.5)
dev.off()

# combined pdf and log - survival

pdf('mcmc_comparison_two_plots_log.pdf',width =11, height =11)
par(cex.axis=1.5,cex.lab=1.5, mai=c(1,1,1,1), mfrow=c(2,1))
plot(density(discharge_df[,1]),lwd=3,col="red",type='l',
     main="",xlab=expression("Discharge (m"^3*"/s)"), xaxt="n")
axis(side=1, at=axTicks(1),
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(density(Q_samp),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)

legend(25000,0.0002, legend=c("Full MCMC Results","MCMC Sample","Probable Maximum Flood"),
       col=c("red", "blue","black"), lty=c(1,2,1), lwd=3, cex=1.5)

plot(mcmc_50k_order,1-cdf_50k(mcmc_50k_order),lwd=3,col="red",type='l',
     log="y",
     main="",xlab=expression("Discharge (m"^3*"/s)"),ylab='Survival (1-CDF)', xaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
lines(q_samp_order,1-cdf_q_samp(q_samp_order),lwd=3,col="blue",type='l',lty=2)
abline(v=PMF,lwd=5)
points(annu_Q_order,1-cdf_Q(annu_Q_order),col="black",pch=16)

legend(25000,1, legend=c("Full MCMC Results", "MCMC Sample","Probable Maximum Flood","Observed"),
       col=c("red", "blue","black","black"), lty=c(1,2,1,NA_integer_),
       pch=c(NA_integer_,NA_integer_,NA_integer_,16), lwd=3, cex=1.5)
dev.off()
