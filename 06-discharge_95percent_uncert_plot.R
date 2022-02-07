##==============================================================================
##
## Script to plot upstream discharge pdfs and survival functions for the
## sampled and observed data as well as the MLE and MCMC results.
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
##==============================================================================

# Global variables
wd<-getwd()
setwd(wd)
set.seed(0)
if("fExtremes" %in% (.packages())){
  detach("package:fExtremes", unload=TRUE)
}
if("evd" %in% (.packages())){
  detach("package:evd", unload=TRUE)
}

library(evir)

load("./GEV_Parameters.RData") # Maximum likelihood GEV parameter set
load("./GEV_Parameters_MCMC.RData") # MCMC parameter sets
load("./annual_maxima_cms.RData") # Annual maxima of discharge in cfs
load('./MCMC_Q_samp_random.RData')
# load('./PMF.RData')
load('./Q500.RData')

cdf_fun<- function(x,mu,sigma,xi) {
  if(xi!=0) t<-(1+xi*((x-mu)/sigma))^(-1/xi) else t<-exp(-(x-mu)/sigma)
  cdf_value<-exp(-t)
  if (is.nan(cdf_value)==T) cdf_value=0
  return(cdf_value)
}

plot_Qs <- c(seq(1000,ceiling(Q500/1000)*1000,500)) # discharge values of axis x 

#Annual maxima of discharge
histogram<-hist(annu_max_Q$peak_va, freq=F) # observed maxima histogram
histogram$density<-histogram$density/(sum(histogram$density)*2000)
annu_Q_order<-annu_max_Q[order(annu_max_Q[,2]),2]
cdf_annu_Q<-ecdf(annu_Q_order)
surv_annu_Q<-1-cdf_annu_Q(annu_Q_order)

#Maximum Likelihood GEV 
MLE_Q<-dgev(plot_Qs[1]:Q500,mu=GEV_params[1],sigma=GEV_params[2],
            xi=GEV_params[3]) # discharge densities based on the maximum likelihood GEV
surv_mle<- sapply(plot_Qs[1]:Q500, function (x){1-cdf_fun(x,
                  mu=GEV_params[1],sigma=GEV_params[2],xi=GEV_params[3])})

#Bayesian maximum a posteriori GEV 
Bayes_Q<-dgev(plot_Qs[1]:Q500,mu=mean(mu_chain),sigma=mean(sigma_chain),
            xi=mean(xi_chain))
surv_bayes<- sapply(plot_Qs[1]:Q500, function (x){1-cdf_fun(x,
                    mu=mean(mu_chain),sigma=mean(sigma_chain),xi=mean(xi_chain))})


#Sampled discharge
sample_Q<-sort(Q_samp) # 2000 samples from MCMC results
pdf_sample_Q<-density(sample_Q)
cdf_sample_Q<-ecdf(sample_Q)
surv_sample_Q<-1-cdf_sample_Q(sample_Q)

#MCMC Discharge
uncert_Q<-data.frame() # data frame of densities of discharge values based on MCMC GEVs
for (i in 1:length(plot_Qs)) { 
  print(plot_Qs[i])
  for (j in 1:50000) {
    uncert_Q[i,j]<-dgev(plot_Qs[i],xi=xi_chain[j+50000],
                        mu=mu_chain[j+50000],sigma=sigma_chain[j+50000])
    if(is.nan(uncert_Q[i,j])==T) uncert_Q[i,j]=1e-10
  }
}
rownames(uncert_Q)<-plot_Qs
save(uncert_Q,file='./discharge_uncertainty_pdf.RData')
  #load('./discharge_uncertainty_pdf.RData')

# Matrix of MCMC CDF values
CDF<-data.frame()
for (k in 1:ncol(uncert_Q)) {
  # for (k in 1:1000) {
    print(k)
    for (l in 1:nrow(uncert_Q)) {
    CDF[l,k]<-cdf_fun(x=plot_Qs[l],mu=mu_chain[k+50000],
                      sigma=sigma_chain[k+50000],xi=xi_chain[k+50000])
      }
  }
rownames(CDF)<-plot_Qs
save(CDF,file='./discharge_CDF.RData')
#load('./discharge_CDF.RData')

 # Find upper and lower density limits of MCMC discharge data
lower_2.5 <- sapply(1:length(plot_Qs), function (x){max(1e-20,quantile(uncert_Q[x,],0.025)$'2.5%')})
upper_97.5 <- sapply(1:length(plot_Qs), function (x) {quantile(uncert_Q[x,],0.975)$'97.5%'})
surv_lower_2.5 <- sapply(1:length(plot_Qs), function (x){max(1e-20,(1-quantile(CDF[x,],0.975)$'97.5%'))})
surv_upper_97.5 <- sapply(1:length(plot_Qs), function (x) {1-quantile(CDF[x,],0.025)$'2.5%'})

pdf("Discharge_95percent_Uncertainty_PDF.pdf",width =4.86,height =7.88)
################################################
### Panel A
################################################

# plot pdfs
par(mfrow=c(2,1))
par(cex=0.5,mai=c(0.5,0.4,0.1,0.3)) # mai   c(bottom, left, top, right)

ymin=0
ymax=1.1*max(upper_97.5)
xmin=0
xmax=ceiling(Q500/1000)*1000

# The base plot
plot(histogram, freq=F,xlim = c(xmin,xmax),xaxt="n",xaxs="i",yaxs="i",
     ylim = c(ymin,ymax),yaxt="n",main="",xlab = "",ylab="")

# Axes 
axis(1,pos=ymin, at=seq(0,xmax,5000),labels=formatC(seq(0,xmax,5000), format="d", big.mark=','),cex.axis=0.8,lwd=0.5)
axis(2, pos=xmin,at = c(0,round(max(upper_97.5),5)),labels=c(0,"0.0004"),lwd=0.5,cex.axis=0.8)

# x and y axis labels 
mtext(expression("Discharge (m"^3*"/s)"),side=1,line=2.5,cex=0.5)
mtext("Density",side=2,line=2.5,cex=0.5)

# Box around the plot 
lines(x=c(xmin,xmax),y=c(ymin,ymin))
lines(x=c(xmin,xmin),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymax,ymax))
lines(x=c(xmax,xmax),y=c(ymin,ymax))

# Uncertainty boundaries
polygon(x = c(plot_Qs,rev(plot_Qs)), 
        y = c(upper_97.5, rev(lower_2.5)),
        border = NA , col=rgb(0.68, 0.85, 0.90,0.5))

# MLE PDF and samples PDF

lines(plot_Qs[1]:Q500,MLE_Q,lty=1,col="red",lwd=2)
lines(plot_Qs[1]:Q500,Bayes_Q,lty=1,col="blue",lwd=2)
lines(pdf_sample_Q$x[-(1:19)],pdf_sample_Q$y[-(1:19)] ,lty=2,col="green",lwd=2)
# Legend
legend(15000,ymax-ymax*0.01,
       c("Maximum Likelihood GEV","Maximum a Posteriori GEV","Sampled Discharge","95% C.I. Uncertainty Range","Observed Annual Maxima of Discharge"),
       col = c('red','blue','green',rgb(0.68, 0.85, 0.90,0.5),'black'),
       pt.bg = c(NA,NA,NA, rgb(0.68, 0.85, 0.90,0.5),"white"),
       pch = c(NA,NA,NA, 22,22),
       lty = c(1,1,2,NA,NA),
       lwd = c(2,2,1,NA,0.5),
       bty = 'n',
       pt.cex = c(NA,NA,NA,2,2),
       cex=0.9)
# Panel indicator 
text(xmin+1000,ymax-ymax*0.05,"a)",cex=1.5)

################################################
### Panel B
################################################
# plot survival functions
#pdf("Discharge_Uncertainty_PDF1.pdf",width =3.94,height =2.43)

par(cex=0.5,mai=c(0.5,0.4,0.1,0.3)) # mai   c(bottom, left, top, right)

ymin=1e-4
ymax=2
xmin=0
xmax=ceiling(Q500/1000)*1000

# The base plot
plot(MLE_Q,surv_mle, xlim = c(xmin,xmax),pch=NA,xaxt="n",xaxs="i",yaxs="i",
     ylim = c(ymin,ymax),yaxt="n",main="",xlab = "",ylab="",log = "y")

# Axes 
axis(1,pos=ymin, at=seq(0,xmax,5000),labels=formatC(seq(0,xmax,5000), format="d", big.mark=','),cex.axis=0.8,lwd=0.5)
axis(2, pos=xmin,at = c(ymin,1),lwd=0.5,cex.axis=0.8)

# x and y axis labels 
mtext(expression("Discharge (m"^3*"/s)"),side=1,line=2.5,cex=0.5)
mtext("Survival (1-CDF)",side=2,line=2.5,cex=0.5)

# Box around the plot 
lines(x=c(xmin,xmin),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymax,ymax))
lines(x=c(xmax,xmax),y=c(ymin,ymax))

# Uncertainty boundaries
polygon(x = c(plot_Qs,rev(plot_Qs)), 
        y = c(surv_upper_97.5,
              rev(surv_lower_2.5)),
        border = NA , col=rgb(0.68, 0.85, 0.90,0.5))

# MLE PDF
lines(plot_Qs[1]:Q500,surv_mle,lty=1,col="red",lwd=2)
lines(plot_Qs[1]:Q500,surv_bayes,lty=1,col="blue",lwd=2)
lines(sample_Q,surv_sample_Q,lty=2,col="green",lwd=2)
points(annu_Q_order,surv_annu_Q)

# Legend
legend(15000,1,
       c("Maximum Likelihood GEV","Maximum a Posteriori GEV","Sampled Discharge","95% C.I. Uncertainty Range","Observed Annual Maxima of Discharge"),
       col = c('red','blue','green',rgb(0.68, 0.85, 0.90,0.5),'black'),
       pt.bg = c(NA,NA,NA, rgb(0.68, 0.85, 0.90,0.5),"white"),
       pch = c(NA,NA,NA, 22,1),
       lty = c(1,1,2,NA,NA),
       lwd = c(2,2,1,NA,0.5),
       bty = 'n',
       pt.cex = c(NA,NA,NA,2,2),
       cex=0.9)

text(xmin+1000,0.5,"b)",cex=1.5)

dev.off()

