##==============================================================================
##
## Script to compare return levels with and without uncertainty quantification
##
## Authors: Mahkameh Zarekarizi (mahkameh.zare@gmail.com)
##==============================================================================
## Copyright 2019 Mahkameh Zarekarizi
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
##    folder (Zarekarizi-Home-Elevation)
##    To do so:
##      1. If on RStudio, open the README.md file. Then on the menu bar, go to
##         Session-->Set Working Directory-->To Source File Location
##      2. If on RStudio, on the lower right box, open "Zarekarizi-Home-Elevation"
##         Then, click on More --> Set as Working Directory
##      3. On console, type the following command:
##         setwd("~/.../../../../Zarekarizi-Home-Elevation")
## 2. To Run:
##      1. Click on Source button or, on console, type: Source("../../....R")
## 3. Outputs:
##      1. output includes a plot that compares return levels with and without uncertainty
##==============================================================================

# Global variables
wd<-getwd()
setwd(wd)
myblue <- rgb(0.68, 0.85, 0.90,0.5)
if("fExtremes" %in% (.packages())){
  detach("package:fExtremes", unload=TRUE)
}
if("evd" %in% (.packages())){
  detach("package:evd", unload=TRUE)
}
library('DEoptim')
#--------------------------------------------------------------
# Functions----------------------------------------------------

## function to estimate the return level from GEV distribution
myreturnlevel <- function(t,mu,sigma,xi){
  library(evir)
  x<-qgev(p=1-(1/t),xi=xi,mu=mu,sigma=sigma)
  return(x)
}

# define empirical probability searching function
median.auxiliary.func <- function(p, e, n){
  out <- abs((1-pbinom(e-1, n, p))-0.5)
  return(out)
}

# Numerical median probability return period formula
median.rt <- function(obs){
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
  out <- sort(1/prob, decreasing = FALSE)
  return(out)
}
#---------------------------------------------------------------------
# Load libraries and data required to run this code
load("./GEV_Parameters_MCMC.RData")
load("./annual_maxima_cms.RData")
load('./GEV_Max_a_Posteriori.RData')
load('./GEV_Posterior_Mean.RData')

library(evir)
plot_rps <- c(seq(1,2,0.1),seq(3,9,1),seq(10,90,10),seq(100,500,100))

# Find return levels for each parameter set
MC_rl <- sapply(1:nrow(mcmcSamples), function(x) {
  myreturnlevel(plot_rps, mu=mcmcSamples[x,1], sigma=mcmcSamples[x,2], xi=mcmcSamples[x,3])
})
rownames(MC_rl)<-plot_rps

# Find upper and lower limits for 90% CI bounds
lower_5 <- sapply(1:length(plot_rps), function (x){quantile(MC_rl[x,],0.05)})
upper_95 <- sapply(1:length(plot_rps), function (x) {quantile(MC_rl[x,],0.95)})

# We need a second panel to show the density at return level of 500
rl_500 <- MC_rl[32,]
rl_500_lb=quantile(rl_500,0.05)
rl_500_ub=quantile(rl_500,0.95)
rl_500[rl_500<rl_500_lb]=NA
rl_500[rl_500>rl_500_ub]=NA
h500<-hist(rl_500,25,plot = T,probability = T)

# maximum a posteriori return levels
MC_rl_MAP<-myreturnlevel(plot_rps,GEV_est_MAP[1],GEV_est_MAP[2],GEV_est_MAP[3])
# posterior mean return levels
MC_rl_mean<-myreturnlevel(plot_rps,bayesEstimator[1,1],bayesEstimator[1,2],bayesEstimator[1,3])
########### PLOT ##############
pdf("Return_Period_90percent_Uncertainty_Plot.pdf",width =3.94, height =2.43)

# plot high-level variables
par(cex=0.5,mai=c(0.4,0.4,0.2,0.15)) #c(bottom, left, top, right)
par(cex=0.5,fig=c(0,0.7,0.05,1))

ymin=0
ymax=50000
xmin=1
xmax=500

# The base plot
plot(plot_rps,MC_rl_mean, log="x", xlim = c(xmin,xmax),type="n",bty="n",xaxt="n",xaxs="i",yaxs="i",
     ylim = c(ymin,ymax),yaxt="n",xlab = "",ylab="")

# Axes 
axis(1,pos=ymin, at=c(1,10,100,1000,500),cex.axis=0.8,lwd=0.5)
axis(2,pos=xmin, at = c(seq(ymin,ymax,by=20000)),labels=formatC(seq(ymin,ymax,by=20000), format="d", big.mark=','),cex.axis=0.8,lwd=0.5)

# x and y axis labels 
mtext("Return period (years)",side=1,line=2.5,cex=0.5)
mtext(expression("Discharge (m"^3*"/s)"),side=2,line=2.5,cex=0.5)

# Box around the plot 
lines(x=c(xmin,xmin),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymax,ymax))
lines(x=c(xmax,xmax),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymin,ymin))


# Uncertainty boundaries
polygon(x = c(plot_rps[2:length(plot_rps)],rev(plot_rps[2:length(plot_rps)])), 
        y = c(upper_95[2:length(plot_rps)], rev(lower_5[2:length(plot_rps)])), border = NA , col = myblue)

# With and without uncertainty lines
lines(plot_rps[2:length(plot_rps)],MC_rl_mean[2:length(MC_rl_mean)],lty=1,col="blue")
lines(plot_rps[2:length(plot_rps)],MC_rl_MAP[2:length(MC_rl_MAP)],lty=1,col="red")

# Observation points
points(median.rt(sort(annu_max_Q[,2])), sort(annu_max_Q[,2]), lwd = 1, cex = 0.75, pch = 20, bg = "white")

# Legend
legend(2,ymax-ymax*0.01,
       c("Maximum a posteriori","Posterior mean",
         "90% C.I. uncertainty area","Observed annual maxima of discharge"),
       col = c('red',"blue", myblue,'black'),
       pt.bg = c(NA, NA, myblue,"white"),
       pch = c(NA, NA, 22,20),
       lty = c(1, 1,NA,NA),
       lwd = c(1, 1, NA,1.5),
       bty = 'n',
       pt.cex = c(NA, NA, 2,1),
       inset = c(0.01, -0.01),
       cex=0.9)
text(xmin+xmin*0.4,ymax-ymax*0.05,"a)",cex=1.5)

## Second panel 
par(cex=0.5,fig=c(0.62,1,0.05,1),new=TRUE)
xmin=0
xmax=max(h500$density)
xmax=xmax+0.1*xmax
ymean=MC_rl_mean[32]
yMAP=MC_rl_MAP[32]

plot(h500$density,h500$mids,ylim=c(ymin,ymax),xlim=c(xmin,xmax),xlab="",ylab="",xaxs="i",yaxs="i",bty="n",xaxt="n",yaxt="n",type="n")
axis(1,pos=ymin,at=c(xmin,max(h500$density)),labels=c(0,signif(max(h500$density),0)),cex.axis=0.8,lwd=0.5)
mtext("\nProbability density of\n500-year return period",side=1,line=3,cex=0.5)

# box around the plot 
lines(x=c(xmin,xmin),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymax,ymax))
lines(x=c(xmax,xmax),y=c(ymin,ymax))
lines(x=c(xmin,xmax),y=c(ymin,ymin))

# main line
mids<-h500$mids
 polygon(x=c(smooth(h500$density[1:length(mids)],"3R"),rep(0,length(mids))),y=c(mids,rev(mids)),col=myblue,border = NA)
#polygon(x=c(h500$density[1:length(mids)],rep(0,length(mids))),y=c(mids,rev(mids)),col=myblue,border = NA)
lines(x=c(xmin,xmax),y=c(ymean,ymean),col="blue")
lines(x=c(xmin,xmax),y=c(yMAP,yMAP),col="red")
# Panel indicator 
text(xmax*0.1,ymax-ymax*0.05,"b)",cex=1.5)

dev.off()

