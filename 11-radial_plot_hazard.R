wd<-getwd() #Gets the current directory 
setwd(wd)
set.seed(0)
mygreen <- rgb(44/255, 185/255, 95/255, 1) 
myblue <- rgb(0/255, 128/255, 1, 1)
myred <- rgb(1, 102/255, 102/255, 1)

# Loading libraries 
library(RColorBrewer)
library(graphics)
library(plotrix)

source('sobol_functions.R')

# input files that contain sobol indices
sobol_file_1 <- "radial_plot_table_1_haz.csv"
sobol_file_2 <- "radial_plot_table_2_haz.csv"

n_params <- 8 # set number of parameters
names=c('Discharge','Rriver Bed\nElevation','River\nWidth','Channel\nRoughness',
        'Floodplain\nRoughness','DEM\nResolution')
cols=c('darkgreen','darkgreen','darkgreen','darkgreen',
       'darkgreen','darkgreen')

## Import data from sensitivity analysis
# First- and total-order indices
s1st <- read.csv(sobol_file_1)[,-1]
parnames <- s1st[,1]

# Import second-order indices
s2_table <- read.csv(sobol_file_2)

# Convert second-order to upper-triangular matrix
s2 <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2)
s2 <- as.data.frame(s2)
colnames(s2) <- rownames(s2) <- s1st$Parameter

# Convert confidence intervals to upper-triangular matrix
s2_conf_low <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_high <- matrix(nrow=n_params, ncol=n_params, byrow=FALSE)
s2_conf_low[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_low)
s2_conf_high[1:(n_params-1), 2:n_params] = upper.diag(s2_table$S2_conf_high)

s2_conf_low <- as.data.frame(s2_conf_low)
s2_conf_high <- as.data.frame(s2_conf_high)
colnames(s2_conf_low) <- rownames(s2_conf_low) <- s1st$Parameter
colnames(s2_conf_high) <- rownames(s2_conf_high) <- s1st$Parameter

# Determine which indices are statistically significant
dummy<-read.csv('./dummy_haz.csv')
sig.cutoff_S1 <- dummy$high.ci[1] #
sig.cutoff_ST <- dummy$high.ci[2]

# S1 & ST: using the confidence intervals
s1st1<-s1st

for (i in 1:nrow(s1st)) {
s1st1$s1_sig[i]<-if(s1st1$S1_conf_low[i]-sig.cutoff_S1>=0) 1 else(0)
s1st1$st_sig[i]<-if(s1st1$ST_conf_low[i]-sig.cutoff_ST>=0) 1 else(0)
s1st1$sig[i]<-max(s1st1$s1_sig,s1st1$st_sig)
}

# S2: using the confidence intervals
s2_sig1 <- stat_sig_s2(s2,s2_conf_low,s2_conf_high,method='gtr',greater=0)


# Settings for the radial plot
cent_x=0
cent_y=0.2
radi=0.6
alph=360/(n_params)


pdf('radial_plot_hazard.pdf',width =3.94, height =3.94)

par(cex=0.5,mai=c(0.1,0.1,0.1,0.1))
plot(c(-1,1),c(-1,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",type="n")
draw.circle(0,.2,0.5,border = NA,col="gray90")

for(j in 1:(n_params)){
  i=j-1
  cosa=cospi(alph*i/180)
  sina=sinpi(alph*i/180)
  text(cent_x+cosa*(radi+radi*.15),cent_y+sina*(radi+radi*.15),names[j],srt=0,cex=1,col=cols[j])
  
  myX=cent_x+cosa*(radi-0.2*radi)
  myY=cent_y+sina*(radi-0.2*radi)
  for (z in j:n_params){ #Second-order interactions 
    if(s2_sig1[j,z]==1){
      g=z-1
      cosaa=cospi(alph*g/180)
      sinaa=sinpi(alph*g/180)
      EndX=cent_x+cosaa*(radi-0.2*radi)
      EndY=cent_y+sinaa*(radi-0.2*radi)
      lines(c(myX,EndX),c(myY,EndY),col='darkblue',
            lwd=qunif((s2[j,z]*s2_sig1[j,z]-min(s2*s2_sig1,na.rm=T))/(max(s2*s2_sig1,na.rm=T)-min(s2*s2_sig1,na.rm=T)),0.5,5))
    }
  }
  
  if(s1st1[j,9]>=1){ #Total-order nodes 
    draw.circle(cent_x+cosa*(radi-0.2*radi),cent_y+sina*(radi-0.2*radi),
                radius = qunif((s1st1[j,4]-min(s1st1[,4]))/(max(s1st1[,4])-min(s1st1[,4])),0.03,0.1),
                col="black")}
  
  if(s1st1[j,8]>=1){ #First-order nodes 
    draw.circle(cent_x+cosa*(radi-0.2*radi),cent_y+sina*(radi-0.2*radi),
                radius = qunif((s1st1[j,2]-min(s1st1[,2]))/(max(s1st1[,2])-min(s1st1[,2])),0.01,0.08),
                col=rgb(1, 102/255, 102/255,1),border = NA)}
}

# Plot the box below the plot 
x1=0.3
y1=0
draw.circle(x1+-0.9,y1+-0.97,0.08,border = NA,col=rgb(1, 102/255, 102/255,1))
draw.circle(x1+-0.7,y1+-0.97,0.01,border = NA,col=rgb(1, 102/255, 102/255,1))
text(x1+-0.9,y1+-0.83,paste(round(100*max(s1st1[s1st1[,"s1_sig"]>=1,2])),'%',sep=""))
text(x1+-0.7,y1+-0.83,paste(round(100*min(s1st1[s1st1[,"s1_sig"]>=1,2])),'%',sep=""))
text(x1+-0.8,y1+-0.75,'First-order')

draw.circle(x1+-0.4,y1+-0.97,0.1,col="black")
draw.circle(x1+-0.2,y1+-0.97,0.03,col="black")
text(x1+-0.4,y1+-0.83,paste(round(100*max(s1st1[s1st1[,"st_sig"]>=1,4])),'%',sep=""))
text(x1+-0.2,y1+-0.83,paste(round(100*min(s1st1[s1st1[,"st_sig"]>=1,4])),'%',sep=""))
text(x1+-0.3,y1+-0.75,'Total-order')

lines(c(x1+0.1,x1+0.2),c(y1+-0.97,y1+-0.97),lwd=5,col="darkblue")
#lines(c(x1+0.3,x1+0.4),c(y1+-0.97,y1+-0.97),lwd=0.5,col="darkblue")
text(x1+0.15,y1+-0.83,paste(round(100*max(s2[s2_sig1>=1],na.rm=T)),'%',sep=""))
#text(x1+0.35,y1+-0.83,paste(round(100*min(s2[s2_sig1>=1],na.rm=T)),'%',sep=""))
#text(x1+0.25,y1+-0.75,'Second-order')
text(x1+0.15,y1+-0.75,'Second-order')

#par(fig=c(0,1,0,1),new=T)
#plot(c(0,1),c(0,1),type="n",bty="n",xaxt="n",yaxt="n")
#text(0.9,0.9,'Earth sciences',cex=1.5,col="darkgreen", font=2)#,srt=-50)
#text(0.5,0.24,'Social sciences',cex=1.5,col="purple", font=2)#,srt=20)
#text(0.08,0.66,'Engineering',cex=1.5,col="darkred", font=2)#,srt=-90)

dev.off()


