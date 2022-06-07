# This script plots the radial plots for the sensitivity analysis.
# This script is modified based on "12-radial_plot_risk.R"
rm(list = ls())
graphics.off()
library(plotrix)

folder <- "/storage/work/h/hxy46/Sensitivity/Iman_test"
setwd(folder)

# selected sample sizes
selected_size<-c(1000,3000,5000,10000,20000,30000)

# Parameters and names
n_params <- 8
names <- c("Q","Z","W","n_ch","n_fp","DEM","V","X")
cols=c('darkgreen','darkgreen','darkgreen','darkgreen',
       'darkgreen','darkgreen','darkred','purple')

# A threshold that determines which indices are statistical significant
thres <- 0.01

# Settings for the radial plot
cent_x=0
cent_y=0.2
radi=0.6
alph=360/(n_params)

# A function to create an upper diagnal matrix for 2nd-order indices
upper.diag <- function(x){
  m<-(-1+sqrt(1+8*length(x)))/2+1
  X<-lower.tri(matrix(NA,m,m),diag=TRUE)
  X[X==FALSE]<-x
  X[X==1] <- 0
  X
}

# Generate a radial plot for each seed and each sample size
for (seed_num in 1:3){
  for (size_num in 1:6){
    load(paste(folder,"/sensitivity_Risk_",seed_num,"_",size_num,sep=""))
    
    # Each index is the mean of the last 100 saved MCMC iterations
    # Total-order indices:
    st <- colMeans(S_Risk$T[c(400:499), ])
    
    # First-order indices:
    s1 <- colMeans(S_Risk$S[c(400:499),c(1:n_params)])
    
    # second-order indices:
    s2 <- S_Risk$S[c(400:499),c((n_params+1):(n_params+n_params*(n_params-1)/2))]
    s2 <- upper.diag(colMeans(s2))
    
    pdfname <- paste("radial_plot_Risk_",seed_num,"_",size_num,".pdf",sep="")
    pdf(pdfname,width =3.94, height =3.94)
    
    par(cex=0.5,mai=c(0.1,0.1,0.1,0.1))
    plot(c(-1,1),c(-1,1),bty="n",xlab="",ylab="",xaxt="n",yaxt="n",type="n")
    draw.circle(0,.2,0.5,border = NA,col="gray90")
    
    for(j in 1:(n_params)){
      i=j-1
      cosa=cospi(alph*i/180)
      sina=sinpi(alph*i/180)
      text(cent_x+cosa*(radi+radi*.15),cent_y+sina*(radi+radi*.15),
           names[j],srt=0,cex=1,col=cols[j])
      myX=cent_x+cosa*(radi-0.2*radi)
      myY=cent_y+sina*(radi-0.2*radi)
      if (j < n_params){
        for (z in (j+1):n_params){ #Second-order interactions 
          if (s2[j,z]>=thres){
            g <- z-1
            cosaa=cospi(alph*g/180)
            sinaa=sinpi(alph*g/180)
            EndX=cent_x+cosaa*(radi-0.2*radi)
            EndY=cent_y+sinaa*(radi-0.2*radi)
            lines(c(myX,EndX),c(myY,EndY),col='darkblue',
                  lwd=qunif(s2[j,z]/max(s2),0,5))
          }
        }
      }
      
      if(st[j]>=thres){ #Total-order nodes 
        draw.circle(cent_x+cosa*(radi-0.2*radi),cent_y+sina*(radi-0.2*radi),
                    radius = qunif(st[j]/max(st),0.03,0.1),
                    col="black")
      }
      
      if(s1[j]>=thres){ #First-order nodes 
        draw.circle(cent_x+cosa*(radi-0.2*radi),cent_y+sina*(radi-0.2*radi),
                    radius = qunif(s1[j]/max(s1),0.01,0.08),
                    col=rgb(1, 102/255, 102/255,1),border = NA)
      }
    }
    
    # Plot the box below the plot 
    x1=0.3
    y1=0
    draw.circle(x1+-0.9,y1+-0.97,0.08,border = NA,col=rgb(1, 102/255, 102/255,1))
    draw.circle(x1+-0.7,y1+-0.97,0.01,border = NA,col=rgb(1, 102/255, 102/255,1))
    text(x1+-0.9,y1+-0.83,paste(round(100*max(s1)),'%',sep=""))
    text(x1+-0.7,y1+-0.83,paste(round(100*max(min(s1),thres)),'%',sep=""))
    text(x1+-0.8,y1+-0.75,'First-order')
    
    draw.circle(x1+-0.4,y1+-0.97,0.1,col="black")
    draw.circle(x1+-0.2,y1+-0.97,0.03,col="black")
    text(x1+-0.4,y1+-0.83,paste(round(100*max(st)),'%',sep=""))
    text(x1+-0.2,y1+-0.83,paste(round(100*max(min(s1),thres)),'%',sep=""))
    text(x1+-0.3,y1+-0.75,'Total-order')
    
    lines(c(x1+0.1,x1+0.2),c(y1+-0.97,y1+-0.97),lwd=5,col="darkblue")
    text(x1+0.15,y1+-0.83,paste(round(100*max(s2,na.rm=T)),'%',sep=""))
    text(x1+0.15,y1+-0.75,'Second-order')
    
    #title of used samples and seed
    x2=0
    y2=1
    text(x2,y2,paste('Sample size: ',selected_size[size_num], 
                     "  seed ",seed_num,sep=""))
    dev.off()
  }
}