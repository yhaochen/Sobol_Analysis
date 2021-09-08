
wd<-getwd() #Gets the current directory 
setwd(wd)
load('./parameter_set.RData')
table<-read.csv('./table.csv')[,-(1:2)]
params<-c('Discharge','River Bed Elevation','River Width','Channel Roughness',
          'Floodplain Roughness','DEM Resolution','Vulnerability','Exposure')

pdf('ensemble_pdfs.pdf',width =8.5, height =11)

par(cex=0.2,mai=c(0.25,0.25,0.25,0.25),mfrow=c(3,3))

for(i in 1:8) plot(density(table[,i]),main = params[i])

dev.off()

pdf('parameters_pdfs.pdf',width =8.5, height =11)

par(cex=0.2,mai=c(0.25,0.25,0.25,0.25),mfrow=c(3,3))

for(i in 1:8) plot(density(para[,i]),main = params[i])

dev.off()

pdf('response_pdfs.pdf',width =8.5, height =11)

par(cex=0.2,mai=c(1,1,1,1),mfrow=c(2,1))

plot(density(table$n.damaged.houses*table$mean.hazard..m.),xlab = "",main = 'Total Hazard (m)')
plot(density(table$n.damaged.houses*table$mean.risk..USD.),xlab = "",main = 'Total Risk (USD)')

dev.off()
