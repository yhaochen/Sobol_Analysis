
wd<-getwd() #Gets the current directory 
setwd(wd)
load('./parameter_set.RData')
table<-read.csv('./table.csv')[,-(1:2)]
params<-c('Discharge (CMS)','River Bed Elevation Error (m)',
          'River Width Error','Channel Roughness',
          'Floodplain Roughness','DEM Resolution (m)',
          'Vulnerability Error','Exposure Error')

pdf('ensemble_pdfs.pdf',width =11, height =11)

par(cex=0.2,mai=c(0.25,0.25,0.25,0.25),mfrow=c(3,3))

for(i in 1:8) plot(density(table[,i]),main = params[i])

dev.off()

pdf('parameters_pdfs.pdf',width =11, height =8.5)

par(cex=0.2,mai=c(0.25,0.25,0.25,0.25),mfrow=c(3,3))

for(i in 1:8) hist(para[,i],main = params[i])

dev.off()

pdf('response_pdfs.pdf',width =8.5, height =11)

par(cex=0.2,mai=c(1,1,1,1),mfrow=c(2,1))

hist(table$n.damaged.houses,xlab = "Number of Damaged Houses",main = 'Total Hazard', xaxt="n", yaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
axis(side=2, at=axTicks(2), 
     labels=formatC(axTicks(2), format="d", big.mark=','))
hist(table$n.damaged.houses*table$mean.risk..USD.,xlab = "Damage (USD)",main = 'Total Risk', xaxt="n", yaxt="n")
axis(side=1, at=axTicks(1), 
     labels=formatC(axTicks(1), format="d", big.mark=','))
axis(side=2, at=axTicks(2), 
     labels=formatC(axTicks(2), format="d", big.mark=','))
dev.off()

