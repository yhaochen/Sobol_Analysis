library(dataRetrieval)
library(extRemes)

# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)

# Download annual maximum instantaneous river flow data of the USGS gauge Susquehanna River at Sunbury
annu_max_Q <- readNWISpeak(
  "01554000",
  startDate = "",
  endDate = "",
  asDateTime = TRUE,
  convertType = TRUE
)
annu_max_Q<-data.frame(annu_max_Q$peak_dt,annu_max_Q$peak_va)[-(1:3),]
colnames(annu_max_Q)<-c("dates","peak_va")
save(annu_max_Q,file = 'annual_maxima_cfs.RData')

# Pick above flood stage flow data
#flood_stage=24 # flood stage in feet based on NOAA (https://water.weather.gov/ahps2/hydrograph.php?gage=sbyp1&wfo=ctp)
#rating<-readNWISrating("01554000")
#min_flood<-approx(rating$INDEP,rating$DEP,flood_stage)
#min_flood=min_flood$y #minimum flood in cfs equivalent of flood stage
#floods<-annu_max_Q[annu_max_Q$peak_va>=min_flood,]
#floods$peak_va<-floods$peak_va*0.3048^3 # conversion to cms
#floods$peak_va<-floods$peak_va/max(floods$peak_va)
#save(floods,file = 'floods.RData')

fit<-fevd(log10(annu_max_Q[,2]*0.3048^3),type = 'GEV', method = 'MLE')

# GEV location, shape, and scale parameters
location<-fit$results$par[1]
scale<-fit$results$par[2]
shape<-fit$results$par[3]

# Print the results
GEV_params=c(location,scale,shape)

#plot(density(floods$peak_va),col='red')
#points(floods$peak_va,devd(floods$peak_va,loc=location,scale=scale,shape=shape,type = 'GEV'))



# Save the parameters
save(GEV_params,file="GEV_Parameters.RData")

jpeg(file="GEV_plot.jpeg",pointsize = 24,width = 1024, height = 1024)
plot(fit)
dev.off()
rm(list=ls())
