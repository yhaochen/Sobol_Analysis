#BASS method
rm(list = ls())
graphics.off()
library(BASS)
folder <- "/storage/work/h/hxy46/Sensitivity/Iman_test"
setwd(folder)

# Different seeds
seed_num <- 3
# Different sample size
size_num <- 6
selected_size <- c(1000,3000,5000,10000,20000,30000)

# MCMC parameter setting
mcmc_length <- 5000000
mcmc_thin <- 10000

#Load data
load("BASS_param_set.RData")

for (seed in 1:seed_num){
  set.seed(seed)
  for (size_index in 1:size_num){
    # Risk model and sensitivity
    Risk <- BASS_param_set[c(1:selected_size[size_index]),10]
    Param_Risk <- BASS_param_set[c(1:selected_size[size_index]),c(1:8)]
    BASS_model_Risk <- bass(Param_Risk,Risk,nmcmc=mcmc_length,thin=mcmc_thin)
    S_Risk<-sobol(BASS_model_Risk,verbose = TRUE)
    save(S_Risk,file = paste("sensitivity_Risk_",seed,"_",size_index,sep=""))
    
    # Hazard model and sensitivity
    Hazard<-BASS_param_set[c(1:selected_size[size_index]),9]
    Param_Hazard<-BASS_param_set[c(1:selected_size[size_index]),c(1:6)]
    BASS_model_Hazard<-bass(Param_Hazard,Hazard,nmcmc=mcmc_length,thin=mcmc_thin)
    S_Hazard<-sobol(BASS_model_Hazard,verbose = TRUE)
    save(S_Hazard,file = paste("sensitivity_Hazard_",seed,"_",size_index,sep=""))
  }
}




