#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library("sensobol")
library("data.table")
library("ggplot2")
library("extRemes")
library("TruncatedDistributions")
library("raster")
library("foreach")
library("doParallel")

# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)



params <- c("Q", "Z", "W", "n_ch", "n_fp", "DEM","V","X") # name of parameters
N<-2000 # number of samples
k <- length(params) # number of parameters
R <- 1e4 # number of bootstrap samples
order<-"second"
type <- "norm" # method to compute the confidence intervals of sobol indices
conf <- 0.95 # confidence interval
matrices=c("A", "B", "AB")
#ensemble of parameters between 0 and 1
set.seed(1)
mat <- sobol_matrices(N = N, params = params,order=order,type="LHS"
                      ,matrices = matrices)
#                      ,matrices = c("A", "B", "AB", "BA"))
save(mat,file = './ensemble_uniform.RData')

load('./GEV_Parameters.RData')
load('./Q100_MCMC.RData')

# from Maggie: function map between [0,1] and a bounded parameter range
map_range <- function(x, bdin, bdout) {
  bdout[1] + (bdout[2] - bdout[1]) * ((x - bdin[1]) / (bdin[2] - bdin[1]))
}

# parameter set
set.seed(1)
#mat[, "Q"] <- qevd(mat[, "Q"], loc=GEV_params[1],
#                          scale=GEV_params[2],shape=GEV_params[3],type='GEV')
      # systematic sampling from Q_100 MCMC data
Q_samp <- 10^Q100_MCMC[seq(2,length(Q100_MCMC),by=floor(length(Q100_MCMC)/N))]
#mat[, "Q"] <- 10^sample(Q100_MCMC,N)
Qrange<-range(Q_samp)
mat[, "Q"] <- map_range(Q_samp, c(Qrange),c(0,1))

set.seed(2)
mat[, "Z"] <- qtbeta(mat[, "Z"], alpha=5, beta=5, a=0, b=1)

set.seed(3)
mat[, "W"] <- qtbeta(mat[, "W"], alpha = 5, beta = 5, a=0, b=1)

set.seed(4)
mat[, "n_ch"] <- qtnorm(mat[, "n_ch"], mean=(0.03-0.02)/(0.2-0.02), sd=0.5,a=0, b=1)

set.seed(5)
mat[, "n_fp"] <- qtnorm(mat[, "n_fp"], mean=(0.12-0.02)/(0.4-0.02), sd=0.5,a=0, b=1)

set.seed(7)
mat[, "V"] <- qtbeta(mat[, "V"], alpha = 5, beta = 5, a=0, b=1)

set.seed(8)
mat[, "X"] <- qtbeta(mat[, "X"], alpha = 5, beta = 5, a=0, b=1)

save(mat,file = './ensemble.RData')

#parameter sets for model run
para<-mat
para[,'Q']<-Q_samp
para[,'Z']<-map_range(para[,'Z'],c(0,1),c(-5,+5))
para[,'W']<-map_range(para[,'W'],c(0,1),c(-0.1,+0.1))
para[,'n_ch']<-map_range(para[,'n_ch'],c(0,1),c(0.02,0.2))
para[,'n_fp']<-map_range(para[,'n_fp'],c(0,1),c(0.02,0.4))
10->para[, "DEM"][para[, "DEM"]<1/3]
30->para[, "DEM"][para[, "DEM"]<2/3]
50->para[, "DEM"][para[, "DEM"]<3/3]
para[,'V']<-map_range(para[,'V'],c(0,1),c(-0.4,+0.4))
para[,'X']<-map_range(para[,'X'],c(0,1),c(-0.4,+0.4))

save(para,file = './parameter_set.RData')

####################################################################################
####################################################################################
# model setup and run

#parameter sets and functions
load('parameter_set.RData')
source('./haz_risk_function.R')

# Read the depth-damage table and hypothetical houses unit price estimate raster
vulner_table <-
  read.csv("./Inputs/vulnerability.csv") # depth-damage table for Selinsgrove for a one story house without basement from (Stuart A. Davis and L. Leigh Skaggs 1992)
house_price <- raster("./Inputs/house_price.asc")

# Read LISFLOOD-FP river and parameter files
sample_river <-
  as.matrix(read.delim2("./LISFLOOD_8/Sample_Selinsgrove.river", header = F)) # sample Susquehanna River file
sample_par <-
  as.matrix(read.delim2("./LISFLOOD_8/Sample_Selinsgrove.par", header = F)) # sample parameter file

# Output folders
if (dir.exists(paste0(wd, "/Outputs")) == F)
  dir.create(paste0(wd, "/Outputs"))
if (dir.exists(paste0(wd, "/Outputs/Extent")) == F)
  dir.create(paste0(wd, "/Outputs/Extent"))
if (dir.exists(paste0(wd, "/Outputs/Hazard")) == F)
  dir.create(paste0(wd, "/Outputs/Hazard"))
if (dir.exists(paste0(wd, "/Outputs/Risk")) == F)
  dir.create(paste0(wd, "/Outputs/Risk"))
if (dir.exists(paste0(wd, "/Outputs/Summary")) == F)
  dir.create(paste0(wd, "/Outputs/Summary"))


run_start = 1 #starting row number of the parameters table to read
run_end = nrow(para) #ending row number of the parameters table to read


#setwd(paste0(wd,'/LISFLOOD_8'))
#system("./chmod a+x lisflood.exe") # activate the .exe file for the first run

#setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) # -1 not to overload system
registerDoParallel(cl)

start<-proc.time()

#foreach (i = run_start:run_end) %dopar% {
#  print(i)
#  haz_risk_run(i)
#  }
for (i in run_start:run_end) haz_risk_run(i)
end<-proc.time()

print(end-start)
stopCluster(cl)

# model run end
####################################################################################
####################################################################################
#Setting up the model response table
setwd(wd)
load('./ensemble.RData')
list <- list.files('./Outputs/Summary', pattern = '.csv')

table <- data.frame()


for (i in seq_along(list)) {
    data = read.csv(paste0('./Outputs/Summary/', list[i]))
  table = rbind(table, data)
}
table <- table[order(table$run.no.),][,-1]
write.csv(table,"table.csv")


# calculation of sobol indices for risk
table<-read.csv("table.csv")
y_risk=table$mean.risk..USD.*table$n.damaged.houses


plot_uncertainty(Y = y_risk, N = N) + labs(y = "Counts", x = "$y$")
plot_scatter(data = mat, N = N, Y = y_risk, params = params)
plot_multiscatter(data = mat, N = N, Y = y_risk, params = params)


ind <- sobol_indices(Y = y_risk, N = N, params = params, boot = T, R = R,
                     type = type, conf = conf,order = order,first = "saltelli")
cols <- colnames(ind$results)[1:5]

ind$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
ind
write.csv(ind$results,"ind_totalDamage.csv")


ind.dummy <- sobol_dummy(Y = y_risk, N = N, params = params, boot = TRUE, R = R)
write.csv(ind.dummy,'./dummy_risk.csv')
plot(ind, dummy = ind.dummy,order = "first")
plot(ind, dummy = ind.dummy,order = "second")




sub.sample <- seq(100, N, 10) # Define sub-samples

convergence<-sobol_convergence(
  matrices,
  Y=y_risk,
  N,
  sub.sample,
  params,
  first="saltelli",
  total = "jansen",
  order = order,
  seed = 666,
  plot.order=order
)
write.csv(convergence[1],"converge_totalDamage.csv")
convergence[2]
convergence[3]

pdf("converge_first_risk.pdf",width =11, height =8.5)
convergence[2]
dev.off()

pdf("converge_second_risk.pdf",width =11, height =8.5)
convergence[3]
dev.off()

######################################################
# calculation of sobol indices for hazard
#params <- c("Q", "Z", "W", "n_ch", "n_fp", "DEM") # name of parameters
#k <- length(params) # number of parameters
y_haz=table$mean.hazard..m.*table$n.damaged.houses
ind <- sobol_indices(Y = y_haz, N = N, params = params, boot = T, R = R,
                     type = type, conf = conf,order = order,first = "saltelli")
cols <- colnames(ind$results)[1:5]

ind$results[, (cols):= round(.SD, 3), .SDcols = (cols)]
ind
write.csv(ind$results,"ind_totalHazard.csv")


ind.dummy <- sobol_dummy(Y = y_haz, N = N, params = params, boot = TRUE, R = R)
write.csv(ind.dummy,'./dummy_haz.csv')
plot(ind, dummy = ind.dummy,order = "first")
plot(ind, dummy = ind.dummy,order = "second")

sub.sample <- seq(100, N, 10) # Define sub-samples

convergence<-sobol_convergence(
  matrices,
  Y=y_haz,
  N,
  sub.sample,
  params,
  first="saltelli",
  total = "jansen",
  order = order,
  seed = 666,
  plot.order="second",
)
write.csv(convergence[1],"converge_totalHazard.csv")
convergence[2]
convergence[3]
pdf("converge_first_haz.pdf",width =11, height =8.5)
convergence[2]
dev.off()

pdf("converge_second_haz.pdf",width =11, height =8.5)
convergence[3]
dev.off()
