# Set working directory
(wd <- getwd())
if (!is.null(wd))
  setwd(wd)

a=read.table("./Sobol-2_mostlikely_scenario.txt",header = T)
#Setting up radial plot tables for first and total order indices
col_names<-c("Parameter", "S1", "S1_conf_low", "S1_conf_high", "ST", "ST_conf_low", "ST_conf_high")
df=NULL
table<-read.csv('./ind_totalDamage.csv')
params <- c("Q", "Z", "W", "n_ch", "n_fp", "DEM","V","X") # name of parameters

S1<-table[table$sensitivity == "Si",]
ST<-table[table$sensitivity == "Ti",]

df<-data.frame(S1$parameters,S1$original,S1$low.ci,S1$high.ci,ST$original,ST$low.ci,ST$high.ci)
colnames(df)<-col_names
write.csv(df,'radial_plot_table_1_risk.csv')

#Setting up radial plot tables for second order indices
col_names<-c("Parameter_1", "Parameter_2",  "S2", "S2_conf_low",  "S2_conf_high")
S2<-table[table$sensitivity == "Sij",]
parameters= t(as.data.frame(strsplit(S2$parameters, split=".",fixed=T)))
df=data.frame(parameters,S2$original,S2$low.ci,S2$high.ci)
rownames(df)<-1:nrow(df)
colnames(df)<-col_names
write.csv(df,'radial_plot_table_2_risk.csv')






