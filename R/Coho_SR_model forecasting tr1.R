###########################################################################
###     Hierarchical SR models w/time variant alphas, PR priors        ####
###########################################################################

#setwd("C:/data/centralcoast")
require(rjags)
require(R2jags)
library(runjags)
require(gsl)

sim_eval <- 20 # forecast model for 16 years, or 4 generations
thinning <- 200
samples <- 1500
n_chains <- 4
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)

ER_values <- read.csv("Data/ER_values.csv", header=T)

SR.dat$er_est <- SR.dat$er_2
SR.dat$er_est[is.na(SR.dat$er_est)] <- SR.dat$er_E[is.na(SR.dat$er_est)]

n.pops<-max(SR.dat$pop_no)
n.years<-length(seq(1980,2016,1))
spawn<-matrix(nrow=n.years,ncol=n.pops,NA)

rec1<-matrix(nrow=n.years,ncol=n.pops,NA)
rec2<-matrix(nrow=n.years,ncol=n.pops,NA)
rec3<-matrix(nrow=n.years,ncol=n.pops,NA)

logRS1<-matrix(nrow=n.years,ncol=n.pops,NA)
logRS2<-matrix(nrow=n.years,ncol=n.pops,NA)
logRS3<-matrix(nrow=n.years,ncol=n.pops,NA)
logRSo<-matrix(nrow=n.years,ncol=n.pops,NA)
er_est <- matrix(nrow=n.years,ncol=n.pops,NA)

for(i in 1:n.pops){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  spawn[,i]<-dat[,"escapement"]
  rec1[,i]<-dat[,"rec_E"]
  rec2[,i]<-dat[,"rec_2"]
  rec3[,i]<-dat[,"rec_3"]
  logRS1[,i]<-dat[,"logRS1"]
  logRS2[,i]<-dat[,"logRS2"]
  logRS3[,i]<-dat[,"logRS3"]
  logRSo[,i]<-ifelse(!is.na(dat[,"logRS2"]),dat[,"logRS2"],dat[,"logRS1"])
  er_est[,i]<-dat$er_est
}


yrs<-matrix(nrow=n.years,ncol=n.pops,NA) 

for(i in 1:n.pops){
  yrs[1:length(which(spawn[,i]>0)),i]<-which(spawn[,i]>0)
}

#### making matrixes JAGS to index
spawn.mat<-matrix(nrow=n.years,ncol=n.pops,NA)
logRS1.mat<-matrix(nrow=n.years,ncol=n.pops,NA)
logRS2.mat<-matrix(nrow=n.years,ncol=n.pops,NA)
logRS3.mat<-matrix(nrow=n.years,ncol=n.pops,NA)
logRSB.mat<-matrix(nrow=n.years,ncol=n.pops,NA)


for(i in 1:n.pops){
  spawn.mat[1:length(na.omit(yrs[,i])),i]<-spawn[na.omit(yrs[,i]),i]
}

for(i in 1:n.pops){
  logRS1.mat[1:length(na.omit(yrs[,i])),i]<-logRS1[na.omit(yrs[,i]),i]
  logRS2.mat[1:length(na.omit(yrs[,i])),i]<-logRS2[na.omit(yrs[,i]),i]
  logRS3.mat[1:length(na.omit(yrs[,i])),i]<-logRS3[na.omit(yrs[,i]),i]
}

pops_E<-unique(SR.dat[which(SR.dat$stat_area<5),1])
pops_2<-seq(1,52,1)
pops_2<-pops_2[-c(pops_E)]

for(i in pops_2){
  logRSB.mat[1:length(na.omit(yrs[,i])),i]<-logRS2[na.omit(yrs[,i]),i]
}

for(i in pops_E){
  logRSB.mat[1:length(na.omit(yrs[,i])),i]<-logRS1[na.omit(yrs[,i]),i]
  
}

## Will eventually need to subset out populations where we want to analyze RS2 instead of RS1

yrs.dat<-rep(NA,n.pops)
for (i in 1:n.pops){
  yrs.dat[i]<-length(na.omit(yrs[,i]))
}

#################################################
##### Model 3 - hierarchical model w/time-variant (recursive) alphas and uninformative priors
###   with reginal hyper priors mu.alpha which is also time varying, estimating overall change in productivity 
###   for each population group

####   The next step would be to add a covariance matrix between annual productivities

## I grouped rivers and smiths inlet with Area 7 & 8

co_pops<-read.table("Data/coho_groups.txt",header=TRUE)
co_pops$mean_total<-NA

for(i in 1:52){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6

### using mean run size (harvest scenario 1) across the timeseries as our prior on capacity
Smax.p<-rep(NA,52)
Smax.tau<-rep(NA,52)
n.pops<-max(SR.dat$pop_no)

for(i in 1:n.pops){
  Smax.p[i]<-log(mean(na.omit(rec1[,i])))
  Smax.tau[i] <- 0.2*Smax.p[i] #scales var to Smax by multiplying log abundance by 0.2 for our SD
}

Kgroups <- 6 # number of groups that we want to see correlated
covarGroups <- matrix(0,nrow=Kgroups,ncol=Kgroups)
diag(covarGroups) <- 1

group <- co_pops$group
Npops <- length(group)

resultSR_B3<-readRDS("coho central coast/COSR3B.cormat_lalpha_MCMC.rds")
mcmc_names <- colnames(resultSR_B3[[1]])
#######################################################################
#### estimating Smsy using lamberts w Scheuerell method
Smsy<-Sgen<-matrix(data=NA,nrow=samples,ncol=n.pops)

samp.chain<-sample(1:n_chains,samples,replace=TRUE)
samp.MCMC<-sample(1:samples, samples, replace=FALSE)
# Sgen; solve Smsy=Sgen*exp(a*(1- Sgen/S(k))) from Holt et al. 2009 Indicators of Status and Benchmarks for Conservation Units in Canada's Wild Salmon Policy 
## just need to adjust the length based on how long the MCMC is
Sgen_find <- function(Smsy,Sgen,a,b)
{
  (Smsy-(Sgen*exp(a*(1- Sgen/(-a/b)))))^2
}
for (j in 1:n.pops){
  samp.chain<-sample(1:n_chains,samples,replace=TRUE)
  samp.MCMC<-sample(1:samples, samples, replace=FALSE)
  
  for (i in 1:samples){
    alpha <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[j]]
    beta <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[j]]
    Smsy[i,j]<-(1 - lambert_W0(exp(1 - (alpha)))) / (beta) #Smsy in 1000 random draws from MCMC
    Sgen[i,j] <- optimize(Sgen_find,c(0,Smsy[i,j]),tol=0.0001,Smsy=Smsy[i,j],a=alpha,b=beta)$minimum
  }
}

Smsy_priors <- apply(Smsy,2,FUN=function(x){c(mean(x),sd(x))})
Smsy_priors <- data.frame("pop"=1:n.pops,"mean"=Smsy_priors[1,],"tau"=1/(Smsy_priors[2,]^2))
Sgen_priors <- apply(Sgen,2,FUN=function(x){c(mean(x),sd(x))})
Sgen_priors <- data.frame("pop"=1:n.pops,"mean"=Sgen_priors[1,],"tau"=1/(Sgen_priors[2,]^2))
sim_year <- rep(1,n.years+sim_eval)
sim_year[min(SR.dat$year):(max(SR.dat$year)+sim_eval)<1997 | (min(SR.dat$year):(max(SR.dat$year)+sim_eval)>2000 & min(SR.dat$year):(max(SR.dat$year)+sim_eval)<2017)] <- 0
future_sim <- rep(0,n.years+sim_eval)
future_sim[min(SR.dat$year):(max(SR.dat$year)+sim_eval)>2017] <- 1
baselines <- sapply(1:n.pops,function(x){mean(SR.dat$total_runE[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines_2 <- sapply(1:n.pops,function(x){mean(SR.dat$total_run2[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines<-ifelse(!is.na(baselines_2),baselines_2,baselines)


ER <- aggregate(mean_ER~group_no+fishery,data=ER_values,mean)
ER_values$cv <- ER_values$sd_ER/ER_values$mean_ER
escape_reg <- matrix(0,nrow=Kgroups,ncol=5,dimnames=list("groups"=ER_values$group_name[ER_values$fishery=="bc_total"],"regs"=c("No harvest","10-year average","50% BC reduction","50% AK reduction","50% AK & BC reduction")))
escape_reg[,1] <- 0
escape_reg[,2] <- aggregate(mean_ER~group_no,data=ER_values,sum)$mean_ER
escape_reg[,3] <- sapply(1:Kgroups,function(x){ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + 0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_reg[,4] <- sapply(1:Kgroups,function(x){0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_reg[,5] <- sapply(1:Kgroups,function(x){0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + 0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})

ER_cv <- aggregate(cv~group_no+fishery,data=ER_values,mean)
escape_sd_reg <- matrix(0,nrow=Kgroups,ncol=5,dimnames=list("groups"=ER_values$group_name[ER_values$fishery=="bc_total"],"regs"=c("No harvest","10-year average","50% BC reduction","50% AK reduction","50% AK & BC reduction")))
escape_sd_reg[,1] <- 1e-6
escape_sd_reg[,2] <- aggregate(sd_ER~group_no,data=ER_values,sum)$sd_ER
escape_sd_reg[,3] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_sd_reg[,4] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_sd_reg[,5] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
prop_mat <- matrix(NA,ncol=3,nrow=n.pops)
for(x in 1:n.pops)
{
  prop_mat[x,] <- c(mean(SR.dat$age_3[SR.dat$pop_no==x],na.rm=TRUE),mean(SR.dat$age_4[SR.dat$pop_no==x],na.rm=TRUE),mean(SR.dat$age_5[SR.dat$pop_no==x],na.rm=TRUE))
}
prop_mat[,1] <- 1-(prop_mat[,2]+prop_mat[,3])
spawners_obs <- rbind(spawn.mat,matrix(NA,nrow=sim_eval,ncol=n.pops))
spawners_obs[is.na(spawners_obs)] <- 0
datamcmc2 <- list(
  "logRS"=logRSo,
  "spawners"=spawn,
  #"yrs.dat"=yrs.dat,
  #"yrs_obs"=yrs,
  "Smax.p"=Smax.p,
  "Smax.tau"=Smax.tau,
  "n.years"=n.years,
  "sim_eval"=sim_eval,
  "Kgroups"=Kgroups,
  "covarGroups"=covarGroups,
  "init_lalpha"=rep(1,Kgroups),
  "group"=group,
  "Npops"=Npops,
  #"yrs.sim"=(1:(n.years+sim_eval))[sim_year==1],
  "Smsy_mean"=Smsy_priors$mean,
  "Smsy_tau"=Smsy_priors$tau,
  "Sgen_mean"=Sgen_priors$mean,
  #"simulate_year"=sim_year,
  #"future_sim"=future_sim,
  "Nregs"=ncol(escape_reg),
  "mu_reg"=escape_reg,
  "sd_reg"=escape_sd_reg,
  "mu_1996"=escape_reg[,2],
  "sd_1996"=escape_sd_reg[,2],
  "baseline_index"=baselines,
  "prop"=prop_mat,
  "init_s0"=colMeans(spawn,na.rm=TRUE)
)
sigma_sd <- apply(logRSo,2,sd,na.rm=TRUE)
spawn_sd <- apply(log(spawn),2,sd,na.rm=TRUE)

variables_1 <- c("msy_ratio","sgen_ratio","run_ratio")#,"beta","tau.alpha","sigma.alpha","tau","sigma") # these are the variables to keep track of
variables_2 <- c("escape_rate","new_run","spawners","new_spawners")#,"beta","tau.alpha","sigma.alpha","tau","sigma") # these are the variables to keep track of
variables_3 <- c("new_recruits","logRS","new_logRS")#,"beta","tau.alpha","sigma.alpha","tau","sigma") # these are the variables to keep track of
init_fx <- function(chain_id)
{
  list("Smax"=Smax.p,
       "tauMVN"=covarGroups,
       "sigma_alpha"=rep(5,Kgroups),
       "mu_lalpha"=matrix(0,nrow=n.years+sim_eval,ncol=Kgroups),
       "sigma_s0"=spawn_sd,
       "sigma"=sigma_sd,
       "spawn_trend"=rep(0,Kgroups),
       "spawn_trend_pop"=rep(0,Npops))
}

temp3B <- jags.model(file="JAGS/forecasting model tr1.jags", data=datamcmc2 ,n.chains = n_chains, n.adapt=50000, quiet=FALSE,inits=init_fx)#, inits=inits)
update(temp3B, n.iter = 75000)  # burnin
result_forecast_metrics <- coda.samples(model=temp3B, variable.names=variables_1, n.iter=thinning*samples, thin = thinning) 
result_forecast_run <- coda.samples(model=temp3B, variable.names=variables_2, n.iter=thinning*samples, thin = thinning) 
result_forecast_RS <- coda.samples(model=temp3B, variable.names=variables_3, n.iter=thinning*samples, thin = thinning) 
saveRDS(result_forecast_metrics,file="Results/coho_forecasting_tr1.rds")
saveRDS(result_forecast_run,file="Results/coho_forecasting_popdyn_tr1.rds")
saveRDS(result_forecast_RS,file="Results/coho_forecasting_rec_tr1.rds")

saveRDS(datamcmc2,file="Data/jags_forecasting_data_tr1.rds")
