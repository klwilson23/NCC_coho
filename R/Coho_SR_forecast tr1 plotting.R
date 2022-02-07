###########################################################################
###     Hierarchical SR models w/time variant alphas, PR priors        ####
###########################################################################

#setwd("C:/data/centralcoast")
require(rjags)
require(R2jags)
library(runjags)
require(gsl)

result_forecast_metrics <- readRDS("Results/coho_forecasting_tr1.rds")
result_forecast_run <- readRDS("Results/coho_forecasting_popdyn_tr1.rds")
result_forecast_RS <- readRDS("Results/coho_forecasting_rec_tr1.rds")

sim_eval <- 20 # forecast model for 16 years, or 4 generations
thinning <- 200
samples <- 1500
n_chains <- 4
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat_full <- SR.dat
SR.dat_predict<-subset(SR.dat,SR.dat$year>=2017)
SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)
pop_names <- unique(SR.dat$population[order(SR.dat$pop_no)])
group_names <- c("Area 7-10","Hecate Lowlands - Area 5/6","Area 6 - Inner Waters","Haida Gwaii - Area 2E","Skeena - Area 4","Nass - Area 3")
group_names <- group_names[c(4,6,5,2,3,1)]

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
pops_2<-seq(1,n.pops,1)
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

for(i in co_pops$pop_no){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6

### using mean run size (harvest scenario 1) across the timeseries as our prior on capacity
Smax.p<-rep(NA,n.pops)
Smax.tau<-rep(NA,n.pops)
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

resultSR_B3<-readRDS("Results/COSR3B.tr1_lalpha_MCMC.rds")
mcmc_names <- colnames(resultSR_B3[[1]])
#######################################################################
#### estimating Smsy using lamberts w Scheuerell method
Smsy<-Sgen<-matrix(data=NA,nrow=samples,ncol=n.pops)

samp.chain<-sample(1:n_chains,samples,replace=TRUE)
samp.MCMC<-sample(1:samples, samples, replace=FALSE)
# Sgen; solve Smsy=Sgen*exp(a*(1- Sgen/S(k))) from Holt et al. 2009 Indicators of Status and Benchmarks for Conservation Units in Canada's Wild Salmon Policy 
## just need to adjust the length based on how long the MCMC is
Smsy<-Sgen<-Umsy<-matrix(data=NA,nrow=samples,ncol=n.pops)

samp.chain<-sample(1:n_chains,samples,replace=TRUE)
samp.MCMC<-sample(1:samples, samples, replace=FALSE)
# Sgen; solve Smsy=Sgen*exp(a*(1- Sgen/S(k))) from Holt et al. 2009 Indicators of Status and Benchmarks for Conservation Units in Canada's Wild Salmon Policy 
# Sgen; solve Smsy=Sgen*exp(a*(1- Sgen/S(k))) from Holt et al. 2009 Indicators of Status and Benchmarks for Conservation Units in Canada's Wild Salmon Policy 
## just need to adjust the length based on how long the MCMC is
Sgen_find <- function(Smsy,Sgen,a,b)
{
  (Smsy-(Sgen*exp(a*(1- Sgen/(a/b)))))^2
}
umsy_find <- function(Smsy,Umsy,a,b)
{
  #Umsy <- Umsy_priors$mean[1]
  #Smsy <- Smsy_priors$mean[1]
  #a <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[1]]
  #b <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[1]]
  (Smsy-(1-Umsy)*(a/b))^2
}
for (j in 1:n.pops){
  samp.chain<-sample(1:n_chains,samples,replace=TRUE)
  samp.MCMC<-sample(1:samples, samples, replace=FALSE)
  
  for (i in 1:samples){
    alpha <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[j]]
    beta <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[j]]
    Smsy[i,j]<-(1 - lambert_W0(exp(1 - (alpha)))) / (beta) #Smsy in 1000 random draws from MCMC
    Sgen[i,j] <- optimize(Sgen_find,c(0,Smsy[i,j]),tol=0.0001,Smsy=Smsy[i,j],a=alpha,b=beta)$minimum
    Umsy[i,j] <- optimize(umsy_find,c(0,1),tol=0.0001,Smsy=Smsy[i,j],a=alpha,b=beta)$minimum
  }
}

Smsy_priors <- apply(Smsy,2,FUN=function(x){c(mean(x),sd(x))})
Smsy_priors <- data.frame("pop"=1:n.pops,"mean"=Smsy_priors[1,],"tau"=1/(Smsy_priors[2,]^2))
Sgen_priors <- apply(Sgen,2,FUN=function(x){c(mean(x),sd(x))})
Sgen_priors <- data.frame("pop"=1:n.pops,"mean"=Sgen_priors[1,],"tau"=1/(Sgen_priors[2,]^2))
Umsy_priors <- apply(Umsy,2,FUN=function(x){c(mean(x),sd(x))})
Umsy_priors <- data.frame("pop"=1:n.pops,"mean"=Umsy_priors[1,],"tau"=1/(Umsy_priors[2,]^2))


sim_year <- rep(1,n.years+sim_eval)
sim_year[min(SR.dat$year):(max(SR.dat$year)+sim_eval)<1997 | (min(SR.dat$year):(max(SR.dat$year)+sim_eval)>2000 & min(SR.dat$year):(max(SR.dat$year)+sim_eval)<2017)] <- 0
future_sim <- rep(0,n.years+sim_eval)
future_sim[min(SR.dat$year):(max(SR.dat$year)+sim_eval)>2017] <- 1
baselines <- sapply(1:n.pops,function(x){mean(SR.dat$total_runE[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines_2 <- sapply(1:n.pops,function(x){mean(SR.dat$total_run2[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines<-ifelse(!is.na(baselines_2),baselines_2,baselines)

baseline_spawn <- sapply(1:n.pops,function(x){mean(SR.dat$escapement[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})

s_over_msy <- sapply(1:n.pops,function(x){sum(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]<=(Smsy_priors$mean[x]),na.rm=TRUE)/sum(!is.na(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]))})
data.frame("pop"=pop_names,"s_over_msy"=s_over_msy)

SR.dat_full$total_run <- ifelse(!is.na(SR.dat_full$total_run2),SR.dat_full$total_run2,SR.dat_full$total_runE)
pop_sub <- aggregate(total_run~pop_no,SR.dat_full,mean,na.rm=TRUE)
run_over_baseline <- 100*(sum(pop_sub$total_run)-sum(baselines[pop_sub$pop_no]))/sum(baselines[pop_sub$pop_no])

SR.dat_full$baseline_run <- (baselines)[SR.dat_full$pop_no]
SR.dat_full$lrp <- (Smsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$ucur <- ifelse(!is.na(SR.dat_full$er_2),SR.dat_full$er_2,SR.dat_full$er_E)
SR.dat_full$umsy <- (Umsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$regime <- ifelse(SR.dat_full$year>=1990 & SR.dat_full$year<=2001,"early",ifelse(SR.dat_full$year<=2011,"recovering","recent"))

pop_sub <- aggregate(cbind(total_run/baseline_run,escapement/lrp,ucur/umsy)~pop_no+regime,SR.dat_full,mean,na.rm=F)
recent <- pop_sub[pop_sub$regime=="recent",]
early <- pop_sub[pop_sub$regime=="early",]
recovering <- pop_sub[pop_sub$regime=="recovering",]
time_series <- merge(merge(early,recovering,by="pop_no"),recent,by="pop_no")
time_series$pop <- pop_names[time_series$pop_no]

SR.dat_2016 <- SR.dat_full[SR.dat_full$year>=2018,]
SR.dat_2016$total_run <- ifelse(!is.na(SR.dat_2016$total_run2),SR.dat_2016$total_run2,SR.dat_2016$total_runE)

par(mar=c(4,4,0.1,0.1))
status <- data.frame("pop"=pop_names[recent$pop_no],"s_over_msy"=recent$V2,"baseline"=recent$V1,"u_over_umsy"=recent$V3)
plot(NA,ylab=expression(u[cur]/u[MSY]),xlab=expression(S[cur]/S[MSY]),xlim=1.1*range(c(0,2,status$s_over_msy)),ylim=1.1*range(c(0,2,status$u_over_umsy)))
polygon(y=c(-5,1,1,-5),x=c(1,1,40,40),col=adjustcolor("seagreen",1),border=NA,xpd=FALSE)
polygon(y=c(1,1,20,20),x=c(-5,1,1,-5),col=adjustcolor("tomato",1),border=NA,xpd=FALSE)
polygon(y=c(1,20,20,1),x=c(1,1,40,40),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
polygon(y=c(-5,1,1,-5),x=c(-5,-5,1,1),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
abline(h=1,v=1,lty=2,lwd=0.5,col="grey10")
points(status$s_over_msy,status$u_over_umsy,pch=21,bg="grey80",col="white")
samp <- sapply(c(0,0.5,1,2),function(x){which.min(abs(status$s_over_msy-x))})
text(status$s_over_msy[samp],status$u_over_umsy[samp],labels=status$pop[samp],adj=c(0,0),offset=1.5,cex=0.7,col="black")
Corner_text("overfishing, not overfished","topright",cex=0.7)
Corner_text("healthy","bottomright",cex=0.7)
Corner_text("overfishing and overfished","topleft",cex=0.7)
Corner_text("no overfishing, but overfished","bottomleft",cex=0.7)

sum(status[,2]>1 & status[,4]<1) # healthy populations
sum(status[,2]>1 & status[,4]>1) # healthy, but overfishing
sum(status[,2]<1 & status[,4]>1) # overfished, overfishing
sum(status[,2]<1 & status[,4]<1) # overfished, but not overfishing

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
escape_sd_reg[,5] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
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
variables_2 <- c("new_recruits","spawners","new_spawners")#,"beta","tau.alpha","sigma.alpha","tau","sigma") # these are the variables to keep track of
init_fx <- function(chain_id)
{
  list("Smax"=Smax.p,
       "tauMVN"=covarGroups,
       "sigma_alpha"=rep(5,Kgroups),
       "mu_lalpha"=matrix(0,nrow=n.years+sim_eval,ncol=Kgroups),
       "sigma_s0"=spawn_sd,
       "sigma"=sigma_sd)
}

############### plotting ################
prob_fore <- c(0.025,0.25,0.5,0.75,0.975)
mcmc_names <- colnames(result_forecast_metrics[[1]])
years<-seq(1980,2016,1)
layout(matrix(1:10,nrow=5,byrow=TRUE))
for(i in 1:datamcmc2$Nregs)
{
  msy_scenario_1 <- result_forecast_metrics[,grepl("\\bmsy_ratio",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  msy_scenario_1 <- as.matrix(msy_scenario_1)
  msy_scenario_1[1,1:datamcmc2$sim_eval]
  msy_scenario_1 <- array(unlist(msy_scenario_1),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  msy_scenario_1[1:10,1:2,1]
  str(msy_scenario_1)
  msy_scen1_ci <- apply(msy_scenario_1,c(2,3),quantile,probs=prob_fore) 
  matplot(msy_scen1_ci[3,,],type="l")
  boxplot(msy_scen1_ci[3,sim_eval,],ylim=c(0,4))
  legend("topleft",colnames(escape_reg)[i],bty="n")
}

prob_fore <- c(0.025,0.25,0.5,0.75,0.975)
mcmc_names <- colnames(result_forecast_metrics[[1]])
years<-seq(1980,2016,1)
layout(matrix(1:10,nrow=5,byrow=TRUE))
for(i in 1:datamcmc2$Nregs)
{
  msy_scenario_1 <- result_forecast_metrics[,grepl("\\bsgen_ratio",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  msy_scenario_1 <- as.matrix(msy_scenario_1)
  msy_scenario_1[1,1:datamcmc2$sim_eval]
  msy_scenario_1 <- array(unlist(msy_scenario_1),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  msy_scenario_1[1:10,1:2,1]
  str(msy_scenario_1)
  msy_scen1_ci <- apply(msy_scenario_1,c(2,3),quantile,probs=prob_fore) 
  matplot(msy_scen1_ci[3,,],type="l")
  boxplot(msy_scen1_ci[3,sim_eval,],ylim=c(0,4))
  legend("topleft",colnames(escape_reg)[i],bty="n")
}

mcmc_names <- colnames(result_forecast_metrics[[1]])
years<-seq(1980,2016,1)
layout(matrix(1:10,nrow=5,byrow=TRUE))
for(i in 1:datamcmc2$Nregs)
{
  msy_scenario_1 <- result_forecast_metrics[,grepl("\\brun_ratio",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  msy_scenario_1 <- as.matrix(msy_scenario_1)
  msy_scenario_1[1,1:datamcmc2$sim_eval]
  msy_scenario_1 <- array(unlist(msy_scenario_1),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  #print(msy_scenario_1[1:10,5:6,1])
  msy_scen1_ci <- apply(msy_scenario_1,c(2,3),quantile,probs=prob_fore) 
  matplot(msy_scen1_ci[3,,],type="l")
  boxplot(msy_scen1_ci[3,sim_eval,],ylim=c(0,4))
  legend("topleft",colnames(escape_reg)[i],bty="n")
}

#### plot population dynamics

mcmc_names <- colnames(result_forecast_run[[1]])

spawn_post <- result_forecast_run[,grepl("\\bspawners",mcmc_names)]
spawn_post[[1]][1:10,1:datamcmc2$sim_eval]
spawn_post <- as.matrix(spawn_post)
spawn_post[1,1:datamcmc2$sim_eval]
spawn_post <- array(unlist((spawn_post)),dim=c(n_chains*samples,datamcmc2$n.years,datamcmc2$Npops))
spawn_ci <- apply(spawn_post,c(2,3),quantile,probs=prob_fore)
rmse <- below_lrp <- below_msy <- matrix(NA,nrow=5,ncol=n.pops,dimnames=list("Scenarios"=colnames(escape_reg),"Populations"=pop_names))

for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  matplot(spawn_forecast,type="l")
  pdf(paste("Figures/coho forecasting tr1 scenario ",i,".pdf",sep=""),height=7,width=10)
  layout(matrix(1:12,ncol=4,byrow=TRUE))
  par(mar=c(5,4,0.1,0.1))
  for(j in 1:n.pops)
  {
    year_seq <- c(years,max(years)+1:sim_eval)
    plot(year_seq,spawn_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato",xlab="",ylab="Spawners",ylim=range(spawn_forecast_arr[,j,]))
    polygon(c(year_seq,rev(year_seq)),c(spawn_forecast_arr[,j,2],rev(spawn_forecast_arr[,j,3])),col=adjustcolor("tomato",0.5),border=NA)
    lines(year_seq,spawn_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato")
    points(years,spawn[,j],pch=21,bg="grey90")
    subdat_predict <- subset(SR.dat_predict,SR.dat_predict$pop_no==j)
    points(subdat_predict$year,subdat_predict$escapement,pch=22,bg="dodgerblue",cex=1.2)
    abline(h=0.8*Smsy_priors$mean[j],lwd=1,lty=1,col="black")
    abline(h=Sgen_priors$mean[j],lwd=1,lty=2,col="red")
    legend("topleft",paste(pop_names[j],",",colnames(escape_reg)[i]),bty="n")
    yr_match <- match(subdat_predict$year[!is.na(subdat_predict$escapement)],c(max(years)+1:sim_eval))
    rmse[i,j] <- ifelse(sum(!is.na(subdat_predict$escapement))>0,median(abs(spawn_forecast_1[,yr_match,j]-subdat_predict$escapement[!is.na(subdat_predict$escapement)])),NA)
    below_lrp[i,j] <- ifelse(spawn_forecast[year_seq==2025,j]<=(Sgen_priors$mean[j]),1,0)
    below_msy[i,j] <- ifelse(spawn_forecast[year_seq==2025,j]<=(0.8*Smsy_priors$mean[j]),1,0)
  }
  dev.off()
}

saveRDS(rmse,file="Results/predictive error tr1.rds")
below_run <- matrix(NA,nrow=5,ncol=n.pops,dimnames=list("Scenarios"=colnames(escape_reg),"Populations"=pop_names))
for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_run",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  matplot(spawn_forecast,type="l")
  layout(matrix(1:12,ncol=4,byrow=TRUE))
  par(mar=c(5,4,0.1,0.1))
  for(j in 1:n.pops)
  {
    year_seq <- c(years,max(years)+1:sim_eval)
    plot(year_seq,spawn_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato",xlab="",ylab="Run size",ylim=range(spawn_forecast_arr[,j,]))
    polygon(c(year_seq,rev(year_seq)),c(spawn_forecast_arr[,j,2],rev(spawn_forecast_arr[,j,3])),col=adjustcolor("tomato",0.5),border=NA)
    lines(year_seq,spawn_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato")
    points(years,spawn[,j],pch=21,bg="grey90")
    subdat_predict <- subset(SR.dat_predict,SR.dat_predict$pop_no==j)
    points(subdat_predict$year,ifelse(!is.na(subdat_predict[,"total_run2"]),subdat_predict[,"total_run2"],subdat_predict[,"total_runE"]),pch=22,bg="dodgerblue",cex=1.2)
    abline(h=0.8*Smsy_priors$mean[j],lwd=1,lty=1,col="black")
    abline(h=Sgen_priors$mean[j],lwd=1,lty=2,col="red")
    legend("topleft",paste(pop_names[j],",",colnames(escape_reg)[i]),bty="n")
    below_run[i,j] <- ifelse(spawn_forecast[year_seq==2025,j]<=(baselines[j]),1,0)
  }
}

rowSums(below_run)/n.pops
rowSums(below_lrp)/n.pops
rowSums(below_msy)/n.pops

layout(1)
plot(spawn_forecast[,30],pch=21,bg="grey90")
lines(spawn[,30],lwd=2,col="red")

mcmc_names <- colnames(result_forecast_RS[[1]])

rec_post <- result_forecast_RS[,grepl("\\bnew_recruits",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
rec_post <- as.matrix(rec_post)
rec_post <- array(unlist((rec_post)),dim=c(n_chains*samples,datamcmc2$n.years+datamcmc2$sim_eval,datamcmc2$Npops))
rec_ci <- apply(rec_post,c(2,3),quantile,probs=prob_fore,na.rm=TRUE)
plot(rec_ci[3,,1],type="l")
lines(rec1[,1],lwd=2,col="red")

for(i in n.pops)
{
  plot(spawn_forecast[,i],log(rec_ci[3,,i]/spawn_forecast[,i]),ylim=range(log(rec_ci[3,,i]/spawn_forecast[,i]),log(rec1[,i]/spawn_forecast[1:n.years,i]),na.rm=TRUE),xlab="spawners",ylab="productivity")
  points(spawn[,i],log(rec1[,i]/spawn[,i]),pch=21,bg="red")
  points(spawn[,i],logRSo[,i],pch=21,bg="tomato")
}

plot(spawn_forecast[,1],ylim=range(rec_ci[3,,1],spawn_forecast[,1]))
lines(rec_ci[3,,1])
lines(rec1[,1],lwd=2,col="red")
plot(rec1[,1],rec_ci[3,1:n.years,1])


# recruitment dynamics:
mcmc_names <- colnames(result_forecast_RS[[1]])
rec_post <- result_forecast_RS[,grepl("\\blogRS",mcmc_names)]
rec_post[[1]][1:10,1:datamcmc2$sim_eval]
rec_post <- as.matrix(rec_post)
rec_post[1,1:datamcmc2$sim_eval]
rec_post <- array(unlist((rec_post)),dim=c(n_chains*samples,datamcmc2$n.years,datamcmc2$Npops))
rec_ci <- apply(rec_post,c(2,3),quantile,probs=prob_fore)
for(i in 1:datamcmc2$Nregs)
{
  rec_forecast_1 <- result_forecast_RS[,grepl("\\bnew_logRS",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  rec_forecast_1[[1]][1:10,1:datamcmc2$sim_eval]
  rec_forecast_1 <- as.matrix(rec_forecast_1)
  rec_forecast_1[1,1:datamcmc2$sim_eval]
  rec_forecast_1 <- array(unlist((rec_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  rec_fore_ci <- apply(rec_forecast_1,c(2,3),quantile,probs=prob_fore)
  rec_forecast_mdn <- sapply(1:n.pops,function(x){rec_forecast_1[which.min(abs(rec_fore_ci[3,datamcmc2$sim_eval,x]-rec_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  rec_forecast <- rbind(rec_ci[3,,],rec_fore_ci[3,,])
  rec_forecast_ui <- rbind(rec_ci[4,,],rec_fore_ci[4,,])
  rec_forecast_li <- rbind(rec_ci[2,,],rec_fore_ci[2,,])
  rec_forecast_arr <- simplify2array(list(rec_forecast,rec_forecast_ui,rec_forecast_li))
  matplot(rec_forecast,type="l")
  pdf(paste("Figures/coho productivity tr1 forecast scenario ",i,".pdf",sep=""),height=7,width=10)
  layout(matrix(1:12,ncol=4,byrow=TRUE))
  par(mar=c(5,4,0.1,0.1))
  for(j in 1:n.pops)
  {
    year_seq <- c(years,max(years)+1:sim_eval)
    plot(year_seq,rec_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato",xlab="",ylab="ln (R/S)",ylim=range(rec_forecast_arr[,j,]))
    polygon(c(year_seq,rev(year_seq)),c(rec_forecast_arr[,j,2],rev(rec_forecast_arr[,j,3])),col=adjustcolor("tomato",0.5),border=NA)
    lines(year_seq,rec_forecast_arr[,j,1],type="l",lty=1,lwd=0.5,col="tomato")
    points(years,logRSo[,j],pch=21,bg="grey90")
    subdat_predict <- subset(SR.dat_predict,SR.dat_predict$pop_no==j)
    points(subdat_predict$year,log(ifelse(!is.na(subdat_predict[,"RS_2"]),subdat_predict[,"RS_2"],subdat_predict[,"RS_E"])),pch=22,bg="dodgerblue",cex=1.2)
    abline(h=mean(logRSo[,j],na.rm=TRUE),lwd=1,lty=1,col="black")
    legend("topleft",paste(pop_names[j],",",colnames(escape_reg)[i]),bty="n")
  }
dev.off()
}


#### result figures:
jpeg(filename="Figures/risk posterior mass results tr1.jpeg", width=8,height=8, units="in", res=600)
mcmc_names <- colnames(result_forecast_run[[1]])
prob_fore <- c(0.025,0.4,0.5,0.6,0.975)

matLayout <- matrix(0,nrow=8,ncol=4,byrow=TRUE)
matLayout[1,2] <- 1
matLayout[1,3] <- 7

matLayout[2,2] <- 2
matLayout[2,3] <- 8

matLayout[3,2] <- 3
matLayout[3,3] <- 9

matLayout[4,2] <- 4
matLayout[4,3] <- 10

matLayout[5,2] <- 5
matLayout[5,3] <- 11

matLayout[6,2] <- 6
matLayout[6,3] <- 12

matLayout[7:8,1:2] <- 13
matLayout[7:8,3:4] <- 14
layout(matLayout)

for(j in 1:Kgroups)
{
  year_seq <- c(years,max(years)+1:sim_eval)
  year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
  par(mar=c(0,0,0,2.5))
  if(j==3)
  {
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,100),ylab="",xaxt="n",pch=0,yaxt="n",xpd=NA,xlab="")
    mtext("% run below baseline",side=2,line=2.5)
  }else{
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,100),ylab="",xaxt="n",pch=0,yaxt="n",xlab="")
  }
  text(median(year_seq[year_clip]),95,group_names[j],font=2)
  axis(2,at=c(20,40,60,80))
  for(i in 1:datamcmc2$Nregs)
  {
    spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_run",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
    spawn_forecast_1 <- as.matrix(spawn_forecast_1)
    spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
    spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
    spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<baselines[x])}))
    spawn_fore_risk <- colSums(spawn_fore_risk[group==j,])/(n_chains*samples*sum(group==j))
    spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
    spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
    spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
    spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
    spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
    group_mdn <- rowSums(spawn_forecast_arr[,group==j,1])/sum(baselines[group==j])
    group_ui <- rowSums(spawn_forecast_arr[,group==j,2])/sum(baselines[group==j])
    group_li <- rowSums(spawn_forecast_arr[,group==j,3])/sum(baselines[group==j])
    ylims <- range(group_mdn,group_ui,group_li)
    #polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
    lines(year_seq[year_clip],100*spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
    #abline(h=1,lwd=1,lty=1,col="black")
  }
}
axis(1)

for(j in 1:Kgroups)
{
  year_seq <- c(years,max(years)+1:sim_eval)
  year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
  par(mar=c(0,2.5,0,0)) 
  if(j==3)
  {
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,100),ylab="",xaxt="n",pch=0,yaxt="n",xpd=NA,xlab="")
    mtext("% spawners below 0.8S(msy)",side=2,line=2.5)
  }else{
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,100),ylab="",xaxt="n",pch=0,yaxt="n",xlab="")
  }
  text(median(year_seq[year_clip]),95,group_names[j],font=2)
  axis(2,at=c(20,40,60,80))
  for(i in 1:datamcmc2$Nregs)
  {
    spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
    spawn_forecast_1 <- as.matrix(spawn_forecast_1)
    spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
    spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
    spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<(0.8*Smsy_priors$mean[x]))}))
    spawn_fore_risk <- colSums(spawn_fore_risk[group==j,])/(n_chains*samples*sum(group==j))
    
    spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
    spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
    spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
    spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
    spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
    group_mdn <- rowSums(spawn_forecast_arr[,group==j,1])/sum(Smsy_priors$mean[group==j])
    group_ui <- rowSums(spawn_forecast_arr[,group==j,2])/sum(Smsy_priors$mean[group==j])
    group_li <- rowSums(spawn_forecast_arr[,group==j,3])/sum(Smsy_priors$mean[group==j])
    ylims <- range(group_mdn,group_ui,group_li)
    #polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
    lines(year_seq[year_clip],100*spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
    #abline(h=1,lwd=1,lty=1,col="black")
  }
}
axis(1)

year_seq <- c(years,max(years)+1:sim_eval)
year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
par(mar=c(4.5,4,3,2))
plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,100),ylab="% run below baseline",pch=0,yaxt="n",xlab="Year")
axis(2,at=c(20,40,60,80))

for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_run",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  
  spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<baselines[x])}))
  spawn_fore_risk <- colSums(spawn_fore_risk)/(n_chains*samples*n.pops)
  
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  group_mdn <- rowSums(spawn_forecast_arr[,,1])/sum(baselines)
  group_ui <- rowSums(spawn_forecast_arr[,,2])/sum(baselines)
  group_li <- rowSums(spawn_forecast_arr[,,3])/sum(baselines)
  ylims <- range(group_mdn,group_ui,group_li)
  #polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
  lines(year_seq[year_clip],100*spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
  #abline(h=1,lwd=1,lty=1,col="black")
  text(x=year_seq[length(year_seq)]-1,y=100*spawn_fore_risk[length(year_seq)],LETTERS[i],font=2)
}

year_seq <- c(years,max(years)+1:sim_eval)
year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
par(mar=c(4.5,4,3,2))
boxplot(NA,xlim=c(0.5,5.5),ylim=c(0,105),ylab="% spawners below 0.8S(msy)",xaxt="n",col=0,yaxt="n",xlab="Management scenario",xpd=NA)
for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  
  spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<(0.8*Smsy_priors$mean[x]))}))
  risk_2030 <- sum(spawn_fore_risk[,year_clip==which(year_seq==2030)])/(n_chains*samples*n.pops)
  spawn_fore_risk <- rowSums(spawn_fore_risk)/(n_chains*samples*sim_eval)
  
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  group_mdn <- rowSums(spawn_forecast_arr[year_clip,,1])/sum(Smsy_priors$mean)
  group_ui <- rowSums(spawn_forecast_arr[year_clip,,2])/sum(Smsy_priors$mean)
  group_li <- rowSums(spawn_forecast_arr[year_clip,,3])/sum(Smsy_priors$mean)
  ylims <- range(group_mdn,group_ui,group_li)
  #polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
  boxplot(100*spawn_fore_risk,at=i,add=TRUE,col=adjustcolor(i,0.5))
  #abline(h=1,lwd=1,lty=1,col="black")
  text(x=i,y=-20,LETTERS[i],font=2,xpd=NA)
  text(x=i,y=1.15*max(100*spawn_fore_risk),round(mean(100*spawn_fore_risk),0),font=2,xpd=NA)
  
}

legend(x=3,y=225,legend=colnames(escape_reg),title="Management scenario",bty="n",pch=22,pt.bg=adjustcolor(1:4,0.9),xpd=NA)

dev.off()


#### result figures:
jpeg(filename="Figures/posterior median results tr1.jpeg", width=8,height=8, units="in", res=600)
mcmc_names <- colnames(result_forecast_run[[1]])
prob_fore <- c(0.025,0.4,0.5,0.6,0.975)

matLayout <- matrix(0,nrow=8,ncol=4,byrow=TRUE)
matLayout[1,2] <- 1
matLayout[1,3] <- 7

matLayout[2,2] <- 2
matLayout[2,3] <- 8

matLayout[3,2] <- 3
matLayout[3,3] <- 9

matLayout[4,2] <- 4
matLayout[4,3] <- 10

matLayout[5,2] <- 5
matLayout[5,3] <- 11

matLayout[6,2] <- 6
matLayout[6,3] <- 12

matLayout[7:8,1:2] <- 13
matLayout[7:8,3:4] <- 14
layout(matLayout)

for(j in 1:Kgroups)
{
  year_seq <- c(years,max(years)+1:sim_eval)
  year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
  par(mar=c(0,0,0,2))
  if(j==3)
  {
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,1.25),ylab="",xaxt="n",pch=0,yaxt="n",xpd=NA,xlab="")
    mtext("Relative run size (to 2000-15 baseline)",side=2,line=2.5)
  }else{
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,1.25),ylab="",xaxt="n",pch=0,yaxt="n",xlab="")
  }
  text(median(year_seq[year_clip]),1.15,group_names[j],font=2)
  axis(2,at=c(0.25,0.5,0.75,1))
  for(i in 1:datamcmc2$Nregs)
  {
    spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_run",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
    spawn_forecast_1 <- as.matrix(spawn_forecast_1)
    spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
    spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
    spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<baselines[x])}))
    spawn_fore_risk <- colSums(spawn_fore_risk[group==j,])/(n_chains*samples*sum(group==j))
    spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
    spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
    spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
    spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
    spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
    group_mdn <- rowSums(spawn_forecast_arr[,group==j,1])/sum(baselines[group==j])
    group_ui <- rowSums(spawn_forecast_arr[,group==j,2])/sum(baselines[group==j])
    group_li <- rowSums(spawn_forecast_arr[,group==j,3])/sum(baselines[group==j])
    ylims <- range(group_mdn,group_ui,group_li)
    polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
    #lines(year_seq[year_clip],spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
    abline(h=1,lwd=1,lty=2,col="red")
  }
}
axis(1)

for(j in 1:Kgroups)
{
  year_seq <- c(years,max(years)+1:sim_eval)
  year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
  par(mar=c(0,2,0,0)) 
  if(j==3)
  {
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,3),ylab="",xaxt="n",pch=0,yaxt="n",xpd=NA,xlab="")
    mtext("Spawners / 0.8S(msy)",side=2,line=2.5)
  }else{
    plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,3),ylab="",xaxt="n",pch=0,yaxt="n",xlab="")
  }
  text(median(year_seq[year_clip]),2.85,group_names[j],font=2)
  axis(2,at=c(0.5,1,1.5,2,2.5))
  for(i in 1:datamcmc2$Nregs)
  {
    spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
    spawn_forecast_1 <- as.matrix(spawn_forecast_1)
    spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
    spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
    spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<(0.8*Smsy_priors$mean[x]))}))
    spawn_fore_risk <- colSums(spawn_fore_risk[group==j,])/(n_chains*samples*sum(group==j))
    
    spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
    spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
    spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
    spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
    spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
    group_mdn <- rowSums(spawn_forecast_arr[,group==j,1])/(0.8*sum(Smsy_priors$mean[group==j]))
    group_ui <- rowSums(spawn_forecast_arr[,group==j,2])/(0.8*sum(Smsy_priors$mean[group==j]))
    group_li <- rowSums(spawn_forecast_arr[,group==j,3])/(0.8*sum(Smsy_priors$mean[group==j]))
    ylims <- range(group_mdn,group_ui,group_li)
    polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
    #lines(year_seq[year_clip],spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
    abline(h=1,lwd=1,lty=2,col="red")
  }
}
axis(1)

year_seq <- c(years,max(years)+1:sim_eval)
year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
par(mar=c(4.5,4,3,2))
plot(0,xlim=range(year_seq[year_clip]),ylim=c(0,1.25),ylab="Relative run size (to 2000-15 baseline)",pch=0,yaxt="n",xlab="Year")
axis(2,at=c(0.25,0.5,0.75,1))

for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_run",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  
  spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<baselines[x])}))
  spawn_fore_risk <- colSums(spawn_fore_risk)/(n_chains*samples*n.pops)
  
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  group_mdn <- rowSums(spawn_forecast_arr[,,1])/sum(baselines)
  group_ui <- rowSums(spawn_forecast_arr[,,2])/sum(baselines)
  group_li <- rowSums(spawn_forecast_arr[,,3])/sum(baselines)
  ylims <- range(group_mdn,group_ui,group_li)
  polygon(c(year_seq[year_clip],rev(year_seq[year_clip])),c(group_ui[year_clip],rev(group_li[year_clip])),col=adjustcolor(i,0.25),border=NA)
  #lines(year_seq[year_clip],spawn_fore_risk,type="l",lty=1,lwd=2,col=i)
  abline(h=1,lwd=1,lty=2,col="red")
  text(x=year_seq[length(year_seq)]-1,y=spawn_forecast_ui[length(year_seq)],LETTERS[i],font=2)
}

year_seq <- c(years,max(years)+1:sim_eval)
year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
par(mar=c(4.5,4,3,2))
boxplot(NA,xlim=c(0.5,5.5),ylim=c(0,3),ylab="Spawners / 0.80 S(msy)",xaxt="n",col=0,yaxt="n",xlab="Management scenario")
for(i in 1:datamcmc2$Nregs)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Npops))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  
  spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<(0.8*Smsy_priors$mean[x]))}))
  spawn_fore_risk <- rowSums(spawn_fore_risk)/(n_chains*samples*sim_eval)
  
  spawn_forecast_mdn <- sapply(1:n.pops,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  spawn_forecast <- rbind(spawn_ci[3,,],spawn_fore_ci[3,,])
  spawn_forecast_ui <- rbind(spawn_ci[4,,],spawn_fore_ci[4,,])
  spawn_forecast_li <- rbind(spawn_ci[2,,],spawn_fore_ci[2,,])
  spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
  group_mdn <- rowSums(spawn_forecast_arr[year_clip,,1])/(0.8*sum(Smsy_priors$mean))
  group_ui <- rowSums(spawn_forecast_arr[year_clip,,2])/(0.8*sum(Smsy_priors$mean))
  group_li <- rowSums(spawn_forecast_arr[year_clip,,3])/(0.8*sum(Smsy_priors$mean))
  ylims <- range(group_mdn,group_ui,group_li)
  boxplot(c(group_ui,group_li,group_mdn),at=i,add=TRUE,col=adjustcolor(i,0.5))
  abline(h=1,lwd=1,lty=2,col="red")
  text(x=i,y=-0.9,LETTERS[i],font=2,xpd=NA)
  text(x=i,y=1.15*max(group_ui),round(mean(group_mdn),1),font=2,xpd=NA)
}

legend(x=3,y=6,legend=colnames(escape_reg),title="Management scenario",bty="n",pch=22,pt.bg=adjustcolor(1:4,0.9),xpd=NA)

dev.off()

100*(risk_2030-sum(status$s_over_msy<1)/nrow(status))/(sum(status$s_over_msy<1)/nrow(status))
