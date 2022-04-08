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

escape_reg <- readRDS("Data/harvest scenarios.rds")
escape_sd_reg <- readRDS("Data/harvest var scenarios.rds")

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

#### plot population dynamics

mcmc_names <- colnames(result_forecast_run[[1]])
prob_fore <- c(0.025,0.25,0.5,0.75,0.975)

spawn_post <- result_forecast_run[,grepl("\\bspawners",mcmc_names)]
spawn_post[[1]][1:10,1:datamcmc2$sim_eval]
spawn_post <- as.matrix(spawn_post)
spawn_post[1,1:datamcmc2$sim_eval]
spawn_post <- array(unlist((spawn_post)),dim=c(n_chains*samples,datamcmc2$n.years,datamcmc2$Npops))
spawn_ci <- apply(spawn_post,c(2,3),quantile,probs=prob_fore)
rmse <- below_lrp <- below_msy <- matrix(NA,nrow=5,ncol=n.pops,dimnames=list("Scenarios"=colnames(escape_reg),"Populations"=pop_names))
pop_plot <- which(pop_names=="kasiks")
jpeg(filename="Figures/example population dynamics.jpeg", width=7,height=5, units="in", res=600)

layout(matrix(1,ncol=1,byrow=TRUE))
par(mar=c(5,4,0.5,0.5))
years<-seq(1980,2016,1)
year_seq <- c(years,max(years)+1:sim_eval)
manage_cols <- c("saddlebrown","orange","purple","green4","dodgerblue")
for(j in pop_plot)
{
  spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",j,",",sep=""),mcmc_names)]
  spawn_forecast_1 <- as.matrix(spawn_forecast_1)
  spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,datamcmc2$sim_eval,datamcmc2$Nregs))
  spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
  spawn_forecast_mdn <- sapply(1:datamcmc2$Nregs,function(x){spawn_forecast_1[which.min(abs(spawn_fore_ci[3,datamcmc2$sim_eval,x]-spawn_forecast_1[,datamcmc2$sim_eval,x])),,x]})
  ylims <- range(c(spawn_ci[,,j], spawn_fore_ci[2:4,,]))
  plot(years,spawn_ci[3,,j],type="l",xlim=range(year_seq),ylim=ylims,lty=1,lwd=0.5,col="black",ylab="Spawners",xlab="Years")
  polygon(c(years,rev(years)),c(spawn_ci[4,,j],rev(spawn_ci[2,,j])),col=adjustcolor("black",0.5),border=NA)
  lines(years,spawn_ci[3,,j],type="l",lty=1,lwd=0.5,col="black")
  subdat_predict <- subset(SR.dat_predict,SR.dat_predict$pop_no==j)
  for(i in c(1,2,5))
  {
    year_sim <- ((max(years)-1)+1:(sim_eval+1))
    year_sel <- (length(years)-1)+1:(sim_eval+1)
    spawn_forecast <- c(spawn_ci[3,,j],spawn_fore_ci[3,,i])
    spawn_forecast_ui <- c(spawn_ci[4,,j],spawn_fore_ci[4,,i])
    spawn_forecast_li <- c(spawn_ci[2,,j],spawn_fore_ci[2,,i])
    spawn_forecast_arr <- simplify2array(list(spawn_forecast,spawn_forecast_ui,spawn_forecast_li))
    polygon(c(year_sim,rev(year_sim)),c(spawn_forecast_arr[year_sel,2],rev(spawn_forecast_arr[year_sel,3])),col=adjustcolor(manage_cols[i],0.5),border=NA)
    lines(year_sim,spawn_forecast_arr[year_sel,1],type="l",lty=1,lwd=0.5,col=manage_cols[i])
  }
  points(years,spawn[,j],pch=21,bg="grey90")
  points(subdat_predict$year,subdat_predict$escapement,pch=22,bg="dodgerblue",cex=1.2)
  legend("topright",colnames(escape_reg)[c(1,2,5)],pch=22,pt.bg=adjustcolor(manage_cols[c(1,2,5)],0.5),bty="n",title=paste("Management scenarios, ",pop_names[j]," river",sep=""))
  abline(h=0.8*Smsy_priors$mean[j],lwd=1,lty=1,col="black")
  abline(h=Sgen_priors$mean[j],lwd=1,lty=2,col="red")
}
dev.off()