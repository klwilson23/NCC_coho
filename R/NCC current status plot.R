###########################################################################
###     Hierarchical SR models w/time variant alphas, PR priors        ####
###########################################################################
Corner_text <- function(text, location="topright",...)
{
  legend(location,legend=text, bty ="n", pch=NA,...)
}
#setwd("C:/data/centralcoast")
require(rjags)
require(R2jags)
library(runjags)
require(gsl)
library(tidyverse)
library(ggplot2)

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

ER_values <- readRDS("Data/harvest scenarios.rds")

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

for(i in co_pops$pop_no){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6

### using mean run size (harvest scenario 1) across the timeseries as our prior on capacity
n.pops<-max(SR.dat$pop_no)

Smax.p<-rep(NA,n.pops)
Smax.tau<-rep(NA,n.pops)

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

Smsy_priors <- readRDS("Results/Smsy.rds")
Sgen_priors <- readRDS("Results/Sgen.rds")
Umsy_priors <- readRDS("Results/Umsy.rds")

alpha <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[j]]
beta <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[j]]
curve(exp(alpha)*x*exp(-beta*x),from=0,to=2*alpha/beta,xlab="Spawners",ylab="Recruits")
abline(v=(1-Umsy_priors$mean[j])*alpha/beta)
abline(a=0,b=1,lty=2)
abline(v=Smsy_priors$mean[j])

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
SR.dat_full$usr <- 0.8*(Smsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$lrp <- (Sgen_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$ucur <- ifelse(!is.na(SR.dat_full$er_2),SR.dat_full$er_2,SR.dat_full$er_E)
SR.dat_full$umsy <- (Umsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$recent <- ifelse(SR.dat_full$year>=1980 & SR.dat_full$year<=1996,"early",ifelse(SR.dat_full$year<=2016,"recovering","recent"))
SR.dat_full$regime <- ifelse(SR.dat_full$year>=1980 & SR.dat_full$year<=1996,"early",ifelse(SR.dat_full$year<=2016,"recovering","recent"))

nsamps <- sapply(1:n.pops,function(x){sum(!is.na(SR.dat_full$escapement[SR.dat_full$pop_no==x]))})
n_rec <- sapply(1:n.pops,function(x){sum(!is.na(SR.dat_full$escapement[SR.dat_full$pop_no==x & SR.dat_full$regime=='recent']))})

SR.dat_full$n <- nsamps[SR.dat_full$pop_no]
SR.dat_full$n_rec <- n_rec[SR.dat_full$pop_no]

# make stacked bar plot

pop_sub <- aggregate(cbind(total_run/baseline_run,escapement/lrp,escapement/usr,ucur/umsy)~population+pop_no+regime,SR.dat_full,mean,na.rm=F)
pop_sub$regime <- factor(pop_sub$regime,levels=c("early","recovering","recent"),labels=c("1980-1996","1997-2017","Since 2017"))
pop_sub$mean_escapement <- co_pops$mean_total[match(pop_sub$population,co_pops$population)]
pop_sub$region_no <- co_pops$group[match(pop_sub$population,co_pops$population)]
pop_sub$Region[pop_sub$region_no==1] <- "Haida Gwaii"
pop_sub$Region[pop_sub$region_no==2] <- "Nass"
pop_sub$Region[pop_sub$region_no==3] <- "Skeena"
pop_sub$Region[pop_sub$region_no==4] <- "Hecate Lowlands"
pop_sub$Region[pop_sub$region_no==5] <- "Inner Waters"
pop_sub$Region[pop_sub$region_no==6] <- "Central Coast (South)"
pop_sub$Year <- pop_sub$year

pop_sub$Region <- factor(pop_sub$Region,levels=c("Haida Gwaii","Nass","Skeena","Hecate Lowlands","Inner Waters","Central Coast (South)"))
pop_sub$status <- ifelse(pop_sub$V2<1,"Below LRP",ifelse(pop_sub$V3<1,"Below USR",ifelse(pop_sub$V1<1,"Below Historical Baseline","Above Historical Baseline")))
pop_sub$status <- factor(pop_sub$status,levels=c("Above Historical Baseline","Below Historical Baseline","Below USR","Below LRP"))

co_pops$Region[co_pops$group==1] <- "Haida Gwaii"
co_pops$Region[co_pops$group==2] <- "Nass"
co_pops$Region[co_pops$group==3] <- "Skeena"
co_pops$Region[co_pops$group==4] <- "Hecate Lowlands"
co_pops$Region[co_pops$group==5] <- "Inner Waters"
co_pops$Region[co_pops$group==6] <- "Central Coast (South)"
co_pops$Region <- factor(co_pops$Region,levels=c("Haida Gwaii","Nass","Skeena","Hecate Lowlands","Inner Waters","Central Coast (South)"))
co_pops2 <- expand.grid("pop_no"=co_pops$pop_no,"regime"=c("1980-1996","1997-2017","Since 2017"))
co_pop2 <- full_join(co_pops,co_pops2,by=c("pop_no"))

co_pop2 <- full_join(co_pop2,pop_sub,by=c("pop_no","Region","regime"))
co_pop2$status <- as.character(co_pop2$status)
co_pop2$status[is.na(co_pop2$status)] <- "Data Deficient"
co_pop2$status <- factor(co_pop2$status,levels=c("Above Historical Baseline","Below Historical Baseline","Below USR","Below LRP","Data Deficient"),labels=c("Healthy","Cautious (Baseline)","Cautious (USR)","Critical (LRP)","Data Deficient"))

co_pop2$value <- 1


ggplot(co_pop2, aes(fill=status,y=value,x=Region)) + 
  geom_bar(position="stack",stat="identity")+
  ylab(expression(N[t]/N[RP])) +
  facet_wrap(~regime,ncol=1) +
  scale_fill_manual(name="Status",values=c("darkgreen","orange4","orange","red4","grey60")) +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(strip.text = element_text(hjust = 0),legend.text = element_text(size=7),legend.title = element_text(size=8))
ggsave("Figures/Relative status across regions and time.jpeg", width = 7, height=6,units="in", dpi=600)


recent <- pop_sub[pop_sub$recent=="recent",]
early <- pop_sub[pop_sub$recent=="early",]
recovering <- pop_sub[pop_sub$recent=="recovering",]
time_series <- merge(merge(early,recovering,by="pop_no"),recent,by="pop_no")
time_series$pop <- pop_names[time_series$pop_no]
jpeg(filename="Figures/current status tr1.jpeg", width=4.5,height=4.5, units="in", res=600)
par(mar=c(4,4,0.1,0.1))
status <- data.frame("pop"=pop_names[recent$pop_no],"s_over_msy"=recent$V3,"baseline"=recent$V1,"u_over_umsy"=recent$V4)
plot(NA,ylab=expression(u[cur]/u[MSY]),xlab=expression(S[cur]/(0.8*S[MSY])),xlim=1.1*range(c(0,2,status$s_over_msy)),ylim=1.1*range(c(0,2,status$u_over_umsy)))
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
dev.off()

jpeg(filename="Figures/time-series of status tr1.jpeg", width=4.5,height=4.5, units="in", res=600)
sub_series <- time_series[c(which.max(time_series$V4.x),which.max(time_series$V3.x-time_series$V3),which.max(time_series$V4.y-time_series$V4),which.max(time_series$V3-time_series$V3.y)),]
par(mar=c(4,4,0.1,0.1))
plot(NA,ylab=expression(u[t]/u[MSY]),xlab=expression(S[t]/S[MSY]),xlim=1.1*range(c(0,2,sub_series$V3.x,sub_series$V3.y,sub_series$V3)),ylim=1.1*range(c(0,2,sub_series$V4.x,sub_series$V4.y,sub_series$V4)))
polygon(y=c(-5,1,1,-5),x=c(1,1,40,40),col=adjustcolor("seagreen",1),border=NA,xpd=FALSE)
polygon(y=c(1,1,20,20),x=c(-5,1,1,-5),col=adjustcolor("tomato",1),border=NA,xpd=FALSE)
polygon(y=c(1,20,20,1),x=c(1,1,40,40),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
polygon(y=c(-5,1,1,-5),x=c(-5,-5,1,1),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
abline(h=1,v=1,lty=2,lwd=0.5,col="grey10")
shape::Arrows(y0=sub_series$V4.x,y1=sub_series$V4.y,x1=sub_series$V3.y,x0=sub_series$V3.x,arr.adj=1,arr.length = 0.3)
shape::Arrows(y0=sub_series$V4.y,y1=sub_series$V4,x1=sub_series$V3,x0=sub_series$V3.y,arr.adj=1,arr.length = 0.3)
points(sub_series$V3.x,sub_series$V4.x,pch=21,bg="black",col="white")
points(sub_series$V3.y,sub_series$V4.y,pch=21,bg="dodgerblue",col="white")
points(sub_series$V3,sub_series$V4,pch=21,bg="grey80",col="white")
samp <- 1:4
text(sub_series$V3[samp],sub_series$V4[samp],labels=sub_series$pop[samp],adj=c(1,-1),offset=1.5,cex=0.7,col="black")
legend("topright",c("'90s","'00s","'10s"),pch=21,col="white",pt.bg=c("black","dodgerblue","grey80"),bty="n",bg=NA,title="Time period",cex=0.8)
dev.off()

sub_series <- time_series
par(mar=c(4,4,0.1,0.1))
plot(NA,ylab=expression(u[t]/u[MSY]),xlab=expression(S[t]/S[MSY]),xlim=1.1*range(c(0,2,sub_series$V3.x,sub_series$V3.y,sub_series$V3)),ylim=1.1*range(c(0,2,sub_series$V4.x,sub_series$V4.y,sub_series$V4)))
polygon(y=c(-5,1,1,-5),x=c(1,1,40,40),col=adjustcolor("seagreen",1),border=NA,xpd=FALSE)
polygon(y=c(1,1,20,20),x=c(-5,1,1,-5),col=adjustcolor("tomato",1),border=NA,xpd=FALSE)
polygon(y=c(1,20,20,1),x=c(1,1,40,40),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
polygon(y=c(-5,1,1,-5),x=c(-5,-5,1,1),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
abline(h=1,v=1,lty=2,lwd=0.5,col="grey10")
shape::Arrows(y0=sub_series$V4.x,y1=sub_series$V4.y,x1=sub_series$V3.y,x0=sub_series$V3.x,arr.adj=1,arr.length = 0.3)
shape::Arrows(y0=sub_series$V4.y,y1=sub_series$V4,x1=sub_series$V3,x0=sub_series$V3.y,arr.adj=1,arr.length = 0.3)
points(sub_series$V3.x,sub_series$V4.x,pch=21,bg="black",col="white")
points(sub_series$V3.y,sub_series$V4.y,pch=21,bg="dodgerblue",col="white")
points(sub_series$V3,sub_series$V4,pch=21,bg="grey80",col="white")
#text(sub_series$V3[samp],sub_series$V1[samp],labels=sub_series$pop[samp],adj=c(-0.1,1),offset=1.5,cex=0.7,col="black")
legend("topright",c("'90s","'00s","'10s"),pch=21,col="white",pt.bg=c("black","dodgerblue","grey80"),bty="n",bg=NA,title="Time period",cex=0.8)

SR.dat_subby <- SR.dat_full[SR.dat_full$n_rec>=2 & SR.dat_full$n>=5,]
pop_sub <- aggregate(cbind(total_run/baseline_run,escapement/lrp,escapement/usr,ucur/umsy)~pop_no+regime,SR.dat_subby,mean,na.rm=F)
colnames(pop_sub) <- c("pop_no","regime","run_ratio","lrp","msy","u_ratio")
recent <- pop_sub[pop_sub$regime=="recent",]
early <- pop_sub[pop_sub$regime=="early",]
recovering <- pop_sub[pop_sub$regime=="recovering",]
time_series <- merge(merge(early,recovering,by="pop_no"),recent,by="pop_no")
time_series$pop <- pop_names[time_series$pop_no]
pop_sub <- aggregate(cbind(total_run/baseline_run,escapement/lrp,escapement/usr,ucur/umsy)~pop_no+recent,SR.dat_subby,mean,na.rm=F)
colnames(pop_sub) <- c("pop_no","recent","run_ratio","lrp","msy","u_ratio")
recent <- pop_sub[pop_sub$recent=="recent",]

jpeg(filename="Figures/current status and time-series tr1.jpeg", width=4.5,height=9, units="in", res=600)
layout(matrix(1:2,nrow=2,ncol=1))
par(mar=c(4,4,0.1,0.1))
status <- data.frame("pop"=pop_names[recent$pop_no],"s_over_msy"=recent$msy,"baseline"=recent$run_ratio,"u_over_umsy"=recent$u_ratio)
plot(NA,ylab=expression(u[cur]/u[MSY]),xlab=expression(S[cur]/0.8*S[MSY]),xlim=1.1*range(c(0,2,status$s_over_msy,3.25)),ylim=1.1*range(c(0,2,status$u_over_umsy)))
polygon(y=c(-5,1,1,-5),x=c(1,1,40,40),col=adjustcolor("seagreen",1),border=NA,xpd=FALSE)
polygon(y=c(1,1,20,20),x=c(-5,1,1,-5),col=adjustcolor("tomato",1),border=NA,xpd=FALSE)
polygon(y=c(1,20,20,1),x=c(1,1,40,40),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
polygon(y=c(-5,1,1,-5),x=c(-5,-5,1,1),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
abline(h=1,v=1,lty=2,lwd=0.5,col="grey10")
points(status$s_over_msy,status$u_over_umsy,pch=21,bg="grey80",col="white")
samp <- sapply(c(0,0.5,1,2),function(x){which.min(abs(status$s_over_msy-x))})
text(status$s_over_msy[samp],status$u_over_umsy[samp],labels=status$pop[samp],adj=c(0,0),offset=1.5,cex=0.7,col="black")
Corner_text("(a)","topleft",cex=1)

sub_series <- time_series[c(which.max(time_series$msy),which.min(time_series$msy),which.max(time_series$u_ratio),which.max(time_series$u_ratio.x)),]
plot(NA,ylab=expression(u[t]/u[MSY]),xlab=expression(S[t]/0.8*S[MSY]),xlim=1.1*range(c(0,2,sub_series$msy.x,sub_series$msy.y,sub_series$msy,3.25)),ylim=1.1*range(c(0,2,sub_series$u_ratio.x,sub_series$u_ratio.y,sub_series$u_ratio)))
polygon(y=c(-5,1,1,-5),x=c(1,1,40,40),col=adjustcolor("seagreen",1),border=NA,xpd=FALSE)
polygon(y=c(1,1,20,20),x=c(-5,1,1,-5),col=adjustcolor("tomato",1),border=NA,xpd=FALSE)
polygon(y=c(1,20,20,1),x=c(1,1,40,40),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
polygon(y=c(-5,1,1,-5),x=c(-5,-5,1,1),col=adjustcolor("orange",1),border=NA,xpd=FALSE)
abline(h=1,v=1,lty=2,lwd=0.5,col="grey10")
shape::Arrows(y0=sub_series$u_ratio.x,y1=sub_series$u_ratio.y,x1=sub_series$msy.y,x0=sub_series$msy.x,arr.adj=1,arr.length = 0.3)
shape::Arrows(y0=sub_series$u_ratio.y,y1=sub_series$u_ratio,x1=sub_series$msy,x0=sub_series$msy.y,arr.adj=1,arr.length = 0.3)
points(sub_series$msy.x,sub_series$u_ratio.x,pch=21,bg="black",col="white")
points(sub_series$msy.y,sub_series$u_ratio.y,pch=21,bg="dodgerblue",col="white")
points(sub_series$msy,sub_series$u_ratio,pch=21,bg="grey80",col="white")
samp <- 1:4
text(sub_series$msy[samp],sub_series$u_ratio[samp],labels=sub_series$pop[samp],adj=c(0,1),offset=1.5,cex=0.7,col="black")
legend("topright",c("'90s","'00s","'10s"),pch=21,col="white",pt.bg=c("black","dodgerblue","grey80"),bty="n",bg=NA,title="Time period",cex=0.8)
Corner_text("(b)","topleft",cex=1)

dev.off()

pop_sub <- aggregate(cbind(escapement/baseline_spawn,escapement/lrp,escapement/usr,ucur/umsy)~pop_no+recent,SR.dat_full,mean,na.rm=F)
colnames(pop_sub) <- c("pop_no","recent","baseline","lrp","msy","u_ratio")
pop_sub$pop <- pop_names[pop_sub$pop_no]
pop_sub <- pop_sub[pop_sub$recent=='recent',]
sum(pop_sub$baseline<=1 & pop_sub$msy<=1)
sum((pop_sub$baseline<=1 & pop_sub$msy>1) | (pop_sub$baseline>1 & pop_sub$msy<=1))
sum(pop_sub$baseline>1 & pop_sub$msy>1)
nrow(pop_sub)
