Smsy_priors <- readRDS("Results/Smsy.rds")
Sgen_priors <- readRDS("Results/Sgen.rds")
Umsy_priors <- readRDS("Results/Umsy.rds")
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


co_pops<-read.table("Data/coho_groups.txt",header=TRUE)
co_pops$mean_total<-NA

for(i in co_pops$pop_no){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6
group <- co_pops$group
Npops <- n.pops <- length(group)
SR.dat_full$group <- group[match(SR.dat_full$pop_no,co_pops$pop_no)]
baselines <- sapply(1:n.pops,function(x){mean(SR.dat$total_runE[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines_2 <- sapply(1:n.pops,function(x){mean(SR.dat$total_run2[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines<-ifelse(!is.na(baselines_2),baselines_2,baselines)

baseline_spawn <- sapply(1:n.pops,function(x){mean(SR.dat$escapement[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})

s_over_msy <- sapply(1:n.pops,function(x){sum(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]<=(Smsy_priors$mean[x]),na.rm=TRUE)/sum(!is.na(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]))})

SR.dat_full$total_run <- ifelse(!is.na(SR.dat_full$total_run2),SR.dat_full$total_run2,SR.dat_full$total_runE)
pop_sub <- aggregate(total_run~pop_no,SR.dat_full[SR.dat_full$year>=2017,],mean,na.rm=TRUE)
run_over_baseline <- 100*((sum(pop_sub$total_run)-sum(baselines[pop_sub$pop_no]))/sum(baselines[pop_sub$pop_no])) # CJFAS revision - declines in total run size (Discussion)

pop_sub <- aggregate(escapement~pop_no,SR.dat_predict,mean,na.rm=TRUE)
escapement_over_baseline <- 100*((sum(pop_sub$escapement)-sum(baseline_spawn[pop_sub$pop_no]))/sum(baseline_spawn[pop_sub$pop_no])) # CJFAS revision - declines in total escapement (Discussion)


SR.dat_full$baseline_run <- (baselines)[SR.dat_full$pop_no]
SR.dat_full$baseline_spawn <- (baseline_spawn)[SR.dat_full$pop_no]
SR.dat_full$lrp <- (Sgen_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$usr <- 0.8*(Smsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$ucur <- ifelse(!is.na(SR.dat_full$er_2),SR.dat_full$er_2,SR.dat_full$er_E)
SR.dat_full$umsy <- (Umsy_priors$mean)[SR.dat_full$pop_no]
SR.dat_full$regime <- ifelse(SR.dat_full$year>=1980 & SR.dat_full$year<=1996,"early",ifelse(SR.dat_full$year<=2016,"recovering","recent"))
pop_sub <- aggregate(cbind(total_run/baseline_run,escapement/lrp,escapement/usr,escapement/baseline_spawn,ucur/umsy)~pop_no+regime,SR.dat_full,mean,na.rm=F)
colnames(pop_sub) <- c("pop_no","regime","run_ratio","lrp_ratio","usr_ratio","baseline_ratio","umsy_ratio")
recent <- pop_sub[pop_sub$regime=="recent",]
early <- pop_sub[pop_sub$regime=="early",]
recovering <- pop_sub[pop_sub$regime=="recovering",]
time_series <- merge(merge(merge(early,recovering,by="pop_no"),recent,by="pop_no"),current,by="pop_no")
time_series$pop <- pop_names[time_series$pop_no]

risk_fn <- function(x)
  {
  return(list("run_ratio"=sum(x$run_ratio<=1)/nrow(x),"lrp_ratio"=sum(x$lrp_ratio<=1)/nrow(x),"usr_ratio"=sum(x$usr_ratio<=1)/nrow(x),"baseline_ratio"=sum(x$baseline_ratio<=1)/nrow(x),"umsy_ratio"=sum(x$umsy_ratio>=1)/nrow(x),"pops"=nrow(x)))
  }

risk_fn(recent)

pop_sub_full <- aggregate(cbind(total_run/baseline_run,escapement/lrp,escapement/usr,escapement/baseline_spawn,ucur/umsy)~pop_no+regime+year,SR.dat_full,mean,na.rm=T)
colnames(pop_sub_full) <- c("pop_no","regime","year","run_ratio","lrp_ratio","usr_ratio","baseline_ratio","umsy_ratio")

risk_fn(pop_sub_full[pop_sub_full$year>=2017,])
risk_fn(recovering)
risk_fn(early)

rel_risk_fn <- function(x)
{
  return(list("run_ratio"=mean(x$run_ratio,na.rm=TRUE),"lrp_ratio"=mean(x$lrp_ratio,na.rm=TRUE),"usr_ratio"=mean(x$usr_ratio,na.rm=TRUE),"baseline_ratio"=mean(x$baseline_ratio,na.rm=TRUE),"umsy_ratio"=mean(x$umsy_ratio,na.rm=TRUE),"pops"=nrow(x)))
}
rel_risk_fn(pop_sub_full[pop_sub_full$year>=2017,])
rel_risk_fn(recent)
rel_risk_fn(recovering)
rel_risk_fn(early)
