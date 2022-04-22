library(ggplot2)
Smsy_priors <- readRDS("Results/Smsy.rds")
Sgen_priors <- readRDS("Results/Sgen.rds")
Umsy_priors <- readRDS("Results/Umsy.rds")
escape_reg <- readRDS("Data/harvest scenarios.rds")
escape_sd_reg <- readRDS("Data/harvest var scenarios.rds")
group_names <- c("Central Coast (South)","Hecate Lowlands","Inner Waters","Haida Gwaii","Skeena","Nass")
sim_eval <- 20 # forecast model for 16 years, or 4 generations
thinning <- 200
samples <- 1500
n_chains <- 4
group_names <- group_names[c(4,6,5,2,3,1)]
group_names <- c(group_names,"Coast-wide")
# plot boxplot
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat_full <- SR.dat
SR.dat_predict<-subset(SR.dat,SR.dat$year>=2017)
SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)

pop_names <- unique(SR.dat$population[order(SR.dat$pop_no)])

co_pops<-read.table("Data/coho_groups.txt",header=TRUE)
co_pops$mean_total<-NA
years<-seq(1980,2016,1)

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
baseline_spawn <- sapply(1:n.pops,function(x){mean(SR.dat$escapement[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})

prob_fore <- c(0.025,0.25,0.5,0.75,0.975)

box_fn <- function(mod="random",year_sim=sim_eval,ref_point_label="S(gen)",ref_point=Sgen_priors$mean)
{
  mod_lab <- ifelse(mod=="random","",ifelse(mod=="tr1","_tr1","_ar1"))
  result_forecast_metrics <- readRDS(paste("Results/coho_forecasting",mod_lab,".rds",sep=""))
  result_forecast_run <- readRDS(paste("Results/coho_forecasting_popdyn",mod_lab,".rds",sep=""))
  result_forecast_RS <- readRDS(paste("Results/coho_forecasting_rec",mod_lab,".rds",sep=""))
  mcmc_names <- colnames(result_forecast_run[[1]])
  df_full <- NULL
  matLayout <- matrix(1:7,nrow=7,ncol=1,byrow=FALSE)
  layout(matLayout)
  par(mar=c(4.5,4,0,0))
  max_box <- 110
  for(j in 1:(Kgroups+1))
  {
    year_seq <- c(years,max(years)+1:sim_eval)
    year_clip <- (length(year_seq)-(sim_eval-1)):length(year_seq)
    if(j==4 & mod=="random")
    {
      boxplot(NA,xlim=c(0.5,5.5),ylim=c(0,max_box),ylab="",xaxt="n",col=0,yaxt="n",xlab="",xpd=NA)
      mtext(paste("% spawners below",ref_point_label,sep=""),side=2,line=2.5)
    }else{
      if(j==7)
      {
        boxplot(NA,xlim=c(0.5,5.5),ylim=c(0,max_box),ylab="",xaxt="n",col=0,yaxt="n",xlab="Management scenario",xpd=NA)
      }else{
        boxplot(NA,xlim=c(0.5,5.5),ylim=c(0,max_box),ylab="",xaxt="n",col=0,yaxt="n",xpd=NA)
      }
    }
    text(3,95,group_names[j],font=2)
    axis(2,at=c(20,40,60,80))
    if(j<7)
    {
      for(i in 1:5)
        {
        spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
        spawn_forecast_1 <- as.matrix(spawn_forecast_1)
        spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,sim_eval,n.pops))
        spawn_forecast_1 <- spawn_forecast_1[,1:year_sim,]
        spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
        spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<ref_point[x])}))
        #spawn_fore_risk <- colSums(spawn_fore_risk[group==j,])/(n_chains*samples*sum(group==j))
        spawn_fore_risk <- rowSums(spawn_fore_risk[group==j,])/(n_chains*samples*year_sim)
        
        #spawn_fore_ratio <- apply(t(apply(spawn_forecast_1[,,group==j],c(1,2),sum))/sum(ref_point[group==j]),1,median)
        spawn_fore_ratio <- apply(t(sapply(1:n.pops,function(x){spawn_forecast_1[,,x]/ref_point[x]}))[group==j,],1,median)
        
       # spawn_fore_ratio <- median(apply(spawn_forecast_1[,,group==j],c(1,2,3),FUN=function(x){sum(x)/sum(ref_point[group==j])}))
        
        
        df <- data.frame("Risk"=spawn_fore_risk,"Region"=group_names[j],"Regulation"=colnames(escape_reg)[i],"Population"=co_pops$population[co_pops$group==j],"Ratio"=spawn_fore_ratio,"Reference Point"=ref_point_label)
        df_full <- rbind(df_full,df)
        
        boxplot(100*spawn_fore_risk,at=i,add=TRUE,col=adjustcolor(i,0.5),yaxt="n")
        #abline(h=1,lwd=1,lty=1,col="black")
        text(x=i,y=1.15*max(100*spawn_fore_risk),round(mean(100*spawn_fore_risk),0),font=2,xpd=NA)
      }
    }
    if(j==7)
      {
      for(i in 1:5)
        {
        spawn_forecast_1 <- result_forecast_run[,grepl("\\bnew_spawners",mcmc_names) & grepl(paste("\\b,",i,"]",sep=""),mcmc_names)]
        spawn_forecast_1 <- as.matrix(spawn_forecast_1)
        spawn_forecast_1 <- array(unlist((spawn_forecast_1)),dim=c(n_chains*samples,sim_eval,n.pops))
        spawn_forecast_1 <- spawn_forecast_1[,1:year_sim,]
        
        spawn_fore_ci <- apply(spawn_forecast_1,c(2,3),quantile,probs=prob_fore)
          
        spawn_fore_risk <- t(sapply(1:n.pops,function(x){colSums(spawn_forecast_1[,,x]<(ref_point[x]))}))
        spawn_fore_risk <- rowSums(spawn_fore_risk)/(n_chains*samples*year_sim)
        spawn_fore_ratio <- apply(t(sapply(1:n.pops,function(x){spawn_forecast_1[,,x]/ref_point[x]})),1,median)
        
        df <- data.frame("Risk"=spawn_fore_risk,"Region"=group_names[j],"Regulation"=colnames(escape_reg)[i],"Population"=pop_names,"Ratio"=spawn_fore_ratio,"Reference Point"=ref_point_label)
        df_full <- rbind(df_full,df)
        
        boxplot(100*spawn_fore_risk,at=i,add=TRUE,col=adjustcolor(i,0.5),yaxt="n")
        #abline(h=1,lwd=1,lty=1,col="black")
        text(x=i,y=-20,LETTERS[i],font=2,xpd=NA)
        text(x=i,y=1.15*max(100*spawn_fore_risk),round(mean(100*spawn_fore_risk),0),font=2,xpd=NA)
      }
    }
  }
  return(df_full)
}

ran_df <- box_fn(mod="random",year_sim=sim_eval,ref_point_label="S(gen)",ref_point=Sgen_priors$mean)
ran_df$Model <- "Random walk"
ar_df <- box_fn(mod="ar1",year_sim=sim_eval,ref_point_label="S(gen)",ref_point=Sgen_priors$mean)
ar_df$Model <- "Mean reverting"
tr_df <- box_fn(mod="tr1",year_sim=sim_eval,ref_point_label="S(gen)",ref_point=Sgen_priors$mean)
tr_df$Model <- "Trending"

df <- rbind(ran_df,ar_df,tr_df)
df$Region <- factor(df$Region,levels=group_names[c(7,1:6)])
df$Regulation <- factor(df$Regulation,levels=colnames(escape_reg)[c(2,4,3,5,1)])

ran_df <- box_fn(mod="random",year_sim=sim_eval,ref_point_label="S(MSY)",ref_point=0.8*Smsy_priors$mean)
ran_df$Model <- "Random walk"
ar_df <- box_fn(mod="ar1",year_sim=sim_eval,ref_point_label="S(MSY)",ref_point=0.8*Smsy_priors$mean)
ar_df$Model <- "Mean reverting"
tr_df <- box_fn(mod="tr1",year_sim=sim_eval,ref_point_label="S(MSY)",ref_point=0.8*Smsy_priors$mean)
tr_df$Model <- "Trending"

df_msy <- rbind(ran_df,ar_df,tr_df)
df_msy$Region <- factor(df_msy$Region,levels=group_names[c(7,1:6)])
df_msy$Regulation <- factor(df_msy$Regulation,levels=colnames(escape_reg)[c(2,4,3,5,1)])

ran_df <- box_fn(mod="random",year_sim=sim_eval,ref_point_label="S(baseline)",ref_point=baseline_spawn)
ran_df$Model <- "Random walk"
ar_df <- box_fn(mod="ar1",year_sim=sim_eval,ref_point_label="S(baseline)",ref_point=baseline_spawn)
ar_df$Model <- "Mean reverting"
tr_df <- box_fn(mod="tr1",year_sim=sim_eval,ref_point_label="S(baseline)",ref_point=baseline_spawn)
tr_df$Model <- "Trending"

df_base <- rbind(ran_df,ar_df,tr_df)
df_base$Region <- factor(df_base$Region,levels=group_names[c(7,1:6)])
df_base$Regulation <- factor(df_base$Regulation,levels=colnames(escape_reg)[c(2,4,3,5,1)])

df_total <- rbind(df,df_msy,df_base)

aggregate(Risk*100~Reference.Point+Regulation+Model,df_total,mean)
aggregate(Risk*100~Reference.Point+Region+Model,df_total[df_total$Regulation=='10-year average',],mean)
aggregate(Risk*100~Reference.Point+Regulation+Model,df_total[df_total$Regulation%in%c('10-year average',"50% AK & BC reduction","No harvest"),],mean)
aggregate(Ratio~Reference.Point+Regulation+Model,df_total[df_total$Regulation%in%c('10-year average',"50% AK & BC reduction","No harvest"),],mean)
dodge <- position_dodge(width = 0.9)
margins <- c(0.25,0.25,0.25,0.5)
ggplot(df,aes(x=Regulation,y=Risk*100,fill=Regulation)) +
  geom_violin(position=dodge,scale="width") +
  geom_jitter(pch=21,color="white",alpha=0.5,fill="grey20",width=0.25) +
  geom_boxplot(position=dodge,width=0.1,color="white",alpha=0.5) +
  ylab(expression("% risk below limit reference point ( "*S["20 years forward"]<=S[GEN]*" )"))+
  facet_grid(rows=vars(Region),cols=vars(Model),labeller=label_wrap_gen(width=15,multi_line = TRUE)) +
  theme_minimal() +
  #coord_flip(clip = "off") +
  #guides(fill = guide_legend(override.aes = list(size=2)))+
  scale_fill_brewer(type="qual",palette=3,direction = -1) +
  theme(legend.position="top",strip.text.y = element_text(size=6.5),strip.text.x = element_text(size=8),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),legend.text=element_text(size=6),legend.title=element_text(size=7),axis.title=element_text(size=8),legend.key.size = unit(0.9, "line"),panel.spacing.y = unit(0.75, "lines"))
ggsave("Figures/scenario comparisons.jpeg",width=5.5,height=7,units="in")


ggplot(df_msy,aes(x=Regulation,y=Ratio,fill=Regulation)) +
  geom_violin(position=dodge,scale="width") +
  geom_jitter(pch=21,color="white",alpha=0.5,fill="grey20",width=0.25) +
  geom_boxplot(position=dodge,width=0.1,color="white",alpha=0.5) +
  geom_hline(yintercept=1,lty=2,colour="red")+
  ylab(expression("Relative population status ( "*S["20 years forward"]/0.8*S[MSY]*" )"))+
  facet_grid(rows=vars(Region),cols=vars(Model),labeller=label_wrap_gen(width=15,multi_line = TRUE)) +
  theme_minimal() +
  #coord_flip(clip = "off") +
  #guides(fill = guide_legend(override.aes = list(size=2)))+
  scale_fill_brewer(type="qual",palette=3,direction = -1) +
  theme(legend.position="top",strip.text.y = element_text(size=6.5),strip.text.x = element_text(size=8),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),legend.text=element_text(size=6),legend.title=element_text(size=7),axis.title=element_text(size=8),legend.key.size = unit(0.9, "line"),panel.spacing.y = unit(0.75, "lines"))
ggsave("Figures/scenario comparisons ratio.jpeg",width=5.5,height=7,units="in")


df_total$Reference.Point <- factor(df_total$Reference.Point,levels=c("S(gen)","S(MSY)","S(baseline)"))
ggplot(df_total[df_total$Region=="Coast-wide",],aes(x=Regulation,y=Ratio,fill=Regulation)) +
  geom_violin(position=dodge,scale="width") +
  geom_jitter(pch=21,color="white",alpha=0.5,fill="grey20",width=0.25) +
  geom_boxplot(position=dodge,width=0.1,color="white",alpha=0.5) +
  geom_hline(yintercept=1,lty=2,colour="red")+
  ylab(expression("Relative population status ( "*S["20 years forward"]/"Reference Point)"))+
  facet_grid(rows=vars(Reference.Point),cols=vars(Model),labeller=label_wrap_gen(width=15,multi_line = TRUE),scales="free") +
  theme_minimal() +
  #coord_flip(clip = "off") +
  #guides(fill = guide_legend(override.aes = list(size=2)))+
  scale_fill_brewer(type="qual",palette=3,direction = -1) +
  theme(legend.position="top",strip.text.y = element_text(size=6.5),strip.text.x = element_text(size=8),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),legend.text=element_text(size=6),legend.title=element_text(size=7),axis.title=element_text(size=8),legend.key.size = unit(0.9, "line"),panel.spacing.y = unit(0.75, "lines"))
ggsave("Figures/scenario comparisons ratios coastal.jpeg",width=5.5,height=7,units="in")

ggplot(df_total[df_total$Region=="Coast-wide",],aes(x=Regulation,y=100*Risk,fill=Regulation)) +
  geom_violin(position=dodge,scale="width") +
  geom_jitter(pch=21,color="white",alpha=0.5,fill="grey20",width=0.25) +
  geom_boxplot(position=dodge,width=0.1,color="white",alpha=0.5) +
  ylab(expression("% risk below reference point ( "*S["20 years forward"]<="Reference Point)"))+
  facet_grid(rows=vars(Reference.Point),cols=vars(Model),labeller=label_wrap_gen(width=15,multi_line = TRUE)) +
  theme_minimal() +
  #coord_flip(clip = "off") +
  #guides(fill = guide_legend(override.aes = list(size=2)))+
  scale_fill_brewer(type="qual",palette=3,direction = -1) +
  theme(legend.position="top",strip.text.y = element_text(size=6.5),strip.text.x = element_text(size=8),axis.text.x=element_text(size=7,angle=45,hjust=1),axis.text.y=element_text(size=7),legend.text=element_text(size=6),legend.title=element_text(size=7),axis.title=element_text(size=8),legend.key.size = unit(0.9, "line"),panel.spacing.y = unit(0.75, "lines"))
ggsave("Figures/scenario comparisons risk coastal.jpeg",width=5.5,height=7,units="in")

library(dplyr)
SR.dat_temp <- SR.dat_full[!is.na(SR.dat_full$escapement),] %>%
  group_by(pop_no) %>%
  slice(which.max(year))
SR.dat_temp$lrp <- SR.dat_temp$escapement/Sgen_priors$mean
SR.dat_temp$usr <- SR.dat_temp$escapement/(0.8*Smsy_priors$mean)
SR.dat_temp$bratio <- SR.dat_temp$escapement/baseline_spawn
SR.dat_temp$umsy <- ifelse(!is.na(SR.dat_temp$er_2),SR.dat_temp$er_2,SR.dat_temp$er_E)/Umsy_priors$mean
SR.dat_temp <- SR.dat_temp[SR.dat_temp$year>=2017,]
sum(SR.dat_temp$lrp<=1)/nrow(SR.dat_temp)
sum(SR.dat_temp$usr<=1 & SR.dat_temp$lrp >1)/nrow(SR.dat_temp)
sum(SR.dat_temp$bratio<=1)/nrow(SR.dat_temp)
sum(SR.dat_temp$umsy>=1)/nrow(SR.dat_temp)

s_over_msy <- sapply(1:n.pops,function(x){sum(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]<=(0.8*Smsy_priors$mean[x]),na.rm=TRUE)/sum(!is.na(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]))})
s_over_gen <- sapply(1:n.pops,function(x){sum(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]<=(Sgen_priors$mean[x]),na.rm=TRUE)/sum(!is.na(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]))})
s_over_base <- sapply(1:n.pops,function(x){sum(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]<=(baseline_spawn[x]),na.rm=TRUE)/sum(!is.na(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]))})
s_msy_ratio <- sapply(1:n.pops,function(x){mean(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]/(0.8*Smsy_priors$mean[x]),na.rm=TRUE)})
s_gen_ratio <- sapply(1:n.pops,function(x){mean(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]/(Sgen_priors$mean[x]),na.rm=TRUE)})
s_base_ratio <- sapply(1:n.pops,function(x){mean(SR.dat_predict[SR.dat_predict$pop_no==x,"escapement"]/(baseline_spawn[x]),na.rm=TRUE)})

pop_risks <- data.frame("pop"=pop_names,"region"=group_names[group],"s_over_gen"=s_over_gen,"s_over_msy"=s_over_msy,"s_over_base"=s_over_base,s_gen_ratio,s_msy_ratio,s_base_ratio)

aggregate(100*(1-s_base_ratio)~region,pop_risks,mean)
aggregate(100*(1-s_base_ratio)~region,pop_risks,range)

colSums(pop_risks[,3:5]>0,na.rm=TRUE)
sum((pop_risks[,4]*pop_risks[,5])>0,na.rm=TRUE)
sum((pop_risks[,4]>0 & pop_risks[,5]<0) | (pop_risks[,4]<0 & pop_risks[,5]>0),na.rm=TRUE)

colSums(pop_risks[,3:5]==1,na.rm=TRUE)
sum((pop_risks[,4]*pop_risks[,5])==1,na.rm=TRUE)
colSums(pop_risks[,3:5]==0,na.rm=TRUE)

SR.dat_predict[SR.dat_predict$population=="east_arm",]
SR.dat_predict[SR.dat_predict$population=="nias",]
pop_risks[pop_risks$s_msy_ratio<1,]
