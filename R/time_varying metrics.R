library(tidyr)
library(dplyr)
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat_full <- SR.dat
SR.dat_predict<-subset(SR.dat,SR.dat$year>=2017)
SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)
time_varying <- readRDS("Results/time_varying_metrics.rds")
time_varying <- tidyr::pivot_wider(time_varying,names_from="Metric",values_from="Value")

Smsy_priors <- readRDS("Results/Smsy.rds")
Sgen_priors <- readRDS("Results/Sgen.rds")
Umsy_priors <- readRDS("Results/Umsy.rds")

df_full <-merge(SR.dat,time_varying,by.x=c("year","population"),by.y=c("Year","Population"))
df_full$Ut <- ifelse(!is.na(df_full$er_2),df_full$er_2,df_full$er_E)

df_full <- tidyr::pivot_longer(df_full,cols=c("Smsy","Sgen","Umsy"),names_to="Metric",values_to = "Time-varying")
df_full$Smsy <- Smsy_priors$mean[match(df_full$pop_no,Smsy_priors$pop)]
df_full$Sgen <- Sgen_priors$mean[match(df_full$pop_no,Sgen_priors$pop)]
df_full$Umsy <- Umsy_priors$mean[match(df_full$pop_no,Umsy_priors$pop)]

df_full <- tidyr::pivot_longer(df_full,cols=c("Smsy","Sgen","Umsy"),names_to="Equilibrium",values_to = "Mean")
df_full <- tidyr::pivot_longer(df_full,cols=c("Equilibrium","Metric"),names_to="Value",values_to = "Metric")
co_pops<-read.table("Data/coho_groups.txt",header=TRUE)
co_pops$mean_total<-NA

for(i in co_pops$pop_no){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6
group_names <- c("Central Coast (South)","Hecate Lowlands","Inner Waters","Haida Gwaii","Skeena","Nass")
group_names <- group_names[c(4,6,5,2,3,1)]
df_full$group <- co_pops$group[match(df_full$population,co_pops$population)]
df_full$region <- group_names[df_full$group]

df_summ <- df_full %>%
  group_by(year,region,Metric) %>%
  summarise(value=mean(Value),mean=mean(Mean))

ggplot(df_summ[df_summ$Metric%in%c("Sgen","Smsy"),],aes(x=year,y=value,colour=Metric,fill=Metric))+
  geom_line()+
  geom_point()+
  geom_line(df_summ[df_summ$Metric%in%c("Sgen","Smsy"),],aes(x=year,y=mean))+
  facet_grid(rows=vars(region),cols=vars(Metric),labeller=label_wrap_gen(width=15,multi_line = TRUE),scales='free') +
  theme_minimal() +
  #coord_flip(clip = "off") +
  #guides(fill = guide_legend(override.aes = list(size=2)))+
  scale_fill_brewer(type="qual",palette=3,direction = -1) +
  theme(legend.position="top")
