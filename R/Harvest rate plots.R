library(tidyverse)
ER_values <- read.csv("Data/ER_values.csv", header=T) # average harvest by sector and region
harvest_rates <- openxlsx::read.xlsx("Data/CohoCU_harvest_by country_NCC_1980-2020.xlsx") # time-series modified from English et al. 2017
harvest_rates <- harvest_rates[complete.cases(harvest_rates),]
harvest_rates <- pivot_longer(harvest_rates,cols=-1,names_to=c("CU_region","sector"),names_pattern = "(.*)_([^_]+)$",values_to="harvest_rate")

harvest_rates$region <- factor(harvest_rates$CU_region,levels=unique(harvest_rates$CU_region),labels=c("Haida Gwaii","Nass","Skeena","Coast-wide","Hecate Lowlands","Inner Waters (Douglas/Kitimat CU)","Inner Waters","Central Coast (South)","Central Coast (Rivers/Smith CUs)"))
harvest_rates$region <- relevel(harvest_rates$region,ref="Coast-wide")
harvest_rates$sector <- factor(harvest_rates$sector,levels=unique(harvest_rates$sector),labels=c("Alaska","Canadian","Total","Canadian"))

harvest_rates$harvest_rate[harvest_rates$region=="Hecate Lowlands"&harvest_rates$sector=="Canadian"] <- harvest_rates$harvest_rate[harvest_rates$region=="Coast-wide"&harvest_rates$sector=="Canadian"]

harvest_rates$harvest_rate[harvest_rates$region=="Inner Waters (Douglas/Kitimat CU)"&harvest_rates$sector=="Canadian"] <- harvest_rates$harvest_rate[harvest_rates$region=="Coast-wide"&harvest_rates$sector=="Canadian"]

harvest_rates$harvest_rate[harvest_rates$region=="Inner Waters"&harvest_rates$sector=="Canadian"] <- harvest_rates$harvest_rate[harvest_rates$region=="Coast-wide"&harvest_rates$sector=="Canadian"]

harvest_rates$harvest_rate[harvest_rates$region=="Central Coast (South)"&harvest_rates$sector=="Canadian"] <- harvest_rates$harvest_rate[harvest_rates$region=="Coast-wide"&harvest_rates$sector=="Canadian"]

harvest_rates$harvest_rate[harvest_rates$region=="Central Coast (Rivers/Smith CUs)"&harvest_rates$sector=="Canadian"] <- harvest_rates$harvest_rate[harvest_rates$region=="Coast-wide"&harvest_rates$sector=="Canadian"]

table(harvest_rates$sector,harvest_rates$region)
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
# SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)

SR.dat$er_est <- SR.dat$er_2
SR.dat$er_est[is.na(SR.dat$er_est)] <- SR.dat$er_E[is.na(SR.dat$er_est)]
SR.dat$Ut <- ifelse(!is.na(SR.dat$er_2),SR.dat$er_2,SR.dat$er_E)
SR.dat$total_run_final <- ifelse(!is.na(SR.dat$total_run2),SR.dat$total_run2,SR.dat$total_runE)

## I grouped rivers and smiths inlet with Area 7 & 8
co_pops<-read.table("Data/coho_groups.txt",header=TRUE)

co_pops[which(co_pops$group==7),4]<-6
group_names <- c("Central Coast (South)","Hecate Lowlands","Inner Waters","Haida Gwaii","Skeena","Nass")
group_names <- group_names[c(4,6,5,2,3,1)]
SR.dat$group <- co_pops$group[match(SR.dat$population,co_pops$population)]
SR.dat$region <- group_names[SR.dat$group]
SR.dat$region <- factor(SR.dat$region,levels=c("Haida Gwaii","Nass","Skeena","Hecate Lowlands","Inner Waters","Central Coast (South)"))
table(SR.dat[SR.dat$CU=="northern_coastal",c("stat_area","region")])
table(SR.dat[SR.dat$CU=="douglas_channel_kitimat",c("stat_area","region")])
table(SR.dat[,c("stat_area","region","CU")])

regional_average <- SR.dat %>%
  group_by(region,year) %>%
  summarise("Ut"=mean(Ut))
SR.dat[SR.dat$year==2004,]
ggplot(data=regional_average,aes(x=year,y=Ut,colour=region))+
  geom_line(data=regional_average,aes(x=year,y=Ut,colour=region))+
  geom_point(data=regional_average,aes(x=year,y=Ut,colour=region))+
  facet_wrap(~region,labeller=label_wrap_gen(width=15,multi_line = TRUE),ncol=3) +
  scale_colour_brewer(name="Region",palette = "Paired",direction=1) +
  theme_minimal() +
  theme(legend.position="top") +
  theme(text = element_text(size=10),axis.text = element_text(size=10),legend.title = element_text(size = 9),legend.text = element_text(size=7)) +
  theme(legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(0.5, "lines"),legend.margin = margin(c(0,1,0,0),unit="lines"))

regional_average_harvest <- harvest_rates %>%
  filter(region!="Coast-wide",sector!="Total") %>%
  group_by(year,region,sector) %>%
  summarise("Ut"=mean(harvest_rate))
regional_average_harvest$region <- factor(regional_average_harvest$region,levels=c("Haida Gwaii","Nass","Skeena","Hecate Lowlands","Inner Waters (Douglas/Kitimat CU)","Inner Waters","Central Coast (South)","Central Coast (Rivers/Smith CUs)"))
regional_average_harvest$year <- as.numeric(regional_average_harvest$year)
ggplot(data=regional_average_harvest,aes(x=year,y=Ut,colour=region,group=interaction(region, sector),shape=sector))+
  geom_line(data=regional_average_harvest,aes(x=year,y=Ut,colour=region,group=interaction(region, sector),linetype=sector))+
  geom_point(data=regional_average_harvest,aes(x=year,y=Ut,colour=region,fill=region,group=interaction(region, sector),shape=sector))+
  facet_wrap(~region,labeller=label_wrap_gen(width=40,multi_line = TRUE),ncol=3) +
  theme_minimal() +
  ylab("Exploitation rate")+xlab("Year")+
  scale_shape_manual(name="Sector",values=c(16,2)) +
  scale_linetype_manual(name="Sector",values=c(1,2)) +
  scale_colour_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#fb9a99","#e31a1c","#e31a1c")) +
  scale_fill_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#fb9a99","#e31a1c","#e31a1c")) +
  theme(legend.position="top") +
  theme(text = element_text(size=10),axis.text = element_text(size=10),legend.title = element_text(size = 8),legend.text = element_text(size=6)) +
  theme(legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(0.5, "lines"),legend.margin = margin(c(0,1,0,-1),unit="lines"))
ggsave("Figures/harvest rates over time.jpeg", width = 7, height=6,units="in", dpi=600)


regional_average_run <- SR.dat %>%
  group_by(region,year) %>%
  summarise("total_run"=sum(total_run_final,na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(region) %>%
  mutate("total_run_scale"=as.numeric(scale(total_run))) %>%
  ungroup()

ggplot(data=regional_average_run,aes(x=year,y=total_run_scale,colour=region,group=region))+
  geom_line(data=regional_average_run,aes(x=year,y=total_run_scale,colour=region,group=region))+
  geom_point(data=regional_average_run,aes(x=year,y=total_run_scale,colour=region,fill=region,group=region))+
  facet_wrap(~region,labeller=label_wrap_gen(width=40,multi_line = TRUE),ncol=3) +
  theme_minimal() +
  ylab("Scaled total abundance (escapement + harvest)")+xlab("Year")+
  scale_colour_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  scale_fill_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  theme(legend.position="top") +
  theme(text = element_text(size=10),axis.text = element_text(size=10),legend.title = element_text(size = 8),legend.text = element_text(size=6)) +
  theme(legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(0.5, "lines"),legend.margin = margin(c(0,1,0,-1),unit="lines"))
ggsave("Figures/regional run sizes over time.jpeg", width = 7, height=6,units="in", dpi=600)

SR.dat <- SR.dat %>%
  group_by(population) %>%
  mutate("total_run_scale"=as.numeric(scale(total_run_final))) %>%
  ungroup()

ggplot(data=SR.dat,aes(x=year,y=total_run_scale,colour=region,group=interaction(region, population)))+
  geom_line(data=SR.dat,aes(x=year,y=total_run_scale,colour=region,group=interaction(region, population)),alpha=0.25)+
  geom_point(data=SR.dat,aes(x=year,y=total_run_scale,colour=region,fill=region,group=interaction(region, population)),alpha=0.25)+
  # geom_smooth(data=SR.dat,aes(x=year,y=total_run_scale,colour=region,fill=region,group=region),method = "gam", formula = y ~ s(x, bs = "cs", fx = FALSE, k = length(unique(SR.dat$year))))+
  geom_line(data=regional_average_run,aes(x=year,y=total_run_scale,colour=region,group=region),lwd=0.9)+
  geom_point(data=regional_average_run,aes(x=year,y=total_run_scale,colour=region,group=region),alpha=0.5)+
  facet_wrap(~region,labeller=label_wrap_gen(width=40,multi_line = TRUE),ncol=3) +
  theme_minimal() +
  ylab("Scaled total abundance (escapement + harvest)")+xlab("Year")+
  scale_colour_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  scale_fill_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  theme(legend.position="top",strip.text.x = element_text(hjust = -0)) #+
  # theme(legend.box.just="center",legend.box="horizontal",legend.justification = "center",legend.key.size=unit(0.5, "lines"),legend.margin = margin(c(0,1,0,-1),unit="lines"))
ggsave("Figures/regional and local run sizes over time.jpeg", width = 7, height=5,units="in", dpi=600)

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s, 1, 1)),
                           {s <- substring(s, 2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

SR_dotchart <- SR.dat
SR_dotchart$population <- gsub("kspx","Kispiox",SR_dotchart$population)
SR_dotchart$population <- gsub("zolzap","Ksi Ts'oohl Ts'ap",SR_dotchart$population)
SR_dotchart$population <- gsub("_"," ",SR_dotchart$population)
SR_dotchart$population <- gsub("ck","",SR_dotchart$population)
SR_dotchart$population <- stringr::str_trim(SR_dotchart$population)
SR_dotchart$population <- capwords(SR_dotchart$population)
SR_dotchart$population <- factor(SR_dotchart$population,levels=rev(sort(unique(SR_dotchart$population))))
SR_dotchart$escapement_logic <- factor(ifelse(!is.na(SR_dotchart$escapement),"Observed","Missing"),levels=c("Observed","Missing"))

ggplot(data=SR_dotchart,aes(x=year,y=population,colour=region,shape=escapement_logic))+
  geom_point(data=SR_dotchart,aes(x=year,y=population,colour=region,shape=escapement_logic),alpha=0.9)+
  # facet_wrap(~region,labeller=label_wrap_gen(width=40,multi_line = TRUE),ncol=1) +
  theme_classic() +
  ylab("Population")+xlab("Year")+
  scale_colour_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  scale_fill_manual(name="Region",values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c")) +
  scale_shape_manual(name="Survey",values=c(16,4))+
  guides(colour=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  theme(legend.position="top",strip.text.x = element_text(hjust = -0),legend.justification = "left",legend.box.margin = margin(c(0,0,0,-4),unit="lines"),legend.box = "horizontal",legend.text = element_text(size=10),legend.title = element_text(size=10))
ggsave("Figures/escapement dotchart.jpeg", width = 6, height=8,units="in", dpi=600)
