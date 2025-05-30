Smsy_priors <- readRDS("Results/Smsy.rds")
Sgen_priors <- readRDS("Results/Sgen.rds")
Umsy_priors <- readRDS("Results/Umsy.rds")
SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)

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
table(SR.dat[sapply(1:max(SR.dat$pop_no),function(x){which(SR.dat$pop_no==x)[1]}),c("CU","region")])

baselines <- sapply(1:n.pops,function(x){mean(SR.dat$total_runE[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines_2 <- sapply(1:n.pops,function(x){mean(SR.dat$total_run2[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})
baselines<-ifelse(!is.na(baselines_2),baselines_2,baselines)

baseline_spawn <- sapply(1:n.pops,function(x){mean(SR.dat$escapement[SR.dat$pop_no==x & (SR.dat$year>=(2000)) & (SR.dat$year<=2015)],na.rm=TRUE)})


SR.dat$baseline_run <- (baselines)[SR.dat$pop_no]
SR.dat$baseline_escapement <- (baseline_spawn)[SR.dat$pop_no]
SR.dat$lrp <- (Sgen_priors$mean)[SR.dat$pop_no]
SR.dat$usr <- 0.8*(Smsy_priors$mean)[SR.dat$pop_no]
SR.dat$ucur <- ifelse(!is.na(SR.dat$er_2),SR.dat$er_2,SR.dat$er_E)
SR.dat$umsy <- (Umsy_priors$mean)[SR.dat$pop_no]
SR.dat$regime <- ifelse(SR.dat$year>=1990 & SR.dat$year<=2001,"early",ifelse(SR.dat$year<=2011,"recovering",ifelse(SR.dat$year>=2017,"current","recent")))

df <- SR.dat[SR.dat$population=="martin",]

df_long <- tidyr::pivot_longer(df,cols=c(baseline_escapement,lrp,usr),values_to = "value",names_to = "metric")
sgen_frac <- unique(df$lrp/max(df$escapement,na.rm = TRUE))
smsy_frac <- unique(df$usr/max(df$escapement,na.rm = TRUE))- sgen_frac
rem_frac <- 1-(sgen_frac+smsy_frac)
sbase_frac <- unique(df$baseline_escapement/max(df$escapement,na.rm = TRUE)) - sgen_frac
rem_frac2 <- 1-(sgen_frac+sbase_frac)

sgen <- unique(df$lrp)
smsy <- unique(df$usr)
sbase <- unique(df$baseline_escapement)


margins <- c(2,0.1,0.1,0.1)
ylabel <- "Relative value"
legend_pos<-c(0.9,0.9)
label_text <- c(as.character(expression(italic(0.8*S[MSY]))),
                as.character(expression(italic(S[GEN]))),
                as.character(expression(italic(S[Baseline]))))
label_df <- data.frame("Year"=1975,"scaled_metric"=200+c(smsy,sgen,sbase),"label"=label_text)
stock_status_df <- data.frame("Reference_point"=c("WSP","WSP","WSP",
                                                  "CCFN","CCFN","CCFN"),
                              "Metric"=c(sgen_frac,0.105,1-(0.105+sgen_frac),sgen_frac,0.39,1-(0.39+sgen_frac)),
                              "Status"=c("Critical","Cautious","Healthy",
                                         "Critical","Cautious","Healthy"))
stock_status_df$Reference_point <- relevel(as.factor(stock_status_df$Reference_point),"WSP")
stock_status_df$Status <- factor(stock_status_df$Status,levels=c("Healthy","Cautious","Critical"))
margin_pBar <- c(0.1,2,0.75,0.1)
pBar <- ggplot()+
  geom_bar(data=stock_status_df,aes(fill=Status, y=Metric, x=Reference_point),position="stack", stat="identity")+
  theme_void()+
  ylim(0,1) +
  scale_fill_manual("Stock status",values=c("darkgreen","goldenrod","red4"),labels=c("Healthy","Cautious","Critical")) +
  ggtitle("Stock status")+
  theme(legend.position="none",legend.key = element_rect(fill = "white",color=NA),legend.text=element_text(size=8),plot.margin=unit(margin_pBar,"line"),text=element_text(size=10),axis.text.x = element_text(size=8),plot.title = element_text(hjust = 0.5))
pBar
p1 <- ggplot(data=df,aes(x=year,y=escapement)) +
  geom_hline(yintercept=sgen,lwd=0.5,lty=3,colour="grey20")+
  geom_hline(yintercept=smsy,lwd=0.5,lty=3,colour="grey20")+
  geom_hline(yintercept=sbase,lwd=0.5,lty=5,colour="darkblue")+
  ylab("Spawner abundance") + xlab("Year") +
  xlim(1975,2020) +
  annotate("text",parse=T,label=expression(italic(0.8*S[MSY])*'; Wild Salmon Policy USR'),x=label_df[1,"Year"],y=label_df[1,"scaled_metric"],hjust=0,size=2.5)+
  annotate("text",parse=T,label=expression(italic(S[GEN])*'; Wild Salmon Policy LRP'),x=label_df[2,"Year"],y=label_df[2,"scaled_metric"],hjust=0,size=2.5)+
  annotate("text",parse=T,label=expression(bar(italic(S[2000-15]))*'; CCFN-alternative USR'),x=label_df[3,"Year"],y=label_df[3,"scaled_metric"],hjust=0,size=2.5,colour="darkblue")+
  geom_line(data=df[!is.na(df$escapement),],aes(x=year,y=escapement),colour="grey30",alpha=0.5) +
  theme_minimal() +
  #theme(legend.position="top",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1),plot.margin=unit(margins,"line"),legend.key = element_rect(fill = "white",color=NA),legend.text=element_text(size=7))
  theme(legend.position="none",axis.line = element_line(colour = "black"),axis.ticks=element_line(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.key = element_rect(fill = "white",color=NA),legend.text=element_text(size=10),plot.margin=unit(margins,"line"),text = element_text(size=10))

megaP <- egg::ggarrange(p1,pBar,ncol=2,nrow=1,widths=c(8,1.5),heights=c(1))
megaP

megaP <- ggpubr::ggarrange(p1,pBar,ncol=2,nrow=1,widths=c(10,3),heights=c(1))

ggsave(plot=megaP,filename="Figures/coho wild salmon policy.jpeg",units="in",height=4,width=5,dpi=800)
knitr::plot_crop("Figures/coho wild salmon policy.jpeg")
