wrapper <- function(x, ...) paste(strwrap(x, ...), collapse = "\n")
color_more = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
new_pall <- function(...){palette(color_more,...)}



library(ggplot2)
library(ggridges)
library(merTools)
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(grid)
library(tidyverse)
library(gganimate)
library(transformr)
library(zoo)
library(randomcoloR)
MPA <- read.csv("Data/MPA NCC trends.csv",stringsAsFactors = TRUE)
MPA <- subset(MPA,MPA$Year>=1900)
MPA <- MPA[complete.cases(MPA),]
MPA$Numbers <- pmax(0,MPA$Numbers)
MPA$Functional_Group <- (ifelse(MPA$Group=="Salmon",MPA$Group,
                               ifelse(MPA$Metric=='Body Size',
                                      ifelse(MPA$Source=="McGreer & Frid","Rockfish","Rockfish"),as.character(MPA$Group))))
MPA$Numbers[MPA$Group=='Salmon'] <- exp(MPA$Numbers[MPA$Group=='Salmon'])
MPA$Taxa_Group <- (ifelse(MPA$Group=="Salmon",as.character(MPA$Group),
                                ifelse(MPA$Metric=='Body Size',
                                       ifelse(MPA$Source=="McGreer & Frid","Rockfish (FI)","Rockfish"),as.character(MPA$Group))))
MPA <- subset(MPA,MPA$Source!='McGreer & Frid')
MPA <- subset(MPA,MPA$Species!="Chinook")
MPA_mean <-  MPA %>%
  group_by(Species,Metric,Group,Population,Functional_Group,Taxa_Group,Year) %>%
  summarise("summed_metric"=mean(Numbers,na.rm=TRUE))
MPA_species <- MPA_mean %>%
  group_by(Species,Metric,Group,Functional_Group,Taxa_Group,Year) %>%
  summarise("summed_metric"=sum(summed_metric,na.rm=TRUE)) %>%
  group_by(Group,Functional_Group,Metric,Taxa_Group,Year) %>%
  summarise("mean_metric"=mean(summed_metric,na.rm=TRUE)) %>%
  group_by(Metric,Functional_Group) %>%
  summarise("Year"=Year,"Group"=Group,"Taxa_Group"=Taxa_Group,"scaled_metric"=if_else(Functional_Group=='Kelp',mean_metric,mean_metric/max(mean_metric,na.rm=TRUE)))

MPA_scaled <- MPA_mean %>%
  group_by(Species,Metric,Group,Year) %>%
  summarise("summed_metric"=sum(summed_metric,na.rm=TRUE)) %>%
  group_by(Species,Metric,Group,Year) %>%
  summarise("mean_metric"=mean(summed_metric,na.rm=TRUE)) %>%
  group_by(Metric,Group) %>%
  summarise("Year"=Year,"scaled_metric"=mean_metric/max(mean_metric,na.rm=TRUE))

MPA_scaled$Taxa <- interaction(MPA_scaled$Group,MPA_scaled$Metric,sep=" ")
MPA_plot <- ggplot(MPA_scaled,aes(x=Year,y=scaled_metric,colour=Taxa)) +
  geom_smooth(data=filter(MPA_scaled,Group%in%c("Salmon")),aes(x=Year,y=scaled_metric,colour=Taxa),lwd=0.5,se=FALSE) +
  geom_smooth(data=filter(MPA_scaled,Metric%in%c("Body Size")),aes(x=Year,y=scaled_metric,colour=Taxa),lwd=0.5,se=FALSE) +
  geom_line(data=filter(MPA_scaled,!Group%in%c("Salmon"),!Metric%in%c("Body Size")),aes(x=Year,y=scaled_metric,colour=Taxa))+
  ylab("Scaled metric") +
  ylim(c(0,1)) + xlim(c(1925,2021)) +
  scale_colour_brewer("Ecosystem metrics",type="qual",palette = "Set2") +
  scale_fill_brewer("Ecosystem metrics",type="qual",palette = "Set2") +
  theme_minimal() +
  theme(legend.position="right",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1)) +
  labs(title="Rationale for Marine Protected Areas",
       subtitle="Ecosystem trends along the Central Coast",
       caption="Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.")
ggsave(filename="Figures/CC MPA trends.jpeg",units="in",height=4,width=8,dpi=800)

MPA_plot <- ggplot(MPA_scaled,aes(x=Year,y=scaled_metric,colour=Taxa)) +
  geom_line(aes(x=Year,y=scaled_metric,colour=Taxa))+
  ylab("Scaled metric") +
  ylim(c(0,1)) + xlim(c(1925,2021)) +
  scale_colour_brewer("Ecosystem metrics",type="qual",palette = "Set2") +
  scale_fill_brewer("Ecosystem metrics",type="qual",palette = "Set2") +
  theme_minimal() +
  theme(legend.position="right",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1)) +
  labs(title="Rationale for Marine Protected Areas",
       subtitle="Ecosystem trends along the Central Coast",
       caption="Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.")

for(i in 1:length(unique(MPA_scaled$Taxa)))
{
  MPA_new <- MPA_scaled
  MPA_new$transparency <- 0.5
  MPA_new$transparency[MPA_new$Taxa==unique(MPA_scaled$Taxa)[i]] <- 1
  MPA_new$line_width <- 0.33
  MPA_new$line_width[MPA_new$Taxa==unique(MPA_scaled$Taxa)[i]] <- 1
  MPA_plot <- ggplot(MPA_new,aes(x=Year,y=scaled_metric,colour=Taxa)) +
    geom_line(stat="smooth",data=filter(MPA_new,Group%in%c("Salmon")),aes(x=Year,y=scaled_metric,colour=Taxa,alpha=transparency,size=line_width),se=FALSE) +
    geom_line(stat="smooth",data=filter(MPA_new,Metric%in%c("Body Size")),aes(x=Year,y=scaled_metric,colour=Taxa,alpha=transparency,size=line_width),se=FALSE) +
    geom_line(data=filter(MPA_new,!Group%in%c("Salmon"),!Metric%in%c("Body Size")),aes(x=Year,y=scaled_metric,colour=Taxa,alpha=transparency,size=line_width))+
    ylab("Scaled metric") +
    ylim(c(0,1)) + xlim(c(1925,2021)) +
    scale_colour_brewer("Ecosystem metrics",type="qual",palette = "Set2") +
    scale_alpha_continuous(guide=FALSE,range=c(0.33,1)) +
    scale_size_continuous(guide=FALSE,range=c(0.5,1.25)) +
    theme_minimal() +
    theme(legend.position="right",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1)) +
    labs(title="Rationale for Marine Protected Areas",
         subtitle="Ecosystem trends along the Central Coast",
         caption="Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.")
  MPA_plot
  fig_name <- paste("Figures/CC MPA trends ",unique(MPA_scaled$Taxa)[i],".jpeg",sep="")
  ggsave(filename=fig_name,units="in",height=4,width=8,dpi=800)
}

set.seed(1)
color_pal <- distinctColorPalette(12)
MPA_species$Taxa <- interaction(MPA_species$Taxa_Group,MPA_species$Metric,sep=" ")
MPA_species$Group_new <- ifelse(MPA_species$Group=="Salmon","Salmon (except Chinook)",
                                ifelse(MPA_species$Group=="Rockfish","Yelloweye rockfish",
                                       ifelse(MPA_species$Group=="Herring","Herring",
                                              ifelse(MPA_species$Group=="Kelp","Sea surface temperatures in April","Crab & eulachon"))))

MPA_species$Group_new <- factor(MPA_species$Group_new,levels=c("Crab & eulachon","Yelloweye rockfish","Salmon (except Chinook)","Herring","Sea surface temperatures in April"))

MPA_species$Taxa <- factor(MPA_species$Taxa,levels=c("Crab Catch per trip","Eulachon Catch","Rockfish Abundance","Rockfish Landings","Rockfish Longline harvest rate","Rockfish Body Size","Salmon Abundance","Herring Abundance","Herring Harvest rate","Herring Natural mortality","Kelp Temperature"),labels=c("Crab catch per trap","Eulachon catch","Biomass","Landings","Harvest rate - longline","Body size","Salmon Abundance","Abundance","Harvest rate","Natural mortality","Kelp Temperature"))

colors <- c("Crab catch per trap"="#a6cee3",
            "Eulachon catch"="#1f78b4",
            "Biomass"="black",
            "Landings"="#33a02c",
            "Longline HR"="#fb9a99",
            "Body size"="#e31a1c",
            "Salmon Abundance"="#fdbf6f",
            "Abundance"="#ff7f00",
            "Harvest rate"="#cab2d6",
            "Natural mortality"="#b2df8a",
            "Kelp Temperature"="#6a3d9a")

axis_names <- c(rep("Relative change",4),"Temperature")
margins <- c(0.1,0.1,0.1,0.1)

MPA_plot <- ggplot(MPA_species[MPA_species$Group!='Kelp',],aes(x=Year,y=scaled_metric,colour=Taxa)) +
  geom_line(aes(x=Year,y=scaled_metric,colour=Taxa),lwd=1) +
  geom_hline(yintercept=c(0.5),lwd=0.5,lty=3,colour="black") +
  facet_wrap(~Group_new,scales='fixed',nrow=5) +
  ylab("Relative change") + xlab("") +
  xlim(1950,2020) +
  scale_colour_manual("",values=colors,labels=names(colors)) +
  scale_fill_manual("",values=colors,labels=names(colors)) +
  theme_minimal() +
  theme(legend.position="top",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1),plot.margin=unit(margins,"line"))
  #labs(title="Rationale for Marine Protected Areas",subtitle="Ecosystem trends along the Central Coast",caption="Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.")
ggsave(filename="Figures/CC MPA trends - multipanel.jpeg",units="in",height=10,width=8,dpi=800)
MPA_plot2 <- ggplot(MPA_species[MPA_species$Group=='Kelp',],aes(x=Year,y=scaled_metric,colour=Taxa)) +
  geom_line(data=MPA_species[MPA_species$Group=='Kelp',],aes(x=Year,y=scaled_metric,colour=Taxa),lwd=1) +
  facet_wrap(~Group_new,scales='free_y',nrow=5) +
  ylab("Temperature") +
  xlim(1950,2020) +
  scale_colour_manual("",values=colors) +
  scale_fill_manual("",values=colors) +
  theme_minimal() +
  theme(legend.position="top",strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1),plot.margin=unit(margins,"line"))
  #labs(title="Rationale for Marine Protected Areas",subtitle="Ecosystem trends along the Central Coast",caption="Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.")

MPA_plot2
megaP <- ggarrange(MPA_plot+ rremove("xlab"),MPA_plot2,labels=NULL,ncol=1,nrow=2,heights=c(4,1),legend="top",common.legend=TRUE)
pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.",width=125),color="black",hjust=0,x=0.01,face="italic",size=10),top=text_grob(wrapper("Rationale for Marine Protected Areas",125),color="black",hjust=0,x=0.01,face="italic",size=10))
pAnnotated
ggsave(pAnnotated,filename="Figures/CC MPA trends - annotated d2.jpeg",units="in",height=10,width=6,dpi=800)

#MPA_species$Taxa <- factor(MPA_species$Taxa,levels=c("Crab Catch per trip","Eulachon Catch","Rockfish Abundance","Rockfish Landings","Rockfish Body Size","Salmon Abundance","Herring Abundance","Kelp Temperature"),labels=c("Crab catch per trap","Eulachon Catch","Abundance","Landings","Body Size","Salmon Abundance","Herring Abundance","Kelp Temperature"))
colors <- c("Crab catch per trap"="#a6cee3",
            "Eulachon catch"="#1f78b4",
            "Biomass"="black",
            "Landings"="#33a02c",
            "Harvest rate - longline"="#fb9a99",
            "Body size"="#e31a1c",
            "Salmon Abundance"="#fdbf6f",
            "Abundance"="#ff7f00",
            "Harvest rate"="#cab2d6",
            "Natural mortality"="#b2df8a",
            "Kelp Temperature"="#6a3d9a")
col_names <- c("Crab CPUE",
               "Eulachon Catch",
               "Abundance",
               "Landings",
               "Longline HR",
               "Body size",
               "Abundance",
               "Biomass",
               "Harvest",
               "SST")
make_plot <- function(data,ylabel,hline=TRUE,legend_logic=TRUE)
{
  sub_col <- colors[match(sort(unique(data$Taxa)),names(colors))]
  ifelse(legend_logic,legend_pos<-c(0.9,0.9),legend_pos<-"none")
  ggplot(data,aes(x=Year,y=scaled_metric,colour=Taxa)) +
    geom_line(aes(x=Year,y=scaled_metric,colour=Taxa),lwd=1) +
    {if(hline)geom_hline(yintercept=c(0.5),lwd=0.5,lty=3,colour="black")}+
    facet_wrap(~Group_new)+
    ylab(ylabel) + xlab("Year") +
    xlim(1950,2020) +
    scale_colour_manual("",values=sub_col,labels=names(sub_col)) +
    scale_fill_manual("",values=sub_col,labels=names(sub_col)) +
    theme_minimal() +
    theme(legend.position=legend_pos,strip.text.x=element_text(hjust=0),axis.text.x=element_text(hjust=1),plot.margin=unit(margins,"line"),legend.key = element_rect(fill = "white",color=NA),legend.text=element_text(size=7))
}
margins <- c(0.1,1.1,0.1,1)

p1 <- make_plot(data=MPA_species[MPA_species$Group=='Crab' | MPA_species$Group=='Eulachon',],ylabel="")
p2 <- make_plot(data=MPA_species[MPA_species$Group=='Rockfish' & MPA_species$Metric!='Landings',],ylabel="")
p3 <- make_plot(data=MPA_species[MPA_species$Group=='Salmon',],ylabel="Relative change",legend_logic=FALSE)
p4 <- make_plot(data=MPA_species[MPA_species$Group=='Herring'& MPA_species$Metric!='Natural mortality',],ylabel="")
p5 <- make_plot(data=MPA_species[MPA_species$Group=='Kelp',],ylabel="Temperature",hline=FALSE,legend_logic=FALSE)

megaP <- ggarrange(p1+ rremove("xlab"),
                   p2+ rremove("xlab"),
                   p3+ rremove("xlab"),
                   p4+ rremove("xlab"),
                   p5,labels=NULL,ncol=1,nrow=5,common.legend=FALSE)

pAnnotated <- annotate_figure(megaP,bottom=text_grob(wrapper("Data sources: DFO; CCIRA; Wild Salmon Center; Ban et al.; Moody et al.",width=125),color="black",face="italic",size=9,hjust=0,x=0.01),top=text_grob("Ecosystem changes along the Central Coast",color="black",hjust=0,x=0.01,face="bold",size=14))
pAnnotated
ggsave(pAnnotated,filename="Figures/CC MPA trends - annotated.jpeg",units="in",height=10,width=6,dpi=800)

p1 <- make_plot(MPA_species[MPA_species$Group=='Crab' | MPA_species$Group=='Eulachon',],ylabel="Relative change")
ggsave(p1,filename="Figures/CC MPA trends - Crab & Eulachon.jpeg",units="in",height=4,width=6,dpi=800)
p2 <- make_plot(MPA_species[MPA_species$Group=='Rockfish'& MPA_species$Metric!='Landings',],ylabel="Relative change")
ggsave(p2,filename="Figures/CC MPA trends - Rockfish.jpeg",units="in",height=4,width=6,dpi=800)
p3 <- make_plot(MPA_species[MPA_species$Group=='Salmon',],ylabel="Relative change",legend_logic=FALSE)
ggsave(p3,filename="Figures/CC MPA trends - Salmon.jpeg",units="in",height=4,width=6,dpi=800)
p4 <- make_plot(MPA_species[MPA_species$Group=='Herring' & MPA_species$Metric!='Natural mortality',],ylabel="Relative change")
ggsave(p4,filename="Figures/CC MPA trends - Herring.jpeg",units="in",height=4,width=6,dpi=800)
p5 <- make_plot(MPA_species[MPA_species$Group=='Kelp',],ylabel="Temperature",hline=FALSE,legend_logic=FALSE)
ggsave(p5,filename="Figures/CC MPA trends - SST.jpeg",units="in",height=4,width=6,dpi=800)
