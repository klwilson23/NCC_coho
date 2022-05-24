SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat<-subset(SR.dat,SR.dat$year<2017)
group_names <- c("Area 7-10","Hecate Lowlands - Area 5/6","Area 6 - Inner Waters","Haida Gwaii - Area 2E","Skeena - Area 4","Nass - Area 3")

resultSR_B3<-readRDS("Results/COSR3B.tr1_lalpha_MCMC.rds")
mcmc_names <- colnames(resultSR_B3[[1]])
resultSR_B3 <- as.matrix(resultSR_B3)
group_names <- c("Central Coast (South)","Hecate Lowlands","Inner Waters","Haida Gwaii","Skeena","Nass")
group_names <- group_names[c(4,6,5,2,3,1)]
region_names <- c("Area 7-10","Area 5-6","Area 6 - Inner","Area 2E","Area 4","Area 3")
region_names <- region_names[c(4,6,5,2,3,1)]
mu_alphaHG<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,1]",mcmc_names)]) # grab Group 1
mu_alphaNass<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,2]",mcmc_names)]) # grab Group 2
mu_alphaSkeena<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,3]",mcmc_names)]) # grab Group 3
mu_alphaHec<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,4]",mcmc_names)]) # grab Group 4
mu_alphaNC<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,5]",mcmc_names)]) # grab Group 5
mu_alphaCC<-colMeans(resultSR_B3[,grepl("\\bmu_lalpha",mcmc_names) & grepl("\\b,6]",mcmc_names)]) # grab Group 6

log_alpha <- data.frame("year"=unique(SR.dat$year),"ln_alpha"=c(mu_alphaHG,mu_alphaNass,mu_alphaSkeena,mu_alphaHec,mu_alphaNC,mu_alphaCC),"region"=rep(group_names,each=length(mu_alphaCC)),"pfma"=rep(region_names,each=length(mu_alphaCC)))
head(log_alpha)
library(ggplot2)
ggplot(log_alpha,aes(x=year,y=ln_alpha,fill=region)) +
  geom_point() +
  geom_line() +
  facet_wrap(~region) +
  scale_fill_brewer()

write.csv(log_alpha,"Data/Coho alpha trends.csv",row.names=FALSE)
