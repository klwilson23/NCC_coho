rmse_rw <- readRDS(file="Results/predictive error random walk.rds")
rmse_ar1 <- readRDS(file="Results/predictive error ar1.rds")
rmse_tr1 <- readRDS(file="Results/predictive error tr1.rds")

rmse <- data.frame("model"=factor(c(rep("Random walk",ncol(rmse_rw)),rep("Mean recursive",ncol(rmse_rw)),rep("Trending",ncol(rmse_tr1)))),"error"=c(rmse_rw[2,],rmse_ar1[2,],rmse_tr1[2,]),"population"=rep(colnames(rmse_rw),3),"scale"="Population")
mean_error <- aggregate(error~model,rmse,sum,na.rm=TRUE)

rmse <- rbind(rmse,data.frame(mean_error,"population"="aggregate","scale"="Aggregate"))

boxplot(error~model,data=rmse)
plot(error~model,data=rmse,type="p")
sub_rmse <- rmse[rmse$model!="Trending",]
segments(x0=sub_rmse$model,x1=sub_rmse$model,y0=sub_rmse$error,y1=sub_rmse$error)
sub_rmse <- rmse[rmse$model!="Random walk",]
segments(x0=sub_rmse$model,x1=sub_rmse$model,y0=sub_rmse$error,y1=sub_rmse$error)
aggregate(error~model,rmse[rmse$scale=="Aggregate",],sum)
100*((33894-33152)/33152)
100*((33894-32712)/32712)
rmse[rmse$scale=="Aggregate",]

library(ggplot2)

ggplot(rmse,aes(x=model,y=error,fill=model,group=population)) +
  geom_line() +
  geom_point(pch=21) +
  facet_wrap(~scale,scales="free") +
  ylab("Median absolute error") +
  xlab("Forecasted productivity model") +
  scale_fill_brewer(name="Model",palette='Set2') +
  labs(title="North & Central Coast coho salmon",
       subtitle="Predictive accuracy for 2017-2020 escapement",
       caption="Data source: DFO") +
  theme_minimal() +
  theme(legend.position="top",strip.text.x=element_text(hjust=0))

ggsave("Figures/predictive accuracy.jpeg",height=5,width=8,units="in")
