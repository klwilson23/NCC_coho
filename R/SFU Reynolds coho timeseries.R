library(brms)

salmon <- read.csv("Data/SFUReachCountData.csv")
salmon$date <- as.Date(with(salmon, paste(year, month, day,sep="-")), "%Y-%m-%d")
salmon$yday <- lubridate::yday(salmon$date)
salmon$coho <- rowSums(cbind(salmon$co.dead.total,salmon$co.live.total),na.rm=TRUE)

df <- subset(salmon,subset=salmon$section=="total",select=c(site,yday,year,section,coho
))
df$year_num <- as.numeric(as.factor(df$year))
df$yday <- as.numeric(as.factor(df$yday))
str(df)

nsamps <- aggregate(year_num~site,df,FUN=function(x){length(unique(x))})
sites_n <- nsamps[nsamps$year_num>=3,]
df$n <- nsamps$year_num[match(df$site,nsamps$site)]
df_sub <- df[df$n>=5,]
m1 <- brm(coho~year_num+(1|site),data=df_sub,family=zero_inflated_poisson(),cores=4)
summary(m1)
m2 <- brm(coho~year_num+(1|site),data=df_sub,family=poisson,cores=4)
summary(m2)
m3 <- brm(bf(coho ~ year_num+(1|site), zi ~ year_num),data=df_sub,family=zero_inflated_poisson(),cores=4)
m4 <- brm(bf(coho ~ year_num+(1|site), zi ~ year_num+(1|site)),data=df_sub,family=zero_inflated_poisson(),cores=4)

m1 <- add_criterion(m1, "loo")
m2 <- add_criterion(m2, "loo")
m3 <- add_criterion(m3, "loo")
m4 <- add_criterion(m4, "loo")

loo_compare(m1,m2,m3,m4, criterion = "loo")
plot(conditional_effects(m4), points = TRUE, ask = FALSE)
brms::fixef(m1)
nyears <- max(df$n)
100*(exp(1.78-0.07*nyears)-exp(1.78-0.07*1))/exp(1.78-0.07*1)
plot(m1)
brms::stancode(m1)
predictions <- predict(m4,type="r")
df_sub$predict <- predictions[,1]
head(df_sub)
library(ggplot2)
ggplot(data = df_sub, aes(y = coho, x = year)) +
  geom_point(aes(y = coho, x = year)) +
  geom_line(aes(x = year, y = predict)) +
  facet_wrap(~site)


zeros <- aggregate(coho~site,data=df,FUN=function(x){sum(x==0)/length(x)})
sum(zeros$coho==1);nrow(zeros)

aggregate(coho~site,df_sub[df_sub$year!=2021,],median)
aggregate(coho~site,df_sub[df_sub$year==2021,],sum)
