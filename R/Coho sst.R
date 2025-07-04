########################################################################################
# Coho sst.R
# quick and dirty exploration of sst-alpha relationships 
# code for sst data scraping and manipulations largely compliments of M. Malick (NOAA) 
# see: https://github.com/michaelmalick/r-ersst
#########################################################################################

library(plyr)
library(tidyverse)
library(broom)
source("R/sst functions.R")

## download and process SST data UNCOMMENT IF YOU NEED TO DOWNLOAD RAW SST AND PROCESS IT

#ersst::sst_download(years = 1950:2018,
#                    months = 1:12,
#                    save.dir = "./data/sst_raw/",
#                    version = 5)
#
#sst.raw.full <- ersst::sst_load(years = 1950:2018,
#                                months = 1:12,
#                                read.dir = "./data/sst_raw/",
#                                version = 5)
#
#sst.raw.np <- ersst::sst_subset_space(sst.raw.full,
#                                      lat.min = 36,
#                                      lat.max = 80,
#                                      lon.min = 170,
#                                      lon.max = 250)
#
#sst.raw.df <- ersst::sst_dataframe(sst.raw.np)
#
#write.csv(sst.raw.df, "Data/sst_raw.csv", row.names = FALSE)

# sst.raw <- read.csv("Data/sst_raw.csv")
# head(sst.raw)
# tail(sst.raw)
# sapply(sst.raw, class)
# summary(sst.raw)

## calculate SST anomalies and average across specified period and region 

# sst.anom <- sst.anomaly(sst.raw, ref.years = 1950:2010)
# head(sst.anom)
# tail(sst.anom)
# summary(sst.anom)
# sapply(sst.anom, class)
# 
# ## Convert longitude to match salmon data
# sst.anom$lon2 <- ifelse(sst.anom$lon > 180, sst.anom$lon - 360, sst.anom$lon)
# 
# write.csv(sst.anom, "Data/sst_raw_anomalies.csv", row.names=F)

## Read in static ocean entry points by region and sst anomalies
ocean_entry <- read.csv("Data/static_ocean_entry.csv")
sst_anom <- read.csv("Data/sst_raw_anomalies.csv")

## Calculate average SST anomaly within area where stock spends few months of marine life 
summer_sst_stock_anomalies <- sst.averager(ocean_entry, sst_anom, distance = 200, which.months = c(5:9))
summer_sst_stock_anomalies <- summer_sst_stock_anomalies %>%
  mutate(region = case_when(stock.id == 1 ~ "Central Coast (South)",
                          stock.id == 2 ~ "Haida Gwaii",
                          stock.id == 3 ~ "Hecate Lowlands",
                          stock.id == 4 ~ "Inner Waters",
                          stock.id == 5 ~ "Nass",
                          stock.id == 6 ~ "Skeena"))

## Calculate average SST anomaly within area where stock spends first year of marine life 
annual_sst_stock_anomalies <- sst.averager(ocean_entry, sst_anom, distance = 200, which.months = c(1:12))
annual_sst_stock_anomalies <- annual_sst_stock_anomalies %>%
  mutate(region = case_when(stock.id == 1 ~ "Central Coast (South)",
                            stock.id == 2 ~ "Haida Gwaii",
                            stock.id == 3 ~ "Hecate Lowlands",
                            stock.id == 4 ~ "Inner Waters",
                            stock.id == 5 ~ "Nass",
                            stock.id == 6 ~ "Skeena"))

## merge with time-varying alphas 
tv_alpha <- read.csv("Data/Coho alpha trends.csv")
tv_alpha$brood_yr <- tv_alpha$year
tv_alpha$year <- tv_alpha$year+2 #align year with ocean entry year assuming one year in freshwater

alpha_summer_sst <- left_join(tv_alpha,summer_sst_stock_anomalies)
alpha_summer_sst$region <- factor(alpha_summer_sst$region, levels = c("Central Coast (South)","Hecate Lowlands","Inner Waters",
                                                                       "Skeena","Haida Gwaii","Nass"))
alpha_annual_sst <- left_join(tv_alpha,annual_sst_stock_anomalies)
alpha_annual_sst$region <- factor(alpha_annual_sst$region, levels = c("Central Coast (South)","Hecate Lowlands","Inner Waters",
                                                                      "Skeena","Haida Gwaii","Nass"))

## what do ssts and anomalies look like by region over time? 
ggplot() +
  geom_line(data = alpha_summer_sst, aes(x=year, y = sst, color = region), size=1.5) +
  scale_color_viridis_d() +
  xlab("Year") +
  ylab("Summer SST") +
  theme_bw() 
ggsave("Figures/SST/summer_sst.jpeg")

ggplot() +
  geom_line(data = alpha_annual_sst, aes(x=year, y = sst, color = region), size=1.5) +
  scale_color_viridis_d() +
  xlab("Year") +
  ylab("Annual SST") +
  theme_bw() 
ggsave("Figures/SST/annual_sst.jpeg")

ggplot() +
  geom_line(data = alpha_summer_sst, aes(x=year, y = sst.anom, color = region), size=1.5) +
  scale_color_viridis_d() +
  xlab("Year") +
  ylab("Summer SST anomaly") +
  theme_bw() 

ggplot() +
  geom_line(data = alpha_annual_sst, aes(x=year, y = sst.anom, color = region), size=1.5) +
  scale_color_viridis_d() +
  xlab("Year") +
  ylab("Annual SST anomaly") +
  theme_bw() 

## what do alpha~sst relationships look like by region over time? 
ggplot() +
  geom_point(data = alpha_summer_sst, aes(x= sst.anom, y = ln_alpha, color=year), size=2) +
  geom_smooth(data = alpha_summer_sst, aes(x= sst.anom, y = ln_alpha), method=lm, color="dark grey") +
# geom_smooth(data = alpha_summer_sst, aes(x= sst.anom, y = ln_alpha), method=loess) +
  scale_color_viridis_c() +
  xlab("Summer SST anomaly") +
  ylab("ln alpha") +
  facet_wrap(~region, ncol=2) +
  theme_bw() 
ggsave("Figures/SST/summer_sst_prod.jpeg")

ggplot() +
  geom_point(data = alpha_annual_sst, aes(x= sst.anom, y = ln_alpha, color=year), size=2) +
  geom_smooth(data = alpha_annual_sst, aes(x= sst.anom, y = ln_alpha), method=lm, color="dark grey") +
# geom_smooth(data = alpha_annual_sst, aes(x= sst.anom, y = ln_alpha), method=loess) +
  scale_color_viridis_c() +
  xlab("Annual SST anomaly") +
  ylab("ln alpha") +
  facet_wrap(~region, ncol=2) +
  theme_bw() 
ggsave("Figures/SST/annual_sst_prod.jpeg")

## and how about some stats?
alpha_summer_sst %>% 
  group_by(region) %>%
  dplyr::summarize(correlation = cor(sst.anom, ln_alpha))

summer_sst_lms <- alpha_summer_sst %>%
  nest(data = -region) %>% 
  mutate(
    fit = map(data, ~ lm(ln_alpha ~ sst.anom, data = .x)),
    tidied = map(fit, tidy,conf.int=TRUE)
  ) %>% 
  unnest(tidied)

newdat <- alpha_summer_sst[alpha_summer_sst$brood_yr%in%c(1980,2016),]
summer_decline <- rep(NA,length(unique(summer_sst_lms$region)))
names(summer_decline) <- unique(summer_sst_lms$region)
for(i in 1:length(summer_decline))
{
  ii <- unique(summer_sst_lms$region)[i]
  sub_df <- summer_sst_lms[summer_sst_lms$region==ii,]
  pred <- sub_df$estimate[sub_df$term=="(Intercept)"] + sub_df$estimate[sub_df$term=="sst.anom"]*newdat$sst.anom[newdat$region==ii]
  summer_decline[i] <- 100*((pred[2]-pred[1])/pred[1])
}

summer_sst_lms %>%
  filter(term=="sst.anom")

alpha_annual_sst %>% 
  group_by(region) %>%
  dplyr::summarize(correlation = cor(sst.anom, ln_alpha))

newdat <- alpha_annual_sst[alpha_annual_sst$brood_yr%in%c(1980,2016),]

annual_sst_lms <- alpha_annual_sst %>%
  nest(data = -region) %>% 
  mutate(
    fit = map(data, ~ lm(ln_alpha ~ sst.anom, data = .x)),
    tidied = map(fit, tidy,conf.int=TRUE),
  ) %>% 
  unnest(tidied)

annual_decline <- rep(NA,length(unique(annual_sst_lms$region)))
names(annual_decline) <- unique(annual_sst_lms$region)
for(i in 1:length(annual_decline))
{
  ii <- unique(annual_sst_lms$region)[i]
  sub_df <- annual_sst_lms[annual_sst_lms$region==ii,]
  pred <- sub_df$estimate[sub_df$term=="(Intercept)"] + sub_df$estimate[sub_df$term=="sst.anom"]*newdat$sst.anom[newdat$region==ii]
  annual_decline[i] <- 100*((pred[2]-pred[1])/pred[1])
}
# pred <- function(x,  ...){
#   z <- predict.lm(x, se.fit = TRUE, ...)
#   as.data.frame(z[1:2])
# }
# 
# annual_sst_pred <- alpha_annual_sst %>%
#   group_by(region) %>%
#   nest() %>% 
#   mutate(m1 = purrr::map(.x = data, .f = ~ lm(ln_alpha ~ sst.anom, data = .))) %>%
#   mutate(decline = purrr::map(.x = m1, ~ pred(.))) %>% 
#   select(data,region, decline) %>%
#   unnest(decline)
# alpha_annual_df <- left_join(alpha_annual_sst,annual_sst_pred,by=c("year","region"))
# 
# alpha_annual_df %>%
#   filter(brood_yr==1982 | brood_yr == 2018)

annual_sst_lm_stats <- annual_sst_lms %>%
  filter(term=="sst.anom")

annual_sst_cors <- alpha_annual_sst %>% 
  group_by(region) %>%
  dplyr::summarize(correlation = cor(sst.anom, ln_alpha))

summer_sst_lm_stats <- summer_sst_lms %>%
  filter(term=="sst.anom")

summer_sst_cors <- alpha_summer_sst %>% 
  group_by(region) %>%
  dplyr::summarise(correlation = cor(sst.anom, ln_alpha))

bind_rows(left_join(annual_sst_lm_stats,annual_sst_cors),
          left_join(summer_sst_lm_stats,summer_sst_cors))
