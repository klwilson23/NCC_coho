SR.dat<-read.table("Data/Coho_Brood_MASTER.txt", header=T)
SR.dat_full <- SR.dat
SR.dat_predict<-subset(SR.dat,SR.dat$year>=2017)
#SR.dat<-subset(SR.dat,SR.dat$year<2017)

SR.dat$logRS1<-log(SR.dat$RS_E)
SR.dat$logRS2<-log(SR.dat$RS_2)
SR.dat$logRS3<-log(SR.dat$RS_3)
pop_names <- unique(SR.dat$population[order(SR.dat$pop_no)])
group_names <- c("Area 7-10","Hecate Lowlands - Area 5/6","Area 6 - Inner Waters","Haida Gwaii - Area 2E","Skeena - Area 4","Nass - Area 3")
group_names <- group_names[c(4,6,5,2,3,1)]

ER_values <- readRDS("Data/harvest scenarios.rds")

SR.dat$er_est <- SR.dat$er_2
SR.dat$er_est[is.na(SR.dat$er_est)] <- SR.dat$er_E[is.na(SR.dat$er_est)]

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
coho_ts <- left_join(SR.dat,alpha_annual_sst)

coho_ts$run_size <-ifelse(!is.na(coho_ts$RS_2),coho_ts$RS_2,coho_ts$RS_E)