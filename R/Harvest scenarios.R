ER_values <- read.csv("Data/ER_values.csv", header=T)
Kgroups <- 6 # number of groups that we want to see correlated
ER <- aggregate(mean_ER~group_no+fishery,data=ER_values,mean)
ER_values$cv <- ER_values$sd_ER/ER_values$mean_ER
escape_reg <- matrix(0,nrow=Kgroups,ncol=5,dimnames=list("groups"=ER_values$group_name[ER_values$fishery=="bc_total"],"regs"=c("No harvest","10-year average","50% BC reduction","50% AK reduction","50% AK & BC reduction")))
escape_reg[,1] <- 0
escape_reg[,2] <- aggregate(mean_ER~group_no,data=ER_values,sum)$mean_ER
escape_reg[,3] <- sapply(1:Kgroups,function(x){ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + 0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_reg[,4] <- sapply(1:Kgroups,function(x){0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_reg[,5] <- sapply(1:Kgroups,function(x){0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + 0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})

ER_cv <- aggregate(cv~group_no+fishery,data=ER_values,mean)
escape_sd_reg <- matrix(0,nrow=Kgroups,ncol=5,dimnames=list("groups"=ER_values$group_name[ER_values$fishery=="bc_total"],"regs"=c("No harvest","10-year average","50% BC reduction","50% AK reduction","50% AK & BC reduction")))
escape_sd_reg[,1] <- 1e-6
escape_sd_reg[,2] <- aggregate(sd_ER~group_no,data=ER_values,sum)$sd_ER
escape_sd_reg[,3] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_sd_reg[,4] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})
escape_sd_reg[,5] <- sapply(1:Kgroups,function(x){ER_cv$cv[ER$group_no==x & ER$fishery=='alaska']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='alaska'] + ER_cv$cv[ER$group_no==x & ER$fishery=='bc_total']*0.5*ER$mean_ER[ER$group_no==x & ER$fishery=='bc_total']})


saveRDS(escape_reg,"Data/harvest scenarios.rds")
saveRDS(escape_sd_reg,"Data/harvest var scenarios.rds")