co_pops<-read.table("Data/coho_groups.txt",header=TRUE)
co_pops$mean_total<-NA

for(i in 1:52){
  dat<-subset(SR.dat,SR.dat$pop_no==i)
  co_pops[i,5]<-mean(na.omit(dat$total_runE))
}

### lumping Rivers Smith Inlet with Area 7-8
co_pops[which(co_pops$group==7),4]<-6
n.pops<-max(SR.dat$pop_no)
group <- co_pops$group
Npops <- length(group)
Kgroups <- 6 # number of groups that we want to see correlated
covarGroups <- matrix(0,nrow=Kgroups,ncol=Kgroups)
diag(covarGroups) <- 1

resultSR_B3<-readRDS("Results/COSR3B.tr1_lalpha_MCMC.rds")
mcmc_names <- colnames(resultSR_B3[[1]])
#######################################################################
#### estimating Smsy using lamberts w Scheuerell method
Smsy<-Sgen<-Umsy<-matrix(data=NA,nrow=samples,ncol=n.pops)

samp.chain<-sample(1:n_chains,samples,replace=TRUE)
samp.MCMC<-sample(1:samples, samples, replace=FALSE)
# Sgen; solve Smsy=Sgen*exp(a*(1- Sgen/S(k))) from Holt et al. 2009 Indicators of Status and Benchmarks for Conservation Units in Canada's Wild Salmon Policy 
## just need to adjust the length based on how long the MCMC is
Sgen_find <- function(Smsy,Sgen,a,b)
{
  (Smsy-(Sgen*exp(a*(1- Sgen/(a/b)))))^2
}
umsy_find <- function(Smsy,Umsy,a,b)
{
  #Umsy <- Umsy_priors$mean[1]
  #Smsy <- Smsy_priors$mean[1]
  #a <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[1]]
  #b <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[1]]
  (Smsy-(1-Umsy)*(a/b))^2
}
for (j in 1:n.pops){
  samp.chain<-sample(1:n_chains,samples,replace=TRUE)
  samp.MCMC<-sample(1:samples, samples, replace=FALSE)
  
  for (i in 1:samples){
    alpha <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("ln_alpha.mu",mcmc_names)[j]]
    beta <- resultSR_B3[[samp.chain[i]]][samp.MCMC[i],grep("beta",mcmc_names)[j]]
    Smsy[i,j]<-(1 - lambert_W0(exp(1 - (alpha)))) / (beta) #Smsy in 1000 random draws from MCMC
    Sgen[i,j] <- optimize(Sgen_find,c(0,Smsy[i,j]),tol=0.0001,Smsy=Smsy[i,j],a=alpha,b=beta)$minimum
    Umsy[i,j] <- optimize(umsy_find,c(0,1),tol=0.0001,Smsy=Smsy[i,j],a=alpha,b=beta)$minimum
  }
}

Smsy_priors <- apply(Smsy,2,FUN=function(x){c(mean(x),sd(x))})
Smsy_priors <- data.frame("pop"=1:n.pops,"mean"=Smsy_priors[1,],"tau"=1/(Smsy_priors[2,]^2))
Sgen_priors <- apply(Sgen,2,FUN=function(x){c(mean(x),sd(x))})
Sgen_priors <- data.frame("pop"=1:n.pops,"mean"=Sgen_priors[1,],"tau"=1/(Sgen_priors[2,]^2))
Umsy_priors <- apply(Umsy,2,FUN=function(x){c(mean(x),sd(x))})
Umsy_priors <- data.frame("pop"=1:n.pops,"mean"=Umsy_priors[1,],"tau"=1/(Umsy_priors[2,]^2))
saveRDS(Smsy_priors,"Results/Smsy.rds")
saveRDS(Sgen_priors,"Results/Sgen.rds")
saveRDS(Umsy_priors,"Results/Umsy.rds")
