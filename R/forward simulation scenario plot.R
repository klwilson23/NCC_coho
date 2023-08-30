productivity <- function(scenario="rw",ln_a,rho=0,trend_1=0,trend=0,mean_alpha=0,nyears,nsims,sigma,seed=1)
{
  set.seed(seed)
  ln_t <- rep(NA,nyears+nsims)
  ln_t[1] <- ln_a
  u_t <- rnorm(nyears+nsims,mean=0,sd=sigma)
  for(i in 2:nyears)
  {
    ln_t[i] <-  (trend_1*i+ln_a)+rho*(ln_t[i-1]-(trend_1*(i-1)+ln_a))+u_t[i]
  }
  if(scenario=="rw")
  {
    for(i in (nyears+1):(nyears+nsims))
    {
      ln_t[i] <- ln_t[i-1]+u_t[i]
    }
  }
  if(scenario=="ar1")
  {
    for(i in (nyears+1):(nyears+nsims))
    {
      ln_t[i] <- rho*(ln_t[i-1])+mean_alpha+u_t[i]
    }
  }

  if(scenario=="trend")
  {
    for(i in (nyears+1):(nyears+nsims))
    {
      ln_t[i] <- (trend*i+ln_a)+rho*(ln_t[i-1]-(trend*(i-1)+ln_a))+u_t[i]
    }
    
  }
  return(ln_t)
}

seed <- rpois(1,2020)
ln_a <- 4
sigma <- 0.2
rho <- 0.9
trend_1 <- -0.05
trend <- -0.05
mu <- 4
mean_alpha <- mu*(1-rho) # if mu is the target, then mean_alpha/(1-rho) is the mean reversion
nyears <- 40
nsims <- 15
rw <- productivity(scenario="rw",ln_a=ln_a,rho=rho,trend_1=trend_1,trend=trend,mean_alpha=mean_alpha,nyears=nyears,nsims=nsims,sigma=sigma,seed=seed)
revert <- productivity(scenario="ar1",ln_a=ln_a,rho=rho,trend_1=trend_1,trend=trend,mean_alpha=mean_alpha,nyears=nyears,nsims=nsims,sigma=sigma,seed=seed)
tr1 <- productivity(scenario="trend",ln_a=ln_a,rho=rho,trend_1=trend_1,trend=trend,mean_alpha=mean_alpha,nyears=nyears,nsims=nsims,sigma=sigma,seed=seed)
future_years <- (nyears):(nyears+nsims)

plot(rw[-(future_years[-1])],type="l",ylim=range(rw,revert,tr1),xlim=range(0,future_years),xlab="Years",ylab="Productivity")
lines(future_years,rw[future_years],col="darkgreen",lwd=2)
lines(future_years,revert[future_years],col="dodgerblue",lwd=2)
lines(future_years,tr1[future_years],col="orange",lwd=2)
legend("bottomleft",title="AR(1) Scenario",c("Mean reverting","Random walk","Trend"),bty="n",lty=1,lwd=2,col=c("darkgreen","dodgerblue","orange"))

