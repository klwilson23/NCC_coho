model{
  
  ## Stock-recruit component
  for(i in 1:Npops){
    sigma_s0[i] ~ dt(0, pow(2.5,-2), 1)T(0,) # vague half-Cauchy prior for standard deviations
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dt(0, pow(2.5,-2), 1)T(0,) # vague half-Cauchy prior for standard deviations
    beta[i] <- 1/Smax[i]
    Smax[i] ~ dlnorm(Smax.p[i],1/(Smax.tau[i]^2))
    #spawn_int[i] ~ dnorm(0,1e-3)
    spawners[1,i] ~ dlnorm(log(init_s0[i]),pow(sigma_s0[i], -2))
    logRS[1,i] ~ dnorm(mu.y[1,i], tau[i])
    mu.y[1,i] <- lalpha[1,i] - beta[i] * spawners[1,i]
    lalpha[1,i] ~ dnorm(mu_lalpha[1,group[i]], pow(sigma_alpha[group[i]],-2)) # alpha must be positive
    alpha[1,i] <- log(exp(lalpha[1,i]+0.5*(sigma_alpha[group[i]]^2))) # bias-corrected log-alpha for population i
    
    for(t in 2:n.years) {
      spawners[t,i] ~ dlnorm(log(spawners[t-1,i]),pow(sigma_s0[i], -2))T(0,)
      logRS[t,i] ~ dnorm(mu.y[t,i], tau[i])
      mu.y[t,i] <- lalpha[t,i] - beta[i] * spawners[t,i]
      lalpha[t,i] ~ dnorm(mu_lalpha[t,group[i]], pow(sigma_alpha[group[i]],-2)) # alpha must be positive
      alpha[t,i] <- log(exp(lalpha[t,i]+0.5*(sigma_alpha[group[i]]^2))) # bias-corrected log-alpha for population i
    }
    
    ln_alpha.mu[i] <- mean(lalpha[1:n.years,i]) # long-term average log-alpha
  }
  
  mu_lalpha[1,1:Kgroups] ~ dmnorm(init_lalpha,tauMVN) # initialize the log-alpha 6x6 matrix
  for(z in 2:n.years)
  {
    for(k in 1:Kgroups)
    {
      mn_lalpha_tv[z,k] <- mu_lalpha[z-1,k] #+mean_alpha[k]
    }
    mu_lalpha[z,1:Kgroups] ~ dmnorm(mn_lalpha_tv[z,1:Kgroups],tauMVN) # make log-alpha recursive, and correlated to each group
  }
  
  for(z in 1:n.years)
  {
    for(k in 1:Kgroups)
    {
      mu_alpha[z,k] <- log(exp(mu_lalpha[z,k] + 0.5*sigma_group[k])) # global bias-corrected mean alphas for each of the 6 groups
    }
  }
  
  #u_alpha ~ dunif(-0.99,0.99)
  #u_spawn ~ dunif(-0.99,0.99)
  
  tauMVN ~ dwish(covarGroups,Kgroups+1)
  sigmaGroups <- inverse(tauMVN)
  for(k in 1:Kgroups)
  {
    #mean_alpha[k] ~ dnorm(0,1e-3*(1-u_alpha*u_alpha))
    sigma_group[k] <- sigmaGroups[k,k]
    sigma_alpha[k] ~ dt(0, pow(2.5,-2), 1)T(0,) # vague half-Cauchy prior on standard deviation
  }
}