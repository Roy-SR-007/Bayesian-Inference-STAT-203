rm(list=ls())


# -------------------------------------------------------------------------

# CREDIBLE INTERVALS - BAYESIAN SET-UP

# Data/Likelihood: X|theta ~ Bin(n=10,theta)

# Prior: theta ~ Beta(a,b)

# Taking choices for the parameters of the
# prior distribution

a = c(0.2,0.3,0.5,1,1.5,2)
b = c(0.7,0.4,0.5,1,1.7,1.5)


credible_intervals_bayesian = function(para_prior,para_data,N)
{
  n = para_data
  
  a = para_prior[,1]
  b = para_prior[,2]

  theta = list(0) # to store the prior

  x_theta = list(0) # to store the data (X|theta)

  for(i in 1:length(a))
  {
    # prior distribution: N observations ~ Beta(a,b)
    theta[[i]] = rbeta(N,a[i],b[i])
  }

  for(i in 1:length(a))
  {
    x = array(0)
    for(j in 1:N)
    {
      x[j] = rbinom(1,n,theta[[i]][j])
    }
    # Data (X|theta) ~ N observations from Bin(n=10,theta)
    x_theta[[i]] = x
  }

  # Credible Intervals, using the posterior theta|x

  u_x = function(n,x,a,b,alpha)
  {
    return(qbeta(1-(alpha/2),x+a,n-x+b))
  }

  l_x = function(n,x,a,b,alpha)
  {
    return(qbeta(alpha/2,x+a,n-x+b))
  }

  # Now to fix the values of x at say, x* and 
  # determine the ratio:

  # (# theta_j in [l(x*),u(x*)])/N(x*)

  # We will be fixing x's to x*'s, for each 
  # choice of the parameters of the prior
  # distribution.

  f = list(0)

  for(i in 1:length(a))
  {
    x_star_choices = sort(unique(x_theta[[i]])) # fixing the x*'s
  
    ff = array(0)
  
    for(j in 1:length(x_star_choices))
    {
      ct = 0
      l_xstar = l_x(n,x_star_choices[j],a[i],b[i],0.05)
      u_xstar = u_x(n,x_star_choices[j],a[i],b[i],0.05)
    
      theta_pos = which(x_theta[[i]] == x_star_choices[j])
    
      for(k in 1:as.vector(table(x_theta[[i]]))[j])
      {
        if(theta[[i]][theta_pos[k]] >= l_xstar & theta[[i]][theta_pos[k]] <= u_xstar)
        {
          ct = ct + 1
        }
      }
      ff[j] = ct/(as.vector(table(x_theta[[i]]))[j])
    }
    f[[i]] = ff # approx. 1-alpha = 0.95
  }
  return(f)
}


# -------------------------------------------------------------------------

# Test Case

n=40

a = c(0.15,0.20,0.25,0.30,0.50,0.55,0.65,1.50,2.00,1)
b = c(3.00,2.00,1.20,1.00,0.50,0.45,0.25,0.20,0.10,1)
x = c(a,b)

set.seed(17)

N = 1000000 # sample size

para = matrix(x,byrow=F,nrow = length(a),ncol=2)

res = credible_intervals_bayesian(para,n,N)

# -------------------------------------------------------------------------


