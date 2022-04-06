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



# Classical Inference - Confidence Intervals ------------------------------

## theta ~ Beta(a,b)
## X|theta ~ Bin(n,theta)




# Wald's Large Sample Intervals -------------------------------------------
rm(list=ls())

set.seed(4)

alpha=.4
beta=2

N=10^6
n=200

theta=rbeta(N,alpha,beta)
x=c()

for(i in 1:length(theta))
{
  x[i]=rbinom(1,n,theta[i])
}
theta=round(theta,3)
theta

theta_star=0.530
x_fix=c()
for(i in 1 : length(x))
{
  if(theta[i]==theta_star)
  {
    x_fix=c(x_fix,x[i])
  }
}
length(x_fix)
  
theta_hat=x_fix/n
theta_hat

l1=theta_hat-(qnorm(0.975))*(sqrt((theta_hat*(1-theta_hat))/n))
u1=theta_hat+(qnorm(0.975))*(sqrt((theta_hat*(1-theta_hat))/n))

k=0
for(i in 1:length(x_fix))
{
  if(l1[i]<theta_star & theta_star<u1[i]){
    k=k+1
  }
}
k
confidence_coeff1=k/length(x_fix)
confidence_coeff1

# Agresti-Coull's Intervals -------------------------------------------
rm(list=ls())

set.seed(4)

alpha=0.4
beta=2

N=10^6
n=50+10+40+100

theta=rbeta(N,alpha,beta)
x=c()

for(i in 1:length(theta))
{
  x[i]=rbinom(1,n,theta[i])
}
theta=round(theta,3)
theta

theta_star=0.530
x_fix=c()
for(i in 1 : length(x))
{
  if(theta[i]==theta_star)
  {
    x_fix=c(x_fix,x[i])
  }
}
length(x_fix)

theta_hat_hat=(x_fix+2)/(n+4)
theta_hat_hat

l1=theta_hat_hat-(2)*(sqrt((theta_hat_hat*(1-theta_hat_hat))/n))
u1=theta_hat_hat+(2)*(sqrt((theta_hat_hat*(1-theta_hat_hat))/n))

k=0
for(i in 1:length(x_fix))
{
  if(l1[i]<theta_star & theta_star<u1[i]){
    k=k+1
  }
}
k
confidence_coeff1=k/length(x_fix)
confidence_coeff1

# Pearson-Clopper's Intervals -------------------------------------------
rm(list=ls())

set.seed(4)

alpha=0.3+0.1
beta=0.1+0.9+1

N=10^6
n=10+40+50+100

theta=rbeta(N,alpha,beta)
x=c()

for(i in 1:length(theta))
{
  x[i]=rbinom(1,n,theta[i])
}
theta=round(theta,3)
theta

theta_star=0.530
x_fix=c()
for(i in 1 : length(x))
{
  if(theta[i]==theta_star)
  {
    x_fix=c(x_fix,x[i])
  }
}
length(x_fix)


l1=qbeta(0.025,x_fix,n-x_fix+1)
u1=qbeta(0.975,x_fix+1,n-x_fix)

k=0
for(i in 1:length(x_fix))
{
  if(l1[i]<theta_star & theta_star<u1[i]){
    k=k+1
  }
}
k
confidence_coeff1=k/length(x_fix)
confidence_coeff1

