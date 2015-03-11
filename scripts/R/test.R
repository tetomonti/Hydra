monotonicity.binomial <- function()
{
    #var Nsamples, Ngroups, mu[Ngroups], p[Nsamples], alpha, tau, diff, count, P;
    
    for ( i in 1:Nsamples ) {
        y[i] ~ dbin(p[i],Nsamples)
        logit(p[i]) <- mu[groups[i]]
    }
    for ( j in 1:Ngroups ) {
        mu[j] ~ dnorm(0,tau)
    }
    tau ~ dgamma( alpha0, beta0 )

    diff <- mu[2:Ngroups] - mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) )    

    for ( j in 1:Ngroups ) {
        P[j] <- exp(mu[j])/(exp(mu[j])+1)
    }
}
binomial <- function()
{
    for ( i in 1:Nsamples )
    {
        y[i] ~ dbin(p[i],Nsize)
        logit(p[i]) <- alpha
    }
    alpha ~ dnorm(0,.01)

    P <- exp(alpha)/(exp(alpha)+1)
}
monotonicity.beta <- function() # UNSTABLE
{    
    for ( i in 1:Nsamples )
    {
        y[i] ~ dbeta(alpha[groups[i]],beta[groups[i]])
    }
    for ( j in 1:Ngroups )
    {
        alpha[j] ~ dgamma(alpha0,beta0)
        beta[j] ~ dgamma(alpha0,beta0)
    }
    mu <- alpha/(alpha+beta)
    diff <- mu[2:Ngroups] - mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) )
}

if ( FALSE )
{
    source("~/Research/CBMgithub/scripts/R/test.R")
    
    set.seed(123)
    Y <- c(rbinom(50,200,.2),
           rbinom(50,200,.4),
           rbinom(50,200,.6),
           rbinom(50,200,.8))
    G <- factor(rep(1:4,each=50),levels=1:4)   

    jags.input <- list(y=Y[G==1],Nsamples=sum(G==1),Nsize=200)
    jags.inits <- list(alpha=logit(.2))                   
    jags.inits <- NULL
    
    OUT <- run.jags(model=binomial,data=jags.input,inits=jags.inits,monitor=c('alpha','P'))
   
    jags.input <- list(y=Y,Nsamples=length(Y),groups=G,Ngroups=nlevels(G),query=c(1,1,1),alpha0=.1,beta0=.1)
    jags.inits <- list(mu=logit(c(.2,.4,.6,.8)),tau=.5)                   
    jags.inits <- NULL
    
    OUT <- run.jags(model=monotonicity.binomial,data=jags.input,inits=jags.inits,monitor=c('count','P','mu','tau'))

    jags.input <- list(y=round(dat*length(dat)),Nsamples=length(dat),groups=groups,Ngroups=nlevels(groups),
                       query=c(1,1,1),alpha0=.1,beta0=.1)
    jags.inits <- list(mu=logit(tapply(dat,groups,mean)),tau=.5)

    OUT <- run.jags(model=monotonicity.binomial,data=jags.input,inits=jags.inits,monitor=c('count','P','mu','tau'))

    jags.input <- list(y=dat,Nsamples=length(dat),groups=groups,Ngroups=nlevels(groups),
                       query=c(1,1,1),alpha0=.1,beta0=.1)
    M1 <- tapply(dat,groups,mean)
    M2 <- tapply(dat,groups,sd)
    A <- M1 * ( M1*(1-M1)/M2 - 1)
    B <- (1-M1) * ( M1*(1-M1)/M2 - 1)
    jags.inits <- list(alpha=A,beta=B)

    OUT <- run.jags(model=monotonicity.beta,data=jags.input,inits=jags.inits,monitor=c('count','mu'))

}

if ( FALSE )
{
beta.estimate <- function()
{
    for ( i in 1:N )
    {
        y[i] ~ dbeta(alpha,beta)
    }
    alpha ~ dlnorm(0,.01)
    beta ~ dlnorm(0,.01)
}
beta.estimate <- function()
{
    for ( i in 1:N )
    {
        y[i] ~ dbeta(alpha[i],beta[i])
        alpha[i] <- mu[i] * phi
        beta[i]  <- (1-mu[i]) * phi
        logit(mu[i]) <- a
    }
    phi ~ dgamma(.1,.1)
    a ~ dnorm(0,.001)
    M <- exp(a)/(1+exp(a))
}

set.seed(123)
Y <- rbeta(100,2,5)

OUT <- run.jags(model=beta.estimate,data=list(y=Y,N=length(Y)),monitor=c('alpha','beta','M','phi'))
OUT <- run.jags(model=beta.estimate,data=list(y=dat[groups==1],N=sum(groups==1)),monitor=c('alpha','beta'))
}
if ( FALSE )
{
    setwd('~/Research/CBMgithub/scripts/R/')
    source('jags.R')
    require(rjags)

    logit.fun <- function(P) log(P/(1-P))
    logit.inv <- function(Y) exp(Y)/(exp(Y)+1)
    
    set.seed(123)
    y <- rnorm(1000)
    y <- logit.inv(y)

    idx <- 1:50
    mOut <- jags.model('test.logit.R',data=list(y=y[idx],Nsamples=length(y[idx])),n.chains=1,n.adapt=1000)

    update(mOut,10000)

    ## actual MCMC sampling
    OUT <- jags.samples(model=mOut,variable.names=c('mu','P','tau'),n.iter=10000)

}
