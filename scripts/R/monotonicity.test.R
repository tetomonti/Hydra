## MONOTONICITY TEST

monotonicity.test <- function
(
    dat,            # n-sized vector
    groups,         # n-sized vector indicating group membership of each of dat's observations
    query,          # vector indicating > (1) or < (-1) relationship between subsequent groups
    model=          # distribution models
        c('normal','lognormal','binomial','beta','logit'),
    model.specs=    # additional parameters for the models' specification
        list(tau0=.001,alpha0=.1,beta0=.1),
    nadapt=500,
    nburnin=1000,
    niter=1000,
    epsilon=1.0e-5, # 'zero' probability for the beta and logit models
    verbose=TRUE,   # output
    no.del=FALSE    # don't delete the JAGS file once completed
)
{
    ## BEGIN checks on inputs
    ##
    model <- match.arg(model)
    ngroups <- length( unique(groups) )
    if ( length(query)!=(ngroups-1) )
        stop( 'query must be vector of length ngroups-1:', length(query) )
    if ( !all(query %in% c(-1,1)) )
        stop( 'query must be a vector of {-1,1} only' )
    if ( any(sort(unique(groups)!=(1:ngroups))) )
        stop( 'groups must be a vector of indices btw 1 and ngroups' )
    if ( any( table(groups)<2 ) )
        stop( 'each group must have at least two observations' )
    if ( any( table(groups)<3 ) )
        warning( 'at least one group with < three observations' )
    ## END checks

    ## model-specific specs
    ##
    model.data <- model.inits <- NULL
    if ( model=='beta' || model=='binomial' ) {
        model.specs <- model.specs[-match('tau0',names(model.specs))]
        dat[dat==0] <- epsilon
        dat[dat==1] <- 1.0-epsilon
    }
    #if ( model=='binomial' && max(dat)<=1 ) { # using the binomial model 
    #    dat <- round(dat*length(dat))
    #}
    if ( model=='logit' ) {
        logit.fun <- function(P) log(P/(1-P))
        model.data <- monotonicity.logit.data # data transformations (see JAGS manual, p.32)
        model.inits <- list(logit.mu=logit.fun(tapply(dat,groups,mean)))
        dat[dat==0] <- epsilon
        dat[dat==1] <- 1.0-epsilon
    }
    ## RUN JAGS
    ##
    OUT <- run.jags(model=switch(model,
                        normal=monotonicity.normal,
                        binomial=monotonicity.binomial,
                        beta=monotonicity.beta,
                        logit=monotonicity.logit,
                        stop('unrecognized model:', model )),
                    model.data=model.data,
                    data=c(list(Nsamples=length(dat),Ngroups=ngroups,y=dat,groups=groups,query=query),
                           model.specs),
                    inits=model.inits,
                    monitor=c('mu','count'),
                    nadapt=nadapt,
                    nburnin=nburnin,
                    niter=niter,
                    no.del=no.del)
    VERBOSE(verbose,"done, variables monitored:\n>>", paste(names(OUT),collapse='\n>> '), '\n')
    OUT    
}

## JAGS models
##
## NORMAL model
##
monotonicity.normal <- function()
{    
    for ( i in 1:Nsamples )
    {
        y[i] ~ dnorm( mu[groups[i]], tau )
    }
    for ( j in 1:Ngroups )
    {
        mu[j] ~ dnorm( 0, tau0 )
    }
    tau ~ dgamma( alpha0, beta0 )

    diff <- mu[2:Ngroups] - mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) ) # P(mu_1 < mu_2 < ... < mu_n)
}
## LOG-NORMAL model
##
monotonicity.lognormal <- function()
{    
    for ( i in 1:Nsamples )
    {
        y[i] ~ dlnorm( mu[groups[i]], tau )
    }
    for ( j in 1:Ngroups )
    {
        mu[j] ~ dnorm( 0, tau0 )
    }
    tau ~ dgamma( alpha0, beta0 )

    diff <- mu[2:Ngroups] - mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) ) # P(mu_1 < mu_2 < ... < mu_n)
}
## BETA model (for proportions)
##
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
    count <- step( sum(step(diff*query)) - (Ngroups-1) ) # P(mu_1 < mu_2 < ... < mu_n)
}
monotonicity.binomial <- function() # this model works by multiplying the proportions by Nsamples
{
    for ( i in 1:Nsamples )
    {
        y[i] ~ dbin(p[i],Nsamples)
        logit(p[i]) <- logit.mu[groups[i]]
    }
    for ( j in 1:Ngroups )
    {
        logit.mu[j] ~ dnorm(0,tau)
    }
    tau ~ dgamma( alpha0, beta0 )

    diff <- logit.mu[2:Ngroups] - logit.mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) ) # P(mu_1 < mu_2 < ... < mu_n)

    for ( j in 1:Ngroups ) {
        mu[j] <- exp(logit.mu[j])/(exp(logit.mu[j])+1)
    }

}
## LOGIT model (for proportions)
##
monotonicity.logit.data <- function()
{
    for ( i in 1:Nsamples ) {
        logit.y[i] <- logit( y[i] )
    }
}
monotonicity.logit <- function()
{
    for ( i in 1:Nsamples )
    {
        logit.y[i] ~ dnorm( logit.mu[groups[i]], tau )
    }
    for ( j in 1:Ngroups )
    {
        logit.mu[j] ~ dnorm(0,tau0)
    }
    tau ~ dgamma( alpha0, beta0 )
    
    diff <- logit.mu[2:Ngroups] - logit.mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) ) # P(mu_1 < mu_2 < ... < mu_n)
    
    for ( j in 1:Ngroups ) {
        mu[j] <- exp(logit.mu[j])/(exp(logit.mu[j])+1)
    }
}
monotonicity.logit.old <- function() # THIS DOESN'T WORK
{
    ## adapted from
    ## http://stats.stackexchange.com/questions/41536/how-can-i-model-a-proportion-with-bugs-jags-stan
    ##
    for ( i in 1:Nsamples )
    {
        y[i] ~ dbeta(alpha[groups[i]],beta[groups[i]])
    }
    for ( j in 1:Ngroups )
    {
        alpha[j] <- mu[j] * phi
        beta[j] <- (1-mu[j]) * phi
        logit(mu[j]) <- a
    }
    phi ~ dgamma(alpha0,beta0)
    a ~ dgamma(0,tau0)
    P <- exp(a)/(exp(a)+1)

    diff <- mu[2:Ngroups] - mu[1:(Ngroups-1)]
    count <- step( sum(step(diff*query)) - (Ngroups-1) )
}

