#####################################################
## FUNCTION: RUN JAGS
#####################################################
##
## a simple wrapper function (so as not to have to remember all the jags-related commands)

run.jags <- function
(
    model,          # a jags model specified as an R function
    model.data=NULL,# data transformations, if any (see JAGS manual, p. 32)
    data,           # list of variables to be passed to jags
    inits=NULL,     # initialization values (optional)
    monitor,        # variables to monitor
    nchains=1,      # number of MCMC chains
    nadapt=100,
    nburnin=1000,
    niter=1000,
    use.coda=FALSE,
    no.del=FALSE    # don't delete the model file
)
{
  ## first, write the jags model to file
  modelFile <- genName(stub=names(model))
  write.jags.model( FUN=model, DAT=model.data, file=modelFile )

  ## compile the model (and check syntax)
  modelOut <- jags.model(modelFile,data=data,inits=inits,n.chains=nchains,n.adapt=nadapt)

  ## initial burn-in
  cat("Burn-in\n" )
  update(modelOut,nburnin)

  ## actual MCMC sampling
  cat("Updating models\n" )
  OUT <- {
    if (use.coda)
      coda.samples(model=modelOut,variable.names=monitor,n.iter=niter)
    else
      jags.samples(model=modelOut,variable.names=monitor,n.iter=niter)
  }
  if ( !no.del )
      system(paste('rm -f',modelFile))
  
  OUT
}
#####################################################
## FUNCTION: WRITE JAGS MODEL
#####################################################
##
## the write.jags.model function is just to avoid having to depend on yet another file (and ..
## ..writing the file on-the-fly instead). I'm looking for ways of eliminating the writing ..
## ..to file completely and instead 'piping' the function to 'jags-model(...)' directly

write.jags.model <- function( FUN, DAT, file, overwrite=TRUE )
{
  if ( !overwrite && file.access(file)!=0 )
    stop( "can't overwrite '",file,"'" )

  if ( !is.null(DAT) ) {
      #writeLines( c('data',capture.output(print(body(DAT)))), con=file )
      cat( c('data',capture.output(print(body(DAT)))), sep='\n', file=file )
  }
  #writeLines( c('model',capture.output(print(body(FUN)))), con=file )
  cat( c('model',capture.output(print(body(FUN)))), sep='\n', file=file, append=!is.null(DAT) )
}
#####################################################
## FUNCTION: JAGS SUMMARY STATS
#####################################################

jags.summary.stats <- function
(
    OBJ,
    names=NULL,
    outfile=NULL,
    header='ID',
    probs=c(.0005,.005,.025,.25,.5,.75,.975,.995,.9995),
    ord=NULL,
    decreasing=TRUE
)
{
  sumStats <-
      data.frame( t(apply(OBJ,1,quantile,probs=probs)), mean=apply(OBJ,1,mean), check.names=FALSE )
  if ( !is.null(ord) ) {
    if ( !(ord %in% colnames(sumStats)) ) {
      stop( 'ord must be one of [', paste(colnames(sumStats),collapse=', '), ']')
    }
    sumStats <- sumStats[order(sumStats[,ord],decreasing=decreasing),,drop=FALSE]
  }
  if (!is.null(names)) rownames(sumStats) <- names
  if ( !is.null(outfile) )
    my.write.table(sumStats,file=outfile,names=header)
  sumStats
}

#####################################################
## HERE'S AN EXAMPLE OF USE 
#####################################################
##
## I like to define the jags model as R functions (so as to take
## advantage of the R editing conventions). However, notice that the
## function is not executable in R. It needs to be written to file (by
## using the 'write.jags.model' function), and then used by the rjags
## functions 'jags.model', 'update', and 'jags.samples' (all wrapped
## in my function 'run.jags')

if (FALSE)
{
## LOGIT BINOMIAL MODEL
##
## hierarchical binomial proportion modeling (comparing two proportions)
## (not sure it's 'kosher')

logitBin <- function()
{
  for ( i in 1:length(y) )
  {
    y[i] ~ dbin(p[i],n[i])
    logit(p[i]) <- z[i]
    z[i] ~ dnorm(mu[i],tau)
    mu[i] <- alpha + beta * x[i] # x[i] indicates group 0 or 1
  }
  alpha ~ dnorm(0,0.001) # very vague priors
  beta ~ dnorm(0,0.001)  # ..
  tau ~ dgamma(0.001,0.001)
  
  theta1 <- exp(alpha) / ( 1+exp(alpha) )       # this is just to monitor the estimate of ..
  theta2 <- exp(alpha+beta)/(1+exp(alpha+beta)) # ..the two proportions
}
## testing it on toy data (randomly generated from two binomial distributions)

N <- 100
set.seed(123)
N1 <- sample(28:32,size=N,replace=TRUE); Y1 <- rbinom(N,N1,0.6); R1 <- Y1/N1 # group 0
N2 <- sample(28:32,size=N,replace=TRUE); Y2 <- rbinom(N,N2,0.7); R2 <- Y2/N2 # group 1

## plot the distribution of the samples generated

hist(R2,xlim=c(min(R1,R2),max(R1,R2)),density=10,angle=135,col='blue')
hist(R1,add=TRUE,density=10,angle=45,col='red')
abline(v=c(median(R1),median(R2)),col=c('red','blue'),lwd=2)

## run the hierarchical binomial proportion model (see 'logitBin' above)

OUTlogit <- run.jags(model=logitBin,
                     data=list(y=c(Y1,Y2),n=c(N1,N2),x=rep(0:1,times=c(length(N1),length(N2)))),
                     monitor=c('alpha','beta','tau','theta1','theta2'))

## beta is the parameter that's tested to be different from 0

round(quantile(OUTlogit$beta,probs=c(0,.025,.05,.5,.95,.975,1)),3)
}
