if ( (CBMGIT <- Sys.getenv('CBMGIT'))=="" )
    stop( "Use 'setenv CBMGIT ..' to set CBMgithub's base directory" )

source(paste(CBMGIT,"scripts/R/CBMRtools/R/misc.R",sep='/'))
source(paste(CBMGIT,"scripts/R/jags.R",sep='/'))
source(paste(CBMGIT,"scripts/R/monotonicity.test.R",sep='/'))

OSCC <- load.var(file='~/Research/Projects/oralcancer/TCGA_HNSC/processed_data/oscc_yap_activity.Rdata')

setwd("~/Research/Projects/oralcancer/taz_yap_dpagt1")

valid.stage <- c("Stage 0","Stage I","Stage II","Stage III","Stage IVA/IVB")
valid.grade <- c("G0","G1","G2","G3/G4")

OSCCg <- OSCC[OSCC$grade %in% valid.grade,]
OSCCg <- OSCCg[order(OSCCg$grade),]

OSCCs <- OSCC[OSCC$stage %in% valid.stage,]
OSCCs <- OSCCs[order(OSCCs$stage),]

grade.col <- col.gradient(c('white','dark green'),length(valid.grade));names(grade.col) <- valid.grade
stage.col <- col.gradient(c('white','dark green'),length(valid.stage));names(stage.col) <- valid.stage

png('results/ASSIGNscore.YAPactivity.bygrade.png')
boxplot(OSCCg$yap_activity_score~OSCCg$grade,col=grade.col,main='YAP activity score (by grade)')
dev.off()

png('results/ASSIGNscore.YAPactivity.bystage.png')
boxplot(OSCCs$yap_activity_score~OSCCs$stage,col=stage.col,main='YAP activity score (by stage)')
dev.off()

require(rjags)
logit <- function(P) log(P/(1-P))

groups <- genIndex(factor(OSCCg$grade,levels=valid.grade),base=1,do.sort=TRUE)
dat <- OSCCg$yap_activity_score

jags.input <- list(y=dat,Nsamples=length(dat),groups=groups,Ngroups=nlevels(groups),query=c(1,1,1),
                   tau0=.01,alpha0=.1,beta0=.1)

Nsize <- 20
jags.input <- list(y=round(dat*Nsize),Nsamples=length(dat),Nsize=Nsize,groups=groups,Ngroups=nlevels(groups),
                   query=c(1,1,1),alpha0=.1,beta0=.1)
jags.inits <- list(mu=logit(tapply(dat,groups,mean)))

OUT <- run.jags(model=monotonicity.binomial,data=jags.input,inits=jags.inits,monitor=c('count','diff','P','mu','tau'))

OUT <- monotonicity.test(dat=dat,
                         groups=groups,
                         query=c(1,1,1),
                         model='binomial',
                         nadapt=1000,
                         nburnin=10000,
                         niter=10000,
                         model.specs=list(tau0=.1,alpha0=.1,beta0=.1),
                         verbose=TRUE,no.del=TRUE)

SS <- jags.summary.stats(OUT$count)

