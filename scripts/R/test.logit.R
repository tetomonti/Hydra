data {
    for ( i in 1:Nsamples ) {
        logit.y[i] <- logit(y[i])
    }
}
model {
    for ( i in 1:Nsamples ) {
        logit.y[i] ~ dnorm(mu,tau)
    }
    mu ~ dnorm(0,.001)
    tau ~ dgamma(.1,.1)

    P <- exp(mu)/(exp(mu)+1)
}
