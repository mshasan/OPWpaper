set.seed(123)
nrep=200
m=100
alpha=0.05

pvals<-matrix(runif(nrep*m),nrow=nrep)


sum(apply(pvals<alpha/m,1,sum))/nrep

sum(pvals<alpha/m)/nrep


rej <- c()
set.seed(123)
for(s in 1:nrep)
{
  pvals = runif(m)
  rej[s] <- sum(pvals <= alpha/m)
}

sum(rej)/nrep
