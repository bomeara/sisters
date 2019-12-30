# This is from Käfer and Mousset, 2014, https://doi.org/10.1093/sysbio/syu024.
# It is copied directly from their CC0 licensed code on Dryad:
# Käfer, Jos; Mousset, Sylvain (2014), Data from: Standard sister clade comparison fails when testing derived character states, v2, Dryad, Dataset, https://doi.org/10.5061/dryad.jd8vg

# I've changed their core function (scc.test) to work on a simper data.frame rather than whatever list structure was being used originally (the example data files are no longer available online).


#################################################
###### Equations (Käfer and Mousset, 2014) ######
#################################################

# Conditional probability of observing K=k, given M=m

PKm <- Vectorize( function(k,m){
  return( 2*(m-k-1)/(m-1)/(m-2) )
})

# Conditional probability of having S>=s given M=m and K=k

PSmk <- Vectorize(function(s, m, k) {
  return( exp( lchoose(m-s, k-1) - lchoose(m-3, k-1) ) )
})

# Conditional Expectation of L' (length of the longer root edge) given M=m and K=k

ELmk <- Vectorize(function(m,k) {
  s <- 3:(m-k+1)
  psmk <- PSmk(s, m, k)
  return(0.5 + sum(psmk/s))
})

# Conditional Expectation of L (length of the root edge with p desendants) given M=m

ELmp <- Vectorize(function(m,p) {
  a <- (m-p-1)/(m-2)
  b <- (p-1)/2/(m-2)
  if (a>0) {
    s <- 3:(m-p+1)
    psmk <- PSmk(s, m, p)
    a <- a*(0.5 + sum(psmk/s))
  }
  return(a+b)
})

# Conditional probability of D (size of the derived clade)
# given M=m and K in {d,m-d}

PDmp <- Vectorize( function(d,m) {
  if (d != m-d) {
    a <- ELmp(m,d)
    b <- ELmp(m,m-d)
    return( a/(a+b) )
  } else {
    return(1)
  }
})

# Conditional distribution of D (size of the derived clade) given M=m
# thus taking into account the sampling bias (based on the length of root edges)

PDm <- Vectorize(function(d, m) {
  kp <- 1:(m-2)
  p <- PKm(k=kp, m=m)
  l <- c(ELmk(k=kp, m=m)*p,0) + 0.5*c(0,rev(p))
  l <- l/sum(l)
  return(l[d])
})

# Cumulative conditional probability of D (size of the derived clade)
# given M=m (default=lower tail)

PDmc <- Vectorize(function(d,m, lower.tail=TRUE) {
  kp <- 1:(m-2)
  p <- PKm(k=kp, m=m)
  l <- c(ELmk(k=kp, m=m)*p,0) + 0.5*c(0,rev(p))
  l <- l/sum(l)
  if (lower.tail) {
    return(sum(l[1:d]))
  } else {
    return(sum(l[d:(m-1)]))
  }
})


##########################################################################
## Sister Clade Comparison Test suggested in Slowinski and Guyer, 1993. ##
##########################################################################

## Original Slowinski and Guyer Test
## Using correct=TRUE uses correction for the intrinsic selection bias
## of sister clade analyses

SG.test <- function(dataset, lower.tail=FALSE, correct=FALSE) {
  pvalorg  <- (dataset$m-dataset$d)/(dataset$m-1)
  x2org <- -2*sum(log(pvalorg))
  porg  <- pchisq(x2org, lower.tail=FALSE, df=2*length(pvalorg))
  if (correct) {
    pvalcorr <- PDmc(d=dataset$d, m=dataset$m, lower.tail=lower.tail)
    x2corr <- -2*sum(log(pvalcorr))
    pcorr <- pchisq(x2corr, lower.tail=FALSE, df=2*length(pvalcorr))
    result <- matrix(data=c(x2org, x2corr, porg, pcorr), nrow=2, byrow=FALSE)
    colnames(result) <- c("X2","p.value")
    rownames(result) <- c("Cladogenetic (Slowinski and Guyer)","Anagenetic (Käfer and Mousset)")
  } else {
    result <- matrix(data=c(x2org, porg), ncol=2, byrow=FALSE)
    colnames(result) <- c("X2","p.value")
    rownames(result) <- c("Cladogenetic (Slowinski and Guyer)")
  }
  return(result)
}

########################################################################
## Sister Clade Comparison Test suggested in Käfer and Mousset, 2014. ##
########################################################################

## Käfer and Mousset Test for Sister Clade Comparison
## iter  : number of resamplings
scc.test <- function(dataset, iter=10000, alternative="two.sided") {

  # Select pairs of clades with m > 2
  dataset <- dataset[dataset$m>2,]

  # resampling function to generate n simulated d/m values
  resample <- function(y, n=1) {
    f <- function() {
           p <- runif(dim(y)[1]) < y[,2]
           return(mean(c(y[p,1],1-y[!p,1])))
         }
    return(replicate(n=n, expr=f()))
  }

  # Generate the matrix for the resampling function
  dpmatrix <- matrix(data=c(dataset$d/dataset$m,
                            PDmp(d=dataset$d, m=dataset$m)),
                     ncol=2, byrow=FALSE)

  # Observed and simulated mean d/m values
  dmobs <- mean(dpmatrix[,1])
  dsim <- NULL
  nqueue <- NULL

  # Iterate the resampling method
  for (it in c(iter[1],diff(iter))) {
    # Sample simulated d/m
    dsim <- c(dsim, resample(dpmatrix, n=it))
    ntotal <- length(dsim)
    # Compute the length of the queue of the distribution
    if (alternative == "two.sided") {
      nupqueue <- sum(dsim >= dmobs)
      nlowqueue <- sum(dsim <= dmobs)
      nqueue <- min(nupqueue, nlowqueue)
    } else if (alternative == "greater") {
      nqueue <- sum(dsim >= dmobs)
    } else if (alternative == "less") {
      nqueue <- sum(dsim <= dmobs)
    }
  }
  # compute the p-value of the test
  if (alternative == "two.sided") {
    p.value <- 2*nqueue/length(dsim)
  } else {
    p.value <- nqueue/length(dsim)
  }
  # return the result of the test
  scc <- list(statistic=dmobs, p.value=p.value, alternative=alternative,  iterations=iter )
  return(scc)
}
