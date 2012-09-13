\name{BayesCslogistic}
\alias{BayesCslogistic}

\title{Perform a Bayesian Analysis of a conditionally specified logistic regression model}
\description{
    This function generates a posterior density sample
    from a conditionally specified logistic regression model 
    for multivariate binary data using a random walk Metropolis
    algorithm. The user supplies data and priors, 
    and a sample from the posterior density is returned as a
    object, which can be subsequently analyzed with functions 
    provided in the coda package.
}
  
  
\usage{

BayesCslogistic(formula, type = TRUE, intercept = TRUE, 
           burnin = 1000, mcmc = 10000, thin=1, 
           tune=1.1, beta.start = NA, b0 = 0, B0 = 0, ...)   
}

\arguments{
    \item{formula}{Model formula.}

   \item{type}{logical variable indicating if covariates have the same effect 'TRUE' or different 
               effect 'FALSE' for each variable.}
               
   \item{intercept}{logical variable indicating if only the intercept 'TRUE' or all the covariates 
                    have different effect 'FALSE' for each variable. The option 'type' must be 'FALSE'.}   

    \item{burnin}{The number of burn-in iterations for the sampler.}

    \item{mcmc}{The number of Metropolis iterations for the sampler.}

    \item{thin}{The thinning interval used in the simulation.  The number of
    mcmc iterations must be divisible by this value.}

    \item{tune}{Metropolis tuning parameter. Make sure that the
      acceptance rate is satisfactory (typically between 0.20 and 0.5)
      before using the posterior density sample for inference.}
    
    \item{beta.start}{The starting value for the \eqn{\beta}{beta} vector.
    This can either be a scalar or a column vector with dimension equal to the number of 
    betas. If this takes a scalar value, then that value will serve as the 
    starting value for all of the betas.  The default value of NA will
    use the maximum likelihood estimate of \eqn{\beta}{beta} as the starting 
    value. Those are obtained using the function Cslogistic}

    \item{b0}{The prior mean of \eqn{\beta}{beta}.  This can either be a 
    scalar or a column      
    vector with dimension equal to the number of betas. If this takes a scalar
    value, then that value will serve as the prior mean for all of the
    betas.}
    
    \item{B0}{The prior precision of \eqn{\beta}{beta}.  This can either be a
    scalar
    or a square matrix with dimensions equal to the number of betas.  If this
    takes a scalar value, then that value times an identity matrix serves
    as the prior precision of \eqn{\beta}{beta}. Default value of 0 is
    equivalent to an improper uniform prior for beta.}
     
    \item{\ldots}{further arguments to be passed.}       
}

\value{
   An mcmc object that contains the posterior density sample.  This 
   object can be summarized by functions provided by the coda package.
}


\seealso{
\code{\link{cslogistic}}, \code{\link{MleCslogistic}}.
}

\references{

Garcia-Zattera, M. J., Jara, A., Lesaffre, E. and Declerck, D. (2007). Conditional 
    independence of multivariate binary data with an application in caries research. 
    Computational Statistics and Data Analysis, 51(6): 3223-3232.

Joe, H. and Liu, Y. (1996). A model for multivariate response with covariates based on 
    compatible conditionally specified logistic regressions. Satistics & Probability Letters 
    31: 113-120.   
}

\examples{
# simulated data set

  library(mvtnorm)

  n <- 400
  mu1 <- c(-1.5,-0.5)
  Sigma1 <- matrix(c(1, -0.175,-0.175,1),ncol=2)
  agev <- as.vector(sample(seq(5,6,0.1),n,replace=TRUE))
  beta1 <- 0.2

  z <- rmvnorm(n,mu1,Sigma1)
  zz <- cbind(z[,1]+beta1*agev,z[,2]+beta1*agev)
  dat <- cbind(zz[,1]>0,zz[,2]>0,agev)
  colnames(dat) <- c("y1","y2","age")
  data0 <- data.frame(dat)
  attach(data0)

# equal effect of age for all the covariates

  y <- cbind(y1,y2)

  fit0 <- BayesCslogistic(y~age)
  fit0
  summary(fit0)
  plot(fit0)

# different effects: only intercept

  fit1 <- BayesCslogistic(y~age,type=FALSE)
  fit1
  summary(fit1)
  plot(fit1)

# different effects: all the covariates 

  fit2 <- BayesCslogistic(y~age,type=FALSE,intercept=FALSE)
  fit2
  summary(fit2)
  plot(fit2)

}

\author{
Alejandro Jara \email{atjara@uc.cl}

Maria Jose Garcia-Zattera \email{mjgarcia@uc.cl}
}

\keyword{regression}
