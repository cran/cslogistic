\name{MleCslogistic}
\alias{MleCslogistic}


\title{Perform a Maximum Likelihood Analysis of a conditionally specified logistic regression model}
\description{
   Fit a conditional specified logistic regression model for multivariate binary responses. 
}
  
\usage{
MleCslogistic(formula,type = TRUE, intercept = TRUE, method = "BFGS",
           maxiter=1000 , data, ...)

}

\arguments{
   \item{formula}{a symbolic description of the model to be fit.}
   \item{type}{logical variable indicating if covariates have the same effect 'TRUE' or different 
               effect 'FALSE' for each variable.}
   \item{intercept}{logical variable indicating if only the intercept 'TRUE' or all the covariates 
                    have different effect 'FALSE' for each variable. The option 'type' must be 'FALSE'.}   
   \item{method}{the optimization method to be used; the default method is
          "BFGS".}
   \item{maxiter}{maximum number of iterations used by the optimization method.}
   \item{data}{an optional data frame containing the variables in the model.
           If not found in 'data', the variables are taken from
          'environment(formula)', typically the environment from which
          'cslogistic' is called..}
   \item{...}{further arguments to be passed.}      
}

\seealso{
\code{\link{cslogistic}}, \code{\link{BayesCslogistic}}.
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

  fit0 <- MleCslogistic(y~age)
  fit0
  summary(fit0)

# different effects: only intercept

  fit1 <- MleCslogistic(y~age,type=FALSE)
  fit1
  summary(fit1)

# different effects: all the covariates

  fit2 <- MleCslogistic(y~age,type=FALSE,intercept=FALSE)
  fit2
  summary(fit2)

}

\author{
Alejandro Jara \email{atjara@uc.cl}

Maria Jose Garcia-Zattera \email{mjgarcia@uc.cl}
}

\keyword{regression}


