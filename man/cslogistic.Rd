\name{cslogistic}
\alias{cslogistic}

\title{Perform an Analysis of a conditionally specified logistic regression model}

\description{
  This package contains functions for likelihood and posterior
  analysis of conditionally specified logistic regression models.
}
  
  
\seealso{
\code{\link{MleCslogistic}}, \code{\link{BayesCslogistic}}.
}


\details{
Assume that for each of \eqn{n} experimental units the values of \eqn{m} binary variables 
\deqn{Y_{i1}, \ldots, Y_{im}}{Yi1, \ldots, Yim} 
are recorded. The 'MleCslogistic' and 'BayesCslogistic' functions fit a 
conditional specified logistic regression model, such that for \eqn{i = 1, \ldots, n} 
and \eqn{j = 1, \ldots, m},

\deqn{
\mbox{logit} P(Y_{ij}=1 | Y_{ik}=y_{k}, k \neq j) = X_{ij}\beta_j + \sum_{k=1, k \neq j}^m \alpha_{jk} y_k }{%
logit P(Yij=1 | Yik=yk, k neq j) = Xij \beta j + \sum_{k=1, k \neq j} \alphajk yk }

where, the parameters \eqn{\alpha_{jk}}{\alphajk} have interpretation as conditional log-odds ratios and 
the parameters \eqn{\beta_j}{\beta j} correspond to the regression coefficients associated to the vector of 
covariates \eqn{X_{ij}}{Xij}. For compatibility of conditional distributions it is assumed that 
\eqn{\alpha_{jk} = \alpha_{kj}}{\alphajk = \alphakj}, \eqn{j \neq k }.

}

\references{

Garcia-Zattera, M. J., Jara, A., Lesaffre, E. and Declerck, D. (2005). On conditional 
    independence for multivariate binary data in caries research. In preparation.

Joe, H. and Liu, Y. (1996). A model for multivariate response with covariates based on 
    compatible conditionally specified logistic regressions. Satistics & Probability Letters 
    31: 113-120.   
}

\author{
Alejandro Jara Vallejos \email{Alejandro.JaraVallejos@med.kuleuven.be}

Maria Jose Garcia-Zattera \email{MariaJose.GarciaZattera@med.kuleuven.be}
}

\keyword{regression}
