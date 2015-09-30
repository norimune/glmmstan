\name{glmmstan}
\alias{glmmstan}
\docType{package}
\title{
Generalized liner mixed model in Rstan
}
\description{
Using standard formula notation from \code{glmer} (\code{lme4}), defines a Stan model (\code{rstan}) and optionally samples from the posterior. Can optionally compute WAIC. Supports model families: "gaussian", "binomial", "poisson", "negative binomial", "beta", "gamma", "lognormal" and "ordered". 
}
\details{
\tabular{ll}{
Package: \tab glmmstan\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-09-30\cr
License: \tab GPL\cr}}
\author{Hiroshi Shimizu}
\references{}
\keyword{ package }
\examples{
data(baseball)
fit <- glmmstan(salary~HR+(1+HR|team),data=baseball,family="lognormal")
output_result(fit)
output_stan(fit)
}
