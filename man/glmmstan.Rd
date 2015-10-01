\name{glmmstan}
\alias{glmmstan}
\title{Generalized liner mixed model in Rstan using glmer-stype formulas}
\description{
  Using standard formula notation from \code{glmer} (\code{lme4}), defines a Stan model (\code{rstan}) and optionally samples from the posterior. Can optionally compute WAIC. Supports model families: "gaussian", "binomial", "poisson", "negative binomial", "beta", "gamma", "lognormal" and "ordered". 
}
\usage{
glmmstan(formula_str, data, family="gaussian",
    center=FALSE, slice = NULL, offset=NULL,
    codeonly=FALSE, dataonly=FALSE, modelonly=FALSE,
    cauchy=2.5, lkj_corr=2,
    stancode=NULL, standata=NULL, stanmodel=NULL, stanfile=NULL, stanfit=NULL,
    parallel=FALSE, cores=NULL,
    iter=2000, warmup = NULL, chains= 2, thin=1)
}
\arguments{
  \item{formula_str}{Model formula or list of formulas, using \code{\link{glmer}} notation for varying effects.}
  \item{data}{Data frame or list.}
  \item{family}{Model family name or list of names for outcome(s). Valid choices are: "gaussian", "binomial" (logit link), "poisson" (log link), "nbinomial" (log link),"gamma" (log link), "beta" (logitlinks), "lognormal" (log link),"ordered" (cumulative logit).}
  \item{center}{If TRUE, center the independent variables(set means=0).}
  \item{offset}{Input offset term in model(only "poisson","nbinomial").}
  \item{codeonly}{If TRUE, output only stan code.}
  \item{dataonly}{If TRUE, output only dataset for stan.}
  \item{modelonly}{If TRUE, output only compiled stanmodel.}
  \item{cauchy}{Scale parameters of cauchy distribution which is prior distribution of variance components.}
  \item{parallel}{If TRUE, run the stanmodel in parallel using doparallel package.}
  \item{cores}{Number of cores to use when executing the chains in parallel.}
  \item{iter}{rstan parameter: total number of samples, including warmup.}
  \item{warmup}{rstan parameter: number of adaptation samples.}
  \item{chains}{rstan parameter: number of chains.}
  \item{thin}{rstan parameter: A positive integer specifying the period for saving sample.}
}
\author{Hiroshi Shimizu}
\seealso{\code{\link{glmer}},\code{\link{stan}}}
\examples{

data(baseball)

#"gaussian"
fit1 <- glmmstan(salary_log~1,data=baseball,family="gaussian")
output_result(fit1) #output glmm result
output_stan(fit1) #output summarized stan result (including rhat index)
print(fit1) #output stan result (same print() in rstan)

#"lognormal" with random effect
fit2 <- glmmstan(salary~HR+(1+HR|team),data=baseball,family="lognormal")
output_result(fit2)$WAIC #output only WAIC

#"negative binomial" with offset term
fit3 <- glmmstan(HR~1,data=baseball,family="nbinomial",offset="ATbats")
output_result(fit3)$beta #output only coefficients and scale parameters

#"ordered" with centering indipendent variables
fit4 <- glmmstan(Cluster~salary,data=baseball,family="ordered",center=TRUE)
output_result(fit4)
output_code(fit4) #confirm the stan code

#output only stan code, datase, and stan model
code1 <- glmmstan(HR~1+(1|player),data=baseball,family="poisson",codeonly=TRUE)
dataset1 <- glmmstan(HR~1+(1|player),data=baseball,family="poisson",dataonly=TRUE)
model1 <- glmmstan(HR~1+(1|player),data=baseball,family="poisson",modelonly=TRUE)

fit5 <- stan(model_code=code1, data=dataset1)
fit6 <- sampling(model1,data=dataset1)
}

\keyword{ }