\name{iccstan}
\alias{glmmstan}
\title{intra-class correlations in rstan}
\description{
  calculating intra-class correlations(ICC) in rstan. 
}
\usage{
iccstan(y, group, model=NULL, iter=2000, warmup = NULL, chains= 2, thin=1, cores=NULL, parallel=TRUE)
}
\arguments{
  \item{y}{list.}
  \item{group}{vector.}
  \item{model}{compiled stanmodel.}
  \item{iter}{rstan parameter: total number of samples, including warmup.}
  \item{warmup}{rstan parameter: number of adaptation samples.}
  \item{chains}{rstan parameter: number of chains.}
  \item{thin}{rstan parameter: A positive integer specifying the period for saving sample.}
  \item{parallel}{If TRUE, run the stanmodel in parallel using doparallel package.}
  \item{cores}{Number of cores to use when executing the chains in parallel.}
}
\author{Hiroshi Shimizu}
\seealso{\code{\link{glmmstan}},\code{\link{stan}}}
\examples{

data(baseball)

#intra-class correlations
model <- iccmodel()
y <- subset(baseball,select=c(HR,HIT))
fit <- iccstan(y,baseball$team,model=model)
output_icc(fit)
}

\keyword{ }