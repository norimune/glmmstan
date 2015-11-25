# glmmstan
Using standard formula notation from \code{glmer} (\code{lme4}), defines a Stan model (\code{rstan}) and optionally samples from the posterior. Can optionally compute WAIC. Supports model families: "gaussian", "binomial", "poisson", "negative binomial", "beta", "gamma", "lognormal", "beta-binomial", "ordered",and "zero-inflated poisson, negative binomial and gamma". 

## Install
If you have not install `RStan`, please read [RStan document](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started) and install it.

And you need `devtools` to install `glmmstan` from GitHub.

```
install.packages("devtools")
```

After install these packages, excute the following in R:

```
library(devtools)
install_github("norimune/glmmstan")
```
## Example
```
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

#intra-class correlations
model <- iccmodel()
y <- subset(baseball,select=c(HR,HIT))
fit7 <- iccstan(y,baseball$team,model=model)
output_icc(fit7)

```