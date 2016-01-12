# glmmstanとは
この関数は `glmer` (`{lme4}`)と同じような関数の書き方をするだけで，Stanのモデル (`{rstan}`) を作り，事後分布を生成するものです。豊富なオプションの指定ができます。また，WAICを出力するオプションもあります。対応している分布族は次の通りです。

+ 正規分布 "gaussian"
+ 二項分布 "binomial"
+ ポアソン分布 "poisson"
+ 負の二項分布 "negative binomial"
+ ベータ分布 "beta"
+ ガンマ分布 "gamma"
+ 対数正規分布 "lognormal"
+ ベータ二項分布 "beta-binomial"
+ 順序変数の分布 "ordered"
+ ゼロ過剰ポアソン分布 "zero-inflated poisson"
+ ゼロ過剰負の二項分布　"zero-inflated negative binomial"
+ ゼロ過剰ガンま分布 "zero-inflated and gamma". 

## インストールの方法

もしまだ `RStan`をインストールしていないのであれば，[RStan document](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)を読んで先にインストールを済ませておいてください。

`glmmstan` をGitHubからインストールするには，`devtools`パッケージも必要です。次のようにして準備してください。

```
install.packages("devtools")
```

これらのパッケージのインストールが終われば，次のコードをRで実行してください。

```
library(devtools)
install_github("norimune/glmmstan")
```

## 使用例
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