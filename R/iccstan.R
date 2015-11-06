iccstan <- function(y,group,model=NULL,iter=2000,warmup=NULL,chains=2,thin=1,cores=NULL,parallel=TRUE){
  library("doParallel")
  
  data <- iccdata(y,group)
  
  if(is.null(warmup)==TRUE) warmup = floor(iter/2)
  
  if(is.null(cores)==TRUE){
    cores <- chains
  }
  if(cores>getOption("mc.cores",detectCores())){
    cores <- getOption("mc.cores",detectCores())
  }
  if(is.null(model)==TRUE){
    model <- iccmodel()
  }
  if(parallel==TRUE && chains > 1){
    cat("\nPreparing parallel process.\n")
    if(.Platform$OS.type != "windows") {
      cat("\nMCMC sampling start.\n")
      fit.list <- mclapply(1:chains, mc.cores=cores,function(i)
        rstan::sampling(model, data=data,chain_id = i,
                        iter=iter, warmup = warmup,chains=1,thin=thin)
      )
    }else{
      stopImplicitCluster()
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      parallel::clusterExport(cl=cl,c("data","model"),envir=environment())
      cat("\nMCMC sampling start.\n")
      fit.list <- parallel::parLapply(cl,1:chains, fun=function(i){
        require(rstan)
        rstan::sampling(model, data=data,chain_id = i,
                        iter=iter, warmup = warmup,chains=1,thin=thin)
      })
    }
    fit <- rstan::sflist2stanfit(fit.list)
  }else{
    fit <- rstan::sampling(model, data=data, iter=iter, warmup = warmup,chains=chains,thin=thin,cores=1)
  }
  cat("\nDone.\n")
  return(fit)
}

iccdata <- function(y,group){
  G <- length(unique(group))
  if(is.list(y)==FALSE){
    y <- matrix(y,length(y),1)
  }
  data <- list(N=nrow(y),G=G,P=ncol(y),y=y,group=group)
  return(data)
}

iccmodel <- function(){
  code <- '
    data{
    int N;
    int G;
    int P;
    int group[N];
    vector[P] y[N];
  }
  parameters{
    vector[P] mu;
    vector<lower=0>[P] sdb;
    vector<lower=0>[P] sdw;
    vector[P] r[G];
  }
  model{
    mu ~ normal(0,100);
    sdb ~ cauchy(0,2.5);
    sdw ~ cauchy(0,2.5);
    for(g in 1:G) r[g] ~ normal(mu,sdb);
    for(n in 1:N) y[n] ~ normal(r[group[n]],sdw);
  }
  generated quantities{
    vector<lower=0,upper=1>[P] icc;
    vector<lower=0>[P] varb;
    vector<lower=0>[P] varw;
    for(p in 1:P){
      varb[p] <- sdb[p]^2;
      varw[p] <- sdw[p]^2;
      icc[p] <- varb[p]/(varb[p]+varw[p]);
    }
  }
  '
  model <- stan_model(model_code=code)
  return(model)
}

output_icc<- function(fit,digits = 3){
  return(print(fit,pars=c("varb","varw","icc"),probs=c(0.025,0.5,0.975),digits=digits))
}