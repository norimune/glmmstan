library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

glmmstan <- function(formula_str,data,family="gaussian",center = FALSE,slice = NULL,offset=NULL,                     
                     codeonly=FALSE,dataonly=FALSE,modelonly=FALSE,cauchy = 2.5,lkj_corr = 2,
                     stancode=NULL,standata=NULL,stanmodel=NULL,stanfile=NULL,stanfit=NULL,
                     parallel=FALSE,cores=NULL,iter=2000,warmup = NULL,chains= 2,thin=1){
  
  require("rstan")
  library("doParallel")
  
  #formula...Model formula. Using "glmer" notation.
  #data...Data.frame or list.
  #family...Model family name for outcome.
  #         Valid choices:gaussian,bernoulli,binomial,poisson,nbinomial,gamma,lognormal,beta,ordered.
  #center...If TRUE, xvalue are centered from grand means.
  #slice...Slice variable name.Simple slope effects are estimated.
  #offset...Offset vaiable name.
  
  #######################
  ###private functions###
  #######################
  
  interaction_data <- function(data,xname){
    for(i in 1:length(xname)){
      temp <- regexpr(":",xname[i],fixed=TRUE)
      if(temp[1]>0){
        temp1 <- unlist(strsplit(xname[i],":",fixed=TRUE))
        data[xname[i]] <- 1
        for(j in 1:length(temp1)){
          if(center==T){
            data[xname[i]] <- data[xname[i]] * (data[temp1[j]]-mean(as.matrix(data[temp1[j]]),na.rm=TRUE))
          }else{
            data[xname[i]] <- data[xname[i]] * data[temp1[j]]
          }
        }
      }
    }
    return(data)
  }
  
  interaction2 <- function(data,idname){
    for(i in 1:length(idname)){
      temp <- regexpr(":",idname[i],fixed=TRUE)
      if(temp[1]>0){
        temp1 <- unlist(strsplit(idname[i],":",fixed=TRUE))
        for(j in 1:length(temp1)-1){
          tempdataa <- paste0(data[temp1[j]][,] ,":", data[temp1[j+1]][,])
          data[idname[i]] <- as.numeric(as.factor(tempdataa))
        }
      }
    }
    return(data)
  }
  
  centering <- function(data,xname){
    for(i in 1:length(xname)){
      temp <- regexpr(":",xname[i],fixed=TRUE)
      if(temp[1]>0){
        for(j in 1:length(data)){
          if(names(data[j])==xname[i] && names(data[j]) != "(intercept)"){
            data[j] <- data[j]-mean(as.matrix(data[j]),na.rm=T)
          }
        }
      }
    }
    return(data)
  }
    
  nobars <- function(term) {
    if (!('|' %in% all.names(term))) return(term)
    if (is.call(term) && term[[1]] == as.name('|')) return(NULL)
    if (length(term) == 2) {
      nb <- nobars(term[[2]])
      if (is.null(nb)) return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobars(term[[2]])
    nb3 <- nobars(term[[3]])
    if (is.null(nb2)) return(nb3)
    if (is.null(nb3)) return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }
  
  findbars <- function(term) {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
  }
  
  subbars <- function(term){
    if (is.name(term) || !is.language(term)) return(term)
    if (length(term) == 2) {
      term[[2]] <- subbars(term[[2]])
      return(term)
    }
    stopifnot(length(term) >= 3)
    if (is.call(term) && term[[1]] == as.name('|'))
      term[[1]] <- as.name('+')
    for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
    term
  }
  dvcalc <- function(data,yname){
    yname2 <- unlist(strsplit(yname,"+",fixed=TRUE))
    yname2 <- unlist(strsplit(yname2,"-",fixed=TRUE))
    yname2 <- unlist(strsplit(yname,"*",fixed=TRUE))
    yname2 <- unlist(strsplit(yname2,"/",fixed=TRUE))
    if(length(yname2)>1){
      lyn2 <- length(yname2)-1
      count <- 0
      for(i in 1:lyn2){
        temp <- substr(yname,count+nchar(yname2[i])+1,count+nchar(yname2[i])+1)
        if(temp=="+"){
          data[yname2[1]] <- data[yname2[1]] + data[yname2[i+1]]
        }else if(temp=="-"){
          data[yname2[1]] <- data[yname2[1]] - data[yname2[i+1]]
        }else if(temp=="*"){
          data[yname2[1]] <- data[yname2[1]] * data[yname2[i+1]]
        }else if(temp=="/"){
          data[yname2[1]] <- data[yname2[1]] / data[yname2[i+1]]
        }
        count <- count + 1 + nchar(yname2[i])
      }
    }
    return(data)
  }
  
  #################
  ###check input###
  #################
  
  if(family=="normal") family <- "gaussian"
  if(family=="Gamma") family <- "gamma"
  if(family=="negative binomial") family <- "nbinomial"
  
  if(family=="gaussian" || family=="binomial" || family=="poisson" || family=="gamma" || family=="beta"
       || family=="nbinomial" || family=="ordered" || family=="bernoulli" || family=="lognormal"){
    
  }else{
    stop(paste0("Input family type(",family,") is incorrect."))
  }
  
  if(is.null(warmup)==TRUE){
    warmup = floor(iter/2)
  }
  
  if(is.null(stancode)==TRUE){
    nocode <- TRUE
  }else{
    nocode <- FALSE
    codestan <- stancode
  }
  if(is.null(stanmodel)==TRUE){
    nomodel <- TRUE
  }else{
    nomodel <- FALSE
    modelname <- stanmodel@model_name
    temp <- unlist(strsplit(modelname,"[",fixed=TRUE))
    #family <- gsub("]", "", temp[2],fixed=TRUE)
  }
  if(is.null(standata)==TRUE){
    nostandata <- TRUE
  }else{
    nostandata <- FALSE
    datastan <- standata
  }
  if(is.null(stanfile)==TRUE){
    nofile <- TRUE
  }else{
    nofile <- FALSE
  }
  if(is.null(stanfit)==TRUE){
    nofit <- TRUE
  }else{
    nofit <- FALSE
    modelname <- stanfit@stanmodel@model_name
    stanmodel <- stanfit@stanmodel
    nomodel <- FALSE
    temp <- unlist(strsplit(modelname,"[",fixed=TRUE))
    #family <- gsub("]", "", temp[2],fixed=TRUE)
  }
  
  if(is.null(cores)==TRUE){
    cores <- chains
  }
  if(cores>getOption("mc.cores",detectCores())){
    cores <- getOption("mc.cores",detectCores())
  }
  
  #####################
  ###variable names####
  #####################
  
  formula <- as.formula(formula_str)
  formula_str <- paste0(formula[2],formula[1],formula[3])
  
  xformula <- nobars(formula)
  if ( class(xformula)=="name" & length(xformula)==1 ) {
    xformula <- nobars( as.formula( paste( deparse(formula) , "+ 1" ) ) )
  }
  xname <- colnames( model.matrix(xformula, data))
  if(family=="ordered"){
    xname <- xname[-1]
  }
  yname <- deparse(xformula[[2]] )
  y <- model.frame(xformula, data)[,1]
  totalname <- "null"
  bitotalcheck <- 0
  if (class(y)=="matrix" ) {
    bitotalcheck <- 1
    yname <- colnames(y)[1]
    totalname  <- "(bitotal)"
  }    
  zname <- list()
  if ( nchar(formula[3]) != nchar(nobars(formula)[3]) ) {
    zformula <- findbars(formula)
    for ( i in 1:length(zformula) ) {    
      idname <- all.vars( zformula[[i]][[3]])    
      v <- zformula[[i]][[2]]
      if ( class(v)=="numeric" ) {
        # just intercept
        zname[[idname]] <- "(Intercept)"
      } else {
        tempv <- gsub(" ", "", deparse( v ),fixed=TRUE)
        tempv <- unlist(strsplit(tempv,"+",fixed=TRUE))   
        if(tempv[1]=="1"){
          f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
          zname[[idname]] <- colnames( model.matrix( f , data ) )
        }else{
          f <- as.formula( paste( "~" , deparse( v ) , sep="" ) )
          f <- as.data.frame(model.matrix( f , data ))[-1]
          zname[[idname]] <- colnames( f )
        }
      }
    }
  }
  
  idname <- list()
  R <- length(zname)
  P <- length(xname)
  idname <- c()
  if(R>0) for(i in 1:R) idname[i] <- attr(zname,"name")[i]
  
  ###slice variable
  checkslice <- 0
  if(is.null(slice)==FALSE){    
    for(i in 1:length(xname)){
      if(xname[i]==slice){
        checkslice <- 1
      }
    }
    if(checkslice==0){
      stop("Slice variable is not found in the formula")
    }
    checkslice <- 0
    slicetemp <- list()
    for(i in 1:P){
      temp <- regexpr(":",xname[i],fixed=TRUE)
      if(temp[1]>0){
        interaction <- unlist(strsplit(xname[i],":",fixed=TRUE))
        for(j in 1:length(interaction)){
          if(interaction[j]==slice && length(interaction)==2){
            checkslice <- checkslice + 1
            slicetemp[[checkslice]] <- c(i,interaction[3-j])
          }
        }
      }
    }
    if(checkslice==0){
      stop("Slice variable is not found in the interaction terms")
    }
    intrctnum <- c()
    simplenum <- c()
    for(i in 1:checkslice){
      intrctnum[i] <- slicetemp[[i]][1]
      checksimple <- 0
      for(j in 1:P){
        if(slicetemp[[i]][2]==xname[j]){
          checksimple <- j
        }
      }
      if(checksimple==0){
        stop("Simple slople variable is not found in the main-effect terms")
      }else{
        simplenum[i] <- checksimple
      }
    }
  }
  
  ###offset variable
  checkoffset <- 0
  if(is.null(offset)==FALSE){
    for(i in 1:length(names(data))){
      if(names(data)[i]==offset){
        checkoffset <- 1
      }
    }
    if(checkoffset==0){
      stop("Offset variable is not found in the data")
    }
  }
  
  ######################
  ###creating dataset###
  ######################
  if(nostandata==TRUE){
  
    dat2 <- data
    y <- model.frame(xformula, dat2, na.action=NULL)[,1]
    if(class(y)=="matrix"){
      dat2[yname] <- y[,1]
      dat2["(bitotal)"] <- y[,2]
    }else{
      dat2[yname] <- y    
    }
    dat2["(Intercept)"] <- 1
    if(center==TRUE) dat2 <- centering(dat2,xname)
    if(R>0) dat2 <- interaction2(dat2,idname)
    dat2 <- interaction_data(dat2,xname)    
    if(R>0) for(i in 1:R)dat2 <- interaction_data(dat2,zname[[i]])
    datname <- c(yname,xname)
    if(R>0) for(i in 1:R) datname <- c(datname,idname[i],zname[[i]])
    if(family=="binomial" && totalname!="null") datname <- c(datname,totalname)
    if(checkoffset==1) datname <- c(datname,offset)
    
    dat3 <- na.omit(subset(dat2,select=unique(datname)))
    
    y <- as.numeric(dat3[,yname])
    if(family=="ordered"){
      y <- as.numeric(as.factor(dat3[,yname]))
    }
    K <- length(summary(as.factor(dat3[,yname])))
    x <- subset(dat3,select=xname)
    z <- list()
    id <- list()
    Q <- array()
    G <- array()
    if(R>0){
      for(i in 1:R){
        Q[i] <- length(zname[[i]])
        if(Q[i]==1){
          z[[i]] <- subset(dat3,select=zname[[i]])[,]
        }else {
          z[[i]] <- subset(dat3,select=zname[[i]])
        }        
        id[[i]] <- as.numeric(as.factor(as.numeric(as.factor(dat3[,idname[i]]))))       
        G[i] <- length(unique(id[[i]]))
      }
    }else{
      Q <- c(0,0,0,0)
      G <- c(0,0,0,0)
    }
    N <- nrow(dat3)
    if(family=="binomial" && totalname !="null"){
      bitotal <- dat3[yname][,] + dat3[,totalname]
    }else if(family=="binomial"){
      dat3$bitotal <- max(y)
      bitotal <- dat3$bitotal
    }
    if(checkoffset==1){
      offsetdata <- dat3[offset]
    }
    
    z2name <- c()
    id2name <- c()
    for(i in 1:R){
      z2name[i] <- paste("z",i,sep="")
      id2name[i] <- paste("id",i,sep="")
    }
    if(R>0){
      datastan <- list(N=N,P=P,R=R,Q=Q,G=G,y=y,x=x)
    }else{
      datastan <- list(N=N,P=P,y=y,x=x)
    }
    if(R>0){
      for(i in 1:R){
        datastan[z2name[i]] <- list(z[[i]])
        datastan[id2name[i]] <- list(id[[i]])
      }
    }
    if(family=="binomial") datastan$bitotal <- bitotal
    if(family=="ordered") datastan$K <- K
    if(checkoffset==1) datastan$offset <- offsetdata[,]  
    if(checkslice>0) slicesd <- sd(dat3[slice][,])
    if(family=="gaussian"||R==1) datastan$idn1 <- as.vector(table(datastan$id1))
  }
  
  if(dataonly==TRUE){
    return(datastan)
  }
  
  ##########################
  ####creating stan code####
  ##########################
  
  if(nocode==TRUE && nomodel==TRUE){  
    
    ###data
    data_code <- 'data{\n\tint<lower=1> N;\n\tint<lower=1> P;\n\trow_vector[P] x[N];\n'
    if(R>0){
      data_code <- paste0(data_code,"\tint<lower=1> R;\n\tint<lower=1> G[R];\n\tint<lower=1> Q[R];\n" )
    }
    if(family=="binomial") data_code <- paste0(data_code,"\tint<lower=1> bitotal[N];\n" )
    if(family=="ordered") data_code <- paste0(data_code,"\tint<lower=2> K;\n" )
    if(checkoffset==1) data_code <- paste0(data_code,"\treal offset[N];\n" )
      
    temp1 <- ''
    if(family=="gaussian"){
      temp1 <- ("\treal y[N];\n")
    }else if(family=="gamma" || family=="lognormal"){
      temp1 <- ("\treal<lower=0> y[N];\n")
    }else if(family=="bernoulli" ||family=="binomial" || family=="poisson"|| family=="ordered" || family=="nbinomial"){
      temp1 <- ("\tint<lower=0> y[N];\n")
    }else if(family=="beta" ){
      temp1 <- ("\treal<lower=0,upper=1> y[N];\n")
    }
    if(R>0){
      for(i in 1:R){
        temp1 <- paste0(temp1,"\t","int<lower=1> ", "id",i,"[N];","\n")
      }
    }
    temp2 <- ''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp2 <- paste0(temp2,"\t","real", " z",i,"[N];","\n")
        }else{
          temp2 <- paste0(temp2,"\t","row_vector[Q[",i, "]] z",i,"[N];","\n")
        }  
      }
    }
    if(family=="gaussian"||R==1) temp2 <- paste0(temp2,"\t","int idn1[G[1]];\n")
    data_code <- paste0(data_code,temp1,temp2,"}")
    
    ###transformed data    
    td_code <-'\ntransformed data{\n'
    if(R>0){
      temp1 <- ''
      for(i in 1:R){
        temp1 <- paste0(temp1,"\t","vector[Q[",i, "]] mu",i,";\n")
      }
      temp2 <- ''
      for(i in 1:R){
        temp2 <- paste0(temp2,"\t","for(q in 1:Q[",i,"]) mu" ,i,"[q] <- 0;\n")
      }
      td_code <- paste0(td_code,temp1,temp2,"}")
    }else{
      td_code <- paste0(td_code,"}")
    }
    
    ###parameters
    para_code <-'\nparameters{\n\tvector[P] beta;\n'
    temp1 <- ''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp1 <- paste0(temp1,"\t","real r",i,"[G[",i,"]];\n")
        }else{
          temp1 <- paste0(temp1,"\t","vector[Q[",i, "]] r",i,"[G[",i,"]];\n")
        } 
      }
    }
    
    temp2 <- ''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp2 <- paste0(temp2,"\t","real<lower=0> tau_sd",i,";\n")
        }else{
          temp2 <- paste0(temp2,"\t","vector<lower=0>[Q[",i, "]] tau_sd",i,";\n")
        } 
      }
      for(i in 1:R){
        if(Q[i]>1){
          temp2 <- paste0(temp2,"\t","corr_matrix[Q[",i, "]] tau_corr",i,";\n")
        } 
      }
    }
    temp3 <- ''
    if(family == "gaussian" || family == "gamma" || family=="nbinomial" || family=="lognormal" || family=="beta"){
      temp3 <- paste0("\t","real<lower=0> s;\n")
    }
    if(family=="ordered"){
      temp3 <- paste0(temp3,"\t","ordered[K-1] cutpoints;\n")
    }
    para_code <- paste0(para_code,temp1,temp2,temp3,"}")
    
    ###transformed parameters
    tp_code <-'\ntransformed parameters{\n'
    temp1 <-''
    if(family == "gaussian" || family == "gamma" || family=="nbinomial" || family=="lognormal" || family=="beta"){
      temp1 <- paste0(temp1,"\t","real<lower=0> scale;\n")
    }
    temp2 <- ''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp2 <- paste0(temp2,"\t","real<lower=0> tau",i,";\n")
        }else{
          temp2 <- paste0(temp2,"\t","cov_matrix[Q[",i, "]] tau",i,";\n")
        } 
      }
    }
    temp3 <- ''
    if(family == "gaussian"){
      temp3 <- paste0(temp3,"\tscale <- s^2;\n")
    }else if(family == "gamma"){      
      temp3 <- paste0(temp3,"\tscale <- 1/s;\n")
    }else if(family == "lognormal"){      
      temp3 <- paste0(temp3,"\tscale <- s^2;\n")
    }else if(family == "nbinomial"){
      temp3 <- paste0(temp3,"\tscale <- 1/s;\n")
    }else if(family == "beta"){
      temp3 <- paste0(temp3,"\tscale <- s;\n")
    } 
    temp4 <- ''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp4 <- paste0(temp4,"\t","tau",i," <- tau_sd",i,"^2;\n")
        }else{
          temp4 <- paste0(temp4,"\ttau",i," <- quad_form_diag(tau_corr",i,",tau_sd",i,");\n")
        } 
      }
    }
    
    tp_code <- paste0(tp_code,temp1,temp2,temp3,temp4,"}")
    
    ###model
    if(family!="beta"){
      model_code <-'\nmodel{\n\treal predict[N];\n\tbeta ~ normal(0,1000);\n'
    }else{
      model_code <-'\nmodel{\n\treal predict[N];\n\treal A[N];\n\treal B[N];\n\tbeta ~ normal(0,1000);\n'
    }
    
    temp3 <-''
    if(R>0){
      for(i in 1:R){
        if(Q[i]==1){
          temp3 <- paste0(temp3,"\t","for(g in 1:G[",i,"]) r",i,"[g] ~ normal(0, tau_sd",i,");\n")
        }else{
          temp3 <- paste0(temp3,"\t","r",i," ~ multi_normal(mu",i,",tau",i,");\n")
        } 
      }
      for(i in 1:R){
        if(Q[i]==1){
          temp3 <- paste0(temp3,"\t","tau_sd",i," ~ cauchy(0,",cauchy,");\n")
        }else{
          temp3 <- paste0(temp3,"\t","tau_sd",i," ~ cauchy(0,",cauchy,");\n")
        } 
      }
      for(i in 1:R){
        if(Q[i]>1){
          temp3 <- paste0(temp3,"\t","tau_corr",i," ~ lkj_corr(",lkj_corr,");\n")
        } 
      }
    }
    
    if(family == "gaussian" || family=="lognormal"){
      temp3 <- paste0(temp3,"\t","s ~ cauchy(0,",cauchy,");\n")
    }
        
    temp1 <- '\tfor(n in 1:N){\n\t\tpredict[n] <- x[n]*beta'
    if(R>0){
      for(i in 1:R){
        temp1 <- paste0(temp1,"+z",i,"[n]*r",i,"[id",i,"[n]]")
      }
    }
    temp1 <-paste0(temp1,";\n")
    
    if(family=="binomial"){
      temp1 <- paste0(temp1,"\t\tpredict[n] <- inv_logit(predict[n]);\n")
    }else if(family=="poisson"){
      if(checkoffset==0){
        temp1 <- paste0(temp1,"\t\tpredict[n] <- exp(predict[n]);\n")
      }else{
        temp1 <- paste0(temp1,"\t\tpredict[n] <- offset[n]*exp(predict[n]);\n")
      }      
    }else if(family=="gamma"){
      temp1 <- paste0(temp1,"\t\tpredict[n] <- s / exp(predict[n]);\n")
    }else if(family=="lognormal"){
      
    }else if(family=="nbinomial" && checkoffset == 1){
      temp1 <- paste0(temp1,"\t\tpredict[n] <- log(offset[n])+predict[n];\n")
    }else if(family=="beta"){
      temp1 <- paste0(temp1,"\t\tpredict[n] <- inv_logit(predict[n]);\n")
      temp1 <- paste0(temp1,"\t\tA[n] <- predict[n]*s;\n")
      temp1 <- paste0(temp1,"\t\tB[n] <- (1.0-predict[n])*s;\n")
    }
    temp1 <-paste0(temp1,"\t}\n")
    
    temp2 <-'\t'
    if(family == "gaussian"){
      temp2 <- paste0(temp2,"y ~ normal(predict, s)")
    }else if(family=="bernoulli"){
      temp2 <- paste0(temp2,"y ~ bernoulli_logit(predict)") 
    }else if(family=="binomial"){
      temp2 <- paste0(temp2,"y ~ binomial(bitotal,predict)") 
    }else if(family=="poisson"){
      temp2 <- paste0(temp2,"y ~ poisson(predict)")
    }else if(family=="gamma"){
      temp2 <- paste0(temp2,"y ~ gamma(s,predict)")
    }else if(family=="lognormal"){
      temp2 <- paste0(temp2,"y ~ lognormal(predict,s)")
    }else if(family=="ordered"){
      temp2 <- paste0(temp2,"for(n in 1:N) y[n] ~ ordered_logistic(predict[n],cutpoints)")
    }else if(family=="nbinomial"){
      temp2 <- paste0(temp2,"y ~ neg_binomial_2_log(predict,s)")
    }else if(family=="beta"){
      temp2 <- paste0(temp2,"y ~ beta(A, B);\n")
    }
    temp2 <- paste0(temp2,";\n") 
    
    model_code <- paste0(model_code,temp3,temp1,temp2,"}")
    
    ###generated quantities
    gq_code <-'\ngenerated quantities{\n\treal predict[N];\n\treal log_lik[N];\n'
    if(family=="gaussian"||R==1) gq_code <- paste0(gq_code,"real log_lik_g[G[1]];\nint count;\n")
    temp1 <- ''
    if(checkslice>0){
      for(i in 1:checkslice){
        temp1 <- paste0(temp1,"\treal simple",i,"_high;\n")
        temp1 <- paste0(temp1,"\treal simple",i,"_low;\n")
      }
    }  
    temp2 <- '\tfor(n in 1:N){\n\t\tpredict[n] <- x[n]*beta'
    if(R>0){
      for(i in 1:R){
        temp2 <- paste0(temp2,"+z",i,"[n]*r",i,"[id",i,"[n]]")
      }  
    }  
    temp2 <- paste0(temp2,";\n")
    
    temp3 <-paste0("\t\tlog_lik[n] <- ")
    if(family == "gaussian"){
      temp3 <- paste0(temp3,"normal_log(y[n], predict[n], s)") 
    }else if(family=="bernoulli"){
      temp3 <- paste0(temp3,"bernoulli_logit_log(y[n],predict[n])") 
    }else if(family=="binomial"){
      temp3 <- paste0(temp3,"binomial_log(y[n], bitotal[n],inv_logit(predict[n]))")
    }else if(family=="poisson"){
      if(checkoffset==0){
        temp3 <- paste0(temp3,"poisson_log(y[n], exp(predict[n]))")
      }else{
        temp3 <- paste0(temp3,"poisson_log(y[n], offset[n]*exp(predict[n]))")
      }
    }else if(family=="gamma"){
      temp3 <- paste0(temp3,"gamma_log(y[n], s, s/exp(predict[n]))")
    }else if(family=="lognormal"){
      temp3 <- paste0(temp3,"lognormal_log(y[n], predict[n],s)")
    }else if(family=="ordered"){
      temp3 <- paste0(temp3,"ordered_logistic_log(y[n], predict[n],cutpoints)")
    }else if(family=="nbinomial"){
      if(checkoffset==0){
        temp3 <- paste0(temp3,"neg_binomial_2_log_log(y[n], predict[n],s)")
      }else{
        temp3 <- paste0(temp3,"neg_binomial_2_log_log(y[n], log(offset[n])+predict[n],s)")
      }
    }else if(family=="beta"){
      temp3 <- paste0(temp3,"beta_log(y[n], inv_logit(predict[n])*s, (1.0-inv_logit(predict[n]))*s)")
    }
    temp3 <- paste0(temp3,";\n\t}\n")
        
    if(checkslice>0){
      for(i in 1:checkslice){
        temp3 <- paste0(temp3,"\tsimple",i,"_high <- beta[",simplenum[i],"]+beta[",intrctnum[i],"]*",round(slicesd,digits=4),";\n")
        temp3 <- paste0(temp3,"\tsimple",i,"_low <- beta[",simplenum[i],"]-beta[",intrctnum[i],"]*",round(slicesd,digits=4),";\n")
      }
    }
    temp4 <- ''
    if(family=="gaussian"||R==1){
      temp4 <- paste0(temp4,"\tcount <- 0;\n\tfor(g in 1:G[1]){\n\t\t{\n")
      temp4 <- paste0(temp4,"\t\t\tvector[idn1[g]] yn;\n")
      temp4 <- paste0(temp4,"\t\t\tvector[idn1[g]] predictn;\n")
      if(Q[1]==1){
        temp4 <- paste0(temp4,"\t\t\tvector[idn1[g]] zn;\n")
      }else{
        temp4 <- paste0(temp4,"\t\t\tmatrix[Q[1],idn1[g]] zn;\n")
      }
      temp4 <- paste0(temp4,"\t\t\tvector[idn1[g]] sn;\n")
      temp4 <- paste0(temp4,"\t\t\tmatrix[idn1[g],idn1[g]] taun;\n")
      temp4 <- paste0(temp4,"\t\t\tfor(i in 1:idn1[g]){\n")
      temp4 <- paste0(temp4,"\t\t\t\tcount <- count + 1;\n")
      temp4 <- paste0(temp4,"\t\t\t\tyn[i] <- y[count];\n")
      temp4 <- paste0(temp4,"\t\t\t\tpredictn[i] <- x[count]*beta;\n")
      if(Q[1]==1){
        temp4 <- paste0(temp4,"\t\t\t\tzn[i] <- z1[count];\n")
      }else{
        temp4 <- paste0(temp4,"\t\t\t\tfor(j in 1:Q[1]) zn[j,i] <- z1[count][j];\n")
      }
      temp4 <- paste0(temp4,"\t\t\t\tsn[i] <- scale;\n\t\t\t{\n")
      if(Q[1]==1){
        temp4 <- paste0(temp4,"\t\t\ttaun <- zn*tau1*zn'+diag_matrix(sn);\n")
      }else{
        temp4 <- paste0(temp4,"\t\t\ttaun <- zn'*tau1*zn+diag_matrix(sn);\n")
      }
      temp4 <- paste0(temp4,"\t\t\tlog_lik_g[g] <- multi_normal_log(yn,predictn,taun);\n")
      temp4 <- paste0(temp4,"\t\t{\n\t{\n")
    }
    gq_code <- paste0(gq_code,temp1,temp2,temp3,temp4,"}")
    
    codestan <- paste(data_code,td_code,para_code,tp_code,model_code,gq_code,"\n")    
  }else if(nomodel==FALSE){
    codestan <- stanmodel@model_code
  }
  
  ###################
  ###running rstan###
  ###################
  
  #not stan output
  modelname <- paste0(formula_str," [",family,"]")
  if(modelonly==TRUE){
    modelstan <- stan_model(model_name=modelname,model_code=codestan)
    attr(modelstan,"code") <- c("code"=codestan)
    return(modelstan)
  }
  if(codeonly==TRUE){
    attr(codestan,"family") <- c("family"=family)
    return(codestan)
  }
  
  #MCMCsampling
  cat("\nCompiling the stan code.\n")
  if(parallel==TRUE && chains > 1){
    if(nomodel==TRUE){
      stanmodel <- rstan::stan_model(model_name=modelname,model_code=codestan)
      nomodel <- FALSE
    }
    if(nofit==FALSE){
      stanmodel <- stanfit@stanmodel
    }
    cat("\nPreparing parallel process.\n")
    if(.Platform$OS.type != "windows") {
      cat("\nMCMC sampling start.\n")
      fit.list <- mclapply(1:chains, mc.cores=cores,function(i)
        rstan::sampling(stanmodel, data=datastan,chain_id = i,
                        iter=iter, warmup = warmup,chains=1,thin=thin)
      )
    }else{
      stopImplicitCluster()
      cl <- parallel::makePSOCKcluster(cores)
      on.exit(stopCluster(cl))
      parallel::clusterExport(cl=cl,c("datastan","stanmodel"),envir=environment())
      cat("\nMCMC sampling start.\n")
      fit.list <- parallel::parLapply(cl,1:chains, fun=function(i){
        require(rstan)
        rstan::sampling(stanmodel, data=datastan,chain_id = i,
                        iter=iter, warmup = warmup,chains=1,thin=thin)
      })
    }
    if(chains>1){
      fitstan <- rstan::sflist2stanfit(fit.list)
    }else{
      fitstan <- fit.list
    }
    #stopCluster(cl)
  }else if(nocode==TRUE && nomodel==TRUE){
    fitstan <- rstan::stan(model_name=modelname,model_code=codestan, data=datastan,
                           iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores)   
  }else{    
    if(nomodel==FALSE){
      fitstan <- rstan::sampling(stanmodel, data=datastan,
                                 iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores)
    }else if(nocode==FALSE){
      if(is.null(family)==TRUE){
        family <- attr(codestan,"family")
      }
      modelname <- paste0(formula_str," [",family,"]")
      fitstan <- rstan::stan(model_name=modelname, model_code=codestan, data=datastan, 
                             iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores)
    }else if(nofile==FALSE){
      modelname <- paste0(formula_str," [",family,"]")
      fitstan <- rstan::stan(stanfile,model_name=modelname, data=datastan,
                             iter=iter, warmup = warmup,chains=chains,thin=thin,cores=cores) 
    }
  }    
  
  ############################
  ###calculating parameters###
  ############################
  
  ###calculating WAIC
  loglik <- rstan::extract(fitstan,"log_lik")$log_lik
  lppd <- sum(log(colMeans(exp(loglik))))
  p_waic <- sum(apply(loglik,2,var))
  waic <- -lppd/N + p_waic/N
  waic2 <- waic * (2*N)
  
  ###calculating global parameter WAIC
  if(family=="gaussian"||R==1){
    loglik_g <- rstan::extract(fitstan,"log_lik_g")$log_lik_g
    lppd_g <- sum(log(colMeans(exp(loglik_g))))
    p_waic_g <- sum(apply(loglik_g,2,var))
    waic_g <- -lppd_g/G[1] + p_waic_g/G[1]
    waic2_g <- waic_g * (2*G[1])
  }
  
  ###calculating beta
  beta <- rstan::extract(fitstan,"beta")$beta
  beta_com <- matrix(c(apply(beta,2,mean),apply(beta,2,sd),apply(beta,2,quantile,0.025),apply(beta,2,quantile,0.975)),ncol=4)
  rownames(beta_com) <- colnames(x)
  colnames(beta_com) <- c("coefficient","stdev","95%lower","95%upper")
  if(family=="gaussian" ||family=="gamma"|| family=="nbinomial" || family=="lognormal" || family=="beta"){
    s <- rstan::extract(fitstan,"scale")$scale
    scale <- matrix(c(mean(s),sd(s),quantile(s,0.025),quantile(s,0.975)),ncol=4)
    rownames(scale) <- "scale"
    beta <- rbind(beta_com,scale)
  }else if(family=="ordered"){
    s <- rstan::extract(fitstan,"cutpoints")$cutpoints
    cutpoints <- matrix(c(apply(s,2,mean),apply(s,2,sd),apply(s,2,quantile,0.025),apply(s,2,quantile,0.975)),ncol=4)
    colnames(cutpoints) <- c("coefficient","stdev","95%lower","95%upper")
    temp <- c()
    for(i in 1:(K-1)) temp[i] <- paste0("cutpoints",i)
    rownames(cutpoints) <- temp
    beta <- rbind(beta_com,cutpoints)
  }else{
    beta <- beta_com
  }
    
  ###calculating tau
  if(R>0){
    tauname <- c()
    for(i in 1:R) tauname[i] <- paste("tau",i,sep="")
    tau <- list()
    for(i in 1:R){
      tau_com <- rstan::extract(fitstan,tauname[i])[tauname[i]]
      if(Q[i]>1){
        tau[[i]] <- matrix(nrow=Q[i],ncol=Q[i])
        for(j in 1:Q[i])for(k in 1:Q[i]) tau[[i]][j,k] <- mean(tau_com[[1]][,j,k])
        rownames(tau[[i]]) <- zname[[i]]
        colnames(tau[[i]]) <- zname[[i]]
      }else{
        tau[[i]] <- mean(tau_com[[1]])
        names(tau[[i]]) <- zname[[i]]
      }
    }
    names(tau) <- idname
  }
  
  ###calculating tau's SD
  if(R>0){
    tausdname <- c()
    for(i in 1:R) tausdname[i] <- paste("tau_sd",i,sep="")
    tausd <- list()
    for(i in 1:R){
      tausd_com <- rstan::extract(fitstan,tausdname[i])[tausdname[i]]
      if(Q[i]>1){
        tausd[[i]] <- matrix(c(apply(tausd_com[[1]],2,mean),apply(tausd_com[[1]],2,sd),
                          apply(tausd_com[[1]],2,quantile,0.025),
                          apply(tausd_com[[1]],2,quantile,0.975)),ncol=4)
        rownames(tausd[[i]]) <- zname[[i]]
      }else{
        tausd[[i]] <- matrix(c(mean(tausd_com[[1]]),sd(tausd_com[[1]]),quantile(tausd_com[[1]],0.025),
                          quantile(tausd_com[[1]],0.975)),ncol=4)
        rownames(tausd[[i]]) <- zname[[i]]
      }
      colnames(tausd[[i]]) <- c("coefficient","stdev","95%lower","95%upper")
    }
    names(tausd) <- idname
  }
  
  ###calculating tau's correlations
  taucorrcheck <- 0
  if(R>0){
    taucorrname <- c()
    for(i in 1:R) if(Q[i]>1) taucorrname[i] <- paste("tau_corr",i,sep="")
    taucorr <- list()
    for(i in 1:R){
      if(Q[i]>1){
        taucorrcheck <- 1
        taucorr_com <- rstan::extract(fitstan,taucorrname[i])[taucorrname[i]]
        taucorr[[i]] <- matrix(nrow=Q[i],ncol=Q[i])
        for(j in 1:Q[i])for(k in 1:Q[i]) taucorr[[i]][j,k] <- mean(taucorr_com[[1]][,j,k])
        rownames(taucorr[[i]]) <- zname[[i]]
        colnames(taucorr[[i]]) <- zname[[i]]
      }
    }
    if(length(taucorr)>0) names(taucorr) <- idname
  }
  
  ###calculating simple slope
  if(checkslice>0){
    for(i in 1:checkslice){
      slopehighname <- paste0("simple",i,"_high")
      slopelowname <- paste0("simple",i,"_low")
      slopehigh <- rstan::extract(fitstan,slopehighname)[slopehighname][[1]]
      slopelow <- rstan::extract(fitstan,slopelowname)[slopelowname][[1]]
      slopehigh_com <- matrix(c(mean(slopehigh),sd(slopehigh),quantile(slopehigh,0.025),quantile(slopehigh,0.975)),ncol=4)
      slopelow_com <- matrix(c(mean(slopelow),sd(slopelow),quantile(slopelow,0.025),quantile(slopelow,0.975)),ncol=4)
      rownames(slopehigh_com) <- paste0(xname[simplenum[i]], "_high")
      rownames(slopelow_com) <- paste0(xname[simplenum[i]], "_low")
      simple <- rbind(slopehigh_com,slopelow_com)
    }
    colnames(simple) <-  c("coefficient","stdev","95%lower","95%upper")
  }
  
  ############
  ###output###
  ############
  
  paraname <- c("beta")
  if(family=="gaussian" ||family=="gamma"|| family=="nbinomial"||family=="lognormal"){
    paraname <- c(paraname,"scale")
  }
  if(R>0){
    #for(i in 1:R){
    #  paratau <- paste0("tau",i)
    #  paraname <- c(paraname,paratau)
    #}
    for(i in 1:R){
      paratausd <- paste0("tau_sd",i)
      paraname <- c(paraname,paratausd)
    }
    for(i in 1:R){
      if(Q[i]>1){
        parataucorr <- paste0("tau_corr",i)
        paraname <- c(paraname,parataucorr)
      }      
    }
  }
  if(family=="ordered"){
    paraname <- c(paraname,"cutpoints")
  }
  if(checkslice>0){
    for(i in 1:checkslice){
      temp1 <- paste0("simple",i,"_high")
      temp2 <- paste0("simple",i,"_low")
      paraname <- c(paraname,temp1,temp2)
    }
  }
  
  attr(fitstan,"family")<-c("family"=family)
  attr(fitstan,"paraname")<-c("paraname"=paraname)
  attr(fitstan,"WAIC") <-c("WAIC"=waic2,"lppd"=lppd,"p_waic" = p_waic)
  if(family=="gaussian"||R==1){
    attr(fitstan,"WAIC_g") <-c("WAIC_g"=waic2_g,"lppd_g"=lppd_g,"p_waic_g" = p_waic_g)
  }
  attr(fitstan,"dataname") <- list("yname"=yname,"xname"=xname,"zname"=zname,"idname"=idname)
  attr(fitstan,"beta")<-list("beta"=beta)
  if(R>0){
    attr(fitstan,"tau")<-list("tau"=tau)
    attr(fitstan,"tau_sd")<-list("tau_sd"=tausd)
    if(taucorrcheck ==1) attr(fitstan,"taucorr")<-list("taucorr"=taucorr)
  }
  attr(fitstan,"formula") <- c("formula"=formula)
  attr(fitstan,"code") <-  c("code"=stancode)
  if(checkslice>0){
    attr(fitstan,"simple")<-list("simple"=simple) 
  }
  if(family=="gaussian"||R==1){
    attr(fitstan,"global") <- TRUE
  }else{
    attr(fitstan,"global") <- FALSE
  }
  
  outputwaic <- paste0("\nlppd = ",round(lppd,digits=4),"\npWAIC = ",round(p_waic,digits=4),
                       "\nWAIC = ", round(waic2,digits=4))
  cat(outputwaic)
  
  return(fitstan)
}

output_stan<- function(fitstan,digits = 3){
  paraname <- attr(fitstan,"paraname")
  return(print(fitstan,pars=paraname,probs=c(0.025,0.975),digits=digits))
}

output_code <- function(fitstan){
  if(class(fitstan)=="stanfit"){
    return(cat(fitstan@stanmodel@model_code))
  }else if(class(fitstan)=="stanmodel"){
    return(cat(fitstan@model_code))
  }else{
    cat(fitstan)
  }  
}

output_result <- function(fitstan){
  family <- attr(fitstan,"family")
  formula <- attr(fitstan,"formula")
  WAIC <- list("WAIC"=attr(fitstan,"WAIC"))
  if(attr(fitstan,"global")==TRUE){
    WAIC_g <- list("WAIC_g"=attr(fitstan,"WAIC_g"))
  }
  beta <- attr(fitstan,"beta")
  simple <- attr(fitstan,"simple")
  tau_sd <- attr(fitstan,"tau_sd")
  taucorr <- attr(fitstan,"taucorr")
  tau <- attr(fitstan,"tau")
  result <- c(formula,WAIC,beta,simple,tau_sd,taucorr,tau)
  return(result)
}

output_beta <- function(fitstan){
  family <- attr(fitstan,"family")
  beta <- as.data.frame(rstan::extract(fitstan,"beta")$beta)
  colnames(beta) <- attr(fitstan,"dataname")$xname
  for(i in 1:length(colnames(beta))){
    if(colnames(beta)[i]=="(Intercept)"){
      colnames(beta)[i] <- "Intercept"
    }
  }
  if(family=="gaussian" ||family=="gamma"|| family=="nbinomial"||family=="lognormal"){
    beta$scale <- rstan::extract(fitstan,"scale")$scale
  }
  return(beta)
}

output_tausd <- function(fitstan,variname=NULL){
  idname <- attr(fitstan,"dataname")$idname
  if(length(idname)>1){
    if(is.null(variname)) stop(paste0("Please input random effect's name"))
    for(i in 1:length(idname)){
      if(idname[i]==variname) idnum <- i
    }
  }else{
    idnum = 1
  }
  paraname <- paste("tau_sd",idnum,sep="")
  tausd <- as.data.frame(rstan::extract(fitstan,paraname)[paraname])
  colnames(tausd) <- attr(fitstan,"dataname")$zname[[idnum]]
  return(tausd)
}

output_ggmcmc <- function(fitstan,para="null"){
  library(ggmcmc)
  if(para=="beta"){
    Label <- attr(fitstan,"dataname")$xname
    Parameter <- c()
    for(i in 1:length(Label)) Parameter[i] <- paste0(para,"[",i,"]")
    labels <- data.frame(Parameter,Label)
    S <- ggs(fitstan,par_labels=labels,family=paste0("^",para))
  }else if(para=="null"){
    S <- ggs(fitstan)
  }else{
    S <- ggs(fitstan,family=paste0("^",para))
  }
  return(S)
}

Pglmmstan <- function(formula_str,data,family="gaussian",center = FALSE,slice = NULL,offset=NULL,                     
                     codeonly=FALSE,dataonly=FALSE,modelonly=FALSE,cauchy = 2.5,lkj_corr = 2,
                     stancode=NULL,standata=NULL,stanmodel=NULL,stanfile=NULL,stanfit=NULL,
                     parallel=TRUE,cores=NULL,iter=2000,warmup = NULL,chains= 2,thin=1){
  
  glmmstan(formula_str,data,family=family,center = center,slice = slice,offset=offset,                     
           codeonly=codeonly,dataonly=dataonly,modelonly=modelonly,cauchy = cauchy,lkj_corr = lkj_corr,
           stancode=stancode,standata=standata,stanmodel=stanmodel,stanfile=stanfile,stanfit=stanfit,
           parallel=TRUE,cores=cores,iter=iter,warmup = warmup,chains= chains,thin=thin)
}

map_mcmc <- function(z){
  density(z)$x[which.max(density(z)$y)]
}

