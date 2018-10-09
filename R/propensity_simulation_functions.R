#==============================================================================
#
# Functions used for propensity score simulations (see propensity_simulations.R)
#
#==============================================================================

#--------------------------------
# function to produce Synthetic data 
#--------------------------------
sim_data <- function(n = 500, scenario=1, pint=0.5, t_coef=1, sigy=1) {
  # Simulate data, outcome and treatment effects with variables affecting
  # both treatment assignment as well as outcome (confounding). Treatment
  # effects are either homogeneous (scenario 1) or heterogeneous (scenario 2)
  # In the following X3 and X4 are confounders

  X <- MASS::mvrnorm(n, mu = rep(0, 6), Sigma=diag(6))
  qint<- qlogis(pint)  # logit scale Treatment model intercept
  var_names<- paste0("V",1:6)
  pvars<- c(3,4,5,6)  # Index of variables affecting treatment assignment
  yvars<- c(1,2,3,4)  # index of variables affecting outcome
  pcoef<- c(qint,c(-0.2, 0.8, -0.2, 0.8)) # Treatment variable coefficients
  ycoef<- c(1, -0.9, 0.9, -0.9, 0.9)      # Outcome variable coefficients

  tprob <- plogis(cbind(1,X[,pvars]) %*% pcoef)

  # Treatment indicator
  Tr<- rbinom(n, 1, prob=tprob)

  if(scenario==1) {
    #Homogeneous treatment effect
    delta<- t_coef
  }
  else if (scenario==2) {
    # heterogeneous tretament effect symmetrical around t_coef
    delta<- ifelse(X[,4] <= 0, t_coef - X[,4]^2, t_coef + X[,4]^3)
  }
  # Untreated state
  y0 <- cbind(1,X[,yvars]) %*% ycoef + rnorm(n, 0, sigy)
  # Treated state
  y1 <- y0 + delta
  # observed outcome
  y <- y1 * Tr + (1 - Tr) * y0

  Data <- as.data.frame(X)
  names(Data)<- var_names
  Data$Tr <- Tr
  Data$delta<- delta
  Data$y0 <- y0
  Data$y1 <- y1
  Data$y <- y
  SATE<- mean(Data$delta) # Sample average treatment effect
  SATT<- mean(Data$delta[Data$Tr==1]) # sample average treatment effect on treated
  list(Data=Data,SATE=SATE,SATT=SATT)
}
#--------------------------------
# Funtion to perform propensity score 
# analysis on simulated data
#--------------------------------
pscore.sim<- function(r, j, nopt, v_opt, scenario=1, t_coef=1, pint=0.5,
                      est="ATE") {
  
  sdata <- sim_data(n = nopt, scenario=scenario, t_coef=t_coef, pint=pint[j])
  the_data<- sdata$Data 
  
  if(est=="ATE") target<- sdata$SATE
  else target<- sdata$SATT
  
  var_names<- paste0("V",1:6)
  if(!is.null(v_opt)) var_names<- var_names[!var_names %in% v_opt]
  
  form1<- fbuild("Tr", var_names)
  form2<- fbuild("y", var_names, "Tr")
  
  # some structures for holding results
  te<- se<- pv<- rep(NA, 5)
  if(est=="ATE"){
    ps.nam<- c("Trt","Reg","IPW","IPW_dr","FM")
  } else {
    ps.nam<- c("Trt","PM","IPW","IPW_dr","FM")
  }
                  
  #  Initial check for separation or non-convergence and exit if fail
  if((sum(the_data$Tr)==0) | !(glm(form1,data=the_data,family=binomial)$converged)) {
    df<- data.frame(TE=te,SE=se,Cov=pv)
    df$rep<- r
    df$n<- nopt
    df$vars_in<- j
    df$est<- est
    df$target<- round(target,2)
    df$method<- ps.nam 
    return(df)
  }
  else pscore<- glm(form1,data=the_data,family=binomial)
    
  
  #-- Simple t-test of treatment effect
  mod1 <- glm(y ~ Tr, data=the_data)
  te[1] <- coef(mod1)["Tr"] 
  se[1] <- SE(mod1)["Tr"]
  pv[1] <- ifelse((te[1]-1.96*se[1] < target) & (te[1]+1.96*se[1] > target), 1, 0)
  
  if(est == "ATE") { # Do regression for ATE
    #-- ANCOVA regression
    mod2 <- glm(form2, data=the_data)
    te[2] <- coef(mod2)["Tr"] 
    se[2] <- SE(mod2)["Tr"]
    pv[2] <- ifelse((te[2]-1.96*se[2] < target) & (te[2]+1.96*se[2] >target), 1, 0)
  }
  else { # do pair matching for ATT
    #--- Matching - optimal for ATT (no caliper)-------------
    # using the optmatch package
    fmatch <- match_on(pscore, data=the_data) 
    sclass<- pairmatch(fmatch, data=the_data)
    the_data$w<- get.match.weights(sclass, the_data$Tr, estimand="ATT") # Always ATT for pairmatching 
    if(!all(the_data$w == 0)) {
      design.ps<- svydesign(ids= ~1, weights= ~w, data=the_data)
      mod2 <- svyglm(y ~ Tr, design=design.ps)
      te[2] <- coef(mod2)["Tr"] 
      se[2] <- SE(mod2)["Tr"]
      pv[2]<- ifelse((te[2]-1.96*se[2] < target) & (te[2]+1.96*se[2] > target), 1, 0)
      }
  }
  #--- IPTW Inverse Probability of Treatment Weights 
  
  the_data$w<- calc.pswts(pscore$fitted.values, the_data$Tr, estimand=est, trim=FALSE) 
  if(!all(the_data$w == 0)) {
    design.ps<- svydesign(ids= ~1, weights= ~w, data=the_data)
    mod3 <- svyglm(y ~ Tr, design=design.ps)
    te[3] <- coef(mod3)["Tr"] 
    se[3] <- SE(mod3)["Tr"]
    pv[3]<- ifelse((te[3]-1.96*se[3] < target) & (te[3]+1.96*se[3] > target), 1, 0)
  }
  
  #--- IPTW Inverse Probability of Treatment Weights doubly robust)
  the_data$w<- calc.pswts(pscore$fitted.values, the_data$Tr, estimand=est, trim=FALSE) 
  if(!all(the_data$w == 0)) {
    design.ps<- svydesign(ids= ~1, weights= ~w, data=the_data)
    mod4 <- svyglm(form2, design=design.ps)
    te[4] <- coef(mod4)["Tr"] 
    se[4] <- SE(mod4)["Tr"]
    pv[4] <- ifelse((te[4]-1.96*se[4] < target) & (te[4]+1.96*se[4] > target), 1, 0)
  }
  
  #--- Matching - full (no caliper)-------------
  # using the optmatch package
  fmatch <- match_on(pscore, data=the_data) 
  sclass<- fullmatch(fmatch, data=the_data)
  the_data$w<- get.match.weights(sclass, the_data$Tr, estimand=est) 
  if(!all(the_data$w == 0)) {
    design.ps<- svydesign(ids= ~1, weights= ~w, data=the_data)
    mod5 <- svyglm(y ~ Tr, design=design.ps)
    te[5] <- coef(mod5)["Tr"] 
    se[5] <- SE(mod5)["Tr"]
    pv[5] <- ifelse((te[5]-1.96*se[5] < target) & (te[5]+1.96*se[5] > target), 1, 0)
  }
  
  
  #-------------------------------------------------------------------------
  # Gather results
  df<- data.frame(TE=te,SE=se,Cov=pv)
  df$rep<- r
  df$n<- nopt
  df$vars_in<- j
  df$est<- est
  df$target<- round(target,2)
  df$method<- ps.nam
  return(df)
}

#--------------------------------
# IPTW weighting
#--------------------------------
calc.pswts<- function(p, trt, estimand="ATE", trim=FALSE, q=c(0.1,0.9), stabilise=FALSE) {
  # Propensity score IPTW calculation
  # p is predicted probability of treatment from logistic regression
  # trt is the treatment indicator variable (0/1)
  # Estimand is either ATT or ATE
  # stabilised weights (stab=T) add ptrt as numerator
  
  if(trim){
    pt<- ifelse((p < q[1]) | (p > q[2]), 0, p)
  } 
  else pt<- p
  
  if(estimand=="ATE") {
    w<- trt/pt + (1-trt)/(1-pt)
    w[pt==0]<- 0
  }
  else if(estimand=="ATT") {
    w<- trt + (1-trt)*pt/(1-pt)
    w[trt==0 & pt==0]<- 0
    w[trt==1]<- 1
  }
  else if(estimand=="ATO") {
    # no trimming required with ATO
    w<- trt * (1-p) + (1-trt)*p
  }
  else stop(cat(estimand," does not exist"))
  
  if(stabilise) {
    ptrt<- mean(trt) # marginal probability of treatment
    w[trt==1]<- w[trt==1]*ptrt
    w[trt==0]<- w[trt==0]*(1-ptrt)
  }   
  return(w)
}
#--------------------------------
# Weights from full matching
#--------------------------------
get.match.weights<- function(wts, Tr, estimand="ATE") {
  # calculate ATE or ATT weights from a fullmatch object from matchit
  q<- mean(Tr)
  cnts<- table(wts,Tr)
  w<- rep(0,length(Tr))
  classes<- levels(wts)
  n<- length(classes)
  if(n < 1) return(w)
  if(ncol(cnts) != 2) stop("no matches for at least one level")
  if(estimand=="ATT") {
    # ATT weights
    w[Tr==1]<- 1
    for(i in 1:n){
      w[wts==classes[i] & Tr==0]<- cnts[i,2]/cnts[i,1]
    }
    numt<- sum(cnts[,2])
    if(numt > 0) w[Tr==0]<- w[Tr==0] * numt/sum(w[Tr==0]) # reweight untreated 
  }
  else {
    # ATE weights 
    w0<- rowSums(cnts)/cnts[,1]
    w1<- rowSums(cnts)/cnts[,2]
    for(i in 1:n) {
      w[wts==classes[i] & Tr==0]<- w0[i]
      w[wts==classes[i] & Tr==1]<- w1[i]
    }
     numb<- sum(cnts)
     if(numb > 0){
      w[Tr==0]<- w[Tr==0] * numb/sum(w[Tr==0])
      w[Tr==1]<- w[Tr==1] * numb/sum(w[Tr==1])
     }
  }
  return(w)
}  

#----------------------------------------------------------
# Formula builder 
#
fbuild<- function(resp, indp, trt=NULL, interaction=NULL) {
  if(is.null(trt)) {
    if(!is.null(interaction)) {
      int2<- t(combn(interaction,2))
      all.int<- apply(int2, 1, function(x) paste(x,collapse=":",sep=""))
      form<- paste0(resp, "~", paste0(c(indp,all.int), collapse="+"))
    } else form<- paste0(resp, "~", paste0(indp, collapse="+"))
  }
  else {
    if(!is.null(interaction)) {
      int2<- t(combn(interaction,2))
      all.int<- apply(int2, 1, function(x) paste(x,collapse=":",sep=""))
      form<- paste0(resp, "~", paste0(c(indp,all.int,trt), collapse="+"))
    } else form<- paste0(resp, "~", paste0(c(indp,trt), collapse="+"))
  }
  as.formula(form)
}


#-----------------------------------------------------
ess<- function(w, Tr) {
  # Effective sample size of weighted data
  w0<- w[Tr==0]
  w1<- w[Tr==1]
  n<- table(Tr)
  ess0<- sum(w0)^2/sum(w0^2)
  ess1<- sum(w1)^2/sum(w1^2)
  n0<- sum(w0>0)
  n1<- sum(w1>0)
  data.frame(n0=n0/n[1],n1=n1/n[2],ess0=round(ess0),ess1=round(ess1))
}
