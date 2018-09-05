library(readr)
library(tidyverse)
library(forcats)
library(optmatch)

source("R/propensity_simulation_functions.R")

Alldata<- read_csv("Data/tree_cover_data.csv")

spp.name<- c("Melicytus ramiflorus","Metrosideros umbellata",
             "Weinmannia racemosa","Raukaua simplex")
n<- length(spp.name)

#-----------------------------
# P.score distribution
#-----------------------------
Data.list<- list()

for(i in 1:n) {
  Data<- filter(Alldata, Species==spp.name[i])
  var_names<- names(dplyr::select(Data, BA:MAS))
  form<- fbuild("Control", var_names)
  
  # Propensity scores
  pscore<- glm(form, data=as.data.frame(Data), family=binomial) 
  Data$p.score<- pscore$fitted.values
  
  # Pair matching on ATT
  fmatch <- match_on(pscore, data=Data) 
  wts<- pairmatch(as.matrix(fmatch), data=Data)
  Data$PM<- get.match.weights(wts, Data$Control, "ATT")  
  
  # Full matching on ATT
  wts<- fullmatch(as.matrix(fmatch), data=Data)
  Data$FM<- get.match.weights(wts, Data$Control, "ATT")  
  
  #IPTW
  Data$IPW<- calc.pswts(Data$p.score, Data$Control, "ATT", stabilise=TRUE)
  
  Data$UA<- 1  # No adjustment
  
  Data<- Data %>% dplyr::select(Control, p.score, UA, PM, FM, IPW) 
  Data.list[[spp.name[i]]]<- Data
}

for(i in 1:n) Data.list[[i]]$Species<- spp.name[i]

Data<- do.call('rbind',Data.list)
Data<- Data %>% gather(Method, w, -Control, -Species, -p.score)
Data<- Data %>% group_by(Species, Control, Method) %>% mutate(weights=w/sum(w))

Data<- Data %>% ungroup(Control) %>% 
  mutate(Control = lvls_revalue(factor(Control), c("Untreated","Treated")),
         Method = lvls_revalue(factor(Method), c("FM","IPW","PM","Unadjusted")),
         Method = lvls_reorder(Method,c(4,3,1,2)))

win.graph(10,10)
Data %>%  ggplot(aes(p.score, fill=Control, weights=weights)) +
  geom_density(alpha=0.4, adjust=1.5) +
  facet_grid(Method ~ Species) +
  labs(x="Propensity score", y="Density", fill="Site") +
  theme_bw() +
  theme(strip.text.x = element_text(face="italic"))

#--------------------------------------------------
# Covariate Balance using IPTW. 
# requires package 'cobalt'
#--------------------------------------------------
library(cobalt)

bal.list<- list()

for(i in 1:n) {
  Data<- filter(Alldata, Species==spp.name[i])
  var_names<- names(dplyr::select(Data, BA:MAS))
  form<- fbuild("Control", var_names)
  
  # Logistic regression
  
  pscore<- glm(form, data=Data, family=binomial) 
  Data$p.score<- pscore$fitted.values
  
  Data$wts<- calc.pswts(Data$p.score, Data$Control, "ATT", stabilise=TRUE)
  
  Data.list[[spp.name[i]]]<- Data
  
  bal.list[[spp.name[i]]]<- bal.tab(form, weights="wts", method="weighting", 
                                    distance ="p.score", disp.v.ratio=T, 
                                    un=T, m.threshold=0.1, v.threshold=2, estimand="ATT", 
                                    data=as.data.frame(Data))$Balance
  }


vv<- names(bal.list)
n<- length(vv)
btable<- list()
for(i in 1:n) {
  tmp<- bal.list[[i]]
  cvars<- row.names(tmp)
  tmp<- tmp %>% dplyr::select(Diff.Un, Diff.Adj, V.Ratio.Un, V.Ratio.Adj) %>% 
    mutate(Covariates=cvars,Species=vv[i])
  btable[[i]]<- tmp
}

btable<- do.call('rbind',btable)

win.graph(7,10)
btable %>% dplyr::select(Diff.Un, Diff.Adj, Covariates, Species) %>% 
  gather(MSD, value, -Covariates, - Species) %>%
  filter(Covariates != "p.score") %>%
  mutate(MSD = lvls_reorder(MSD, c(2,1)), 
         MSD = lvls_revalue(MSD, c("Unadjusted","Adjusted")),
         value = if_else(value > 1, 1, value),
         value = if_else(value < -1, -1, value)) %>% 
  ggplot(aes(x = value, shape=MSD, y = Covariates)) +
  geom_point(size=2) +
  scale_shape_manual(values=c(1,16)) +
  geom_line(aes(group=Covariates),color="black") +
  coord_cartesian(xlim = c(-1.1,1.1)) +
  facet_wrap( ~ Species, nrow=2) +
  labs(y = "Covariate", x = "Mean standardised difference") +
  geom_vline(xintercept = c(-0.1,0.1), linetype=2) +
  theme_bw() +
  theme(legend.position="none",
        strip.text.x = element_text(face="italic")) 

#-----------------------------
#   Analysis
#-----------------------------
# Three species where balance was reasonable

library(survey)

spp1<- c("Metrosideros umbellata","Weinmannia racemosa","Raukaua simplex")
n<- length(spp1)
# To hold results
results <- list()

for(i in 1:n) {
  Data<- filter(Alldata, Species==spp1[i])
  var_names<- names(dplyr::select(Data, BA:MAS))
  form<- fbuild("Control", var_names)
  
  nplots<- Data %>% group_by(Control) %>% summarise(Nplots=n())
  #==================================================================
  # Calc propensity weighting
  
  cf<- se<- pv<- rep(NA, 3)  # for results
  #---------------------------
  # Logistic regression
  
  pscore<- glm(form, data=Data, family=binomial) 
  Data$p.score<- pscore$fitted.values
  
  Data$wts<- calc.pswts(Data$p.score, Data$Control, "ATT", trim=FALSE)
  
  design.ps<- svydesign(ids= ~1, weights=~wts, data=Data)
  
  form<- fbuild("log(Y)", var_names, "Control")
  
  mod1<- glm(log(Y) ~ Control, data=Data)
  mod2<- svyglm(log(Y) ~ Control, design=design.ps)
  mod3 <- svyglm(form, design=design.ps)
  
  cf[1] <- coef(mod1)["Control"]
  se[1] <- SE(mod1)["Control"]
  pv[1] <- coef(summary(mod1))['Control','Pr(>|t|)']
  
  cf[2] <- coef(mod2)["Control"]
  se[2] <- SE(mod2)["Control"]
  pv[2] <- coef(summary(mod2))['Control','Pr(>|t|)']
  
  cf[3] <- coef(mod3)["Control"]
  se[3] <- SE(mod3)["Control"]
  pv[3] <- coef(summary(mod3))['Control','Pr(>|t|)']
  
  
  results[[i]]<- data.frame(Species=spp1[i],method=c("Trt","IPW","IPW_dr"),coef=cf,se=se,pv=pv)
  
}

results<- do.call('rbind',results)

#-----

new_levels<- c("Trt","IPW","IPW_dr")


res<- results %>% mutate(LCL=coef-1.96*se,UCL=coef+1.96*se)

win.graph(4,8)
res %>%  mutate(method= lvls_reorder(method,c(3,1,2))) %>% 
  ggplot(aes(x = method, y = coef)) +
  geom_crossbar(aes(ymin=LCL, ymax=UCL), width=0.1) +
  geom_point() +
  facet_wrap(~ Species, ncol=1) +
  labs(x="Method", y="Treatment estimate", color="") +
  geom_hline(yintercept = 0, linetype=2) +
  theme_bw() +
  theme(strip.text.x = element_text(face="italic"))

#--------------------------------
#
# Sensitivty analysis
#

library(treatSens)

Data<- filter(Alldata, Species==spp1[1])
vars<- dplyr::select(Data, BA:MAS) # These covariates already standardised
var_names<- names(vars)

X<- as.matrix(vars)
Y<- as.vector(log(Data$Y))
Z<- as.vector(Data$Control)

nsim<- 100 # increase to 200 to reproduce figure from MS
sens <- treatSens(Y ~ Z+X, trt.family = binomial(link="probit"), nsim = nsim, 
                     spy.range = c(0,2), spz.range = c(-2, 2),grid.dim = c(17,9),
                     standardize = FALSE, verbose = F, weights="ATT")

sens$varnames<- c("Y","Z",var_names)

win.graph(8,8)
sensPlot(sens, data.line=F, txtlab=T,which.txtlab = c(2,4,5,6,8,13,9))
points(0.8, 0.5, pch=16)
       





