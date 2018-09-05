
library(tidyverse)
library(forcats)
library(survey)
library(optmatch)
library(parallel)


source("R/propensity_simulation_functions.R")
#=====================loop=================

cl<- makeCluster(25)  # Set to desired number of cores 

clusterEvalQ(cl, {
  library(dplyr)
  library(survey)
  library(optmatch)
})

clusterExport(cl, varlist=c("pscore.sim","sim_data","calc.pswts","fbuild",
                            "get.match.weights","ess"))
seed<- 99
clusterSetRNGStream(cl, seed)  # For reproducability

# Sample size
n_size <- c(50, 100, 200, 500, 1000)
t_coef<- 1 # Treatment coefficient
estimand<- "ATE"

scenario<- 2  # i.e. heterogeneous treatment effect

pint<- c(0.5, 0.2)  # Marginal probability of treatment

vars_out<- c(NULL) # covariates to exclude
#vars_out<- c("V1","V2")

results<- list()
ind<- 1
# Number of reps. Set to 100 for illustration.  Increase for more robust inferences:
reps <- 1000

for(i in 1:length(n_size)) {
  cat("Doing option ",n_size[i],"\n")

  for(j in 1:length(pint)) {
    
  out <- parLapply(cl, 1:reps, pscore.sim, j=j, nopt=n_size[i], 
                   v_opt=vars_out, scenario=scenario, 
                   t_coef=t_coef, pint=pint, est=estimand)
                    
  res<- do.call('rbind',out)
  results[[ind]]<- res
  ind<- ind+1

  }
}

stopCluster(cl)

results<- do.call('rbind', results)


#=================================================
#
# Some plots
#
#=================================================
library(ggplot2)
library(RColorBrewer)
palette <- brewer.pal(9, "Set1")[c(1:5)]
if(estimand=="ATE") {
  ps.nam<- c("Trt","Reg","IPW","IPW_dr","FM")
} else {
  ps.nam<- c("Trt","PM","IPW","IPW_dr","FM")
}

names(palette) <- ps.nam

new_levels<- c("Treatment probability 50%","Treatment probability 25%")

# Bias


win.graph(10,3)
 results %>% mutate(vars_in = lvls_revalue(factor(vars_in),new_levels),
         Bias=(TE - target), 
         method = fct_relevel(method, ps.nam)) %>%
  ggplot(aes(x = as.factor(n), color = method, y = Bias)) +
  geom_boxplot(outlier.colour = NA, position ="dodge") +
  coord_cartesian(ylim = c(-3, 3)) +
  geom_hline(yintercept=0, linetype=1) +
  facet_wrap( ~ vars_in, ncol=2) +
  scale_colour_manual(values = palette, guide = guide_legend(title="Method")) +
  labs(x = "Sample size", y = "Bias", colour = "") +
  theme(legend.title = element_text(colour="blue", size=16, face="bold")) +
  theme_bw() 
 
 # RMSE

 win.graph(10,3)
 results %>% mutate(vars_in = lvls_revalue(factor(vars_in),new_levels),
          RMSE=sqrt((TE - target)^2), 
          method = fct_relevel(method, ps.nam)) %>%
   ggplot(aes(x = as.factor(n), color = method, y = RMSE)) +
   geom_boxplot(outlier.colour = NA) +
   coord_cartesian(ylim = c(0, 3)) +
   facet_wrap( ~ vars_in, ncol=2) +
   scale_colour_manual(values = palette, guide = guide_legend(title="Method")) +
   labs(x = "Sample size", y = "RMSE", colour = "") +
   theme(legend.title = element_text(colour="blue", size=16, face="bold")) +
   theme_bw() 
 
# Coverage

win.graph(8,6)
results %>% mutate(vars_in = lvls_revalue(factor(vars_in), new_levels)) %>%
  group_by(method, vars_in, n) %>% summarise(Cov=mean(Cov, na.rm=T)) %>% 
  ungroup(method) %>% 
  mutate(method = fct_relevel(method, ps.nam)) %>% 
  ggplot(aes(x = factor(n), shape = vars_in, y = Cov)) +
  geom_jitter(size=2, width=0.05, height=0) +
  geom_hline(yintercept = 0.95, linetype=2) +
  facet_wrap(~ method, ncol=3) +
  scale_shape_manual(values = c(1,2,3,0), guide = guide_legend(title="Model")) +
  labs(x = "Sample size", y = "Proportional Coverage of 95% CI", colour = "") +
  theme_bw() + theme(legend.position=c(0.8,0.15)) 




