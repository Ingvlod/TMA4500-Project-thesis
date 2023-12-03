library(lme4)
library(ggplot2)
library(latex2exp)

#CASE A

#Function to simulate data for given number of simulations and regression coefficients where two of the covariates are correlated. 
simulate_data_correlated_covariate <- function(sim, b0, b1, b2, n = 100, m = 10){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- rep(rep(1:m, each = n), sim)
  obs_id <- rep(as.factor(rep(1:n,m)), sim)
  
  #fixed effect covariate
  x <- rnorm(n*m*sim, mean = 10, sd = 2)
  
  #fixed effect covariate correlated with x
  a <- 2+0.5*x + rnorm(n*m*sim) 
  
  #random cluster-specific effect
  u <- rep(rnorm(m*sim, 0, 0.5), each = n)
  
  #random error
  e <-  rnorm(n*m*sim, 0, 0.5)
  
  #response
  y <- b0 + b1*x + b2*a + u + e
  
  #Create data frame
  df <- data.frame(sim_id, cluster_id, obs_id, x, a, u, y)
  return(df)
}

#Function to estimate regression coefficients by using simulated data from 
#simulate_data_correlated_covariate() and omit one covariate
REmodel_estimate_betas_omitted <- function(sim, b0, b1, b2){
  beta0=c()
  beta1=c()
  
  #simulate data
  df = simulate_data_correlated_covariate(sim, b0, b1, b2)
  
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    REmodel_omitted = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
    beta0=append(beta0,fixef(REmodel_omitted)[1])
    beta1=append(beta1,fixef(REmodel_omitted)[2])
  }
  betas=data.frame(beta0, beta1)
  return(betas)
}

set.seed(222)
b0 = 5 
b1 = 0.25
b2 = 1

betas = REmodel_estimate_betas_omitted(1000, b0, b1, b2) 
beta0_ci=quantile(betas$beta0, prob=c(.025, .975))
beta1_ci=quantile(betas$beta1, prob=c(.025, .975))

hist1 <- ggplot(betas, aes(x=beta0))+geom_histogram(binwidth = 0.1, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(6,8)+
  geom_vline(aes(xintercept=mean(betas$beta0), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
  
hist1
hist2 <- ggplot(betas, aes(x=beta1))+geom_histogram(binwidth = 0.008, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(0.675,0.83)+
  geom_vline(aes(xintercept=mean(betas$beta1), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist2
mean(betas$beta0)
beta0_ci
mean(betas$beta1)
beta1_ci

#CASE B

#Function to simulate data for given number of simulations and regression coefficients where two of the covariates are correlated. 
simulate_data_correlated_covariate_constant_within <- function(sim, b0, b1, b2, n = 100, m = 10){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- rep(rep(1:m, each = n), sim)
  obs_id <- rep(as.factor(rep(1:n,m)), sim)
  
  #fixed effect covariate
  x <- rnorm(n*m*sim, mean = 10, sd = 2)
  
  
  #fixed effect covariate correlated with x
  x_mean_cluster <- c()
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  a <- 2 + 4*x_mean_cluster + rep(rnorm(m*sim), each = n)
  a_wrong <- a - 4*x_mean_cluster
  #random cluster-specific effect
  u <- rep(rnorm(m*sim, 0, 0.5), each = n)
  
  #random error
  e <-  rnorm(n*m*sim, 0, 0.5)
  
  #response
  y <- b0 + b1*x + b2*a + u + e
  
  #Create data frame
  df <- data.frame(sim_id, cluster_id, obs_id, x, x_mean_cluster, a, a_wrong, u, y)
  return(df)
}


#Function to estimate regression coefficients by using simulated data from 
#simulate_data_correlated_covariate_constant_within() and omit one covariate
REmodel_estimate_betas_omitted_constant <- function(sim, b0, b1, b2){
  beta0=c()
  beta1=c()
  
  #simulate data
  df = simulate_data_correlated_covariate_constant_within(sim, b0, b1, b2)
  
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    REmodel_omitted = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
    beta0=append(beta0,fixef(REmodel_omitted)[1])
    beta1=append(beta1,fixef(REmodel_omitted)[2])
  }
  betas=data.frame(beta0, beta1)
  return(betas)
}

set.seed(221)
b0 = 5 
b1 = 0.25
b2 = 1
betas = REmodel_estimate_betas_omitted_constant(100, b0, b1, b2) 

beta0_ci=quantile(betas$beta0, prob=c(.025, .975))
beta1_ci=quantile(betas$beta1, prob=c(.025, .975))

hist1 <- ggplot(betas, aes(x=beta0))+geom_histogram(binwidth = 0.2, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(betas$beta0), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))

hist1
hist2 <- ggplot(betas, aes(x=beta1))+geom_histogram(binwidth = 0.005, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(0.22,0.28)+
  geom_vline(aes(xintercept=mean(betas$beta1), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist2
mean(betas$beta0)
beta0_ci
mean(betas$beta1)
beta1_ci

#CASE C

#Function to estimate regression coefficients by using simulated data from 
#simulate_data_correlated_covariate_constant_within() and omit one covariate
REmodel_estimate_betas_omitted_constant_include_mean <- function(sim, b0, b1, b2){
  beta0=c()
  beta1=c()
  
  #simulate data
  df = simulate_data_correlated_covariate_constant_within(sim, b0, b1, b2)
  
  #for each simulation, estimate regression coefficients and store them
  for (i in 1:sim){
    REmodel_omitted = lmer(y~ x + x_mean_cluster + (1|cluster_id), data = df[df$sim_id == i,])
    beta0=append(beta0,fixef(REmodel_omitted)[1])
    beta1=append(beta1,fixef(REmodel_omitted)[2])
  }
  betas=data.frame(beta0, beta1)
  return(betas)
}

set.seed(221)
b0 = 5 
b1 = 0.25
b2 = 1
betas = REmodel_estimate_betas_omitted_constant_include_mean(100, b0, b1, b2) 

beta0_ci=quantile(betas$beta0, prob=c(.025, .975))
beta1_ci=quantile(betas$beta1, prob=c(.025, .975))

hist1 <- ggplot(betas, aes(x=beta0))+geom_histogram(binwidth = 5, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(-40,80)+
  geom_vline(aes(xintercept=mean(betas$beta0), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))

hist1
hist2 <- ggplot(betas, aes(x=beta1))+geom_histogram(binwidth = 0.005, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(betas$beta1), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist2
mean(betas$beta0)
beta0_ci
mean(betas$beta1)
beta1_ci
