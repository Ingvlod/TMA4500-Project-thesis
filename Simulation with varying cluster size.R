library(lme4)
library(ggplot2)
library(latex2exp)
library(plm)

#Function to simulate data for given number of simulations and regression coefficients where two of the covariates are correlated. 
simulate_data<- function(sim, b0, b1, b2, n, m, constantwithin=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- rep(rep(1:m, each = n), sim)
  obs_id <- rep(rep(1:n,m), sim)
  
  #Covariate
  mean_x = runif(m, min = 5, max = 10)
  x <- rnorm(n*m*sim, mean = rep(rep(mean_x, each=n), sim), sd = 0.4)
  
  #x <- rnorm(n*m*sim, mean = 0, sd = 0.5)
  #Covariate correlated with x 

  if(constantwithin){
    #Covariate correlated with x and constant wihtin clusters
    x_mean_cluster <- c()
    for (i in seq(1, sim*n*m, n)){
      x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
    }
    a <- 2 + 0.5*x_mean_cluster + rep(rnorm(m*sim), each = n)
    }else{
    #Covariate correlated with x 
    a <- 2 + 4*x + rnorm(n*m*sim)
    }

  #Cluster-specific effect
  u <- rep(rnorm(m*sim, 0, 0.5), each = n)
  
  #Random error
  e <-  rnorm(n*m*sim, 0, 0.5)
  
  #Response
  y <- b0 + b1*x + b2*a + u + e
  
  #Create data frame
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, a, e)
  return(df)
}

#Function to estimate regression coefficients by using simulated data from 
#simulate_data_correlated_covariate() and omit one covariate
model_and_estimate_betas_omitted <- function(sim, b0, b1, b2, n, m, constantwithin=FALSE){
  cluster_size=rep(n, each=sim)
  beta0=c()
  beta1=c()
  beta1FE=c()
  #For each cluster size
  for (i in n){
    df = simulate_data(sim, b0, b1, b2, i, m, constantwithin)
    pdata = pdata.frame(df, index="cluster_id")
    #for each simulation, estimate regression coefficients and store them
    for (i in 1:sim){
      FEmodel_omitted = plm(y~ x, data = pdata[pdata$sim_id == i,], model = "within")
      REmodel_omitted = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
      beta0=append(beta0,fixef(REmodel_omitted)[1])
      beta1=append(beta1,fixef(REmodel_omitted)[2])
      beta1FE=append(beta1FE,coefficients(FEmodel_omitted)[1])
    }
  }
  betas=data.frame(cluster_size,beta0,beta1,beta1FE)
  return(betas)
}
b0=5
b1=0.25
b2=1


#Case A
set.seed(222)
b0 = 5 
b1 = 0.25
b2 = 1
n=100
m=10
sim=100
betas = model_and_estimate_betas_omitted(sim, b0, b1, b2,n ,m) 
beta0_ci=quantile(betas$beta0, prob=c(.025, .975))
beta1_ci=quantile(betas$beta1, prob=c(.025, .975))
beta1FE_ci=quantile(betas$beta1FE, prob=c(.025, .975))

hist1 <- ggplot(betas, aes(x=beta0))+geom_histogram(binwidth = 0.1, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(6,8)+
  geom_vline(aes(xintercept=mean(beta0), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
  
hist1
hist2 <- ggplot(betas, aes(x=beta1))+geom_histogram(binwidth = 0.008, fill="cyan3", color="darkcyan", alpha=0.3)+
  xlim(0.675,0.83)+
  geom_vline(aes(xintercept=mean(beta1), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist2
hist3 <- ggplot(betas, aes(x=beta1FE))+geom_histogram(binwidth = 0.008, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(beta1FE), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1FE_ci[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1FE_ci[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(FE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist3
mean(betas$beta0)
beta0_ci
mean(betas$beta1)
beta1_ci

#CASE B


set.seed(234)
b0 = 5
b1 = 3
b2 = 1
n=c(2,5,50)
m=100
sim=100

betas = model_and_estimate_betas_omitted(sim, b0, b1, b2, n, m, constantwithin = TRUE)
hist <- ggplot(betas, aes(x=beta1, fill=as.factor(cluster_size)))+
  geom_histogram(binwidth = 0.03, color="darkslategray", alpha=0.5, position='identity')+
  geom_vline(aes(xintercept=b1), color= "coral", linewidth=1.3, linetype="solid")+
  scale_fill_manual(name="Cluster size, n", values=c("brown1","cyan3","darkolivegreen2"))+ 
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))
hist

betas[betas$cluster_size == n[1],]
beta0_ci_1=quantile(betas[betas$cluster_size == n[1],]$beta0, prob=c(.025, .975))
beta1_ci_1=quantile(betas[betas$cluster_size == n[1],]$beta1, prob=c(.025, .975))
beta1FE_ci_1=quantile(betas[betas$cluster_size == n[1],]$beta1FE, prob=c(.025, .975))


hist1 <- ggplot(betas[betas$cluster_size == n[1],], aes(x=beta0))+geom_histogram(binwidth = 0.2, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(betas[betas$cluster_size == n[1],]$beta0), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci_1[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta0_ci_1[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))

hist1
hist2 <- ggplot(betas[betas$cluster_size == n[1],], aes(x=beta1))+geom_histogram(binwidth = 0.005, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(betas[betas$cluster_size == n[1],]$beta1), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci_1[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1_ci_1[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(Estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist2
hist3 <- ggplot(betas[betas$cluster_size == n[1],], aes(x=beta1FE))+geom_histogram(binwidth = 0.008, fill="cyan3", color="darkcyan", alpha=0.3)+
  geom_vline(aes(xintercept=mean(betas[betas$cluster_size == n[1],]$beta1FE), color="mean", linetype="mean"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1FE_ci_1[1], color="CI", linetype="CI"), linewidth=1.3)+
  geom_vline(aes(xintercept=beta1FE_ci_1[2], color="CI", linetype="CI"), linewidth=1.3)+
  scale_color_manual(name = "Values", values = c("mean" = "coral", "CI"="coral"), labels=c("95% quantile", "Mean"))+
  scale_linetype_manual(name = "Values", values = c("mean" = "solid", "CI"="dotted"), labels=c("95% quantile", "Mean"))+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(FE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  guides(color = guide_legend(override.aes = list(size = 8)))
hist3
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
b2 = 4
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
