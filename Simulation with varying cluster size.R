library(lme4)
library(ggplot2)
library(ggridges)
library(latex2exp)
library(plm)
library(RColorBrewer)

#Function to simulate data for given number of simulations and 
#regression coefficients where two of the covariates are correlated. 
simulate_data<- function(sim, b0, b1, b2, n, m, constantwithin=FALSE){
  #simulation id, cluster id and observation within cluster id
  sim_id <- rep(1:sim, each = n*m)
  cluster_id <- rep(rep(1:m, each = n), sim)
  obs_id <- rep(rep(1:n,m), sim)

  x <- rnorm(n*m*sim, mean = 10, sd = 0.4)
  #Covariate correlated with x 
  x_mean_cluster <- c()
  for (i in seq(1, sim*n*m, n)){
    x_mean_cluster <- append(x_mean_cluster, rep(mean(x[i:(i+n-1)]), n))
  }
  if(constantwithin){
    #Covariate correlated with x and constant within clusters
    a <- 2 + 4*x_mean_cluster + rep(rnorm(m*sim), each = n)
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
  df <- data.frame(y, sim_id, cluster_id, obs_id, x, x_mean_cluster, a, e, u)
  return(df)
}

#Function to estimate regression coefficients by using simulated data from 
#simulate_data_correlated_covariate() and omit one covariate
model_and_estimate_betas_omitted <- function(sim, b0, b1, b2, n, m, constantwithin=FALSE, include_mean=FALSE){
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
      if (include_mean){
        REmodel_omitted = lmer(y~ x + x_mean_cluster + (1|cluster_id), data = df[df$sim_id == i,])
      }else{
        REmodel_omitted = lmer(y~ x + (1|cluster_id), data = df[df$sim_id == i,])
      }
      beta0=append(beta0,fixef(REmodel_omitted)[1])
      beta1=append(beta1,fixef(REmodel_omitted)[2])
      beta1FE=append(beta1FE,coefficients(FEmodel_omitted)[1])
    }
  }
  betas=data.frame(cluster_size,beta0,beta1,beta1FE)
  return(betas)
}


set.seed(235)
b0 = 5
b1 = 3
b2 = 1
n= 2
m= 100
sim = 1000

#Case A
betas_A=model_and_estimate_betas_omitted(sim, b0, b1, b2, n, m, constantwithin=FALSE, include_mean=FALSE)
betas_A$case=rep(as.factor("A"), sim)

mean(betas_A$beta0)
quantile(betas_A$beta0, prob=c(.025, .975))
mean(betas_A$beta1)
quantile(betas_A$beta1, prob=c(.025, .975))
mean(betas_A$beta1FE)
quantile(betas_A$beta1FE, prob=c(.025, .975))

#Case B
betas_B=model_and_estimate_betas_omitted(sim, b0, b1, b2, n, m, constantwithin=TRUE, include_mean=TRUE)
betas_B$case=rep(as.factor("B"), sim)

mean(betas_B$beta0)
quantile(betas_B$beta0, prob=c(.025, .975))
mean(betas_B$beta1)
quantile(betas_B$beta1, prob=c(.025, .975))
mean(betas_B$beta1FE)
quantile(betas_B$beta1FE, prob=c(.025, .975))

#Case C
betas_C=model_and_estimate_betas_omitted(sim, b0, b1, b2, n, m, constantwithin=TRUE, include_mean=FALSE)
betas_C$case=rep(as.factor("C"), sim)

mean(betas_C$beta0)
quantile(betas_C$beta0, prob=c(.025, .975))
mean(betas_C$beta1)
quantile(betas_C$beta1, prob=c(.025, .975))
mean(betas_C$beta1FE)
quantile(betas_C$beta1FE, prob=c(.025, .975))

betas_all_cases = rbind(betas_A, betas_B, betas_C)

histogram_beta0<-ggplot(betas_all_cases, aes(x = beta0, y = case, group = case, fill= case))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  #scale_fill_manual(name="", values=c("darkolivegreen2","brown1","cyan3"))+
  scale_fill_brewer(name="Case", palette="Set1")+
  geom_vline(aes(xintercept=b0, color="True beta"), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(RE estimated $\beta_0$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=8, y=7, label=TeX(r'($\beta_0=5$)'), size=5, color="#bf5252")
histogram_beta0

histogram_beta1<-ggplot(betas_all_cases, aes(x = beta1, y = case, group = case, fill= case))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  #scale_fill_manual(name="", values=c("darkolivegreen2","brown1","cyan3"))+
  scale_fill_brewer(name="Case", palette="Set1")+
  geom_vline(aes(xintercept=b1, color="True beta"), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(RE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=3.3, y=1.5, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
histogram_beta1

histogram_beta1FE<-ggplot(betas_all_cases, aes(x = beta1FE, y = case, group = case, fill= case))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  #scale_fill_manual(name="", values=c("darkolivegreen2","brown1","cyan3"))+
  scale_fill_brewer(name="Case", palette="Set1")+
  geom_vline(aes(xintercept=b1, color="True beta"), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  labs(x=TeX(r'(FE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=15),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=3.3, y=1.5, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
histogram_beta1FE




#Cluster size dependence in case C
set.seed(234)
b0 = 5
b1 = 3
b2 = 1
n=c(2,5,50)
m=100
sim=500

betas = model_and_estimate_betas_omitted(sim, b0, b1, b2, n, m, constantwithin = TRUE, include_mean=FALSE)
display.brewer.all()
hist_RE <- 
  ggplot(betas, aes(x=beta1, fill=as.factor(cluster_size)))+
  geom_histogram(binwidth = 0.03, color="darkslategray", alpha=0.5, position='identity')+
  geom_vline(aes(xintercept=b1), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  scale_fill_brewer(name="Cluster size, n", palette="Set2")+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_1$)'),y=TeX(r'(Count)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=3.05, y=200, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
hist_RE

hist_RE <- 
  ggplot(betas, aes(x=beta1, y = as.factor(cluster_size), group = cluster_size, fill=as.factor(cluster_size)))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  geom_vline(aes(xintercept=b1), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  scale_fill_brewer(name="Cluster size, n", palette="Set2")+
  theme(legend.position = c(0.88,0.80),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(RE estimated $\beta_1$)'),y=TeX(r'(Cluster size)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=2.95, y=7, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
hist_RE
hist_FE <- ggplot(betas, aes(x=beta1FE, fill=as.factor(cluster_size)))+
  geom_histogram(binwidth = 0.03, color="darkslategray", alpha=0.5, position='identity')+
  geom_vline(aes(xintercept=b1), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  scale_fill_brewer(name="Cluster size, n", palette="Set2")+
  theme(legend.position = c(0.88,0.85),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(FE estimated $\beta_1$)'),y=TeX(r'(Cluster size)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=3.05, y=200, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
hist_FE


hist_FE <- 
  ggplot(betas, aes(x=beta1FE, y = as.factor(cluster_size), group = cluster_size, fill=as.factor(cluster_size)))+
  geom_density_ridges2(stat="binline", bins=50, alpha=0.6, scale=5, color = "azure4")+
  geom_vline(aes(xintercept=b1), color= "#bf5252", linewidth=0.9, linetype="dashed", alpha=0.8)+
  scale_fill_brewer(name="Cluster size, n", palette="Set2")+
  theme(legend.position = c(0.88,0.80),legend.key = element_rect(colour = "white"))+
  labs(x=TeX(r'(FE estimated $\beta_1$)'),y=TeX(r'(Cluster size)'))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))+
  annotate(geom = "text", x=2.95, y=7, label=TeX(r'($\beta_1=3$)'), size=5, color="#bf5252")
hist_FE

