setwd("C:/Users/.../Hurdle model")
library(readxl)
library(data.table)
library(dplyr)
library(readr)
library(ordinal)

# D_det<-read_csv("LT_length_data.csv") #<-cannot yet provide

### recast variables as factors ###
D_det$SID<-as.factor(D_det$SID)
D_det$Private<-as.factor(D_det$Private)
D_det$Class<-as.factor(D_det$Class)
D_det$Sex<-as.factor(D_det$Sex)
D_det$Condition<-as.factor(D_det$Condition)
D_det$Item<-as.factor(D_det$Item)
D_det$Soph_post<-as.ordered(D_det$Soph_post)

#** Frequentist model on observed subset of responses **#
m0_det<-clmm(as.ordered(Soph_post) ~ Theta + Condition + Sex + Private +
            (1|SID) + (1|Item),
            data=D_det %>% filter(Soph_post!="H"))
summary(m0_det)

#** Ordinal Hurdle Logit **#
source("oclhm_function_08022023.R")
library(brms)
library(rstan)
gc()
# PER DOCUMENTATION: ** non-coded/non-detected should be 0 **
D_det$Response_H<-if_else(D_det$Soph_post != 'H',
                          as.numeric(D_det$Soph_post) + 1, 0)

# run model #
m_h1<-oclhm(
  formula=as.ordered(Response_H) ~ Theta + Condition + Sex + Private + 
          (1|SID) + (1|Item), 
  data=D_det,
  intercept_prior = "normal(0.0, 2.0)",
  coeff_prior = "normal(0.0, 1.8)",
  sd_prior = "half_normal", 
  chains=3,cores=3,
  corr_variant = TRUE,
  warmup = 500, iter = 2500) #was 500, 2500

# get (poorly-formatted) Stan code that ran #
cat(get_stancode(m_h1$Stan_Obj))

# Posterior summary, ICCs, convergence #
m_h1$Post_Summary
m_h1$Intra_Class_Corr
m_h1$Convergence

# Fitted values #
fitted_1<-oclhm.fitted(m_h1$Stan_Obj)
str(fitted_1)
library(ggplot2)
# For ex, no difference in prob of detection between public and private, boys and girls
ggplot(data=fitted_1,aes(x=fit_detect)) + geom_histogram() +
  facet_wrap(~data_used.Private1)
ggplot(data=fitted_1,aes(x=fit_detect)) + geom_histogram() +
  facet_wrap(~data_used.SexM)

# get and compare LOOIC #
library(loo)
log_Lik1<-extract_log_lik(m_h1$Stan_Obj)
loo_1<-loo(log_Lik1, 
           r_eff = relative_eff(exp(log_Lik1),chain_id = rep(1:3,each=2000))) #<-need to match MCMC iterations post adaptation
#log_Lik2<-extract_log_lik(m_h2$Stan_Obj)
#loo_2<-loo(log_Lik2, 
#           r_eff = relative_eff(exp(log_Lik2),chain_id = rep(1:3,each=500)))
#loo_compare(loo_1,loo_2)
# PSIS diagnostic plot 
extract_log_lik(m_h1$Stan_Obj) %>% loo() %>% plot()

# Various MCMC diagnostics (trace plot) #
library(bayesplot)
out_array<-as.array(m_h1$Stan_Obj)
all_pars<-dimnames(out_array)$parameters
plot_pars<-all_pars[!(all_pars %like% "log_lik")&
                    !(all_pars %like% "mu_fit")&
                    !(all_pars %like% "detect_fit")&
                    !(all_pars %like% "Z_Mat")&
                    !(all_pars %like% "L_")&
                    !(all_pars %like% "R_")&
                    !(all_pars %in% "Corr_1[2,2]")&
                    !(all_pars %in% "Corr_2[2,2]")& 
                    !(all_pars %in% "Corr_1[1,1]")&
                    !(all_pars %in% "Corr_2[1,1]")& 
                    !(all_pars %in% "Corr_1[2,1]")&
                    !(all_pars %in% "Corr_2[2,1]")&  
                    !(all_pars %like% "lp")] #<-select parameters
# make trace plots, rhat plots, sampling efficiency plots
mcmc_trace(out_array,pars=plot_pars)
mcmc_rhat(rhat(m_h1$Stan_Obj,pars=plot_pars))+ yaxis_text(hjust = 1)
mcmc_neff(neff_ratio(m_h1$Stan_Obj,pars=plot_pars))+ yaxis_text(hjust = 1)

# Classroom effects (16 classrooms) #
R_Class<-extract(m_h2$Stan_Obj,pars="R_1")$R_1
dim(R_Class) #returns an array
Class_soph<-apply(R_Class,c(2,3),mean)[1,]
Class_det<-apply(R_Class,c(2,3),mean)[2,]
plot(Class_soph ~ Class_det, 
     ylim=c(-0.8,0.8),xlim=c(-0.8,0.8),
     main='Class random effects',
     xlab='Latent detection', 
     ylab='Latent sophistication after detection')
text(x=Class_det,y=Class_soph,levels(D_det$Class),pos=1,offset=0.25)
abline(v=0,h=0,lty=2,col='green',lwd=2)

# Item effects (26 Items) #
R_Item<-extract(m_h2$Stan_Obj,pars="R_2")$R_2
dim(R_Item) #returns an array
# 1 = ordinal response after detection
Item_soph<-apply(R_Item,c(2,3),mean)[1,]
# 2 = detection
Item_det<-apply(R_Item,c(2,3),mean)[2,]
plot(Item_soph ~ Item_det, 
     ylim=c(-4.5,4.5),xlim=c(-4.5,4.5),
     main='Item random effects',
     xlab='Latent detection', 
     ylab='Latent sophistication after detection')
text(x=Item_det,y=Item_soph,levels(D_det$Item),pos=1,offset=0.25)
abline(v=0,h=0,lty=2,col='green',lwd=2)

# Student effects (186 Students) #
# *AFTER* including: pre-Soph score, intervention, sex, public or private
R_Student<-extract(m_h2$Stan_Obj,pars="R_3")$R_3
dim(R_Student) #returns an array

# 1 = ordinal response after detection
Student_soph<-apply(R_Student,c(2,3),mean)[1,]
# 2 = detection
Student_det<-apply(R_Student,c(2,3),mean)[2,]

plot(Student_soph ~ Student_det, 
     xlim=c(-1.8,1.8),ylim=c(-1.8,1.8),
     main='Student effects',
     xlab='Latent detection', 
     ylab='Latent sophistication after detection')
abline(v=0,h=0,lty=2,col='green',lwd=2)
