oclhm <- function(formula, data, 
                 chains = 3, cores = chains, iter = 2000, warmup = floor(iter/2), 
                 intercept_prior = "normal(0, 2.5)", 
                 coeff_prior = "normal(0, 2.0)", 
                 sd_prior = "half_normal",
                 stan_control = list(adapt_delta = 0.8, max_treedepth=10),
                 corr_variant = TRUE){

# Make fixed effects labels for sophistication & detection#
# Recast as brms formula
formula_Model <- bf(formula = as.formula(formula),family = cumulative())
fe_labels <- brmsterms(formula_Model$formula)$dpars$mu$fe %>%
             model.matrix(data=data) %>% as.data.frame() %>%
             select(-1) %>% colnames() #%>% 
fe_det_labels<-paste(fe_labels, "detection", sep= "_")

# Make random effects labels #
re_labels <- brmsterms(formula_Model$formula)$dpars$mu$re[1] 
fix.eff <- fe_labels
rand.eff <- re_labels

#define random effect block of Stan CODE
r.eff.block <- NA
for(i in 1:nrow(rand.eff)){
  r.eff.block[i] <- paste(
    "int<lower=1> N_", i, "; ", 
    "int<lower=1> M_", i, "; ",      
    "int<lower=1> J_", i, "[N]; ",
    "vector[N] Z_", i, "_1; ", 
    sep = "",
    collapse = ""
  )
}
r.eff.block <- paste(r.eff.block, collapse = " ")

#define fixed effect block of Stan CODE
fixed.eff.block <- "
    data {                      
      int<lower=1> N;          
      int Y[N];               
      int<lower=1> K;        
      matrix[N, K] X;
    "

# define data block of Stan CODE & define number of Thresholds
data_block <- paste(fixed.eff.block, " int<lower=2> nthres; ", 
                    r.eff.block, "int prior_only;}")

# define parameters block of Stan CODE
parms.bl <- NA
#** independent effects variant **#
for(i in 1:nrow(rand.eff)){
  parms.bl[i] <- paste(  
    "vector<lower=0> [M_", i, "]", "sd_", i, ";",         
    "vector[N_",i, "]", "z_",i,"[M_", i,"];",              
    "vector<lower=0>[M_", i,"]", "sd_", i, "_d;", 
    "vector[N_",i,"]"," z_",i, "_d[M_", i,"];",   
    sep = ""
  )
}
#** correlated effects variant **#
#NOTE: Z_Mat_ contains 2 rows: 1 for sophistication; 2 for detection; 2 hard-coded b/c always 2 linear predictors
if(corr_variant){
  for(i in 1:nrow(rand.eff)){
    parms.bl[i] <- paste(  
      "vector<lower=0> [M_", i, "]", "sd_", i, ";",         
      "vector<lower=0>[M_", i,"]", "sd_", i, "_d;", 
      
      "matrix[2, N_",i, "]", "Z_Mat_",i,";",              
      "cholesky_factor_corr[2]","L_",i, ";",
      sep = ""
    )
  }
}

parms.bl <- paste(parms.bl, collapse = " ")
parms.bl <- paste(parms.bl, 
                  "vector[K] b;                
                   ordered[nthres] Intercept; 
                   real Intercept_detect; 
                   vector[K] b_detect;"
)
parms.bl <- paste("parameters {", parms.bl, "}")

# define transformed parameters block of Stan CODE & initialize vectors
transf.parms.bl.v <- NA
#** independent effects variant **#
for(i in 1:nrow(rand.eff)){
  transf.parms.bl.v[i] <- paste(
    "vector[N_", i, "] r_", i, "_1;",
    "vector[N_", i," ] r_", i, "_d;",     
    sep = ""
  )
}
#** correlated effects variant **#
#NOTE: the matrix R_ contains 2 rows: 1 goes into sophistication; 2 goes into detection
if(corr_variant){
  for(i in 1:nrow(rand.eff)){
  transf.parms.bl.v[i] <- paste(
    "matrix [2, N_", i, "] R_", i, ";",
    sep = ""
  )
  }
}
transf.parms.bl.v <- paste(transf.parms.bl.v, collapse = " ")

# transformed parameters block of Stan CODE continued & Cholesky transform of ran effs
transf.parms.bl.r <- NA 
#** independent effects variant **#
for(i in 1:nrow(rand.eff)){ 
  transf.parms.bl.r[i] <- paste(
    "r_", i, "_1 = sd_", i, "[1]*(z_", i, "[1]);",    # ordinal outcome random effs
    "r_", i, "_d = sd_", i,"_d[1]*(z_", i, "_d[1]);", # detection random effs
    sep = ""
  )
}
#** correlated effects variant **#
#NOTE: the matrix R_ contains 2 rows: 1 goes into sophistication; 2 goes into detection
if(corr_variant){
  for(i in 1:nrow(rand.eff)){ 
  transf.parms.bl.r[i] <- paste(
    "R_", i, " = diag_pre_multiply(append_row(sd_",i,",sd_",i,"_d),L_", i, ")*Z_Mat_",i,";",    
    sep = ""
  )
  }
}  
transf.parms.bl.r <- paste(transf.parms.bl.r, collapse = " ")
transf.parms.bl <- paste("transformed parameters {", transf.parms.bl.v,transf.parms.bl.r, "}")

#*********************************************************************************************************#
# begin model block of Stan code (initialize ordinal var mean, detection mean, disc par removed here)
# start with fixed effects only
model.predictors <- "
    vector[N] mu = X*b;
    vector[N] detect = Intercept_detect + X*b_detect; 
    "
# initialize mean formulas to be updated with random effects
lin.pred.re <- NA 
mu.form <- NA 
detect.form <- NA 
# creates the code based on random effect structure to add to ordinal and detection means
#** independent effects variant **#
for(i in 1:nrow(rand.eff)){
  mu.form[i] = paste("r_", i, "_1[J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
  detect.form[i] = paste("r_", i, "_d[J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
}
#** correlated effects variant **#
if(corr_variant){
  for(i in 1:nrow(rand.eff)){
      mu.form[i] =     paste("R_", i, "[1, J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
      detect.form[i] = paste("R_", i, "[2, J_", i,"[n]] * Z_", i, "_1[n] +", sep = "")
  }
}  

# process characters to make a linear formula for random effs to add to ordinal mean
mu.form[nrow(rand.eff)] <- gsub("+", ";", mu.form[nrow(rand.eff)],fixed = T)
mu.form <- paste(mu.form, collapse = " ")
detect.form[nrow(rand.eff)] <- gsub("+", ";", detect.form[nrow(rand.eff)],fixed = T)
detect.form <- paste(detect.form, collapse = " ")

# make complete set of linear predictors with fixed and random effects
# uses mu and detect lin preds from beginning of Model block
lin.pred.re <- paste(
  "for (n in 1:N) {
       mu[n] +=", mu.form, 
  "detect[n] +=", detect.form, "}",
  sep = "")

# inverse-logit transform detection probs from linear scale to probabilities scale
detect.loop <- "for (n in 1:N) { 
        detect[n] = inv_logit(detect[n]);
      }"

#******************************************************************************************************#
# begin priors Stan code block (all priors assigned here)
# note that both detection and ordinal variable on logit scale, same priors appropriate
coeff.priors <- paste(
  "b_detect ~", coeff_prior, ";",                #<-detection coeffs
  "Intercept_detect ~", intercept_prior, ";",    #<-detection Intercept
  "b ~",  coeff_prior, ";",                      #<-ordinal variable coeffs
  "Intercept ~", intercept_prior, ";", sep = ""  #<-ordinal variable thresholds
)

# random effect SDs for ordinal variable and detection prob
re.priors <- NA 
if (sd_prior == "half_normal") {
  for (i in 1:nrow(rand.eff)) {
    re.priors[i] <- paste(
      "target += normal_lpdf(sd_",
      i,
      " | 0, 2.5) - normal_lccdf(0|0, 2.5);",   #<-ordinal
      "target += normal_lpdf(sd_",
      i,
      "_d | 0, 2.5) - normal_lccdf(0|0, 2.5);", #<-detection (note _d syntax)
      sep = ""
    )
  }}
re.priors <- paste(re.priors, collapse = " ")

# indepedent random effect Z-scores to be transformed via Cholesky
z.vars <- NA 
#** independent effects variant **#
for(i in 1:nrow(rand.eff)){
  z.vars[i] <- paste(
    "target += std_normal_lpdf(z_", i,"[1]);", 
    "target += std_normal_lpdf(z_", i, "_d[1]);", 
    sep = ""
  )
}
#** correlated effects variant **#
if(corr_variant){
  for(i in 1:nrow(rand.eff)){
    z.vars[i] <- paste(
      "L_", i," ~ lkj_corr_cholesky(2);",
      "to_vector(Z_Mat_",i,") ~ std_normal();",
      sep = ""
    )
  }
}
z.vars <- paste(z.vars, collapse = " ")

#*********************************************************************************************************#
# define model likelihood block of Stan CODE
# NOTE: Y=0 indicates "hurdle" condition ie non-detectable or non-codable response
likelihood.bl <- "
  for (n in 1:N) {
    if(Y[n] == 0)                     
      target += log1m(detect[n]); 
    else   
      target += log(detect[n]) + ordered_logistic_lpmf(Y[n]|mu[n], Intercept);
  }
"

#**********************************************************************************************************#
# write out the entire model Stan code prior to generated quantities
model_bl <- paste(
  "model { ", 
   model.predictors, 
   lin.pred.re, 
   detect.loop,  
   coeff.priors, 
   re.priors, 
   z.vars, 
   likelihood.bl, 
   "}"
)

#****************************************************************************************************************#
# generated quantities block Stan CODE
# fitted ordinal mean (mu_fit) and fitted detection probabilities (detect_fit)
# compute log-likelihood for diagnostics and Info Criteria
#** independent effects variant **#
gen_quan_bl.p <- paste(
  "generated quantities {
    vector[N] mu_fit = X*b;
    vector[N] detect_fit = Intercept_detect + X*b_detect;
    real log_lik[N];   
    for(n in 1:N){
          mu_fit[n] +=", mu.form, 
     "detect_fit[n] +=", detect.form, 
  "if(Y[n] == 0) 
      log_lik[n] = log1m(inv_logit(detect_fit[n])); 
      else   
      log_lik[n] = log(inv_logit(detect_fit[n])) + ordered_logistic_lpmf(Y[n] | mu_fit[n], Intercept);
      };
    }"
)

#** correlated effects variant **#
if(corr_variant){
gen_quan.cor<-NA  
 for(i in 1:nrow(rand.eff)){
   #make all the between-random effect correlation matrices
   gen_quan.cor[i] <- paste(
     "matrix[2,2] Corr_",i," = multiply_lower_tri_self_transpose(L_",i,");",
     sep=""
    )  
}
gen_quan_bl.p <- paste(
  "generated quantities {",
   paste(gen_quan.cor,collapse = ""),
   "vector[N] mu_fit = X*b;
    vector[N] detect_fit = Intercept_detect + X*b_detect;
    real log_lik[N];   
    for(n in 1:N){
          mu_fit[n] +=", mu.form, 
  "detect_fit[n] +=", detect.form, 
  "if(Y[n] == 0) 
      log_lik[n] = log1m(inv_logit(detect_fit[n])); 
      else   
      log_lik[n] = log(inv_logit(detect_fit[n])) + ordered_logistic_lpmf(Y[n] | mu_fit[n], Intercept);
      };
    }"
)
gen_quan_bl.p <- paste(gen_quan_bl.p, collapse = " ")
}

#*********************************************************************************************#
# write out the ENTIRE MODEL Stan CODE including generated quantities
stan_model_code <- paste(
  data_block,
  parms.bl, 
  transf.parms.bl,
  model_bl,
  gen_quan_bl.p,
  "     ")
stan_model_code <- paste(stan_model_code, collapse = "\n")

# make Standata
standat_Model <- make_standata(
  formula = formula_Model, 
  data = data,
  family = cumulative()
)
# make hurdle condition (where outcome = 0)
standat_Model$Y <- standat_Model$Y - 1
# find thresholds based on unique number of ordinal outcomes
standat_Model$nthres <- sum(unique(standat_Model$Y) != 0) - 1 

# make neat threshold labels
threshold_lab <- NA
for(i in 1:max(standat_Model$Y-1)){
  threshold_lab[i] <- paste(i, "|", i+1, sep = "")
}

# random effect SD labels for PARS to sample (must match names above)
sd.r.d <- NA 
for(i in 1:nrow(rand.eff)){
  sd.r.d[i] <- paste(
    'sd_', i, " ",
    'sd_', i, '_d', "",
    sep = ""
  )
}
sd.r.d <- paste(sd.r.d, collapse = " ")
sd.r.d <- strsplit(sd.r.d, " ")[[1]]

# random effect vector labels for PARS to sample (must match names above)
r.d <- NA 
for(i in 1:nrow(rand.eff)){
  r.d[i] <- paste(
    'r_', i,"_1", " ",
    'r_', i, '_d', "",
    sep = ""
  )
}
r.d <- paste(r.d, collapse = " ")
r.d <- strsplit(r.d, " ")[[1]]

# correlation matrices
cor.m <- NA 
for(i in 1:nrow(rand.eff)){
  cor.m[i] <- paste(
    "Corr_", i,"[1,2]",
#    "Corr_", i,
    sep = ""
  )
}
cor.m <- paste(cor.m, collapse = " ")
cor.m <- strsplit(cor.m, " ")[[1]]

# set which parameters to sample using Stan
pars <- c(
  'Intercept','b',
  'Intercept_detect','b_detect',
  'log_lik', 'mu_fit', 'detect_fit',
  cor.m,
#  r.d, #<-can comment in if needed
  sd.r.d
)
gc()

#****************************************************************************************************#
# *CREATE* and *RUN* Stan MODEL OBJECT
Stan_Model <- stan(
  model_code = stan_model_code,
  data = standat_Model,
#  pars = pars,
  chains = chains,
  warmup = warmup,
  iter = iter,
  cores = cores,
  control = stan_control
) %>% base::suppressWarnings()
# all_samps<-rstan::extract(Stan_Model) #<-full output 
# rstan::get_stancode(Stan_Model) #<-code used


#***********************************#
# summary function to process output 
stan_sum <- function(model, pars){
  
  # gather parameters and generate initial output #
  model_sum <- summary(
    model, 
    pars = pars[!pars %in% c('log_lik', 'mu_fit', 'detect_fit', r.d)],
    probs = c(0.025, 0.5, 0.975)
    )$summary
  sd_vec<-as.data.frame(rstan::extract(Stan_Model, unlist(sd.r.d)))
  
  # re-label output
  to<- length(row.names(model_sum))
  from<- to - length(sd.r.d) - nrow(rand.eff)
  row.names(model_sum)[-(from+1):-to] <- c(threshold_lab, fe_labels, 
                                           "Intercept_detection", fe_det_labels)
  # make random effect SD labels #
  for(j in 1:nrow(rand.eff)){
    sd.r.d[sd.r.d %like% paste("_",j,sep="")]<-c(paste("SD_",rand.eff$group[j],sep=""),
                                                 paste("SD_",rand.eff$group[j],"_d",sep=""))
  }
  row.names(model_sum)[(from+nrow(rand.eff)+1):to]<-sd.r.d
  
  if(corr_variant){
  corr.labs<-NA
  # make correlation labels #
  for(j in 1:nrow(rand.eff)){
    corr.labs[j]<-paste("Ord_Det_Cor_",rand.eff$group[j],sep="")
  }
  row.names(model_sum)[which(row.names(model_sum) %>% like("Corr"))]<-corr.labs
  }  
  
  # pick off Rhats > 1.05 and indicate possibe non-covergence to user
  rhat <- data.frame(rhat=data.frame(model_sum)$Rhat)
  rhat$parm_num <- row.names(rhat)  
  rhat <- rhat %>% dplyr::filter(rhat > 1.05)
  conv_fun <- function(x){
    if(nrow(rhat) == 0){
      converge <- paste("Convergence: all Rhats less than 1.05")
    }
    else{
      converge <- paste("Possible non-convergence in ", 
                        row.names(model_sum)[as.numeric(rhat$parm_num)], sep = "")
    }}
  
  # compute ICCs using Logit ICC formula and raneff SD estimates
  logit_icc<-function(s){ s^2 / (s^2 + (pi^2/3)) }  
  icc_sum <- sd_vec %>%
             logit_icc %>% 
             apply(2, function(x) quantile(x, probs = c(0.025, 0.50, 0.975)))
  # assign labels to ICCs
  colnames(icc_sum)<-row.names(model_sum)[(from+nrow(rand.eff)+1):to] %>% substring(4)

  data_used<- data.frame(standat_Model$X,Y=standat_Model$Y)
  output <- list(
    Post_Summary = round(model_sum, 3), 
    Intra_Class_Corr = round(icc_sum, 3),
    Convergence = conv_fun(rhat),
    data_used = data_used,
    Stan_Obj = Stan_Model
  )
  return(output)
}

invisible(stan_sum(Stan_Model, pars = pars))

} #<- END estimation function


########################################################################################################
##** OCLHM FITTED FUNCTION BELOW HERE **##
oclhm.fitted = function(Stan_Model){

  # extract necessary estimates from Stan model
  inv_logit_fn = function(x){exp(x)/(1 + exp(x))}
  detect_fit <- rstan::extract(Stan_Model)$detect_fit
  mu_fit <- rstan::extract(Stan_Model)$mu_fit
  thres <- rstan::extract(Stan_Model)$Intercept
  p_detect <- inv_logit_fn(detect_fit)   
  
  # matrix of fitted strategy probabilities assuming ordinal logit
  p.mat <- matrix(nrow = ncol(p_detect), ncol = ncol(thres)+1)  # Unconditional Probs (i.e. not conditional on detect)
  # P(Y = first category)
  p.mat[, 1] <- colMeans(p_detect*inv_logit_fn(thres[, 1] - mu_fit))  
  # P(Y = last category)
  p.mat[, ncol(thres)+1] <- colMeans(p_detect*
                                    (1 - inv_logit_fn(thres[, ncol(thres)] - mu_fit)))
  # P(Y = all middle categories)
  for(i in 2:ncol(thres)){
    p.mat[, i] = colMeans(p_detect*(inv_logit_fn(thres[, i] - mu_fit) - 
                                    inv_logit_fn(thres[, i-1] - mu_fit)))
  }
  ### post-process to assign labels (p.mat is N x number categories)
  p.mat <- as.data.frame(p.mat) 
  for(i in 1:ncol(p.mat)){ 
    names(p.mat)[i] = paste("Y = ", i, sep = "")
  }
  
  ### Conditional probabilities
  cp.mat <- matrix(nrow = ncol(p_detect), ncol = ncol(thres)+1)   # Conditional Probs (i.e. conditional on detect)
  ### First
  cp.mat[, 1] <- colMeans(inv_logit_fn(thres[, 1] - mu_fit))  
  ### Last
  cp.mat[, ncol(thres)+1] <- colMeans((1 - inv_logit_fn(thres[, ncol(thres)] - mu_fit)))
  ### All middle categories
  for(i in 2:ncol(thres)){
    cp.mat[, i] = colMeans((inv_logit_fn(thres[, i] - mu_fit) - inv_logit_fn(thres[, i-1] - mu_fit)))
  }
  cp.mat <- as.data.frame(cp.mat)
  ### Post process to assign labels 
  for(i in 1:ncol(cp.mat)){ 
    names(cp.mat)[i] = paste("Y Cond. = ", i, sep = "")
  }
  
  ### Collect unconditional and conditional probabilities into output dataset
  out.mat <- data.frame(
    p.mat,
    fit_detect = colMeans(p_detect),
    cp.mat, 
    data_used = data.frame(standat_Model$X,Y=standat_Model$Y),
    check.names = FALSE
  )
  return(out.mat)
}

