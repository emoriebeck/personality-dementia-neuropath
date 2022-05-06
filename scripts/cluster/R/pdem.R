############## IDA METHOD 1B: Individual Participants Cluster Corrected Linear Regression ##############

setwd("/projects/p31385/dementia")

library(methods)
library(brms)        # bayesian models
library(broom.mixed) # summaries of models
library(rstan)       # bayes underpinnings
library(plyr)        # data wrangling
library(tidyverse)   # data wrangling

sessionInfo()

jobid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))+1
print(jobid)
args <- read.table("scripts/cluster/args/args.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

ipd_mega_mod_fun <- function(trait, outcome, mod, cov){
  ## load the data
  load(sprintf("data/mega-analysis/%s/%s-%s-%s.RData", cov, trait, outcome, mod))
  
  ## compiled Bayesian model to speed up processing and avoid crashing
  if(outcome %in% c("dementia", "hipSclerosis", "lewyBodyDis", "tdp43", "vsclrInfrcts", "vsclrMcrInfrcts")) load("results/bayes_sample_mod_binomial.RData") else load("results/bayes_sample_mod_continuous.RData")
  
  ## formula 
  cv <- c("age", "gender", "education", "smokes", "alcohol")
  #if (cov == "shared") cv <- cv
  if (cov == "butOne")  cv <- c(cv, "CC", "cognition")
  if (cov == "fully")    cv <- c(cv, "CC", "cognition", "BMI", "race", "SRhealth")
  if (cov == "standard") cv <- c("age", "gender", "education")
  rhs <- "p_value"
  rhs <- if(cov != "unadjusted") c(rhs, cv) else rhs
  if(mod != "none") rhs <- c(rhs, paste("p_value", mod, sep = "*"))
  re <- if(mod == "none") "(p_value | study)" else paste(paste("(p_value", mod, sep = " * "), "| study)")
  rhs <- paste(c(rhs, re), collapse = " + ")
  f <- paste("o_value ~ ", rhs, collapse = "")
  
  ## run the models & save
  m <- update(m
              , formula = f
              , newdata = d2
              , iter = 2000
              , warmup = 1000
              , cores = 4
  )
  save(m, file = sprintf("results/models/%s/%s_%s_%s.RData"
                         , cov, outcome, trait, mod))
  
  ## extract model terms and confidence intervals & save
  fx <- tidy(m, conf.int = T) %>%
    select(term, estimate, conf.low, conf.high)
  rx <- std_eff_fun(m)
  save(fx, rx, file = sprintf("results/summary/%s/%s_%s_%s.RData"
                              , cov, outcome, trait, mod))
  
  ## extract heterogeneity estimates
  het <- hetero_fun(m)
  save(het, file = sprintf("results/heterogeneity/%s/%s_%s_%s.RData"
                           , cov, outcome, trait, mod))
  
  ## simple effects for moderators
  if(mod != "none"){
    pred.fx <- fx_pred_fun(m, mod)
    pred.rx <- rx_pred_fun(m, mod)
    save(pred.fx, pred.rx, file = sprintf("results/predicted/%s/%s_%s_%s.RData"
                                          , cov, outcome, trait, mod))
  }
  
  ## clean up the local function environment
  rm(list = c("d", "f", "rhs", "m", "fx", "rx", "het"))
  gc()
}

std_eff_fun <- function(m){
  coef(m, probs = c(0.025, 0.975))[[1]] %>% array_tree(3) %>% 
    tibble(names = names(.), data = .) %>% 
    mutate(data = map(data, ~(.) %>% data.frame %>% 
                        rownames_to_column("study"))) %>% 
    unnest(data) %>% 
    select(names, study, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
}

hetero_fun <- function(m){
  args <- list(x = m, effects = "ran_pars", conf.int = T)
  do.call(tidy, args) %>%
    select(group, term, estimate, conf.low, conf.high) %>%
    separate(term, c("est", "term"), sep = "__") %>%
    mutate_at(vars(estimate:conf.high), ~ifelse(est == "sd", .^2, .)) %>%
    mutate(est = ifelse(est == "sd", "var", est))
}

fx_pred_fun <-function(m, moder){
  d <- m$data 
  d <- d %>% select(-o_value, -p_value, -study)
  cols <- colnames(d)
  md_cl <- class(d[,moder])
  if(any(sapply(d, class) == "numeric")){
    msd <- d %>%
      select_if(is.numeric) %>%
      pivot_longer(everything()
                   , names_to = "item"
                   , values_to = "value") %>%
      group_by(item) %>%
      summarize_at(vars(value), lst(mean, sd)) %>%
      ungroup()
  }
  if(any(sapply(d, class) == "factor")){
    fct_lev <- d %>% 
      select_if(is.factor) %>%
      summarize_all(~list(levels(.)))
  }
  d <- d %>% select(-one_of(moder))
  
  md_levs <- if(md_cl == "numeric"){
    if(moder %in% c("age")) {
      c(-10, 0, 10)
    } else if (moder %in% c("education")) {
      c(-5, 0, 5) 
    } else {
      with(msd, c(mean[item == moder] - sd[item == moder], mean[item == moder], mean[item == moder] + sd[item == moder]))
    }
  } else { 
    sort(unique(fct_lev[,moder][[1]]))
  }
  
  md_fac <- if(moder == "age") c("-10 yrs", "60", "+10 yrs") else if (moder == "education") c("-5 yrs", "12 years", "+5 yrs") else if(moder == "gender") c("Male", "Female") else c("-1 SD", "M", "+1 SD")
  
  mod_frame <- expand.grid(
    p_value = seq(0,10,.1)
    , modvalue = md_levs
    , stringsAsFactors = F
  ) %>% 
    mutate(mod_fac = factor(modvalue, levels = unique(modvalue), labels = md_fac)) %>%
    setNames(c("p_value", moder, "mod_fac")) 
  
  if(ncol(d) > 0){
    if(any(sapply(d, class) == "numeric")){
      mod_frame <- tibble(mod_frame, d %>% select_if(is.numeric) %>% summarize_all(mean))
    } 
    if(any(sapply(d, class) == "factor")){
      fcts <- colnames(d)[sapply(d, class) == "factor"]
      for(i in 1:length(fcts)){
        fct <- fcts[i]
        mod_frame <- crossing(
          mod_frame
          , levels(d[,fct])
        ) %>% setNames(c(colnames(mod_frame), fct))
      }
    }
  }
  
  pred.fx <- 
    bind_cols(
      mod_frame, 
      fitted(m
             , newdata = mod_frame
             , re_formula = NA) %>% data.frame
    ) %>%
    select(one_of(colnames(m$data)), mod_fac, pred = Estimate, lower = Q2.5, upper = Q97.5) %>%
    group_by(p_value, mod_fac) %>% 
    summarize_at(vars(pred, lower, upper), mean) %>%
    ungroup()
  
  rm(list = c("m", "mod_frame", "d", "md_levs"))
  gc()
  return(pred.fx)
}

rx_pred_fun <- function(m, moder){
  d <- m$data
  d <- d %>% select(-o_value, -p_value)
  cols <- colnames(d)
  md_cl <- class(d[,moder])
  if(any(sapply(d, class) == "numeric")){
    msd <- d %>%
      group_by(study) %>%
      select_if(is.numeric) %>%
      pivot_longer(-study
                   , names_to = "item"
                   , values_to = "value") %>%
      group_by(study, item) %>%
      summarize_at(vars(value), lst(mean, sd), na.rm = T) %>%
      ungroup() 
  }
  if(any(sapply(d, class) == "factor")){
    fct_lev <- d %>% 
      group_by(study) %>%
      select_if(is.factor) %>%
      summarize_all(~list((levels(.)))) %>%
      ungroup()
  }
  d <- d %>% select(-one_of(moder))
  
  md_fac <- if(moder == "age") c("-10 yrs", "60", "+10 yrs") else if (moder == "education") c("-5 yrs", "12 years", "+5 yrs") else if(moder == "gender") c("Male", "Female") else c("-1 SD", "M", "+1 SD")
  
  md_levs <- if(moder != "cognition"){
    crossing(
      study = unique(d$study)
      , modvalue = if(moder == "age") c(-10, 0, 10) else if (moder == "education") c(-5,0,5) else c(0,1)
    ) %>%
      mutate(mod_fac = factor(modvalue, levels = unique(modvalue), labels = md_fac)) %>%
      setNames(c("study", moder, "mod_fac"))
  } else {
    msd %>% 
      filter(item == moder) %>%
      mutate(lower = mean - sd, upper = mean + sd) %>%
      select(-sd) %>% 
      pivot_longer(cols = c(mean, lower, upper)
                   , names_to = "meas"
                   , values_to = "modvalue") %>%
      pivot_wider(names_from = "item", values_from = "modvalue") %>%
      select(study, one_of(moder)) %>%
      setNames(c("study", "modvalue")) %>%
      group_by(study) %>%
      mutate(mod_fac = factor(modvalue, levels = unique(modvalue), labels = md_fac)) %>%
      ungroup() %>%
      setNames(c("study", moder, "mod_fac"))
  }
  
  mod_frame <- crossing(
    p_value = seq(0,10,.1)
    , md_levs
  )
  
  if(ncol(d) > 0){
    if(any(sapply(d, class) == "numeric")){
      mod_frame <- d %>% 
        group_by(study) %>% 
        select_if(is.numeric) %>% 
        summarize_all(mean, na.rm = T) %>%
        ungroup() %>%
        full_join(mod_frame)
    } 
    if(any(sapply(d, class) == "factor")){
      fcts <- colnames(d)[sapply(d, class) == "factor"]
      for(i in 1:length(fcts)){
        fct <- fcts[i]
        mod_frame <- crossing(
          mod_frame
          , levels(d[,fct])
        ) %>% setNames(c(colnames(mod_frame), fct))
      }
    }
  }
  
  pred.rx <- 
    bind_cols(
      mod_frame, 
      fitted(m
             , newdata = mod_frame) %>% data.frame
    ) %>%
    select(one_of(colnames(m$data)), mod_fac, pred = Estimate, lower = Q2.5, upper = Q97.5) %>%
    group_by(p_value, mod_fac, study) %>% 
    summarize_at(vars(pred, lower, upper), mean) %>%
    ungroup()
  
  rm(list = c("m", "mod_frame", "d"))
  gc()
  return(pred.rx)
}

ipd_mega_mod_fun(args[,1], args[,2], args[,3], args[,4])