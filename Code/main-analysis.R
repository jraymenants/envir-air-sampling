# Article title: Indoor air surveillance and factors associated with respiratory pathogen detection in community settings in Belgium
# Authors: Joren Raymenants, Caspar Geenen, Lore Budts, Jonathan Thibaut, Marijn Thijssen, Hannelore De Mulder, Sarah Gorissen, Bastiaan Craessaerts, Lies Laenen, Kurt Beuselinck, Sien Ombelet, Els Keyaerts, Emmanuel André

# This script consists of the following sections:
# - Load main dataframe: run this section before running any of the other sections
# - Model outcome A (probability of detection): generalised estimating equations to correct for sample id
# - Model outcome A (probability of detection): mixed effect model to correct for sample id
# - Model outcome A (probability of detection): using regular logit regression
# - Model outcome B: Linear model for CT value with backward elimination
# - Model outcome B: mixed effect model for CT value to correct for sample id
# - Model outcome A (probability of detection): GLM for each pathogen separately
# - Model outcome B: GLM for each pathogen separately, using the retained variables from model 3

library(tidyverse) # general data handling
library(readxl) # read Excel files
library(lubridate) # date handling
library(geepack) # generalised estimating equations
library(lme4) # mixed effect modelling
library(MESS) # for drop1 function on GEE
source("Code/data_cleaning.R") # import functions for data cleaning

# Choose data file location
data_file <- "Input/source_data_v2.xlsx"

########################################################################
# Load main dataframe: run this before running any other section below #
########################################################################

load_source_data <- function(file) {
  
  # Read data file, group pathogens, rename pathogens and remove pathogens where less than 10 samples were positive
  data <- read_excel(data_file, sheet=1, na = "NA") %>%
    recode_pathogens() %>%
    filter(pathogen!="SARS-CoV-2 (respiratory panel)") %>% # Replace "respiratory panel" with "TaqPath" to use respiratory panel results
    group_corona_parainf() %>%
    
    group_by(pathogen) %>%
    filter(sum(detected)>=10) %>% # remove pathogens with less than 10 positive results
    ungroup()
  
  # set purifier to false if data missing
  # and drop observations with missing co2 or humidity values
  # also drop observations with missing inferred values for temperature, ventilation, mask use or vocalisation
  # "babbelen" refers to vocalisation
  data$purifier <- replace_na(data$purifier,F)  # no value means no purifier
  sapply(data, function(x) sum(is.na(x))) # which data is missing?
  data <- data %>% drop_na(co2_avg,humidity_avg,temperature_inferred, ventilation_inferred, mask_inferred, babbelen_inferred)
  
  # Calculate attendee density and rescale variables.
  extra_vars <- data %>%
    rowwise() %>%
    mutate(
      attendee_density = n_attendees_mean/volume,
    ) %>%
    
    # Rescale:
    mutate(
      month=as.factor(month),
      year_week=as.factor(year_week),
      co2_avg=co2_avg/100,
      cases_per_100=cases_per_100/1000 )
  
  return(extra_vars)
}

#################################################
### Model outcome A: probability of detection ###
### Generalised estimating equations          ###
#################################################

# Function to summarise model results (excluding factors pathogen, age group and month)
gee_results <- function(m) {
  results <- as_tibble( summary(m)$coefficients, rownames="variable") %>%
    filter(!grepl("factor", variable, fixed = TRUE)) %>%
    mutate(
      odds_ratio_CI_low = exp( Estimate - Std.err*1.96 ),
      odds_ratio_CI_high = exp( Estimate + Std.err*1.96),
      p = `Pr(>|W|)`,
      significant = if_else(
        odds_ratio_CI_low > 1, "higher",
        if_else(
          odds_ratio_CI_high < 1, "lower",
          "no"
        ))) %>%
    return()
}

# define data
df_for_gee <- load_source_data(data_file) %>% arrange(sample_id)

# Full binomial model with all variables
gee_model <- geeglm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + sampling_time
  + co2_avg + humidity_avg + temperature_inferred + ventilation_inferred  + purifier + hvac
  + n_attendees_mean + attendee_density + mask_inferred + babbelen_inferred
  + cases_per_100
  , data=df_for_gee, id=sample_id, family="binomial", corstr="exchangeable")
gee_model %>% summary
print(gee_results(gee_model))  # Odds ratio 95% confidence interval and p-value for variables which are not factors.
drop1(gee_model, test="Wald")  # Determine which variable has the highest p-value (including factors)
# --> with all variables included, the following have a significant effect:
#       pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Binomial model after backward elimination
# Backward elimination procedure: starting from the full model above, the variable with the highest
# p-value is removed and the model run again in a stepwise process, until all remaining variables
# have a p-value higher than 0.05.
# The resulting model is:
gee_model <- geeglm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + ventilation_inferred
  + babbelen_inferred
  , data=df_for_gee, id=sample_id, family="binomial", corstr="exchangeable")
gee_model %>% summary
print(gee_results(gee_model))
drop1(gee_model, test="Wald")
# --> Retained variables: pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Next, we replace the inferred values of ventilation and babbelen in the above model with only observed values.
# Omit any observations where observed values are missing.
df_for_gee_not_inferred <- df_for_gee %>% drop_na(ventilation_mean, babbelen_mean)
paste0(nrow(df_for_gee_not_inferred)," out of ", nrow(df_for_gee), " observations remain after excluding inferred variables.")
gee_model <- geeglm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + ventilation_mean + babbelen_mean
  , data=df_for_gee_not_inferred, id=c(sample_id), family="binomial", corstr="exchangeable")
gee_model %>% summary
print(gee_results(gee_model))
drop1(gee_model, test="Wald")
# --> babbelen is no longer significant

# Final model after further backward elimination (babbelen removed)
df_for_gee_not_inferred <- df_for_gee %>% drop_na(ventilation_mean)
paste0(nrow(df_for_gee_not_inferred)," out of ", nrow(df_for_gee), " observations remain after excluding inferred variables.")
gee_model <- geeglm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month) # confounders which are not dropped from the model
  + co2_avg + ventilation_mean
  , data=df_for_gee_not_inferred, id=c(sample_id), family="binomial", corstr="exchangeable")
gee_model %>% summary
(res <- gee_results(gee_model))  # Odds ratio 95% confidence interval and p-value
drop1(gee_model, test="Wald")
paste0("Odds ratio of increase of 100 in CO2: ", exp(res$Estimate[res$variable=="co2_avg"]))
paste0("Odds raio of step increase in ventilation: ", exp(res$Estimate[res$variable=="ventilation_mean"]))

#################################################
### Model outcome A: probability of detection ###
### Generalised mixed-effect model            ###
#################################################

# Function to summarise model results
glmer_results <- function(m) {
  results <- as_tibble( summary(m)$coefficients, rownames="variable") %>%
    filter(!grepl("pathogen", variable, fixed = TRUE)) %>%
    filter(!grepl("month", variable, fixed = TRUE)) %>%
    filter(!grepl("age_group", variable, fixed = TRUE)) %>%
    mutate(
      odds_ratio = exp(Estimate),
      odds_ratio_CI_low = exp( Estimate - `Std. Error`*1.96 ),
      odds_ratio_CI_high = exp( Estimate + `Std. Error`*1.96),
      ) %>%
    return()
}

df_for_glmer <- load_source_data(data_file)

# Full binomial model with all variables
glmer_model <- glmer(
  detected ~ 
    pathogen + month + age_group
  + sampling_time
  + co2_avg + humidity_avg + temperature_inferred + ventilation_inferred + purifier + hvac
  + n_attendees_mean  + attendee_density + mask_inferred + babbelen_inferred
  + cases_per_100
  + (1|sample_id)
  , data=df_for_glmer, family="binomial")
glmer_model %>% summary
# --> with all variables included, the following have a significant effect:
#       pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Binomial model after backward elimination
# Backward elimination procedure: starting from the full model above, the variable with the highest
# p-value is removed and the model run again in a stepwise process, until all remaining variables
# have a p-value higher than 0.05.
# The resulting model is:
glmer_model <- glmer(
  detected ~ 
    pathogen + month + age_group
  + co2_avg + ventilation_inferred
  + babbelen_inferred
  + (1|sample_id)
  , data=df_for_glmer, family="binomial")
glmer_model %>% summary
drop1(glmer_model, test="Chisq")
# --> Retained variables: pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Next, we replace the inferred values of ventilation and babbelen in the above model with only observed values.
# Omit any observations where observed values are missing.
df_for_glmer_not_inferred <- df_for_glmer %>% drop_na(ventilation_mean, babbelen_mean)
paste0(nrow(df_for_glmer_not_inferred)," out of ", nrow(df_for_glmer), " observations remain after excluding inferred variables.")
glmer_model <- glmer(
  detected ~ 
    pathogen + month + age_group
  + co2_avg + ventilation_mean + babbelen_mean
  + (1|sample_id)
  ,data=df_for_glmer_not_inferred, family="binomial")
glmer_model %>% summary
drop1(glmer_model, test="Chisq")
# --> babbelen is no longer significant

# Final model after further backward elimination (babbelen removed)
df_for_glmer_not_inferred <- df_for_glmer %>% drop_na(ventilation_mean)
paste0(nrow(df_for_glmer_not_inferred)," out of ", nrow(df_for_glmer), " observations remain after excluding inferred variables.")
glmer_model <- glmer(
  detected ~ 
    pathogen + month + age_group
  + co2_avg + ventilation_mean
  + (1|sample_id)
  ,data=df_for_glmer_not_inferred, family="binomial")
(res <- glmer_model %>% summary)
drop1(glmer_model, test="Chisq")
print(glmer_results(glmer_model))
paste0("Odds ratio of increase of 100 in CO2: ", exp(res$coefficients['co2_avg','Estimate']))
paste0("Odds ratio of one step increase in ventilation: ", exp(res$coefficients['ventilation_mean','Estimate']))

#############################################
# Model outcome A: probability of detection #
# Regular GLM model                         #
#############################################

# function for summarising model results
glm_results <- function(m) {
  results <- as_tibble( summary(m)$coefficients, rownames="variable") %>%
    #filter(!grepl("pathogen", variable, fixed = TRUE)) %>%
    #filter(!grepl("month", variable, fixed = TRUE)) %>%
    #filter(!grepl("age_group", variable, fixed = TRUE)) %>%
    mutate(
      odds_ratio = exp(Estimate),
      odds_ratio_CI_low = exp( Estimate - `Std. Error`*1.96 ),
      odds_ratio_CI_high = exp( Estimate + `Std. Error`*1.96),
    ) %>%
    return()
}

df_for_logit <- load_source_data(data_file)

# Binomial model with all selected initial variables and all pathogens
model <- glm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + sampling_time
  + co2_avg + humidity_avg + temperature_inferred + ventilation_inferred  + purifier + hvac
  + n_attendees_mean + attendee_density + mask_inferred + babbelen_inferred
  + cases_per_100
  ,data=df_for_logit, y=T, x=T, family="binomial")
model %>% summary
# --> with all variables included, the following have a significant effect:
#       pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Manual backward elimination results in:
model <- glm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + ventilation_inferred
  + babbelen_inferred
  ,data=df_for_logit, y=T, x=T, family="binomial")
model %>% summary
drop1(model,test="Chisq")
# --> Retained variables: pathogen, age group, month, co2 (increase), ventilation inferred (decrease), babbelen inferred (decrease)

# Next, we replace the inferred values of ventilation and babbelen in the above model with only observed values.
# Omit any observations where observed values are missing.
df_for_logit_not_inferred <- df_for_logit %>% drop_na(ventilation_mean, babbelen_mean)
paste0(nrow(df_for_logit_not_inferred)," out of ", nrow(df_for_logit), " observations remain after excluding inferred variables.")
model <- glm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + ventilation_mean + babbelen_mean
  ,data=df_for_logit_not_inferred, y=T, x=T, family="binomial")
model %>% summary
drop1(model,test="Chisq")
# --> babbelen no longer significant

# Final model with babbelen removed:
df_for_logit_not_inferred <- df_for_logit %>% drop_na(ventilation_mean)
paste0(nrow(df_for_logit_not_inferred)," out of ", nrow(df_for_logit), " observations remain after excluding inferred variables.")
model <- glm(
  detected ~ 
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + ventilation_mean
  ,data=df_for_logit_not_inferred, y=T, x=T, family="binomial")
(res <- model %>% summary)
drop1(model,test="Chisq")
print(glm_results(model), n=50)
paste0("Odds ratio of increase of 100 in CO2: ", exp(res$coefficients['co2_avg','Estimate']))
paste0("Odds ratio of one step increase in ventilation: ", exp(res$coefficients['ventilation_mean','Estimate']))


##################################
# Model outcome B: CT value      #
# Regular GLM                    #
##################################

# function for summarising results
glm_ct_results <- function(m) {
  results <- as_tibble( summary(m)$coefficients, rownames="variable") %>%
    mutate(
      CI_low = Estimate - `Std. Error`*1.96 ,
      CI_high =  Estimate + `Std. Error`*1.96,
    ) %>%
    return()
}

df_for_ct <- load_source_data(data_file) %>%
  filter(detected, !is.na(ct_value)) # only include positive samples with ct value
  
# Linear model with all selected initial variables and all pathogens
model <- glm(
    ct_value ~
    factor(pathogen) + factor(age_group) + factor(month)
    + sampling_time
    + co2_avg + humidity_avg + temperature_inferred + ventilation_inferred  + purifier + hvac
    + n_attendees_mean + attendee_density + mask_inferred + babbelen_inferred
    + cases_per_100
    ,data=df_for_ct)
model %>% summary
# --> with all variables included, the following have a significant effect:
#       pathogen, age group, month, co2 (more pathogen), purifier (less pathogen)

# Manual backward elimination results in:
model <- glm(
  ct_value ~
    factor(pathogen) + factor(age_group) + factor(month)
  + co2_avg + purifier
  ,data=df_for_ct)
model %>% summary
(res <- model %>% summary)
drop1(model,test="Chisq")
print(glm_ct_results(model),n=50)
# --> Retained variables: pathogen, age group, month, co2 (more pathogen), purifier (less pathogen)
paste0("Fewer PCR cycles per increase of 100 in CO2: ", res$coefficients['co2_avg','Estimate'])
paste0("More PCR cycles in case of air filtration: ", res$coefficients['purifierTRUE','Estimate'])

##################################
# Model outcome B: CT value      #
# Mixed effects model            #
##################################

# function for summarising results
glmer_ct_results <- function(m) {
  results <- as_tibble( summary(m)$coefficients, rownames="variable") %>%
    mutate(
      CI_low = Estimate - `Std. Error`*1.96 ,
      CI_high =  Estimate + `Std. Error`*1.96,
    ) %>%
    return()
}

df_for_ct <- load_source_data(data_file) %>%
  filter(detected, !is.na(ct_value)) # only include positive samples with ct value

# Linear model with all selected initial variables and all pathogens
# Correcting for sample ID with mixed effects model:
model <- lmer(
  ct_value ~
    (1|sample_id) + pathogen + age_group + month
  + sampling_time
  + co2_avg + humidity_avg + temperature_inferred + ventilation_inferred  + purifier + hvac
  + n_attendees_mean + attendee_density + mask_inferred + babbelen_inferred
  + cases_per_100
  ,data=df_for_ct)
as_tibble(summary(model)$coef,rownames="variable") %>%
  mutate(p_value=2 * (1 - pnorm(abs(`t value`)))) %>% print(n=50)
# --> with all variables included, the following have a significant effect:
#       pathogen, age group, month, co2 (more pathogen), purifier (less pathogen)

# Manual backward elimination results in:
model <- lmer(
  ct_value ~
    (1|sample_id) + pathogen + age_group + month
  + co2_avg + purifier
  ,data=df_for_ct)
as_tibble(summary(model)$coef,rownames="variable") %>%
  mutate(p_value=2 * (1 - pnorm(abs(`t value`))))
drop1(model,test="Chisq")
print(glmer_ct_results(model), n=50)
# --> Retained variables: pathogen, age group, month, co2 (more pathogen), purifier (less pathogen)
paste0("Fewer PCR cycles per increase of 100 in CO2: ", summary(model)$coefficients['co2_avg','Estimate'])
paste0("More PCR cycles in case of air filtration: ", summary(model)$coefficients['purifierTRUE','Estimate'])


#############################################
# Model outcome A: probability of detection #
# Single pathogen GLM models                #
#############################################
# Run binomial model using only the retained variables from the models above

df_for_single_path <- load_source_data(data_file)

backw_elim = T # Run further backward elimination on single pathogen models

# Initialise
pathogens <- unique(df_for_single_path$pathogen)
results <- tibble()

# Repeat modelling for each pathogen
for (i in pathogens) {
  print(i)
  df <- df_for_single_path %>% filter(pathogen==i)
  df_not_inferred <- df %>% drop_na(ventilation_mean)
  paste0(nrow(df_not_inferred)," out of ", nrow(df), " observations remain after excluding inferred variables.") %>% print
  
  # Model probability of detection
  model <- glm( detected ~ factor(age_group) + factor(month) + co2_avg + ventilation_mean, data=df_not_inferred, y=T, x=T, family="binomial")
  if (backw_elim) {
    p_co2 <- summary(model)$coefficients["co2_avg","Pr(>|z|)"]
    p_vent <- summary(model)$coefficients["ventilation_mean","Pr(>|z|)"]
    
    if(p_co2 >= 0.05 | p_vent >= 0.05) {
      if( p_co2 >= p_vent ) {
        model <- glm( detected ~ factor(age_group) + factor(month) + ventilation_mean, data=df_not_inferred, family="binomial")
      } else {
        paste0("All ", nrow(df), " observations included after elimination of ventilation.") %>% print
        model <- glm( detected ~ factor(age_group) + factor(month) + co2_avg, data=df, family="binomial")
      }
    }
  }
  coef <- as_tibble(summary(model)$coefficients,rownames="retained_variable") %>%
    filter(!grepl( "factor", retained_variable, fixed = TRUE) , !grepl( "Intercept", retained_variable, fixed = TRUE)) %>% # do not include age_group, month or intercept
    mutate(pathogen=i) %>%
    filter(`Pr(>|z|)` <= 0.05)
    
  results <- bind_rows(results, coef)
}

results <- tibble(pathogen=pathogens) %>%
  full_join(results, multiple="all") %>%
  mutate(
    odds_ratio = exp(Estimate),
    odds_ratio_CI_low = exp(Estimate - `Std. Error`*1.96) ,
    odds_ratio_CI_high = exp(Estimate + `Std. Error`*1.96)
  )

print(results)
write_csv2(results, "Output/single_pathogen_models_outcome_A.csv")


##########################################
# Model outcome 3: CT value              #
# Single pathogen GLM models             #
##########################################
# Run linear model using only the retained variables from the multi-pathogen models above

df_for_single_path <- load_source_data(data_file)

backw_elim = T  # Run further backward elimination on single pathogen models

# Initialise
pathogens <- unique(df_for_single_path$pathogen)
results <- tibble()

# Repeat modelling for each pathogen
for (i in pathogens) {
  print(i)
  df <- df_for_single_path %>% filter(pathogen==i, detected, !is.na(ct_value)) # only include positive samples with ct value

  # Model probability of detection
  model <- glm( ct_value ~ factor(age_group) + factor(month) + co2_avg + purifier, data=df)
  if (backw_elim) {
    p_co2 <- summary(model)$coefficients["co2_avg","Pr(>|t|)"]
    p_purifier <- summary(model)$coefficients["purifierTRUE","Pr(>|t|)"]
    
    if(p_co2 >= 0.05 | p_purifier >= 0.05) {
      if( p_co2 >= p_purifier ) {
        model <- glm( ct_value ~ factor(age_group) + factor(month) + purifier, data=df)
      } else {
        model <- glm( ct_value ~ factor(age_group) + factor(month) + co2_avg, data=df)
      }
    }
  }
  coef <- as_tibble(summary(model)$coefficients,rownames="retained_variable") %>%
    filter(!grepl( "factor", retained_variable, fixed = TRUE) , !grepl( "Intercept", retained_variable, fixed = TRUE)) %>% # do not include age_group, month or intercept
    mutate(pathogen=i) %>%
    filter(`Pr(>|t|)` <= 0.05)
  print(summary(model)$coefficients)
  
  results <- bind_rows(results, coef)
}

results <- tibble(pathogen=pathogens) %>%
  full_join(results, multiple="all") %>%
  mutate(
    CI_low = Estimate - `Std. Error`*1.96 ,
    CI_high = Estimate + `Std. Error`*1.96
  )

print(results, n=50)
write_csv2(results, "Output/single_pathogen_models_outcome_B.csv")


