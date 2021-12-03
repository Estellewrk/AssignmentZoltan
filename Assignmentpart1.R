## Loading packages
library(tidyverse)
library(dplyr)
library(psych)
library(lm.beta)
library(car)
library(gridExtra)
library(lme4)
library(cAIC4)
library(r2glmm)
library(lmerTest)
library(MuMIn)

#### Assignment Part 1: ####

data_sample_1 = read.csv("https://tinyurl.com/ha-dataset1")
DS1 <- data_sample_1

### Function coef_table in order to use it later in the assignment
coef_table = function(model){	
  require(lm.beta)	
  mod_sum = summary(model)	
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,4], 3))		
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"]))		
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"		
  
  
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model), confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])), 2)), mod_sum_p_values)		
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta", "p-value")		
  mod_sum_table["(Intercept)","Std.Beta"] = "0"		
  return(mod_sum_table)	
}

## Overall view of the data 
view(DS1)
describe(DS1)
summary(DS1)
## We can already detected 2 anomalies in pain and STAI_trait: 
## maximum of pain indicate 55 and minimum of STAI indicate 4,20

## Cook's distance and graph for every predictor
Modeldiagnostic <- lm(pain~age+sex+STAI_trait+pain_cat+mindfulness+cortisol_serum+cortisol_saliva, data=DS1)
plot(x = Modeldiagnostic, which = 4)
plot(x = Modeldiagnostic, which = 5)

DS1 %>% 	
  ggplot() +	
  aes(x = age) +	
  geom_histogram()

DS1 %>% 	
  ggplot() +	
  aes(x = STAI_trait) +	
  geom_histogram( bins = 50)

DS1 %>% 	
  ggplot() +	
  aes(x = pain_cat) +	
  geom_histogram( bins = 50)

DS1 %>% 	
  ggplot() +	
  aes(x = mindfulness) +	
  geom_histogram( bins = 50)

DS1 %>% 	
  ggplot() +	
  aes(x = cortisol_serum) +	
  geom_histogram( bins = 50)

DS1 %>% 	
  ggplot() +	
  aes(x = cortisol_saliva) +	
  geom_histogram( bins = 50)

DS1 %>% 	
  ggplot() +	
  aes(x = pain) +	
  geom_histogram( bins = 50)

## Errors detected: ID_43 for STAI_trait and ID_88 for Pain
## I choose to exclude these outliers instead of considering them as typo:
## in research we try not to touch/modify the data obtained. 

DS1 = DS1[-c(34,88),]
view(DS1)
describe(DS1)

################################### 
## My first choice was to considered cases 88 and 34 as typo but after a reflexion about 
# the data in research I choose to exclude them. 
# first anomaly -> max. of pain = 55 (row n°88) (pain is supposed to be rated between 0 to 10)
# second anomaly -> min of STAI = 4,20 (row n°34) is supposed to be 20 (bc STAI is rated from on scale of 20 to 80 )
# so, we assume that those 'anomalies' are considered as typographical error (a mistake in written, a misprint). 
# Then, we are gonna correct those by using the function mutate()
#mutate the first typo (in pain column)
DS1 = DS1 %>%
  mutate(pain = replace(pain, pain == 55, 5))
#mutate the second typo (in STAI_trait column)
DS1 = DS1 %>%
  mutate(STAI_trait = replace(STAI_trait, STAI_trait == 4.20, 42.0))
describe(DS1)
summary(DS1)
#####################################

# Recode the characters female and male as female=1 and male=0
DS1 = DS1 %>%
  mutate(sex_recoded = recode(sex, 
                       "male" = 0, 
                       "female" = 1))
view(DS1)

####Conduct a hierarchical regression
# model1 
model1 <- lm(pain ~ sex_recoded + age, data = DS1)
# Check the assumptions of linear regression
# Normality 
model1 %>%
  plot(which = 2)
# Linearity
model1 %>%
  residualPlots()
--> not violated but almost
# Homoscedasticity
model1 %>%
  plot(which = 3)
--> assumption of homoscedasticity is not violated (but just in case we check that with the NCV test and the Breush-Pagan test)
# NCV test 
model1 %>%
  ncvTest()
=> p value is not significant 
# Breush-Pagan test
model1 %>%
  bptest()
=> p value is not significant 
# Collinearity
model1 %>%
  vif()
#results => sex_recoded: 1.004713 age: 1.004713, means that mulicollinearity is not problematic 


# model2 
model2 <- lm(pain ~ sex_recoded + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, data = DS1)
## Check the assumptions of linear regression
# Normality
model2 %>%
  plot(which = 2)
=> assumption of normality, slightly deviant from normal distribution: 104, 85, 74
   but I choose to not exclude them
# Linearity
model2 %>%
  residualPlots()
# Homoscedacity
model2 %>%
  plot(which = 3)
--> same as with model 1: check with NCV test and Breush-Pagan test
# NCV test 
model2 %>%
  ncvTest()
=> p value not significant
# Breush-Pagan test
model2 %>%
  bptest()
=> p value not significant 

# Multicollinearity
model2 %>%
  vif()
# Multicollinearity detected due to the double cortisol predictors (saliva (VIF=5.07) and serum (VIF=4.79))
# => Removing highly correlated predictors - exclude one cortisol measurement (in the cortisol description
# part of the assignment it is mentioned that "serum cortisol is often regarded in medical research as more 
# reliably related to stress": that is why I choose to remove the predictor cortisol_saliva)
# Then, create a new final model without cortisol_saliva predictor
model2b <- lm(pain ~ sex_recoded + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = DS1)
# Check of multicollinearity again
model2b %>%
  vif()
## Check the assumptions of linear regression again with the new model: model2b
# Normality
model2b %>%
  plot(which = 2)
# Linearity
model2b %>%
  residualPlots()
# Homoscedacity
model2b %>%
  plot(which = 3)
# NCV test
model2b %>%
  ncvTest()
=> p value not significant
# Breush-Pagan test
model2b %>%
  bptest()
=> p value not significant

#### Result part model1:
summary(model1)
AIC(model1)
confint(model1)
lm.beta(model1)
coef_table(model1)
summary(model1)$r.squared 

#### Result part model2:
summary(model2b)
confint(model2b)
lm.beta(model2b)
coef_table(model2b)
summary(model2b)$r.squared

## The predictors used in model 1 are a subset of the predictors used in model 2
# Model comparison in order to assess if substantial new information was gained about pain in model2
# compared to model1
AIC(model1)
-> result 574,13
AIC(model2b)
-> result 479,26
AIC(model1) - AIC(model2b)

anova(model1, model2b)
#-> p-value is significant which in this case means that model2 is better than model1





#### Assignment 2: ####

# As in the first part fo the assignment: check the variables included in the model for coding errors and the model itself for influential outliers
# Case n°34 and 88 are still excluded from this dataframe 
view(DS1)
describe(DS1)
# Creation of a new model with new predictors: weight, IQ and household income
ini_model <- lm(pain ~ sex_recoded + age + STAI_trait + pain_cat + mindfulness + cortisol_serum + weight + household_income + IQ, data = DS1)

# Cook's distance and model diagnostics
plot(x = ini_model, which = 4)
plot(x = ini_model, which = 5)
# Assumptions of linear regression
# Normality
ini_model %>%
  plot(which = 2)
--> normality is not violated
# Linearity 
ini_model %>%
  residualPlots()
--> linearity is not violated 
## homoscedasticity
ini_model %>%
  plot(which = 3)
--> assumption of homoscedasticity is not violated (but just in case we check that with the NCV test and the Breush-Pagan test)
# NCV test
ini_model %>%
  ncvTest()
=> p value is not significant 
# Breush-Pagan test
ini_model %>%
  bptest()
=> p value is not significant 
# Collinearity
ini_model %>%
  vif()
--> multicollinearity is not problematic 

# Backward regression 
backward_model = step(ini_model, direction = "backward")
# Creation of theory_based_model from the end of assignment part1
theory_based_model <- lm(pain ~ sex_recoded + age + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = DS1)


# Compare the initial model and the backward model
AIC(ini_model)
AIC(backward_model)
AIC(ini_model)-AIC(backward_model)
anova(backward_model, ini_model)
# Compare the backward model and the theory-based model on AIC and anova
AIC(backward_model)
AIC(theory_based_model)
AIC(theory_based_model)-AIC(backward_model)



## Results (data collection for the two models)
# backward model
summary(backward_model)
confint(backward_model)
lm.beta(backward_model)
coef_table(backward_model)
# theory_based_model
summary(theory_based_model)
confint(theory_based_model)
lm.beta(theory_based_model)
coef_table(theory_based_model)

## Sample n°2
data_sample_2 = read.csv("https://tinyurl.com/87v6emky")
DS2 <- data_sample_2
view(DS2)
describe(DS2)


DS2 %>% 	
  ggplot() +	
  aes(x = age) +	
  geom_histogram()

DS2 %>% 	
  ggplot() +	
  aes(x = STAI_trait) +	
  geom_histogram( bins = 50)

DS2 %>% 	
  ggplot() +	
  aes(x = pain_cat) +	
  geom_histogram( bins = 50)

DS2 %>% 	
  ggplot() +	
  aes(x = mindfulness) +	
  geom_histogram( bins = 50)

DS2 %>% 	
  ggplot() +	
  aes(x = cortisol_serum) +	
  geom_histogram( bins = 50)

DS2 %>% 	
  ggplot() +	
  aes(x = cortisol_saliva) +	
  geom_histogram( bins = 50)

DS2 %>% 	
  ggplot() +	
  aes(x = pain) +	
  geom_histogram()

# Recode the characters female and male as female=1 and male=0
DS2 = DS2 %>%
  mutate(sex_recoded = recode(sex, 
                              "male" = 0, 
                              "female" = 1))
view(DS2)

# Test of the models on data sample 2 (DS2)
# Compare predicted values and actual pain ratings
# "On data file 2, make predictions on pain using the regression models or equations of the backward model
# and the theory-based model which were "trained" on data file"
predict_test_theory = predict(theory_based_model, DS2)
predict_test_backward = predict(backward_model, DS2)

# "Which model was able to predict the actual pain ratings in DS2 better?"

# First option: calculate the sum of squared differences between predicted and actual pain
RSS_theory = sum((DS2$pain - predict_test_theory)^2)
RSS_theory
RSS_backwards = sum((DS2$pain - predict_test_backward)^2)
RSS_backwards
RSS_backwards - RSS_theory2

# Second option : calculate the sum of absolute differences for each model
RAD_theory = sum(abs(DS2$pain - predict_test_theory))
RAD_theory
RAD_backwards = sum(abs(DS2$pain - predict_test_backward))
RAD_backwards
RAD_backwards - RAD_theory






###### Assignment part3 ######

# Useful function later in the assignment
stdCoef.merMod <- function(object) {	
  sdy <- sd(getME(object,"y"))	
  sdx <- apply(getME(object,"X"), 2, sd)	
  sc <- fixef(object)*sdx/sdy	
  se.fixef <- coef(summary(object))[,"Std. Error"]	
  se <- se.fixef*sdx/sdy	
  return(data.frame(stdcoef=sc, stdse=se))	
}	

DS3 <- read_csv("https://tinyurl.com/b385chpu")
view(DS3)
describe(DS3)
# Error: the minimum for the household_income in negative so I exclude this case
# + one case is describe as woman instead of female, so replace woman as female
DS3 = DS3[-c(2),]
DS3 = DS3 %>%
  mutate(sex = replace(sex, sex=="woman", "female"))
view(DS3)
describe(DS3)

# Building factor for hospital
DS3 %>% mutate(hospital = factor(hospital))
# Build a linear mixed model on data file 3, accounting for the clustering of the data at different hospital sites 
# Random intercept model: including the random intercept of hospital-ID + the fixed effect predictors in assignment 1
Random_int = lmer(pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum + (1|hospital), data = DS3)
Random_int
summary(Random_int)
confint(Random_int)
stdCoef.merMod(Random_int)

## Compute: variance explained by the fixed effect predictors using marginal R2
#          + the variance explained by the fixed and random effect terms combined using conditional R2
# Marginal R2 & Conditional R2
r2beta(Random_int, method = "nsj", data = DS3)
r.squaredGLMM(Random_int)
# Variance for residuals
sum(residuals(Random_int)^2)

# Use the regression equation obtained on data file 3 to predict pain in data file 4
DS4 <- read.csv("https://tinyurl.com/4f8thztv")
describe(DS4)
view(DS4)

pain_pred_random4 <- predict(Random_int, DS4, allow.new.levels = TRUE)

# Compute the variance explained by the model on data file 4
# The formula: 1-(RSS/TSS)
RSS = sum((DS4$pain - pain_pred_random4)^2)
RSS
mod_mean <- lm(pain ~ 1, data = DS4)
mod_mean
TSS=sum((DS4$pain - predict(mod_mean))^2)
TSS
# Formula 1-(RSS/TSS)
R<-1-(RSS/TSS)
R

# Compare this R2 to the marginal and conditional R2 values computed for the model on data file 3
#  Formula 1-(RSS_random/TSS_random)
r.squaredGLMM(Random_int)
# results of marginal and conditional R2: R2m=0.38 and R2c=0.46
1-(RSS/TSS)
# result 0.3806956

# Build a new linear mixed effects model on dataset 3 predicting pain, include the most influential predictor 
# from the previous model, which is cortisol_serum
# Allow for both random intercept and random slope 
Random_int_slope = lmer(pain ~ cortisol_serum + (cortisol_serum|hospital), data = DS3)
Random_int_slope
Random_int_slope_opt = lmer(pain ~ cortisol_serum + (cortisol_serum|hospital), control = lmerControl(optimizer = "Nelder_Mead"), data = DS3)
Random_int_slope_opt

# Visualize the fitted regression lines for each hospital separately
predict_DS3 = DS3 %>% 
  mutate(pred_int = predict(Random_int),
         pred_slope = predict(Random_int_slope_opt))

# Ordered hostpital in order to have it in this order for the graph 
order <- c("hospital_1", "hospital_2", "hospital_3", "hospital_4", "hospital_5", "hospital_6", "hospital_7", "hospital_8", "hospital_9", "hospital_10")
predict_DS3_ordered <- arrange(transform(predict_DS3,
                                         hospital=factor(hospital,levels=order)),hospital)

############ revoir ça: Graph displaying the separate fitted regression lines for the hospitals from the mixed model 
# including only the most influential predictor in the model
predict_DS3_ordered %>% 		
  ggplot() +		
  aes(y = pain, x = cortisol_serum, group = hospital)+		
  geom_point(aes(color = hospital), size = 1) +		
  geom_line(color='red', aes(y=pred_slope, x=cortisol_serum))+		
  facet_wrap( ~ hospital, ncol = 2)
# Discuss whether the random intercept or the random slope model is a better fit for the data, and why
# Compare model fit indices/Compare RSS
sum(residuals(Random_int)^2)
sum(residuals(Random_int_slope_opt)^2)
# "relying on RSS alone in our original dataset when comparing prediction efficiency would be misleading"

# conditional AIC: cAIC 
cAIC(Random_int)$caic
cAIC(Random_int_slope_opt)$caic

# anova
anova(Random_int_slope_opt, Random_int)


