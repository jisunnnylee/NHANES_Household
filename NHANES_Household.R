# ===================================================================================================================
# project: Qol research in cancer patients - household
#
# author: Jiseon Lee, Sooyeon Kim
# created: 2024-07-10
# last edited: 2024-12-31
#
# OBJECTIVE: Household NHANES/K-NHANES
#
#
# ===================================================================================================================
# install.packages("survey")
library(survey)
library(RcmdrMisc)
library(readxl)
library(dplyr)
library(tidyverse)
library(tableone)
library(factoextra)
library(GPArotation)
library(corrplot)
library(ggplot2)
library(semPlot)
library(semTools)
library(lavaan)
library(DACF)
library(lavaanPlot)
library(skimr)
library(knitr)
library(mirt)
library("irr")
library(SimplyAgree)
library(pROC)
library(lordif)
library(jmv)
library(car)
library(reshape2)
library(plyr)
library(tidyr)
library(broom)
library(readstata13) 
library(nhanesA)
library(foreign)
library(data.table)
library(magrittr)
library(stats)
library(survey)
library(survival)
library(car)
library(stringr)
library(readr)
library(ggeffects)
library(effects)
library(openxlsx)
# ====================================================================================================================================================================
# Data 불러오기 ------------------------------------------------------------------------------------------------------------------------------------------------------
load("/Users/sweetpotatogirl/Desktop/SMC/국건영/data/data_NHANES_20240909.Rdata") # N = 6,400

weight <- svydesign(id = ~ p_id,
                    strata = ~ strata,
                    weights = ~ persweight,
                    nest = TRUE,
                    data = df)
options(survey.lonely.psu = "adjust")

##Table 1 
df$unemployed %>% skim
svyCreateTableOne(vars = c("p_age", "age_cat","p_sex", "p_edu", "p_marri2", "house_inc", "smking", "alcohol", "unemployed", "country", "cancer_type", "yr_cancer", "yr_cancer_cat"),
                  factorVars = c("age_cat", "p_sex", "p_edu", "p_marri2", "house_inc", "smking", "alcohol",  "unemployed", "country", "cancer_type", "yr_cancer_cat"),
                  includeNA = F, strata = "alone", data = weight)


svyby(~p_age, by = ~alone, design = weight, FUN = svymean, na.rm = T) 
svyby(~yr_cancer, by = ~alone, design = weight, FUN = svymean, na.rm = T) 


# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Table 2. Adjusted OR with 95% CI 
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

svyCreateTableOne(vars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  factorVars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  includeNA = F, strata = "alone", data = weight) %>% print(showAllLevels = T, noSpaces = T)



#OR function 생성
or1 <- function(outcome){
  y <- outcome
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  res <- data.frame(
    outcome = outcome, 
    hr = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"))
  
  # result <- bind_rows(model1 %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term)),
  #                     model2 %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term)),
  #                     model3 %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term)))
  
  return(res)
}

bind_rows(or1("as.numeric(poor)"),
          or1("as.numeric(smking)"),
          or1("as.numeric(alcohol)"),
          or1("as.numeric(underweight)"),
          or1("as.numeric(hypertension)"),
          or1("as.factor(dm)"),
          or1("as.numeric(dyslipidemia)"),
          or1("as.numeric(pf_lim_e)"),
          or1("as.numeric(suicide)")) #%>% write.csv("tab2_20241231.csv")

or2 <- function(outcome){
  y <- outcome
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer_cat")),
                   design = weight, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  res <- data.frame(
    outcome = outcome, 
    hr = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    p_for_int = round(car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(age_cat) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                                        design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)], 2))
  
  return(res)
}

or2("as.numeric(nofood)")

# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# Table 3. Subgroup analysis #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #
# ------------------------------------------------------------------------------------- #

#OR function 생성
?svyglm
as.svrepdesign(weight)

weight_age1 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(age_cat == 1))
weight_age2 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(age_cat == 2))
weight_age3 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(age_cat == 3))
weight_age4 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(age_cat == 4))


subgroup_age <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_age1, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_age2, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_age3, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model4 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_age4, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  
  res <- data.frame(
    outcome = out2come, 
    age1 = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    age2 = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    age3 = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    age4 = paste0(sprintf("%.2f", model4$estimate), " (", sprintf("%.2f", model4$conf.low), ", ", sprintf("%.2f", model4$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(age_cat) + p_sex + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)] %>% round(3)))
  
  return(res)
}



weight_male <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(p_sex == 1))
weight_female <- svydesign(id = ~ p_id,
                           strata = ~ strata,
                           weights = ~ persweight,
                           nest = TRUE,
                           data = df %>% filter(p_sex == 2))


subgroup_sex <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_male, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                   design = weight_female, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    male = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    female = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(p_sex) + p_age + factor(country) + cancer_type + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}

weight_US <- svydesign(id = ~ p_id,
                       strata = ~ strata,
                       weights = ~ persweight,
                       nest = TRUE,
                       data = df %>% filter(country == "US"))
weight_KOR <- svydesign(id = ~ p_id,
                        strata = ~ strata,
                        weights = ~ persweight,
                        nest = TRUE,
                        data = df %>% filter(country == "KOR"))

subgroup_country <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer_cat")),
                   design = weight_US, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer_cat")),
                   design = weight_KOR, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    us = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    kor = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(country) + p_age + p_sex + cancer_type + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}

weight_dig1 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(yr_cancer_cat == 1))
weight_dig2 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(yr_cancer_cat == 2))
weight_dig3 <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(yr_cancer_cat == 3))


subgroup_yr_cancer <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + cancer_type + yr_cancer")),
                   design = weight_dig1, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + cancer_type + yr_cancer")),
                   design = weight_dig2, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + cancer_type + yr_cancer")),
                   design = weight_dig3, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    dig1 = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    dig2 = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    dig3 = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(yr_cancer_cat) + p_age + p_sex + factor(country) + cancer_type")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}


weight_Gynereproc <- svydesign(id = ~ p_id,
                               strata = ~ strata,
                               weights = ~ persweight,
                               nest = TRUE,
                               data = df %>% filter(cancer_type == "Gynereproc"))
weight_Stomach <- svydesign(id = ~ p_id,
                            strata = ~ strata,
                            weights = ~ persweight,
                            nest = TRUE,
                            data = df %>% filter(cancer_type == "Stomach"))
weight_Breast <- svydesign(id = ~ p_id,
                           strata = ~ strata,
                           weights = ~ persweight,
                           nest = TRUE,
                           data = df %>% filter(cancer_type == "Breast"))
weight_Skin <- svydesign(id = ~ p_id,
                         strata = ~ strata,
                         weights = ~ persweight,
                         nest = TRUE,
                         data = df %>% filter(cancer_type == "Skin"))
weight_Thyroid <- svydesign(id = ~ p_id,
                            strata = ~ strata,
                            weights = ~ persweight,
                            nest = TRUE,
                            data = df %>% filter(cancer_type == "Thyroid"))
weight_Others <- svydesign(id = ~ p_id,
                           strata = ~ strata,
                           weights = ~ persweight,
                           nest = TRUE,
                           data = df %>% filter(cancer_type == "Others"))

subgroup_cancer_type <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Gynereproc, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Stomach, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Breast, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model4 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Skin, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model5 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Thyroid, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model6 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                   design = weight_Others, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    gynerproc = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    stomach = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    breast = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    skin = paste0(sprintf("%.2f", model4$estimate), " (", sprintf("%.2f", model4$conf.low), ", ", sprintf("%.2f", model4$conf.high), ")"),
    thyroid = paste0(sprintf("%.2f", model5$estimate), " (", sprintf("%.2f", model5$conf.low), ", ", sprintf("%.2f", model5$conf.high), ")"),
    others = paste0(sprintf("%.2f", model6$estimate), " (", sprintf("%.2f", model6$conf.low), ", ", sprintf("%.2f", model6$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(cancer_type) + p_age + p_sex + factor(country) + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}


#각 outcome 변수명을 괄호안에 넣으면 OR+95%CI 결과 추출 가능
bind_rows(subgroup_age("as.numeric(poor)"),
          subgroup_age("as.numeric(smking)"),
          subgroup_age("as.numeric(alcohol)"),
          subgroup_age("as.numeric(underweight)"),
          subgroup_age("as.factor(hypertension)"),
          subgroup_age("as.numeric(dm)"),
          subgroup_age("as.numeric(dyslipidemia)"),
          subgroup_age("as.numeric(pf_lim_e)"),
          subgroup_age("as.numeric(suicide)")) #%>% write.csv("tab3_age_20241231.csv")

bind_rows(subgroup_sex("as.numeric(poor)"),
          subgroup_sex("as.numeric(smking)"),
          subgroup_sex("as.numeric(alcohol)"),
          subgroup_sex("as.numeric(underweight)"),
          subgroup_sex("as.factor(hypertension)"),
          subgroup_sex("as.numeric(dm)"),
          subgroup_sex("as.numeric(dyslipidemia)"),
          subgroup_sex("as.numeric(pf_lim_e)"),
          subgroup_sex("as.numeric(suicide)")) #%>% write.csv("tab3_sex_20241231.csv")

bind_rows(subgroup_country("as.numeric(poor)"),
          subgroup_country("as.numeric(smking)"),
          subgroup_country("as.numeric(alcohol)"),
          subgroup_country("as.numeric(underweight)"),
          subgroup_country("as.factor(hypertension)"),
          subgroup_country("as.numeric(dm)"),
          subgroup_country("as.numeric(dyslipidemia)"),
          subgroup_country("as.numeric(pf_lim_e)"),
          subgroup_country("as.numeric(suicide)")) #%>% write.csv("tab3_country_20241231.csv")

bind_rows(subgroup_yr_cancer("as.numeric(poor)"),
          subgroup_yr_cancer("as.numeric(smking)"),
          subgroup_yr_cancer("as.numeric(alcohol)"),
          subgroup_yr_cancer("as.numeric(underweight)"),
          subgroup_yr_cancer("as.factor(hypertension)"),
          subgroup_yr_cancer("as.numeric(dm)"),
          subgroup_yr_cancer("as.numeric(dyslipidemia)"),
          subgroup_yr_cancer("as.numeric(pf_lim_e)"),
          subgroup_yr_cancer("as.numeric(suicide)")) #%>% write.csv("tab4_yrcancer_20241231.csv")


bind_rows(subgroup_cancer_type("as.numeric(poor)"),
          subgroup_cancer_type("as.numeric(smking)"),
          subgroup_cancer_type("as.numeric(alcohol)"),
          subgroup_cancer_type("as.numeric(underweight)"),
          subgroup_cancer_type("as.factor(hypertension)"),
          subgroup_cancer_type("as.numeric(dm)"),
          subgroup_cancer_type("as.numeric(dyslipidemia)"),
          subgroup_cancer_type("as.numeric(pf_lim_e)"),
          subgroup_cancer_type("as.numeric(suicide)")) # %>% write.csv("tab5_type_cancer_20241231.csv")

with(df, table(primary_c, pf_lim_e))

#한국에만 존재하는 데이터#
subgroup_age2 <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_age1, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_age2, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_age3, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model4 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_age4, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  
  res <- data.frame(
    outcome = out2come, 
    age1 = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    age2 = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    age3 = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    age4 = paste0(sprintf("%.2f", model4$estimate), " (", sprintf("%.2f", model4$conf.low), ", ", sprintf("%.2f", model4$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(age_cat) + p_sex + p_age + cancer_type + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}

subgroup_sex2 <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_male, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + cancer_type + yr_cancer_cat")),
                   design = weight_female, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    male = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    female = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(p_sex) + p_age + cancer_type + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}



subgroup_yr_cancer2 <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer")),
                   design = weight_dig1, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer")),
                   design = weight_dig2, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + cancer_type + yr_cancer")),
                   design = weight_dig3, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    dig1 = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    dig2 = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    dig3 = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(yr_cancer_cat) + p_age + p_sex + cancer_type")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}

subgroup_cancer_type2 <- function(out2come){
  y <- out2come
  model1 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Gynereproc, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model2 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Stomach, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model3 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Breast, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model4 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Skin, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model5 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Thyroid, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  model6 <- svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) + p_age + p_sex + yr_cancer_cat")),
                   design = weight_Others, family = "binomial") %>% tidy(exponentiate = T, conf.int = T) %>% filter(grepl("alone", term))
  
  res <- data.frame(
    outcome = out2come, 
    gynerproc = paste0(sprintf("%.2f", model1$estimate), " (", sprintf("%.2f", model1$conf.low), ", ", sprintf("%.2f", model1$conf.high), ")"),
    stomach = paste0(sprintf("%.2f", model2$estimate), " (", sprintf("%.2f", model2$conf.low), ", ", sprintf("%.2f", model2$conf.high), ")"),
    breast = paste0(sprintf("%.2f", model3$estimate), " (", sprintf("%.2f", model3$conf.low), ", ", sprintf("%.2f", model3$conf.high), ")"),
    skin = paste0(sprintf("%.2f", model4$estimate), " (", sprintf("%.2f", model4$conf.low), ", ", sprintf("%.2f", model4$conf.high), ")"),
    thyroid = paste0(sprintf("%.2f", model5$estimate), " (", sprintf("%.2f", model5$conf.low), ", ", sprintf("%.2f", model5$conf.high), ")"),
    others = paste0(sprintf("%.2f", model6$estimate), " (", sprintf("%.2f", model6$conf.low), ", ", sprintf("%.2f", model6$conf.high), ")"),
    p_for_int = sprintf("%.2f", car::Anova(svyglm(as.formula(paste0(y, "~ C(factor(alone), base = 1) * factor(cancer_type) + p_age + p_sex + yr_cancer_cat")),
                                                  design = weight, family = "binomial")) %>% .$`Pr(>Chisq)` %>% .[length(.)]))
  
  return(res)
}
subgroup_age2("as.numeric(nofood)")
subgroup_sex2("as.numeric(nofood)")
subgroup_yr_cancer2("as.numeric(nofood)")
subgroup_cancer_type2("as.numeric(nofood)")
with(df, table(country))
#Figure 1
svyCreateTableOne(vars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  factorVars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  includeNA = F, strata = "alone", data = weight_age1) %>% print(showAllLevels = T, noSpaces = T)

svyCreateTableOne(vars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  factorVars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  includeNA = F, strata = "alone", data = weight_age2) %>% print(showAllLevels = T, noSpaces = T)

svyCreateTableOne(vars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  factorVars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  includeNA = F, strata = "alone", data = weight_age3) %>% print(showAllLevels = T, noSpaces = T)

svyCreateTableOne(vars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  factorVars = c("poor", "smking", "alcohol", "nofood", "underweight", "hypertension", "dm", "dyslipidemia", "pf_lim_e", "suicide"),
                  includeNA = F, strata = "alone", data = weight_age4) %>% print(showAllLevels = T, noSpaces = T)

# save(df, file = "/Users/sweetpotatogirl/Desktop/SMC/국건영/data/NHANES_20241231.RData")
