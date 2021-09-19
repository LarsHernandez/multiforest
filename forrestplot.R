

library(forestmodel)
library(survival)
library(dplyr)
library(tidyverse)
library(forestmodel)
library(rlang)
library(pec)
library(labelled)



source("forestmodels_new.R")

data(Pbc3,package="pec")

pb <- Pbc3 %>%
  transmute(
    dead1yr = if_else(days < (365 * 5) & status == 2, T, F),
    tment   = if_else(tment == 0, "Placebo", "CyA"),    
    sex     = if_else(sex == 1, "Male", "Female"),
    age = case_when(
      between(age, 0, 50) ~ "18 - 50",
      between(age, 50, 60) ~ "50 - 60",
      between(age, 60, 140) ~ "60 - 75"),
    weight = case_when(
      between(weight, 0, 55) ~ "< 55kg",
      between(weight, 55, 65) ~ "55 - 65kg",
      between(weight, 65, 100) ~ "65 - 100kg"),
    stage = case_when(
      stage < 3 ~ "I-II",
      stage == 3 ~ "III",
      stage == 4 ~ "IV",
      is.na(stage) ~ "IV"),
    unit = case_when(
      unit == 1 ~ "Hvidovre",
      unit == 2 ~ "London",
      unit == 3 ~ "Copenhagen",
      unit == 4 ~ "Barcelona",
      unit == 5 ~ "Munich",
      unit == 6 ~ "Lyon"),
    gibleed = if_else(gibleed == 1, "Yes", "No")) %>%
  mutate(unit  = fct_relevel(unit, "London"),
         tment = fct_relevel(tment, "Placebo"))

var_label(pb) <- list(tment = "Treatment",
                      sex = "Sex",
                      age = "Age group",
                      weight = "Weight group",
                      stage = "Cancer stage",
                      unit = "Hospital",
                      gibleed = "Gastrointestinal bleeding")


p <- mforestmodel(pb, dependent="dead1yr", pala="grey75", palb="grey15", lim=c(-2.4,2.4))
p

ggsave(p, file="multivariate_forrestplot3.png",width=10, height=8, dpi=300)



 

library(survival)

cg <-cgd %>% 
  select(sex, age, treat, hos.cat, status, weight, inherit, height) %>%
  mutate(
    treat = as.character(treat),
    hos.cat = as.character(hos.cat),
    weight = case_when(
      between(weight, 0, 55) ~ "< 55kg",
      between(weight, 55, 65) ~ "55 - 65kg",
      between(weight, 65, 100) ~ "65 - 100kg"),
    height = case_when(
      between(height, 0, 120) ~ "< 1.2m",
      between(height, 120, 150) ~ "1.2 - 1.5m",
      between(height, 150, 220) ~ "1.5 - 2.2m"),
    age = case_when(
      between(age, 0, 10) ~ "0 - 10",
      between(age, 11, 20) ~ "11 - 20",
      between(age, 21, 140) ~ "21 - 40"))
  

p2 <- mforestmodel(cg, dependent="status", pala="grey75", palb="grey15", lim=c(-2.4,2.4))
p2

ggsave(p2, file="multivariate_forrestplota.png",width=10, height=8, dpi=300)
