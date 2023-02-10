# Paper title:
# Authors:
# Figure description: this figure shows all test results over the course of the study period for each pathogen,
# by timing and age group of air sample collection.

library(tidyverse)
library(lubridate)
library(readxl)
library(ggtext)
source("Code/data_cleaning.R")
theme_set(theme_bw())

data_file <- "Input/source_data_v2.xlsx"

data <- read_excel(data_file, sheet=1, na = "NA") %>%
  recode_pathogens() %>%
  group_by(pathogen) %>%
  filter(sum(detected)>=1) %>% # remove pathogens with less than 1 positive results
  ungroup() %>%
  
  mutate(
    # Order age group from old to young
    `Age group` = fct_relevel(age_group, "Elderly", "Sportsbar", "University", "Secondary", "Primary", "Kindergarten", "Child care group"),
    # Rename age groups
    `Age group` = recode_factor(`Age group`,
                                `Elderly`="Elderly care home",
                                `Sportsbar`="University / Sports club",
                                `University`="University / Sports club",
                                `Secondary`="Secondary school",
                                `Primary`="Primary school",
                                `Kindergarten`="Preschool",
                                `Child care group`="Nursery"
    ),
    # Rename date
    Date=as.Date(date_start),
    # Rename results
    Result=if_else(detected,"Positive","Negative"),
  ) %>%
  filter(pathogen!="SARS-CoV-2 (respiratory panel)")


ggplot(data, aes(Date,`Age group`,colour=Result)) +
  scale_color_manual(values=c("black", "red")) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  facet_wrap(~pathogen, ncol=4) +
  geom_jitter(width=0.4, height=0.3, size=0.6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_rect(fill="white")) +
  theme(strip.text = ggtext::element_markdown())

ggsave("Output/supp_figure_1.png", width=11.2, height=10)
ggsave("Output/supp_figure_1.pdf", width=11.2, height=10)
