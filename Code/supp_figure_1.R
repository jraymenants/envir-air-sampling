# Paper title:
# Authors:
# Figure description: this figure shows all test results over the course of the study period for each pathogen,
# by timing and age group of air sample collection.

library(tidyverse)
theme_set(theme_bw())

data <- read_csv2("Input/organised_data_with_date_ungrouped_recoded.csv") %>%
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
    # Rename pathogens
    # pathogen = recode(pathogen,
    #                   `TaqPath COVID-19`="SARS-CoV-2",
    #                   `Adenovirus`="",
    #                   ``
    # ),
    # Rename date
    Date=date_start,
    # Rename results
    Result=if_else(detected,"Positive","Negative"),
  )


ggplot(data, aes(Date,`Age group`,colour=Result)) +
  scale_color_manual(values=c("black", "red")) + 
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  facet_wrap(~pathogen, ncol=4) +
  geom_jitter(width=0.4, height=0.3, size=0.6) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(strip.background = element_rect(fill="white"))

ggsave("Output/supp_figure_1.png", width=11.2, height=10)
