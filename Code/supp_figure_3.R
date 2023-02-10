

library(tidyverse)
library(lubridate)
library(gridExtra)
library(readxl)
library(ggtext)
source("Code/data_cleaning.R")
theme_set(theme_bw())

data_file <- "Input/source_data_v2.xlsx"

data <- read_excel(data_file, sheet=1, na = "NA") %>%
  recode_pathogens() %>%
  group_corona_parainf() %>%
  group_by(pathogen) %>%
  filter(sum(detected)>=10) %>% # remove pathogens with no positive results
  ungroup() %>%
  rename(Date=date_start) %>%
  mutate(
    pathogen = if_else(pathogen=="TaqPath COVID-19", "SARS-CoV-2", pathogen)
  )

# number of aerosol samples included
length(unique(data$sample_id))

# prob by co2 and pathogen
p1 <- ggplot(data, aes(co2_avg, as.numeric(detected))) +
  facet_wrap(~pathogen) +
  geom_jitter(size= 0.3, height=0.05, width=0.05, aes(colour=age_group)) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  labs(y= "Pathogen detection rate", x = "CO2 concentration", colour="Age group", tag="a") +
  theme(strip.text = ggtext::element_markdown())

# prob by ventilation and pathogen
p2 <- ggplot(data, aes(ventilation_mean, as.numeric(detected))) +
  facet_wrap(~pathogen) +
  geom_jitter(size= 0.3, height=0.05, width=0.2, aes(colour=age_group)) +
  geom_smooth(method = glm, method.args= list(family="binomial")) +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  labs(y= "Pathogen detection rate", x = "Ventilation", colour="Age group", tag="b") + 
  theme(strip.text = ggtext::element_markdown())

sup_fig_3 <- grid.arrange(p1, p2)

ggsave("Output/supp_figure_3.pdf", plot=sup_fig_3, width=11, height=16)
ggsave("Output/supp_figure_3.png", plot=sup_fig_3, width=11, height=16)

