
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
  filter(sum(detected)>0) %>% # remove pathogens with no positive results
  ungroup() %>%
  mutate(date=as.Date(date_start)) %>%
  filter(pathogen!="SARS-CoV-2 (respiratory panel)")

# number of aerosol samples included
length(unique(data$sample_id))

respi_data <- read_csv2("Input/hospital-respi-data.csv") %>%
  filter(
    `aanvragende-bron` == "LIS/Laboratoria UZ Leuven",
    # age <= 3,
    resultaat == "positief" | resultaat == "negatief",
    is.na(studie),
    date >= ymd(20211014),
    date <= ymd(20220427)
  ) %>%
  mutate(
    detected = (resultaat=="positief"),
    pathogen = substr(labotest,6,200),
    pathogen = if_else(pathogen=="Coronavirus SARS/Wuhan", "SARS-CoV-2", pathogen),
    date = as.Date(date)
  ) %>%
  recode_pathogens() %>%
  inner_join(tibble(pathogen=unique(data$pathogen)), by="pathogen") # exclude pathogens that are not in the aerosol data (because they were excluded)

respi_data %>% group_by(staalsoort) %>% summarise(n=n())

# number of hospital samples included
length(unique(respi_data$sample_id))

# prob by month and pathogen
ggplot(data, aes(date, as.numeric(detected))) +
  coord_cartesian(ylim=c(0, 1)) +
  facet_wrap(~pathogen, ncol=5) +
  geom_jitter(size= 0.3, height=0.03, width=0.3, aes(colour=age_group)) +
  geom_smooth(method = loess, colour="red", fill="red") +
  geom_smooth(method = loess, data=respi_data, colour="black", fill="black") +
  guides(colour = guide_legend(override.aes = list(size=2))) +
  labs(y= "Pathogen detection rate", x = "Sampling date", colour="Age group") +
  theme(strip.text = ggtext::element_markdown()) +
  scale_x_date(date_labels = "%b")


ggsave("Output/supp_figure_2.pdf", width=13, height=11)
ggsave("Output/supp_figure_2.png", width=13, height=11)
