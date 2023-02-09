
library(tidyverse)
library(lubridate)
theme_set(theme_bw())

# Read source data
data <- read_csv2("Input/organised_data_with_date_ungrouped_recoded.csv") %>%
  group_by(pathogen) %>%
  filter(sum(detected)>0) %>% # remove pathogens with no positive results
  ungroup() %>%
  rename(date=date_start) %>%
  mutate(
    pathogen = if_else(pathogen=="TaqPath COVID-19", "SARS-CoV-2", pathogen)
  )

# number of aerosol samples included
length(unique(data$sample_id))

respi_data <- read_csv2("data/hospital-respi-data.csv") %>%
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
    pathogen = recode(pathogen,
                      `TaqPath COVID-19`="SARS-CoV-2",
                      `Adenovirus`="human adenovirus",
                      `Bocavirus`="human bocavirus",
                      `Chlamydophila pneumoniae`="Chlamydia pneumoniae",
                      `Coronavirus 229E`="Human coronavirus 229E",
                      "Coronavirus HKU-1"="Human coronavirus HKU-1",
                      `Coronavirus NL63`="Human coronavirus NL63",
                      `Coronavirus OC43`="Human coronavirus OC43",
                      `Coxiella burnetti`="Coxiella burnetii",
                      `Cytomegalovirus`="human cytomegalovirus",
                      "Entero-/Rhinovirus"="human enterovirus (incl. rhinovirus)",
                      `Enterovirus D68`="enterovirus D68",
                      `Herpes simplex virus 1`="herpes simplex virus type 1",
                      `Humaan metapneumovirus`="Human metapneumovirus",
                      `Influenza A virus`="influenza A virus",
                      `Mycoplasma pneumoniae`="Mycoplasma pneumoniae",
                      `Parainfluenzavirus 3`="human parainfluenza virus 3",
                      `Parainfluenzavirus 4`="human parainfluenza virus 4",
                      `Pneumocystis jiroveci`="Pneumocystis jirovecii",
                      `Respiratoir syncytieel virus`="respiratory syncytial virus",
                      `Streptococcus pneumoniae`="Streptococcus pneumoniae",
                      "Parainfluenzavirus"="Parainfluenzaviruses",
                      "Other coronavirus"="Other coronaviruses"
    )
  ) %>%
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
  labs(y= "Pathogen detection rate", x = "Sampling date", colour="Age group")


ggsave("Output/supp_figure_2.pdf", width=13, height=11)
ggsave("Output/supp_figure_2.png", width=13, height=11)
