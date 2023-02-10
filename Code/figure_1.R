
library(tidyverse)
library(lubridate)
library(readxl)
source("Code/data_cleaning.R")
theme_set(theme_bw())

data_file <- "Input/source_data_v2.xlsx"

# Read source data
data <- read_excel(data_file, sheet=1, na = "NA") %>%
  recode_pathogens() %>%
  filter(age_group=="Child care group") %>%
  group_by(pathogen) %>%
  filter(sum(detected)>0) %>% # remove pathogens with no positive results
  ungroup() %>%
  rename(date=date_start) %>%
  filter(pathogen!="SARS-CoV-2 (respiratory panel)")
  

# number of aerosol samples included
length(unique(data$sample_id))

respi_data <- read_csv2("Input/hospital-respi-data.csv") %>%
  filter(
    `aanvragende-bron` == "LIS/Laboratoria UZ Leuven",
    age <= 3,
    resultaat == "positief" | resultaat == "negatief",
    is.na(studie),
    date >= ymd(20211117),
    date <= ymd(20220427)
  ) %>%
  mutate(
    detected = (resultaat=="positief"),
    pathogen = substr(labotest,6,200),
    pathogen = if_else(pathogen=="Coronavirus SARS/Wuhan", "SARS-CoV-2", pathogen),
  ) %>%
  recode_pathogens() %>%
  inner_join(tibble(pathogen=unique(data$pathogen)), by="pathogen") # exclude pathogens that are not in the aerosol data (because they were excluded)

respi_data %>% group_by(staalsoort) %>% summarise(n=n())

# number of hospital samples included
length(unique(respi_data$sample_id))

all_dates <- expand_grid(
  pathogen = unique(respi_data$pathogen),
  date = seq(ymd(20211117), ymd(20220427), by="days")
)
respi_data_daily <- respi_data %>%
  filter(detected) %>%
  group_by(date, pathogen) %>%
  summarise(n=n()) %>%
  right_join(all_dates) %>%
  mutate(n=replace_na(n,0))



# prob by month and pathogen (child care only), black is hospital detection
plot_data <- bind_rows(
  data %>%
    select(date, detected, pathogen) %>%
    mutate(`Sample type` = "Air sample (nursery)"),
  respi_data %>%
    select(date, detected, pathogen) %>%
    mutate(`Sample type` = "Clinical sample (hospital)")
)
ggplot(plot_data, aes(x=date, y=as.numeric(detected), colour=`Sample type`)) +
  coord_cartesian(ylim=c(0, 1)) +
  facet_wrap(~pathogen, ncol=5) +
  geom_jitter(size= 0.3, height=0.03, width=0.3) +
  geom_smooth(method = loess, aes(fill=`Sample type`)) +
  scale_fill_manual(values=c("red","black")) +
  scale_color_manual(values=c("red","black")) +
  labs(y= "Pathogen detection rate", x = "Sampling date") +
  scale_y_continuous(
    name = "Pathogen detection rate (air samples)",
    sec.axis = sec_axis(~.x * 1, 
                        name = "Pathogen detection rate (clinical samples)",
                        labels = function(x) { format(round(x, 2), nsmall = 2) }
    )
  ) +
  theme(
    axis.title.y.left = element_text(colour = "red"),
    axis.line.y.left = element_line(colour = "red"),
    axis.ticks.y.left = element_line(colour = "red"),
    axis.text.y.left = element_text(colour = "red"),
    legend.justification = c(1,0),
    legend.position = c(1,0),
    strip.text = ggtext::element_markdown()
  )

ggsave("Output/figure_1.pdf", width=13.2, height=10)
ggsave("Output/figure_1.png", width=13.2, height=10)

# daily numbers instead of positivity rate
ggplot(respi_data_daily, aes(date, n)) +
  coord_cartesian(ylim=c(0, 11)) +
  facet_wrap(~pathogen, ncol=8) +
  geom_smooth(method = loess, colour="black", fill="black") +
  geom_point(colour="black", size=0.5) + 
  labs(y= "Daily detected cases", x = "Sampling date")

