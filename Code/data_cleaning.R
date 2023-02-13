library(tidyverse)

recode_pathogens <- function(data) {
  return (data %>%
            mutate( pathogen = recode(pathogen,
                                      `TaqPath COVID-19`="SARS-CoV-2 (TaqPath)",
                                      `Coronavirus SARS / COVID-19`="SARS-CoV-2 (respiratory panel)",
                                      `Adenovirus`="human adenovirus",
                                      `Bocavirus`="human bocavirus",
                                      `Chlamydophila pneumoniae`="*Chlamydia pneumoniae*",
                                      `Chlamydia psittaci`="*Chlamydia psittaci*",
                                      `Coronavirus 229E`="*Human coronavirus 229E*",
                                      "Coronavirus HKU-1"="*Human coronavirus HKU-1*",
                                      `Coronavirus NL63`="*Human coronavirus NL63*",
                                      `Coronavirus OC43`="*Human coronavirus OC43*",
                                      `Coxiella burnetti`="*Coxiella burnetii*",
                                      `Cytomegalovirus`="human cytomegalovirus",
                                      "Entero-/Rhinovirus"="human enterovirus (incl. rhinovirus)",
                                      `Enterovirus D68`="enterovirus D68",
                                      `Herpes simplex virus 1`="herpes simplex virus type 1",
                                      `Herpes simplex virus 2`="herpes simplex virus type 2",
                                      `Humaan metapneumovirus`="*Human metapneumovirus*",
                                      `Influenza A virus`="influenza A virus",
                                      `Influenza B virus`="influenza B virus",
                                      `Legionella pneumophila`="*Legionalla pneumophila*",
                                      `MERS`="MERS-CoV",
                                      `Mycoplasma pneumoniae`="*Mycoplasma pneumoniae*",
                                      `Parainfluenzavirus 1`="human parainfluenza virus 1",
                                      `Parainfluenzavirus 2`="human parainfluenza virus 2",
                                      `Parainfluenzavirus 3`="human parainfluenza virus 3",
                                      `Parainfluenzavirus 4`="human parainfluenza virus 4",
                                      `Parechovirus`="human parechovirus",
                                      `Pneumocystis jiroveci`="*Pneumocystis jirovecii*",
                                      `Respiratoir syncytieel virus`="respiratory syncytial virus",
                                      `Streptococcus pneumoniae`="*Streptococcus pneumoniae*",
                                      `Parainfluenzavirus`="human parainfluenza virus",
                                      `Other coronavirus`="Other human coronavirus"
            )
        )
  )
}

group_corona_parainf <- function(data) {
  return (
    data %>%
      mutate(
        
        # group non-COVID-19 coronaviruses
        pathogen=if_else(
          grepl( "oronavirus", pathogen, fixed = TRUE), # first letter removed: not case sensitive
          "Other coronavirus",
          pathogen),
        
        # group parainfluenza viruses
        pathogen=if_else(
          grepl( "arainfluenza", pathogen, fixed = TRUE), # first letter removed: not case sensitive
          "human parainfluenza",
          pathogen)
        
      ) %>%
      group_by(pathogen,sample_id) %>%
      slice_min(ct_value, n = 1, with_ties = FALSE) %>% # if multiple pathogens were detected from the same group, keep the lowest ct value
      ungroup
  )
}