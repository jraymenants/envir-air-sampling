library(RVAideMemoire)
library(rcompanion)
library(rstatix)
library(jmv)
library(ggplot)

data <- source.data

week_nr=c("2022W10","2022W08","2022W06","2022W07","2022W09","2022W11","2022W12","2022W13","2022W14",
          "2022W16","2022W17")

df_child_care=data[grepl("Child",data$location),]
group_names=unique(df_child_care$location)

df_child_care$value_preab=ifelse(is.na(df_child_care$ct_value),0,1)
df_child_care$week_pathogen=paste(df_child_care$year_week,df_child_care$pathogen,sep='_')

df_child_care$weekday_nl="maandag"
df_child_care$weekday_nl=ifelse(df_child_care$weekday=="Wednesday","woensdag",df_child_care$weekday_nl)
df_child_care$weekday_nl=ifelse(df_child_care$weekday=="Friday","vrijdag",df_child_care$weekday_nl)

df_child_care_select=df_child_care[df_child_care$year_week %in% week_nr,]

######## Cochran's q-test ########

### Get weeks with data for all three days
completed_weeks=df_child_care_select %>%
  group_by(location) %>%
  count('year_week') %>%
  filter(freq==198) %>% 
  pull(1)

df_completed=df_child_care_select[df_child_care_select$year_week %in% completed_weeks,]

## Filter based on location

sign_loc=list()
for (i in 1:length(group_names)) {
  df_qtest=df_completed %>% filter(location==group_names[i]) 
  res.q=cochran.qtest(value_preab ~ weekday_nl | week_pathogen, data = df_qtest)
  if (res.q$p.value<0.05) {
    sign_loc[i]=group_names[i]
  }
}

sign_loc[sapply(sign_loc, is.null)] <- NULL

## If significant pairwise test
for (k in 1:length(sign_loc)) {
  df_mcnemar=df_completed %>% filter(location==sign_loc[[k]])
  print(sign_loc[[k]])
  print(pairwiseMcnemar(value_preab ~ weekday_nl | week_pathogen,
                  data   = df_mcnemar,
                  test   = "mcnemar",
                  method = "holm",
                  digits = 3))
  }

######## Friedman test ########

### Select data with virus detected every day of the week: select by location variable

weekly_present=list()
for (n in 1:length(group_names)) {
  weekly_present[[n]]=data.frame(df_child_care_select %>%
                                 filter(location==group_names[n]) %>%
                                 with(table(week_pathogen,value_preab))) %>%
    filter(value_preab==1 & Freq==3) %>%
    pull(1)
}

for (t in 1:length(group_names)) {
  df_friedman[t]=paste("df_friedman_loc",t,sep = "_")
  fried_data=df_child_care_select[df_child_care_select$location==group_names[t] &
                                    df_child_care_select$week_pathogen %in% weekly_present[[t]],]
  assign(df_friedman[t], fried_data)
}

## Prepare Friedman table with three columns, one for a single day
list_friemand=list(df_friedman_loc_1,df_friedman_loc_2,df_friedman_loc_3)

for (m in 1:length(list_friemand)) {
  df_friedtest[m]=paste("df_friedman_test",m,sep = "_")
  x=list_friemand[[m]]
  y=data.frame(monday=x$ct_value[x$weekday_nl=="maandag"],
                          wednesday=x$ct_value[x$weekday_nl=="woensdag"],
                          friday=x$ct_value[x$weekday_nl=="vrijdag"])
  assign(df_friedtest[m], y)
}

## Perform test

res.loc1=anovaRMNP(df_friedman_test_1, measures = vars(monday,wednesday,friday))  
res.loc2=anovaRMNP(df_friedman_test_2, measures = vars(monday,wednesday,friday)) 
res.loc3=anovaRMNP(df_friedman_test_3, measures = vars(monday,wednesday,friday)) 

## If significant, post-hoc Dunn test

df_dunn_loc1=melt(df_friedman_test_1)  
df_dunn_loc2=melt(df_friedman_test_2)  
df_dunn_loc3=melt(df_friedman_test_3)  

dunn.test::dunn.test(df_dunn_loc1$value,df_dunn_loc1$variable)
dunn.test::dunn.test(df_dunn_loc2$value,df_dunn_loc2$variable)
dunn.test::dunn.test(df_dunn_loc3$value,df_dunn_loc3$variable)