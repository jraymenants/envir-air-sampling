library(RVAideMemoire)
library(rcompanion)
library(rstatix)
library(jmv)
library(ggplot)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggpubr)

data <- test_source

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
  count(year_week) %>%
  filter(n==66) %>% 
  distinct(year_week) %>%
  pull()

completed_weeks=unique(completed_weeks)[-9]

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

df_friedman=list()
for (t in 1:length(group_names)) {
  df_friedman[[t]]=df_child_care_select[df_child_care_select$location==group_names[t] &
                                    df_child_care_select$week_pathogen %in% weekly_present[[t]],]
  }

## Prepare Friedman table with three columns, one for a single day
df_friedman_test=list()
for (m in 1:length(df_friedman)) {
  x=df_friedman[[m]]
  df_friedman_test[[m]]=data.frame(monday=x$ct_value[x$weekday_nl=="maandag"],
                          wednesday=x$ct_value[x$weekday_nl=="woensdag"],
                          friday=x$ct_value[x$weekday_nl=="vrijdag"])
}

## Perform test

res.loc1=anovaRMNP(df_friedman_test[[1]], measures = vars(monday,wednesday,friday))  
res.loc2=anovaRMNP(df_friedman_test[[2]], measures = vars(monday,wednesday,friday)) 
res.loc3=anovaRMNP(df_friedman_test[[3]], measures = vars(monday,wednesday,friday)) 

## If significant, post-hoc Dunn test

df_dunn_loc1=melt(df_friedman_test[[1]])  
df_dunn_loc2=melt(df_friedman_test[[2]])  
df_dunn_loc3=melt(df_friedman_test[[3]])  

dunn.test::dunn.test(df_dunn_loc1$value,df_dunn_loc1$variable)
dunn.test::dunn.test(df_dunn_loc2$value,df_dunn_loc2$variable)
dunn.test::dunn.test(df_dunn_loc3$value,df_dunn_loc3$variable)

############ Figure ################
figure_names=c("location","value_preab","weekday_nl")

df_figure_average=aggregate(.~location+weekday_nl,df_completed[,colnames(df_completed) %in%
                                                                 figure_names],sum)
df_figure_average$value_weekly=df_figure_average$value_preab/length(completed_weeks)
df_figure_average$day_num=1
df_figure_average$day_num[df_figure_average$weekday_nl=="woensdag"]=2
df_figure_average$day_num[df_figure_average$weekday_nl=="vrijdag"]=3
df_figure_average$day="Monday"
df_figure_average$day[df_figure_average$day_num==2]="Wednesday"
df_figure_average$day[df_figure_average$day_num==3]="Friday"

df_figure_average$day <- factor(df_figure_average$day,
                                         levels = c('Monday','Wednesday','Friday'),ordered = TRUE)

df_figure_average$loc_new="Location 1"
df_figure_average$loc_new[df_figure_average$location=="Child care group 2"]="Location 2"
df_figure_average$loc_new[df_figure_average$location=="Child care group 3"]="Location 3"

## Significance 
stat.test <- tibble::tribble(~group1,~group2,~p.adj,
                             "1","2","*","1","3","*")

fig_barplot=ggplot(df_figure_average) +
  geom_bar(stat="identity",position = position_dodge(width = 0.6),
           aes(factor(day_num),value_weekly,color=loc_new),
           fill="white",width = 0.5) +
  geom_line(position = position_dodge(width = 0.6),
            aes(day_num,value_weekly,color=loc_new)) +
  geom_point(position = position_dodge(width = 0.6),
             aes(day_num,value_weekly,color=loc_new)) +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title.y = element_text(size=9)) +
  labs(x="",y="Average number of detected pathogens\n",color="") +
  scale_color_manual(values=c("#010221","#C43302","#0A7373")) +
  stat_pvalue_manual(stat.test, 
                     y.position = 8.8, step.increase = 0.1,
                     label = "p.adj",
                     color = "#C43302") +
  scale_x_discrete(breaks = c(1,2,3),label = c("Monday","Wednesday","Friday"))

##### Loc1 
plot_box_loc1=ggplot(df_friedman[[1]],aes(factor(weekday_nl,levels = c('maandag','woensdag','vrijdag')),ct_value)) +
  geom_boxplot() +
  ylim(40,25) +
  geom_hline(yintercept = median(df_friedman[[1]]$ct_value[df_friedman[[1]]$weekday_nl=="maandag"]),
             size=1,
             color="#010221") +
  scale_x_discrete(breaks=c('maandag','woensdag','vrijdag'),
                   labels=c("No filtration", "No filtration", "No filtration")) +
  theme_bw() +
  labs(y="Ct-value\n",x="",title="Location 1") +
  theme(legend.position = "none",
        plot.title = element_text(size=11),
        axis.title.y = element_text(size=9))

##### Loc2

plot_box_loc2=ggplot(df_friedman[[2]],aes(factor(weekday_nl,levels = c('maandag','woensdag','vrijdag')),ct_value)) +
  geom_boxplot() +
  ylim(40,25) +
  geom_hline(yintercept = median(df_friedman[[2]]$ct_value[df_friedman[[2]]$weekday_nl=="maandag"]),
             size=1,
             color="#C43302") +
  scale_x_discrete(breaks=c('maandag','woensdag','vrijdag'),
                   labels=c("No filtration", "48h filtration", "96h filtration")) +
  theme_bw() +
  labs(y="Ct-value\n",x="",title="Location 2") +
  theme(legend.position = "none",
        plot.title = element_text(size=11),
        axis.title.y = element_text(size=9))

##### Loc3

plot_box_loc3=ggplot(df_friedman[[3]],aes(factor(weekday_nl,levels = c('maandag','woensdag','vrijdag')),ct_value)) +
  geom_boxplot() +
  ylim(40,25) +
  geom_hline(yintercept = median(df_friedman[[3]]$ct_value[df_friedman[[3]]$weekday_nl=="maandag"]),
             size=1,
             color="#0A7373") +
  scale_x_discrete(breaks=c('maandag','woensdag','vrijdag'),
                   labels=c("No filtration", "48h filtration", "96h filtration")) +
  theme_bw() +
  labs(y="Ct-value\n",x="",title="Location 3") +
  theme(legend.position = "none",
        plot.title = element_text(size=11),
        axis.title.y = element_text(size=9))

### Panels

fig_line_box=ggarrange(plot_box_loc1,plot_box_loc2,plot_box_loc3,ncol=1,nrow=3)

ggarrange(fig_barplot,NULL,fig_line_box,ncol=3,widths = c(2,0.2,1.4),
          labels=c("A","","B"))

