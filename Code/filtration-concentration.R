# Linear regression analysis of the change in Ct value (virus concentration) from Wednesday to Monday and from Friday to Monday, on the three different locations #

library(lme4)
library(pbkrtest)

df_completed <- df_kabouters_selected_melt[which(df_kabouters_selected_melt$value_preab == 1),]

## Select data with virus detected every day of the week: select by location variable

df_weekly_present_3=data.frame(df_completed %>%
                               filter(location=="Child care group 3") %>%
                               with(table(week_pathogen,value_preab)))

df_weekly_present_2=data.frame(df_completed %>%
                                 filter(location=="Child care group 2") %>%
                                 with(table(week_pathogen,value_preab)))

df_weekly_present_1=data.frame(df_completed %>%
                                 filter(location=="Child care group 1") %>%
                                 with(table(week_pathogen,value_preab)))

## Select with frequency 3 and presence 1 per location
week_pathogen_loc3=df_weekly_present_3$week_pathogen[df_weekly_present_3$value_preab==1 &
                                                     df_weekly_present_3$Freq==3]

df_loc3_completed=df_completed[df_completed$location=="Child care group 3" &
                                                        df_completed$week_pathogen %in% week_pathogen_loc3,]

week_pathogen_loc2=df_weekly_present_2$week_pathogen[df_weekly_present_2$value_preab==1 &
                                                       df_weekly_present_3$Freq==3]

df_loc2_completed=df_completed[df_completed$location=="Child care group 2" &
                                 df_completed$week_pathogen %in% week_pathogen_loc2,]


week_pathogen_loc1=df_weekly_present_1$week_pathogen[df_weekly_present_1$value_preab==1 &
                                                       df_weekly_present_1$Freq==3]

df_loc1_completed=df_completed[df_completed$location=="Child care group 1" &
                                 df_completed$week_pathogen %in% week_pathogen_loc1,]


## renaming of variables ##

ct1 <- df_loc1_completed$value
ct2 <- df_loc2_completed$value
ct3 <- df_loc3_completed$value
ct1 <- as.numeric(ct1)
ct2 <- as.numeric(ct2)
ct3 <- as.numeric(ct3)
day1 <- df_loc1_completed$day
day2 <- df_loc2_completed$day
day3 <- df_loc3_completed$day
ID1 <- df_loc1_completed$week_pathogen
ID2 <- df_loc2_completed$week_pathogen
ID3 <- df_loc3_completed$week_pathogen

## random effects linear regression models to model change in Ct value on Wednesday/Friday as compared to Monday ##
### random effects on "clusters" of week&pathogen per location (ID1-ID2-ID3) 

model1 <- lmer(ct1 ~day1 + (1|ID1), REML=FALSE, data=df_loc1_completed)
summary(model1)
df.model1 <- get_Lb_ddf(model1, fixef(model1))
coefs1 <- data.frame(coef(summary(model1)))
coefs1$p.model1 <- 2 * (1-pt(abs(coefs1$t.value), df.model1))
coefs1
confint(model1)

model2a <- lmer(ct2~day2 + (1|ID2), REML=FALSE, data=df_loc2_completed) 
summary(model2a)
df.model2 <- get_Lb_ddf(model2a, fixef(model2a))
coefs2 <- data.frame(coef(summary(model2a)))
coefs2$p.model2a <- 2 * (1-pt(abs(coefs2$t.value), df.model2))
coefs2
confint(model2a)

model3a <- lmer(ct3~day3 + (1|ID3), REML=FALSE, data=df_loc3_completed)
summary(model3a)
df.model3 <- get_Lb_ddf(model3a, fixef(model3a))
coefs3 <- data.frame(coef(summary(model3a)))
coefs3$p.model2a <- 2 * (1-pt(abs(coefs3$t.value), df.model3))
coefs3
confint(model3a)
