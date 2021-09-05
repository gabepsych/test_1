gc(rm(list=ls()))
##### PACKAGES
library(dplyr)
library(jtools)
library(emmeans)
library(interactions)
library(lme4)
library(lmerTest)
library(psycho)
library(psych)
library(labelled)
library(ggplot2)
library(plyr)
library(effects)
library(extraoperators)
library(JWileymisc)
library(multilevelTools)
library(mixedup)
library(haven)
library(lattice)
library(sjPlot)
library(psych)
library(insight)
library(afex)
library(sjPlot)
library(visreg)

##### DATA PREP

#Load full dataset
full <- read_spss("COVID-19_DC ACROSS COUNTRIES_04.28.21.sav")

#Coerse COUNTRY to "factor"
full$COUNTRY <- labelled::to_factor(full$COUNTRY)
class(full$COUNTRY)
#COUNTRY is now a factor w/ 27 levels

#Select only Females and Males
full <- dplyr::filter(full, Gender < 3)
full$Gender <- as.numeric(full$Gender)
full$Gender <- full$Gender - 1
# 0 = male; 1 = female

#Add new items
morph1 <- function(x){
  x <-
    x %>% mutate(
      preDASSTOT = preDASSAnxTOT + preDASSDepTOT + preDASSStrTOT,
      postDASSTOT = postDASSAnxTOT + postDASSDepTOT + postDASSStrTOT,
      diff_DASSTOT = postDASSTOT - preDASSTOT,
      posDC = (Sdcemotp + Sdcprobp + Ddcp)/3
    )
  
  data_means1 <- plyr::ddply(x, "COUNTRY", summarize,
                             preDASSTOT_country = mean(preDASSTOT, na.rm = TRUE),
                             postDASSTOT_country = mean(postDASSTOT, na.rm = TRUE),
                             posDC_country = mean(posDC, na.rm = TRUE),
                             negDC_country = mean(Ndcp, na.rm = TRUE))
  
  x <- merge(x, data_means1, by = "COUNTRY")
  
  x <-
    x %>% mutate(
      preDASSTOT_country_c = scale(preDASSTOT_country, center = TRUE, scale = FALSE),
      postDASSTOT_country_c = scale(postDASSTOT_country, center = TRUE, scale = FALSE),
      posDC_country_c = scale(posDC_country, center = TRUE, scale = FALSE),
      negDC_country_c = scale(negDC_country, center = TRUE, scale = FALSE),
      
      preDASSTOT_dev = preDASSTOT - preDASSTOT_country,
      postDASSTOT_dev = postDASSTOT - postDASSTOT_country,
      posDC_dev = posDC - posDC_country,
      negDC_dev = Ndcp - negDC_country,
      
      preDASSTOT_c = scale(preDASSTOT, center = TRUE, scale = TRUE),
      Stresscoms_c = scale(Stresscoms, center = TRUE, scale = TRUE),
      Gender_c = scale(Gender, center = TRUE, scale = FALSE)
    )
}

full1<- morph1(full)
write.csv(full1,"COVID-19_DC ACROSS COUNTRIES_05.08.21.csv")
write_sav(full1,"COVID-19_DC ACROSS COUNTRIES_05.08.21.sav")
data <- morph1(full)

names(data)  
#select only the variables we need
names <- c("COUNTRY", "ID", "Gender", "Stresscoms","preDASSTOT", "postDASSTOT", "diff_DASSTOT",
           "PRQCTOT", "Stresscoms_c", "preDASSTOT_c", "postDASSTOT_country_c", "postDASSTOT_dev", "posDC_country_c",
           "posDC_dev", "negDC_country_c", "negDC_dev")
data <- dplyr::select(data, all_of(names))

# #Visualize raw Relationship Satisfaction
# ggplot(data=data, aes(x=PRQCTOT)) +
#   geom_histogram(fill="white", color="black",bins=50) +
#   labs(x = "Relationship Satisfaction")
# 
# #Visualize raw pre-COVID 19 DASS
# ggplot(data=data, aes(x=preDASSTOT)) +
#   geom_histogram(fill="white", color="black",bins=50) +
#   labs(x = "DASS prior to COVID-19 Travel Ban")
# 
# #Visualize raw post-COVID 19 DASS
# ggplot(data=data, aes(x=postDASSTOT)) +
#   geom_histogram(fill="white", color="black",bins=50) +
#   labs(x = "DASS after COVID-19 Travel Ban")
# 
# #Visualize raw change in DASS
# ggplot(data=data, aes(x=diff_DASSTOT)) +
#   geom_histogram(fill="white", color="black",bins=30) +
#   labs(x = "Change in DASS pre- to post-COVID-19 Travel Ban")

#######
# Hypothesis #1 - Do self-reported symptoms of psychological distress (i.e., depression, anxiety, and stress) 
  #increase from pre- to post-COVID-19 travel bans.

#Dependent T-test
fit0 <- t.test(data$postDASSTOT, data$preDASSTOT, paired = TRUE, alternative = "two.sided")
fit0

#Dependent T-test in regression framework (no-random intercept)
fit1 <- lm(diff_DASSTOT ~ 1, data = data)
summary(fit1) 
summ(fit1)

#Allowing intercepts to vary across countries
fit2 <- lmer(diff_DASSTOT ~ 1 + (1|COUNTRY), data = data)
fit2.sum <- summary(fit2)
summ(fit2, digits = 4)
confint(fit2)

#On average, distress increased from pre- to post-COVID-19 travel bans
#even after accounting for between-country differences

####Plotting Random Intercepts

#Extract random intercepts
rand1 <- as.data.frame(mixedup::extract_random_effects(fit2))
#Extract fixed intercept
fixedef1 <- summary(fit2)$coefficients[1,1]

#Add fixed intercept to random intercepts and their 95% CIs
rand1 <- within(rand1,{
  diff.pred <- rand1$value+fixedef1
  lower_2.5.pred <- rand1$lower_2.5+fixedef1
  upper_97.5.pred <- rand1$upper_97.5+fixedef1
})

#Reorder effects in descending order
rand1 <- within(rand1, {
  group.2 <- factor(group, 
                      levels = group[order(value)])
})

write.csv(rand1, "h1_table.csv")

#####################################################################################################
#Hypothesis #2 - Does Post-COVID-19 distress predict Relationship Quality
                #after controlling for pre-COVID-19 distress, gender, and stress communication.
                #
#n = 13749
#Omit missing data to aid in model comparisons
data <- na.omit(data)
#n = 11968
write.csv(data, "COVID-19_DC ACROSS COUNTRIES_05.08.21_NO_MISSING.csv")
write_sav(data, "COVID-19_DC ACROSS COUNTRIES_05.08.21_NO_MISSING.sav")

#Unconditional Means model
fit3 <- lmer(PRQCTOT ~ 1 + (1|COUNTRY),
             data = data)
fit3.sum <- summary(fit3)
summ(fit3, digits = 5)

# #Visual exploration
# tmp <- multilevelTools::meanDecompose(PRQCTOT ~ COUNTRY, data = data)
# str(tmp, nchar.max = 30)

# #Test Between-Country Distribution
# plot(testDistribution(tmp[["PRQCTOT by COUNTRY"]]$X,
#                       extremevalues = "theoretical", ev.perc = .001),
#      varlab = "Between-Country Relationship Quality")
# 
# #Test Within-Country Distribution
# plot(testDistribution(tmp[["PRQCTOT by residual"]]$X,
#                       extremevalues = "theoretical", ev.perc = .001),
#      varlab = "Between-Country Relationship Quality")
# #Within-country relationship quality is left-skewed

#####Include fixed effects
fit3.1 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms +
               postDASSTOT_country_c + postDASSTOT_dev +
                (1|COUNTRY), data = data)

fit3.1.summ <- summary(fit3.1)
summ(fit3.1)
modelPerformance(fit3.1)
anova(fit3, fit3.1, refit = TRUE) #fit improves. Refit with ML because I am comparing fixed effects

#####Include random effect for postDASSTOT_dev

# fit3.2 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                (1 + postDASSTOT_dev|COUNTRY), data = data)
# #Initial model failed to converge

##Try different optimizations techniques

NEWUOAcontrol <- lmerControl(optCtrl = list(
  algorithm = "NLOPT_LN_NEWUOA",
  xtol_abs = 1e-12,
  ftol_abs = 1e-12))

PRAXIScontrol <- lmerControl(optCtrl = list(
  algorithm = "NLOPT_LN_PRAXIS",
  xtol_abs = 1e-12,
  ftol_abs = 1e-12))

NELDERMEADcontrol <- lmerControl(optCtrl = list(
  algorithm = "NLOPT_LN_NELDERMEAD",
  xtol_abs = 1e-12,
  ftol_abs = 1e-12))

COBYLAcontrol <- lmerControl(optCtrl = list(
  algorithm = "NLOPT_LN_COBYLA",
  xtol_abs = 1e-12,
  ftol_abs = 1e-12))

SUBPLEXcontrol <- lmerControl(optCtrl = list(
  algorithm = "NLOPT_LN_SBPLX",
  xtol_abs = 1e-12,
  ftol_abs = 1e-12))

fit3.2 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                (1 + postDASSTOT_dev|COUNTRY), data = data,
               control = SUBPLEXcontrol)
#model converged; no warnings

fit3.21 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                 (1 + postDASSTOT_dev|COUNTRY), data = data,
               control = PRAXIScontrol)
#model converged; no warnings

fit3.22 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                 (1 + postDASSTOT_dev|COUNTRY), data = data,
               control = NELDERMEADcontrol)
#model converged; no warnings

summ(fit3.2, digits = 5)
summ(fit3.21, digits = 5)
summ(fit3.22, digits = 5)
#estimates and random effects are the exact same

anova(fit3.1, fit3.2, refit = FALSE) #fit improves. Do not refit using ML

#Include random effect for pre-COVID distress
fit3.3 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                 (1 + postDASSTOT_dev + preDASSTOT_c|COUNTRY), data = data, control = SUBPLEXcontrol)
#model converged; no warnings

anova(fit3.2, fit3.3, refit = FALSE) #fit improves, barely: 0.02019 *

#Include random effect for Gender
fit3.4 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                 (1 + postDASSTOT_dev + preDASSTOT_c + Gender|COUNTRY), data = data, control = SUBPLEXcontrol)
#model converged; singular fit
anova(fit3.3, fit3.4, refit = FALSE) #fit improves

#Include random effect for Stress Communication
fit3.5 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                  (1 + postDASSTOT_dev + Gender + Stresscoms|COUNTRY), data = data, control = SUBPLEXcontrol)
#model converged; no warning
anova(fit3.1, fit3.5, refit = FALSE) #fit improves

# #Include random effect for Stress Communication
# fit3.51 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                  (1 + postDASSTOT_dev + preDASSTOT_c + Gender + Stresscoms|COUNTRY), data = data, control = SUBPLEXcontrol)
# #singular fit


fit3.5.sum <- summary(fit3.5)
summ(fit3.5, digits = 5)
fit3.5.sum
write.csv(fit3.5.sum$coefficients, "fit3.5coef.csv")
write.csv(fit3.5.sum$varcor, "fit3.5_varcorr.csv")
#confint(fit3.5)
confint(fit3.5, method = "Wald")

fit3.5.sum$coefficients
summary(rePCA(fit3.5))
#####Explore final model (fit3.5)

# #Model Diagnostics
# md <- JWileymisc::modelDiagnostics(fit3.5, ev.perc = 0.001)
# plot(md, ask = FALSE, ncol = 2, nrow = 3)
# 
# mvextreme <- subset(md$extremeValues,
#                     EffectType == "Multivariate Random Effect COUNTRY")
# mvextreme
# #no extreme values

fit3.5.sum$vcov
fit3.5.sum$varcor
12.84130^2
0.12419^2
2.29072^2
2.49962^2
14.65842^2
fit3.5.sum$optinfo
fit3.5.sum$sigma
confint(fit3.5, method = "Wald")

insight::get_variance(fit3.5)

#########Plot Random Effects for postDASSTOT_dev

#Extract Random Effects
rand2 <- as.data.frame(mixedup::extract_random_effects(fit3.5))

#Filter for randome effects of postDASSTOT_dev
rand2.dass <- dplyr::filter(rand2, effect == "postDASSTOT_dev")

#Extract fixed effect of postDASSTOT_dev
fixedef2 <- summary(fit3.5)$coefficients[6,1]

#Add fixed effect of postDASSTOT_dev to random effect and 95% CIs
rand2.dass <- within(rand2.dass,{
  postDASSTOT_dev.pred <- rand2.dass$value+fixedef2
  lower_2.5.pred <- rand2.dass$lower_2.5+fixedef2
  upper_97.5.pred <- rand2.dass$upper_97.5+fixedef2
})

#Reorder effects in descending order
rand2.dass <- within(rand2.dass, {
  group.2 <- factor(group, 
                    levels = group[order(value)])
})

write.csv(rand2, "h2_table1.csv")
############################Plot random effects for Gender

#Filter for random effects of postDASSTOT_dev
rand2.gender <- dplyr::filter(rand2, effect == "Gender")

#Extract fixed effect of postDASSTOT_dev
fixedef2.gender <- summary(fit3.5)$coefficients[3,1]

#Add fixed effect of postDASSTOT_dev to random effect and 95% CIs
rand2.gender <- within(rand2.gender,{
  Gender.pred <- rand2.gender$value+fixedef2.gender
  lower_2.5.pred <- rand2.gender$lower_2.5+fixedef2.gender
  upper_97.5.pred <- rand2.gender$upper_97.5+fixedef2.gender
})

#Reorder effects in descending order
rand2.gender <- within(rand2.gender, {
  group.2 <- factor(group, 
                    levels = group[order(value)])
})

##################################Plot random effects for Stresscoms

#Filter for random effects of stresscoms
rand2.stress <- dplyr::filter(rand2, effect == "Stresscoms")

#Extract fixed effect of stresscoms
fixedef2.stress <- summary(fit3.5)$coefficients[4,1]

#Add fixed effect of stresscoms to random effect and 95% CIs
rand2.stress <- within(rand2.stress,{
  Stresscoms.pred <- rand2.stress$value+fixedef2.stress
  lower_2.5.pred <- rand2.stress$lower_2.5+fixedef2.stress
  upper_97.5.pred <- rand2.stress$upper_97.5+fixedef2.stress
})

#Reorder effects in descending order
rand2.stress <- within(rand2.stress, {
  group.2 <- factor(group, 
                    levels = group[order(value)])
})


#####Plot scatterplot matrix of partial reg line for RelQuality on postDASSTOT_dev

lm1 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms,
          data = data)

lm2 <- lm(postDASSTOT ~ 1 + preDASSTOT + Gender + Stresscoms,
          data = data)

resid.lm1 <- residuals(lm1)
resid.lm2 <- residuals(lm2)

# ggplot(data, aes(x=resid.lm2, y=resid.lm1)) + 
#   geom_point(size=0.05) +   
#   geom_jitter(position = position_jitter(w = 0.4, h = 0.4)) +   
#   geom_smooth(method="lm", se=FALSE,size = 1) + 
#   theme_bw() +
#   theme(axis.title = element_text(size = 15, face = "bold"),
#         axis.text = element_text(size = 11, face = "bold"),
#         strip.text = element_text(size = 11, face = "bold")) +
#   theme(axis.line = element_line(color = 'black')) + 
#   labs(y = "Relationship Quality\n",
#        x = "\nPost-COVID-19 Distress") + 
#   facet_wrap( ~COUNTRY)
# 
# ggplot(data, aes(x=diff_DASSTOT, y=PRQCTOT)) + 
#   geom_point(size=0.05) +   
#   geom_jitter(position = position_jitter(w = 0.4, h = 0.4)) +   
#   geom_smooth(method="lm", se=FALSE,size = 1) + 
#   theme_bw() +
#   theme(axis.title = element_text(size = 15, face = "bold"),
#         axis.text = element_text(size = 11, face = "bold"),
#         strip.text = element_text(size = 11, face = "bold")) +
#   theme(axis.line = element_line(color = 'black')) + 
#   labs(y = "Relationship Quality\n",
#        x = "\nPost-COVID-19 Distress") + 
#   facet_wrap( ~COUNTRY)


tiff('plot5_scattermatrix.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
ggplot(data, aes(x=diff_DASSTOT, y=PRQCTOT)) + 
  geom_point(size=0.05) +   
  geom_jitter(position = position_jitter(w = 0.4, h = 0.4)) +   
  geom_smooth(method="lm", se=FALSE,size = 1) + 
  theme_bw() +
  theme(axis.title = element_text(size = 15, face = "bold"),
        axis.text = element_text(size = 11, face = "bold"),
        strip.text = element_text(size = 11, face = "bold")) +
  theme(axis.line = element_line(color = 'black')) + 
  labs(y = "Relationship Quality\n",
       x = "Difference in Pre- and Post-COVID-19 Distress") + 
  facet_wrap( ~COUNTRY)
dev.off()

vis1 <- visreg(fit = fit3.5, "postDASSTOT_dev", by = "COUNTRY", data = data, type = "conditional",
          re.form = ~(1 + postDASSTOT_dev + Gender + Stresscoms|COUNTRY), plot = FALSE)

plot(vis1, ylab = "Relationship Quality", xlab = "Post-COVID-19 Distress\n(between-person)",
     layout=c(4,6), ylim = c(75,150))

###############################################################################################################
#Hypothesis #3 - The effect of country-level distress and person-level distress will be
#moderated by positive dyadic coping

fit4 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
                  (1 + postDASSTOT_dev + Gender + Stresscoms|COUNTRY), data = data, control = SUBPLEXcontrol)
#very large eigenvalue warning
#anova(fit3.5, fit4) #fit improved

# #Add random effect for PosDC_dev
# fit4.1 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                (1 + postDASSTOT_dev + Gender + Stresscoms + posDC_dev|COUNTRY), data = data, control = NELDERMEADcontrol)
# #failed to converge
# 
# #Remove random effect for Stresscoms and add random effect for PosDC_dev
# fit4.2 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                (1 + postDASSTOT_dev + Gender + posDC_dev|COUNTRY), data = data, control = strictControl)
# anova(fit4, fit4.2, refit = FALSE) #fit improves

# #Try adding random effect for interaction term
# fit4.3 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                  posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                  (1 + postDASSTOT_dev + Gender + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data, control = strictControl)
# #failed to converge


# #Try adding random effect for interaction term and removing random effect for Gender
# fit4.4_PRAXIS <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                          posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                          (1 + postDASSTOT_dev + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
#                        control = PRAXIScontrol)
# #failed to converge - negative eigenvalue
# 
# fit4.4_NELDERMEAD <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                         posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                         (1 + postDASSTOT_dev + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
#                       control = NELDERMEADcontrol)
# #failed to converge - negative eigenvalue
# 
# fit4.4_COBYLA <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                             posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                             (1 + postDASSTOT_dev + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
#                           control = COBYLAcontrol)
# #failed to converge - 2 negative eigenvalues
# 
# fit4.4_NEWUOUA <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                          posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                          (1 + postDASSTOT_dev + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
#                        control = NEWUOAcontrol)
# #singular fit

fit4.4_SBPLX <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
                        posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
                        (1 + postDASSTOT_dev + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
                      control = SUBPLEXcontrol)
#warning - very large eigenvalue
anova(fit3.5, fit4.4_SBPLX, refit = FALSE) #fit improves

# fit4.41_SBPLX <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                        posDC_country_c + posDC_dev + posDC_dev*postDASSTOT_dev +
#                        (1 + postDASSTOT_dev + Gender + Stresscoms + posDC_dev + posDC_dev*postDASSTOT_dev|COUNTRY), data = data,
#                      control = SUBPLEXcontrol)
# #warning - did not converge negative eigenvalue

# summ(fit4.4_NELDERMEAD, digits = 5)
# summ(fit4.4_NEWUOUA, digits = 5)
# summ(fit4.4_COBYLA, digits = 5)
# summ(fit4.4_PRAXIS, digits = 5)
summ(fit4.4_SBPLX, digits = 5)
fit4.4_SBPLX.sum <- summary(fit4.4_SBPLX)
confint(fit4.4_SBPLX, method = "Wald")

summary(rePCA(fit4.4_SBPLX))
fit4.4_SBPLX.sum$varcor
5.025329^2
0.070910^2
2.915363^2
0.082199^2
12.753348^2

write.csv(fit4.4_SBPLX.sum$coefficients, "fit4.4_SBPLX_coef.csv")
write.csv(fit4.4_SBPLX.sum$varcor, "fit4.4_SBPLX_varcorr.csv")

# md1 <- JWileymisc::modelDiagnostics(fit4.4_SBPLX, ev.perc = 0.001)
# mvextreme1 <- subset(md1$extremeValues,
#                     EffectType == "Multivariate Random Effect COUNTRY")
# mvextreme1
#no extreme values


#Plot random effects for interaction term
rand3 <- as.data.frame(mixedup::extract_random_effects(fit4.4_SBPLX))
rand3

rand3.interact <- dplyr::filter(rand3, effect == 'postDASSTOT_dev:posDC_dev')

#Extract fixed effect of postDASSTOT_dev:posDC_dev
fixedef3.interact <- summary(fit4.4_SBPLX)$coefficients[9,1]

#Add fixed effect of postDASSTOT_dev to random effect and 95% CIs
rand3.interact <- within(rand3.interact,{
  Interact.pred <- rand3.interact$value+fixedef3.interact
  lower_2.5.pred <- rand3.interact$lower_2.5+fixedef3.interact
  upper_97.5.pred <- rand3.interact$upper_97.5+fixedef3.interact
})

#Reorder effects in descending order
rand3.interact <- within(rand3.interact, {
  group.2 <- factor(group, 
                    levels = group[order(value)])
})

#####Probe interactions using Fit4.4_SBPLX

# posDCa <- mean(data$posDC) + sd(data$posDC)
# posDC <- mean(data$posDC)
# posDCb <- mean(data$posDC) - sd(data$posDC)

# posDC_country_c_A <- mean(data$posDC_country_c) + sd(data$posDC_country_c)
# posDC_country_c_m <- mean(data$posDC_country_c)
# posDC_country_c_B <- mean(data$posDC_country_c) - sd(data$posDC_country_c)
# 
posDC_dev_A <- mean(data$posDC_dev, na.rm = TRUE) + (sd(data$posDC_dev, na.rm = TRUE))
posDC_dev_m <- mean(data$posDC_dev, na.rm = TRUE)
posDC_dev_B <- mean(data$posDC_dev, na.rm = TRUE) - (sd(data$posDC_dev, na.rm = TRUE))

# mylist0 <- list(posDC = c(posDCa, posDC, posDCb))
# mylist1 <- list(posDC_country_c=c(posDC_country_c_A, posDC_country_c_m, posDC_country_c_B))
mylist2 <- list(posDC_dev=c(posDC_dev_A, posDC_dev_m, posDC_dev_B))

simple_slopes1 <- emtrends(fit4.4_SBPLX, ~posDC_dev, var="postDASSTOT_dev", at=mylist2,
                           lmer.df = c("satterthwaite"))
summary(simple_slopes1)

###############################################################################################################
#Hypothesis #4 - The effect of country-level distress and person-level distress will be
#moderated by negative dyadic coping

fit5 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
               negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
               (1 + postDASSTOT_dev + Gender + Stresscoms|COUNTRY), data = data, control = SUBPLEXcontrol)
#no warning
#anova(fit3.5, fit5) #fit improved

# #Add random effect for NegDC_dev
# fit5.1 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
#                (1 + postDASSTOT_dev + Gender + Stresscoms + negDC_dev|COUNTRY), data = data, control = NELDERMEADcontrol)
# #failed to converge
# 
#Remove random effect for Stresscoms and Gender and add random effect for NegDC_dev and interaction term
fit5.2 <- lmer(PRQCTOT ~ 1 + preDASSTOT_c + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
               negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
               (1 + postDASSTOT_dev + negDC_dev + negDC_dev*postDASSTOT_dev|COUNTRY), data = data, control = SUBPLEXcontrol)
#large eigenvalue warning
anova(fit3.5, fit5.2, refit = FALSE) #fit improves

# # #Try adding random effect for interaction term
# fit5.3 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                  negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
#                  (1 + postDASSTOT_dev + Gender + negDC_dev + negDC_dev*postDASSTOT_dev|COUNTRY), data = data, control = NELDERMEADcontrol)
# #failed to converge

# fit5.3 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                  negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
#                  (1 + postDASSTOT_dev + Gender + negDC_dev + negDC_dev*postDASSTOT_dev|COUNTRY), data = data, control = SUBPLEXcontrol)
# #large eigenvalue
# anova(fit5.2, fit5.3, refit = FALSE) #Fit improves
# 
# 
# fit5.4 <- lmer(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT_country_c + postDASSTOT_dev +
#                  negDC_country_c + negDC_dev + negDC_dev*postDASSTOT_dev +
#                  (1 + postDASSTOT_dev + Gender + Stresscoms + negDC_dev + negDC_dev*postDASSTOT_dev|COUNTRY), data = data, control = SUBPLEXcontrol)
# #large eigenvalue
# anova(fit5.3, fit5.4, refit = FALSE) #fit improves

summ(fit5.2, digits = 5)
confint(fit5.2, method = "Wald")
fit5.2.sum <- summary(fit5.2)
write.csv(as.data.frame(fit5.2.sum$varcor), "fit5.2_varcorr.csv")
write.csv(fit5.2.sum$coefficients, "fit5.2_coef.csv")

fit5.2.sum$varcor
4.954254^2
0.077656^2
4.175561^2
0.085733^2
13.634226^2

# md2 <- JWileymisc::modelDiagnostics(fit5.2, ev.perc = 0.001)
# mvextreme2 <- subset(md2$extremeValues,
#                      EffectType == "Multivariate Random Effect COUNTRY")
# mvextreme2
# #no extreme values

#Plot random effects for interaction term
rand4 <- as.data.frame(mixedup::extract_random_effects(fit5.2))
rand4

rand4.interact <- dplyr::filter(rand4, effect == 'postDASSTOT_dev:negDC_dev')

#Extract fixed effect of negDC_dev:postDASSTOT_dev
fixedef4.interact <- summary(fit5.2)$coefficients[9,1]

#Add fixed effect of postDASSTOT_dev to random effect and 95% CIs
rand4.interact <- within(rand4.interact,{
  Interact.pred <- rand4.interact$value+fixedef4.interact
  lower_2.5.pred <- rand4.interact$lower_2.5+fixedef4.interact
  upper_97.5.pred <- rand4.interact$upper_97.5+fixedef4.interact
})

#Reorder effects in descending order
rand4.interact <- within(rand4.interact, {
  group.2 <- factor(group, 
                    levels = group[order(value)])
})


#####Probe interactions using Fit5.2

# negDCa <- mean(data$negDC) + sd(data$negDC)
# negDC <- mean(data$negDC)
# negDCb <- mean(data$negDC) - sd(data$negDC)

# negDC_country_c_A <- mean(data$negDC_country_c) + sd(data$negDC_country_c)
# negDC_country_c_m <- mean(data$negDC_country_c)
# negDC_country_c_B <- mean(data$negDC_country_c) - sd(data$negDC_country_c)
# 
negDC_dev_A <- mean(data$negDC_dev, na.rm = TRUE) + (sd(data$negDC_dev, na.rm = TRUE))
negDC_dev_m <- mean(data$negDC_dev, na.rm = TRUE)
negDC_dev_B <- mean(data$negDC_dev, na.rm = TRUE) - (sd(data$negDC_dev, na.rm = TRUE))

# mylist0 <- list(negDC = c(negDCa, negDC, negDCb))
# mylist1 <- list(negDC_country_c=c(negDC_country_c_A, negDC_country_c_m, negDC_country_c_B))
mylist2 <- list(negDC_dev=c(negDC_dev_A, negDC_dev_m, negDC_dev_B))

simple_slopes2 <- emtrends(fit5.2, ~negDC_dev, var="postDASSTOT_dev", at=mylist2,
                           lmer.df = c("satterthwaite"))
summary(simple_slopes2)


# #Reliabilities
# names2 <- c("COUNTRY", "ID", "Gender", "preDASSDepTOT", "preDASSAnxTOT", "preDASSStrTOT",
#             "postDASSDepTOT", "postDASSAnxTOT", "postDASSStrTOT",
#             "PRQCTOT", "Stresscoms", "Sdcemotp", "Sdcprobp", "Ddcp", "Ndcp",
#             "DCI22","DCI25","DCI26","DCI27")
# test <- dplyr::select(full, all_of(names2))
# test$COUNTRY <- labelled::to_factor(test$COUNTRY)
# test <- dplyr::select(test,COUNTRY,DCI22,DCI25,DCI26,DCI27)
# test1 <- dplyr::select(test,DCI22,DCI25,DCI26,DCI27)
# test <- dplyr::filter(test, COUNTRY == "USA")
# test <- dplyr::select(test, -COUNTRY)
# omega(test)
# alpha(test)
# library(ltm)
# cronbach.alpha(test, na.rm = TRUE)
# cronbach.alpha(test1, na.rm = TRUE)
# install.packages("coefficientalpha")
# library(coefficientalpha)
# x<-coefficientalpha::alpha(test)
# coefficientalpha::omega(test)
# coefficientalpha::plot.alpha(x)
# 
# canada <- dplyr::filter(data, COUNTRY == "CANADA")
# usa <- dplyr::filter(data, COUNTRY == "USA")
# table(canada$preDASSTOT)
# table(data$preDASSTOT)
# ftable(canada$preDASSTOT)
# ftable(usa$preDASSTOT)
# hist(canada$preDASSTOT)
# hist(usa$preDASSTOT)
# 
# hist(canada$preDASSAnxTOT)
# hist(canada$preDASSDepTOT)
# hist(canada$preDASSStrTOT)

#######################################################################################
#FINAL PLOTS


p1 <- ggplot(rand1, aes(x = diff.pred, y = group.2,
                        xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef1, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "\nDifference in Pre- and Post-COVID-19\nDistress Random Intercepts",
       y = "Country")

p2 <- ggplot(rand2.dass, aes(x = postDASSTOT_dev.pred, y = group.2,
                             xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef2, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  # scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "Post-COVID-19 Distress (within-country)\nRandom Effects",
       y = "Country")

p3 <- ggplot(rand2.gender, aes(x = Gender.pred, y = group.2,
                               xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef2.gender, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  # scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "Gender\nRandom Effects",
       y = "Country")

p4<- ggplot(rand2.stress, aes(x = Stresscoms.pred, y = group.2,
                              xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef2.stress, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  # scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "Stress Communication\nRandom Effects",
       y = "Country")

library(ggpubr)

tiff('plot1_randomeffects_1.tiff', units="in", width=10, height=5, res=2400, compression = 'lzw')
ggarrange(p1,p2, labels = c("A","B"))
dev.off()

tiff('plot1_randomeffects_2.tiff', units="in", width=10, height=5, res=2400, compression = 'lzw')
ggarrange(p3,p4, labels = c("C","D"))
dev.off()

##################

p5 <- ggplot(rand3.interact, aes(x = Interact.pred, y = group.2,
                                 xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef3.interact, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  # scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "Post-COVID-19 Distress by Positive DC\nRandom Effects",
       y = "Country")

p6 <- interact_plot(fit4.4_SBPLX, pred = postDASSTOT_dev, modx = posDC_dev,
                    fit4.4_SBPLX = TRUE,
                    #colors = "Greys",
                    plot.points = FALSE,
                    line.thickness = 1.5,
                    interval = TRUE, #95% confidence prediction interval,
                    int.type = "confidence",
                    #partial.residuals = TRUE,
                    #linearity.check = TRUE,
                    main.title = "Post-COVID-19 Distress by Positive DC",
                    legend.main = "Positive DC",
                    x.label = "Post-COVID-19 Distress",
                    y.label = "Relationship Quality"
) + ylim(80,120)

p7 <- ggplot(rand4.interact, aes(x = Interact.pred, y = group.2,
                                 xmin = lower_2.5.pred, xmax = upper_97.5.pred)) +
  geom_vline(xintercept = fixedef4.interact, linetype = "dashed", col = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "solid", col = "blue", size = 0.5) +
  geom_errorbarh() +
  geom_point() +
  # scale_x_continuous(breaks = seq(-3, 15, 1))+
  theme_bw() +
  labs(x = "Post-COVID-19 Distress by Negative DC\nRandom Effects",
       y = "Country")

p8 <- interact_plot(fit5.2, pred = postDASSTOT_dev, modx = negDC_dev,
                    fit5.2 = TRUE,
                    #colors = "Greys",
                    plot.points = FALSE,
                    line.thickness = 1.5,
                    interval = TRUE, #95% confidence prediction interval,
                    int.type = "confidence",
                    #partial.residuals = TRUE,
                    #linearity.check = TRUE,
                    main.title = "Post-COVID-19 Distress by Negative DC",
                    legend.main = "Negative DC",
                    x.label = "Post-COVID-19 Distress",
                    y.label = "Relationship Quality"
) + ylim(80,120)

tiff('plot2_posDC.tiff', units="in", width=10, height=5, res=2400, compression = 'lzw')
ggarrange(p5,p6, labels = c("A", "B"))
dev.off()

tiff('plot2_negDC.tiff', units="in", width=10, height=5, res=2400, compression = 'lzw')
ggarrange(p7,p8, labels = c("C", "D"))
dev.off()

##############################################
#Analyses for Australia, Romania, and Portugal

#Australia
aust <- dplyr::filter(full, COUNTRY == "AUSTRALIA")
#Make sure DC measures are computed correctly
aust$Sdcprobp <- (aust$DCI8 + aust$DCI13)/2
aust$Sdcemotp <- (aust$DCI5 + aust$DCI6)/2
aust$Ddcp <-     (aust$DCI12 +aust$DCI14)/2
aust$Ndcp <-     (aust$DCI7 + aust$DCI10 + aust$DCI11 + aust$DCI15)/4
aust <- morph1(aust)

#Romania
rom <- dplyr::filter(full, COUNTRY == "ROMANIA")
#Compute PRQCTOT
rom <- rom %>% mutate(
  PRQCTOT = rowSums(dplyr::across(PRQC1:PRQC18), na.rm = TRUE)
)
rom <- morph1(rom)

#Portugal
port <- dplyr::filter(full, COUNTRY == "PORTUGAL")
#Compute PRQCTOT
port <- port %>% mutate(
  PRQCTOT = rowSums(dplyr::across(PRQC1:PRQC18), na.rm = TRUE)
)
port <- morph1(port)

names1 <- c("COUNTRY", "ID", "Gender", "Stresscoms","preDASSTOT", "postDASSTOT", "diff_DASSTOT",
           "PRQCTOT","posDC", "Ndcp")

aust <- dplyr::select(aust, all_of(names1))
rom <- dplyr::select(rom, all_of(names1))
port <- dplyr::select(port, all_of(names1))

####Australia
##H2
fit_aust1 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + postDASSTOT, data = aust)
fit_aust1.summ <- summary(fit_aust1)
summ(fit_aust1, digits = 5)
confint(fit_aust1)
##H3
#posDC
fit_aust2 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + postDASSTOT + posDC + postDASSTOT*posDC, data = aust)
fit_aust2.summ <- summary(fit_aust2)
summ(fit_aust2, digits = 5)
confint(fit_aust2)
#negDC
fit_aust3 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + postDASSTOT + Ndcp + postDASSTOT*Ndcp, data = aust)
fit_aust3.summ <- summary(fit_aust3)
summ(fit_aust3, digits = 5)
confint(fit_aust3)

###Romania
##H2
fit_rom1 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT, data = rom)
fit_rom1.summ <- summary(fit_rom1)
summ(fit_rom1, digits = 5)
##H3
#posDC
fit_rom2 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT + posDC + postDASSTOT*posDC, data = rom)
fit_rom2.summ <- summary(fit_rom2)
summ(fit_rom2, digits = 5)
confint(fit_rom2)
#negDC
fit_rom3 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT + Ndcp + postDASSTOT*Ndcp, data = rom)
fit_rom3.summ <- summary(fit_rom3)
summ(fit_rom3, digits = 5)

###Portugal
##H2
fit_port1 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT, data = port)
fit_port1.summ <- summary(fit_port1)
summ(fit_port1, digits = 5)
##H3
#posDC
fit_port2 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT + posDC + postDASSTOT*posDC, data = port)
fit_port2.summ <- summary(fit_port2)
summ(fit_port2, digits = 5)
#negDC
fit_port3 <- lm(PRQCTOT ~ 1 + preDASSTOT + Gender + Stresscoms + postDASSTOT + Ndcp + postDASSTOT*Ndcp, data = port)
fit_port3.summ <- summary(fit_port3)
summ(fit_port3, digits = 5)

##Probe Significant interactions
#Australia - Ndcp
interact_plot(fit_aust3, pred = postDASSTOT, modx = Ndcp)

negDC_A <- mean(aust$Ndcp, na.rm = TRUE) + (sd(aust$Ndcp, na.rm = TRUE))
negDC_m <- mean(aust$Ndcp, na.rm = TRUE)
negDC_B <- mean(aust$Ndcp, na.rm = TRUE) - (sd(aust$Ndcp, na.rm = TRUE))

list1 <- list(Ndcp=c(negDC_A, negDC_m, negDC_B))

simple_slopes_aust <- emtrends(fit_aust3, ~Ndcp, var="postDASSTOT", at=list1)
summary(simple_slopes_aust)

#Romania - posDC
interact_plot(fit_rom2, pred = postDASSTOT, modx = posDC)

posDC_A <- mean(rom$posDC, na.rm = TRUE) + (sd(rom$posDC, na.rm = TRUE))
posDC_m <- mean(rom$posDC, na.rm = TRUE)
posDC_B <- mean(rom$posDC, na.rm = TRUE) - (sd(rom$posDC, na.rm = TRUE))

list2 <- list(posDC=c(posDC_A, posDC_m, posDC_B))

simple_slopes_rom <- emtrends(fit_rom2, ~posDC, var="postDASSTOT", at=list2)
summary(simple_slopes_rom)
