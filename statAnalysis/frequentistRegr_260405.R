library(dplyr)
library(ggplot2)
library(patchwork)
# for methodological guidance please consider
# https://stats.oarc.ucla.edu/r/dae/logit-regression/
setwd("~/Dropbox/Coding/Code/recur/statAnalysis")

# -------- Load and standardize data ----------------------------------------
#mydata = tibble(read.csv('mydata.csv'))
source('dataGatherer.R')
# sd from 0
data_std = mydata |>  
  #filter(iROC_plan <=.6) |>
  mutate(across(c(CLro:ks, AGE, BW), ~c(scale(.))))

# -------. Risk of RNB depending on how the patient is-----------------------
formulaPat = paste0("isRNB ~ 1 + AGE + BW + RAC + ",
                    "CLro + Q2ro + V1ro + V2ro + V1su + ",
                    "V2su + CLsu + Q2su + ke0 + EC50 + Hill + ks")
modelPat = glm(formula=formulaPat,
               family = binomial,
               data = data_std)
print (exp(cbind(OR=coef(modelPat), confint(modelPat, level=.9 ))))

# ------- Risk of RNB depending of what we do to the patient------------------
modelWe = glm(formula="isRNB ~ 1 + catSUG + cumROC",
              family = binomial,
              data = data_std ) 
print (exp(cbind(OR=coef(modelWe), confint(modelWe, level=.9 ))))

# ------- Diff ratio SUG/remROC among isRNB           ------------------
modelRatio = glm(formula="rSUG_remROC ~ 1 + isRNB",
              family = gaussian,
              data = data_std)
print (exp(cbind(OR=coef(modelRatio), confint(modelRatio, level=.9 ))))

# ------- Diff amount roc just before SUG, among isRNB------------------
modelA1 = glm(formula="A1ro ~ 1 + isRNB",
                 family = gaussian,
                 data = data_std)
print (exp(cbind(OR=coef(modelA1), confint(modelA1, level=.9 ))))

# ------- Diff amount roc just before SUG, among isRNB------------------
modelA2 = glm(formula="A2ro ~ 1 + isRNB",
              family = gaussian,
              data = data_std)
print (exp(cbind(OR=coef(modelA2), confint(modelA2, level=.9 ))))

# A try on plot logistic regression for continuous variables------------
# Remember, cumROC in mikroM / Kg
ndata = tibble(cumROC = seq(0.6,4.0,
                            length.out=1000),
               catSUG = rep(c('fixed','adjusted'),500))

pred = predict.glm(modelWe, newdata=ndata,type='link',se.fit=TRUE)
plot(ndata$cumROC,exp(pred$fit),type='l', lwd=2,
     main='OR per unit change.OBS! SomethingWRONG here', 
     ylab='OR', 
     xlab='cumROC, mg/Kg')
lines(ndata$cumROC,exp(pred$fit + 1.96*pred$se.fit), type='l')
lines(ndata$cumROC,exp(pred$fit - 1.96*pred$se.fit), type='l')
# ------- prediction RESULT
pred = predict.glm(modelWe, newdata=ndata,type='response',se.fit=TRUE)
plot(ndata$cumROC,pred$fit,type='l', lwd=2,
     main='Actual probabilities', ylab='P(RNB)', xlab='cumROC, mg/Kg')
lines(ndata$cumROC,pred$fit + 1.96*pred$se.fit, type='l', lwd=.5)
lines(ndata$cumROC,pred$fit - 1.96*pred$se.fit, type='l', lwd=.5)

# as stats.oarc.ucla.edu------------------------------------------------- 
ndataPlus = cbind(ndata,predict(modelWe, newdata=ndata, type='link',se=TRUE))
ndataPlus = within(ndataPlus, {PredictedProb <- plogis(fit)
PredictedProb <- plogis(fit)
LL <- plogis(fit - 1.96*se.fit)
UL <- plogis(fit + 1.96*se.fit) })
print(head(ndataPlus))

gP_cumROC = ggplot(ndataPlus, aes(x = cumROC, y = PredictedProb)) + 
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = catSUG), alpha = 0.2) + 
  geom_line(aes(colour = catSUG),size = 1) 

plot(gP_cumROC)













