# harmonisation prep / data prep for harmonisied test data
# Last change: 27 Jun 2025
# Max Korbmacher, max.korbmacher@gmail.com
#
# prep ----
# wash your hands before eating
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# define the savepath
savepath = "/Users/max/Documents/Local/MS/NormativeModels/results/"
#
pacman::p_load(tidyverse,psych,effsize,ggseg,patchwork,rstatix,ggpubr,
               caret,lme4,lmerTest,haven,reshape2,stats,entropy,
               ggseg3d, longCombat, neuroCombat)
long = read.csv("/Users/max/Documents/Local/MS/NormativeModels/code/long_output.csv")
cross = read.csv("/Users/max/Documents/Local/MS/NormativeModels/code/test.csv")
#
# harmonise longitudinal data ----
LC = function(dat){
  features = c(dat %>% select(starts_with("lh") & ends_with("volume")) %>% names(),
               dat %>% select(starts_with("rh") & ends_with("volume")) %>% names(),
               names(dat)[grepl(c("halamus|allidum|mygdala|
                                campus|utamen|audate|
                                CC|Cerebellum.Cortex"),names(dat))],
               dat %>% select(TotalGrayVol, EstimatedTotalIntraCranialVol,
                              ends_with("Hippocampus")) %>% names())
  dat = dat %>% select(all_of(features),eid,data,sex,session,scanner, age, edss)%>%na.omit
  covars = dat %>% dplyr::select(sex,edss,data,age)
  dat = longCombat(idvar = "eid", timevar = "session", batchvar = "scanner", 
                   features = features,
                   formula = "age + sex", ranef = "(1|eid)", data = dat)
  dat = dat$data_combat
  dat = cbind(covars,dat)
  colnames(dat) = gsub('.combat','',colnames(dat))
  return(dat)
}
long = LC(long)
# harmonise cross-sectional data ----
covars = cross %>% dplyr::select(eid,sex,scanner,age,data)
covars$sex = ifelse(covars$sex == "F" | covars$sex == "Female",0,1)
datasets = covars$data
covars$data = as.numeric(factor(cross$data))
cross = neuroCombat(t(cross%>%dplyr::select(EstimatedTotalIntraCranialVol,starts_with("Left"),starts_with("Right"), starts_with("lh"),starts_with("rh"))),batch=as.numeric(factor(cross$scanner)),mod=model.matrix(~covars$age+covars$sex), mean.only = T)
cross = data.frame(t(cross$dat.combat))
cross = cbind(covars,cross)
# save ----
write.csv(file = paste(savepath,"cross_raw_harmonised.csv",sep=""),cross)
write.csv(file = paste(savepath,"long_raw_harmonised.csv",sep=""),long)
