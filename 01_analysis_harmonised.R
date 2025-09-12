# Application of Brain Reference to cross sectional and longitudinal MS data
# Max Korbmacher, April 2025
#
# ----------------------------------- #
# --------------Structure------------ #
# ----------------------------------- #
# 0. Data wrangling------------------ #
# 1. Case-control checks------------- #
# 1.1 number of deviations----------- #
# 1.2 Z-score comparison------------- #
# 1.3 Diagnostic prediction---------- #
# 2. Longitudinal assessment--------- #
# 2.1 Development of deviations------ #
# 3. Individual-level assessments---- #
# 3.1 At baseline-------------------- #
# 3.2 Longitudinally----------------- #
# ----------------------------------- #
# ----------------------------------- #
#
# wash your hands before eating
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# define the savepath
savepath = "/Users/max/Documents/Local/MS/NormativeModels/Regular/results/"
#
# 0. Data wrangling----------------
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,psych,effsize,ggseg,patchwork,rstatix,ggpubr,
               caret,lme4,lmerTest,haven,reshape2,stats,entropy,
               ggseg3d, ggsegExtra, mice, pwr, longCombat)
# power
pwr.f2.test(u = 3, v = 43, f2 = NULL, sig.level = 0.05, power = .80)

# load data
cross = read.csv("/Users/max/Documents/Local/MS/NormativeModels/Regular/code/Zscores_harmonised.csv")
long = read.csv("/Users/max/Documents/Local/MS/NormativeModels/Regular/code/Zscores_long_harmonised.csv")
demo10 = read.csv('/Users/max/Documents/Local/MS/data/OFAMS88_OFAMS10_lifestylepaper_ updated_beskyttet.csv',sep = ";") # 10 years follow up demo & clin scores
PASAT_OSL = read.csv("/Users/max/Documents/Local/Data/Oslo/Database_cognitive_phenotypes_MS_Oslo.csv")# zscored PASAT
fati = read_sas("/Users/max/Documents/Local/MS/demographics/Statistikk-filer/fatigue.sas7bdat") # fatigue scores
fati_OSL = read.csv("/Users/max/Documents/Local/Data/Oslo/fatigue_Oslo.csv")

# missingness
summary(is.na(long %>% select(c(eid, age, sex, data, edss, session), ends_with("z_score"))))
paste("Total number of available MRI sessions:",long %>% select(c(eid, age, sex, data, session), ends_with("z_score")) %>% na.omit %>% nrow)
paste("EDSS is missing for this proportion of sessions:",
      (1-long %>% select(c(eid, age, sex, data, edss, session), ends_with("z_score")) %>% na.omit %>% nrow /
         long %>% select(c(eid, age, sex, data, session), ends_with("z_score")) %>% na.omit %>% nrow)*100,
      "%, corresponding to this number of sessions:", long %>% select(c(eid, age, sex, data, session), ends_with("z_score")) %>% na.omit %>% nrow-
        long %>% select(c(eid, age, sex, data, edss, session), ends_with("z_score")) %>% na.omit %>% nrow)
# imputation
long_no_na = long %>% select(c(eid, age, sex, data, session), ends_with("z_score")) %>% na.omit
long00 = long[long$eid %in% long_no_na$eid,] %>% select(eid,age,sex,edss,session)
EID = long_no_na$eid
# Convert eid to integer as required by mice for clustering
long00$eid <- as.integer(as.factor(long00$eid))  # ensures grouping is preserved
long00$sex <- as.factor(long00$sex)              # recommended for categorical

# Setup mice
imp0 <- mice(long00, maxit = 0)
pred <- imp0$predictorMatrix
meth <- imp0$method

# Define method and predictor structure
meth["edss"] <- "2l.norm"
pred["edss", ] <- 0
pred["edss", "eid"] <- -2      # cluster ID (must be integer!)
pred["edss", "age"] <- 1       # time
pred["edss", "sex"] <- 1       # covariate

# Perform imputation
imp <- mice(long00, method = meth, predictorMatrix = pred, m = 5, seed = 123)

# Extract all completed datasets as a list
completed_list <- lapply(1:imp$m, function(i) complete(imp, i))

# Bind all imputations together with an imputation index
all_imputations <- bind_rows(
  lapply(seq_along(completed_list), function(i) {
    completed_list[[i]] %>% mutate(.imp = i)
  })
)

# Identify missing rows in original data
missing_idx <- which(is.na(long00$edss))

# Extract identifying keys for those missing rows from original data
keys_missing <- long00[missing_idx, c("eid", "session", "age")]

# Calculate average imputed edss values for missing rows, including keys
average_imputed_edss <- all_imputations %>%
  filter(row_number() %in% missing_idx) %>%    # Filter only missing rows (across all imputations)
  group_by(eid, session, age) %>%               # Group by unique row keys
  summarise(mean_edss = mean(edss), .groups = "drop")

# make sure that the first session is labelled correctly
long <- long %>%
  group_by(eid) %>%
  mutate(session = session - min(session) + 1) %>%
  ungroup()

long_surrogate = long[long$eid %in% long_no_na$eid,]
long_surrogate$eid <- as.integer(as.factor(long_surrogate$eid))  # ensures grouping is preserved

long <- long_surrogate %>%
  left_join(average_imputed_edss, by = c("eid", "session", "age")) %>%
  mutate(edss = if_else(is.na(edss), mean_edss, edss)) %>%
  select(-mean_edss)
long = na.omit(long)
long$eid = EID
cross$diagnosis= c(replicate(nrow(cross)/2,"MS"),replicate(nrow(cross)/2,"HC"))

# descriptives
## rows cross
cross %>% filter(diagnosis == "MS") %>%nrow # matched sample
length(unique(long %>% filter(data == "MS" & session ==1)  %>% na.omit %>% pull(eid))) # Oslo
length(unique(long %>% select(c(eid, age, sex, data, edss, session), ends_with("z_score")) %>% filter(data == "OFAMS" & session ==1)  %>% na.omit %>% pull(eid))) # Bergen

## rows long
length(unique(long  %>% na.omit %>% pull(eid)))
length(long  %>% na.omit %>% pull(eid))

length(unique(long %>% filter(data == "MS")  %>% na.omit %>% pull(eid)))
long %>% filter(data == "MS") %>% na.omit %>% nrow

length(unique(long %>% filter(data == "OFAMS")%>% na.omit %>% pull(eid)))
long %>% filter(data == "OFAMS") %>% na.omit %>% nrow

## age
long %>% filter(session == 1) %>% select(age) %>% na.omit %>% summarize(M = round(mean(age),1), SD = round(sd(age),1))
long %>% filter(data == "OFAMS" & session == 1) %>% select(age) %>% na.omit %>% summarize(M = round(mean(age),1), SD = round(sd(age),1))
long %>% filter(data == "MS" & session == 1) %>% select(age) %>% na.omit %>% summarize(M = round(mean(age),1), SD = round(sd(age),1))

## sex
long %>% filter(session == 1) %>% group_by(sex) %>% summarise(n = n()) %>%mutate(freq = n / sum(n))
long %>% filter(data == "OFAMS" & session == 1) %>% group_by(sex) %>% summarise(n = n()) %>%mutate(freq = n / sum(n))
long %>% filter(data == "MS" & session == 1) %>% group_by(sex) %>% summarise(n = n()) %>%mutate(freq = n / sum(n))

## edss
long %>% filter(session == 1) %>% select(edss) %>% na.omit %>% summarize(M = round(mean(edss),1), SD = round(sd(edss),1))
long %>% filter(data == "OFAMS"& session == 1) %>% select(edss) %>% na.omit %>% summarize(M = round(mean(edss),1), SD = round(sd(edss),1))
long %>% filter(data == "MS"& session == 1) %>% select(edss) %>% na.omit %>% summarize(M = round(mean(edss),1), SD = round(sd(edss),1))

# 1. Case-control checks-----------
# 1.1 number of deviations---------
cross$nb_deviations = cross %>%
  select(ends_with("z_score")) %>%
  mutate_all(~ ifelse((.) <= -1.96, 1, 0)) %>% # mutate_all(~ ifelse(abs(.) >= 1.96, 1, 0)) %>%
  transmute(z_score_sum = rowSums(.)) %>%
  pull(z_score_sum)
cross %>% group_by(diagnosis) %>% summarize(M = mean(nb_deviations), SD = sd(nb_deviations))
psych::cohen.d(cross$nb_deviations,factor(cross$diagnosis))[1]
rstatix::t_test(cross, nb_deviations~diagnosis)

# 1.2 Z-score comparison-----------
# Select z_score columns
z_scores = cross %>% select(ends_with("z_score"))

# Add diagnosis vector as a column if not already part of the data frame
z_scores = z_scores %>% mutate(diagnosis = cross$diagnosis)

# Calculate the group-wise mean for each column
group_means = z_scores %>%
  group_by(diagnosis) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(-diagnosis, names_to = "variable", values_to = "mean_value") %>%
  pivot_wider(names_from = diagnosis, values_from = mean_value)

# Calculate difference (MS - HC or vice versa)
diffs = group_means %>%
  mutate(difference = MS - HC)

# Calculate mean and SD of absolute values
diffs %>% summarise(mean_abs = mean(abs(difference)),sd_abs = sd(abs(difference)))
diffs %>% summarise(mean_abs = mean((difference)),sd_abs = sd((difference)))

# express this as an effect size (difference of means) and add a p-value
effsize::cohen.d(diffs%>% gather(variable, value) %>%filter(variable!="difference")%>%pull(value),diffs%>%gather(variable, value) %>%filter(variable!="difference")%>%pull(variable))
t_test(diffs%>%gather(variable, value) %>%filter(variable!="difference"),value~variable)

diffs[order(diffs$MS),]

# plot group-wise Z-vals and their differences
z_long <- diffs %>%
  rename(region = variable) %>%
  pivot_longer(cols = c(HC, MS, difference), names_to = "group", values_to = "value") %>%
  mutate(
    group = factor(group, levels = c("HC", "MS", "difference")),
    region = str_remove(region, "_z_score"), # Remove suffix
  )

# Split into cortical vs subcortical based on naming patterns
z_cortical <- z_long %>% filter(str_detect(region, "^lh_|^rh_")) 
#%>%mutate(region = str_remove(region, "^lh_|^rh_"))
z_subcortical <- z_long %>% filter(!str_detect(region, "^lh_|^rh_"))



# order the data frame
z_cortical = z_cortical[order(z_cortical$region),]
z_subcortical = z_subcortical[order(z_subcortical$region),] 

# label data correctly
z_cortical = (z_cortical %>% group_by(group)) %>% mutate(label = c(replicate(2,(brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",z_cortical$region)][1:34]))))
z_subcortical = (z_subcortical %>% group_by(group)) %>% mutate(label = (c(brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))],"Left-Cerebellum-Cortex","Right-Cerebellum-Cortex"))[order(c(brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))],"Left-Cerebellum-Cortex","Right-Cerebellum-Cortex"))])
# add hemi
z_cortical$hemi = ifelse(grepl("lh_",z_cortical$region)==T,"left","right")
z_subcortical$hemi = ifelse(grepl("Left.",z_subcortical$region)==T,"left","right")

names(z_cortical) = c("label","group","Z","region","hemi")
names(z_subcortical) = c("region","group","Z","label","hemi")
#
#
#Cortical plots
C1 = ggplot(dk %>% as_tibble() %>% left_join(z_cortical %>% select(region,hemi,Z) %>% filter(group=="HC") %>% as_tibble())) + geom_brain(atlas = dk,aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) +
  theme_void() + theme(legend.position="none")
C2 = ggplot(dk %>% as_tibble() %>% left_join(z_cortical %>% select(region,hemi,Z) %>% filter(group=="MS") %>% as_tibble())) + geom_brain(atlas = dk,aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) +
  theme_void() + theme(legend.position="none")
C3 = ggplot(dk %>% as_tibble() %>% left_join(z_cortical %>% select(region,hemi,Z) %>% filter(group=="difference") %>% as_tibble())) + geom_brain(atlas = dk,aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) +
  theme_void() + theme(legend.position="none")

# Subcortical plots
## This plots exclusively the coronal stats (all others are troublesome)
coronal_brain_aseg = as_tibble(aseg) %>%
  filter(side == "coronal", !grepl("\\d", label))
z_subcortical = merge(z_subcortical,coronal_brain_aseg,by="label")
S1 = ggplot(z_subcortical%>%filter(group=="HC")) + geom_brain(atlas = aseg, side = "coronal",aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) +
  #labs(title="Regional volume loss") + 
  theme_void()
S2 = ggplot(z_subcortical%>%filter(group=="MS")) + geom_brain(atlas = aseg, side = "coronal",aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) +
  #labs(title="Regional volume loss") + 
  theme_void()
S3 = ggplot(z_subcortical%>%filter(group=="difference")) + geom_brain(atlas = aseg, side = "coronal",aes(fill = Z),color="black")+
  scale_fill_gradient2(low = "blue",mid = "white",high="red",limits = c(-1.25,1.75)) + 
  #labs(title="Regional volume loss") + 
  theme_void()

# merge the plots
p1 = ggarrange(C1,S1,nrow=1,widths=c(2,.5))
p2 = ggarrange(C2,S2,nrow=1,widths=c(2,.5))
p3 = ggarrange(C3,S3,nrow=1,widths=c(2,.5))
Zplot = ggarrange(p1,p2,p3,ncol=1,labels = c("HC","MS","Diff",
                                             hjust = c(0,0,0),common.legend = T, legend = "right"))
# save a single figure
ggsave(paste(savepath,"Zplot_harmonised.pdf",sep=""),plot=Zplot, width = 10, height = 6)

# 1.3 Diagnostic prediction--------

# Here, we compare whether the Z scores or regular volumetrics perform better.

# Code to train the SVM 
set.seed(1234) 
# set the 3 fold crossvalidation with AU  
# to pick for us what we call the best model 
control <- trainControl(method="cv",number=3, classProbs = TRUE) 
metric <- "Accuracy"
# model <- train(diagnosis ~., data = cross%>%select(ends_with("z_score"),diagnosis), 
#                method = "svmRadial", 
#                tuneLength = 50,preProc = c("center","scale"),  
#                metric=metric, trControl=control)
# model
# plot(model)
# predict <- predict(model, newdata = cross) 
# confusionMatrix(predict, factor(cross$diagnosis))
# #
# #
# model1 <- train(diagnosis ~., data = cross%>%select(ends_with("volume"),diagnosis), 
#                method = "svmRadial", 
#                tuneLength = 50,preProc = c("center","scale"),  
#                metric=metric, trControl=control)
# model1
# plot(model1)
# predict <- predict(model1, newdata = cross) 
# confusionMatrix(predict, factor(cross$diagnosis))
# #
# #
# # Try again without CV to check whether Z scores are still performing better
# model3 <- train(diagnosis ~., data = cross%>%select(ends_with("z_score"),diagnosis), 
#                method = "svmRadial", 
#                tuneLength = 50,preProc = c("center","scale"),  
#                metric=metric)
# model3
# plot(model3)
# predict <- predict(model3, newdata = cross) 
# confusionMatrix(predict, factor(cross$diagnosis))
# 
# 
# model4 <- train(diagnosis ~., data = cross%>%select(ends_with("volume"),diagnosis), 
#                 method = "svmRadial", 
#                 tuneLength = 50,preProc = c("center","scale"),  
#                 metric=metric)
# model4
# plot(model4)
# predict <- predict(model4, newdata = cross) 
# confusionMatrix(predict, factor(cross$diagnosis))
#
# The comparison is somewhat inconclusive.
# Z-values seem to be better, but since the sample is so small and the group 
# differences so big, it is difficult to say which set of variable is really better.
#
# Hence, we comment this out for now.
#
#
#
#
# 2. Longitudinal assessment------- 
# 2.1 Development of deviations----
# 2.1.1 number of deviations and edss----
long$nb_deviations = long %>%
  select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(. <= -1.96, 1, 0)) %>%
  transmute(z_score_sum = rowSums(.)) %>%
  pull(z_score_sum)
m = lmer(edss~nb_deviations+age+sex+(1|eid),long)
summary(m)
m = lm(edss~nb_deviations+age+sex,long%>%filter(session == 1))
summary(m)
effectsize::standardize_parameters(m)
print("The number of deviations, measured by EDSS, does not tell us something about the disability development but baseline state.")

# 2.1.2 number of deviations and PASAT----
#
# OSL
PASAT_OSL$eid = gsub("MS","MS_",PASAT_OSL$subject_id)
PASAT_OSL$session = PASAT_OSL$tpoint
PASAT_OSL$PASAT = PASAT_OSL$MACFIMS_PASAT3_zscore
PASAT_OSL = PASAT_OSL %>% select(eid,session,PASAT)
PASAT_OSL %>% select(PASAT)%>% na.omit %>% nrow
length(unique((PASAT_OSL %>% na.omit)$eid))

# OFAMS
pasat1 = demo10%>%dplyr::select(Patnr,BL_PASATcorrect,PASAT_24M, PASAT_OFAMS10)
pasat1 = melt(pasat1, id.vars = c("Patnr"))
names(pasat1) = c("eid","session","PASAT")
pasat1$session = ifelse(pasat1$session == "BL_PASATcorrect",1,0)+ifelse(pasat1$session == "PASAT_24M",24,0)+ifelse(pasat1$session == "PASAT_OFAMS10",145,0)
long1 = merge(rbind(PASAT_OSL,pasat1),long,by=c("eid","session"))
m = lmer(PASAT~nb_deviations+age+sex+TotalGrayVol+(1|eid),long1)
summary(m)
m = lm(PASAT~nb_deviations+age+sex+TotalGrayVol,long1%>%filter(session == 1))
summary(m)
effectsize::standardize_parameters(m)
#
# 2.1.3 number of deviations and fatigue----
#
# prep eid and session
fati_OSL$session = substr(fati_OSL$eid,9,11)
fati_OSL$eid = substr(fati_OSL$eid,1,7)
fati_OSL = fati_OSL %>% select(eid,session,fatigue)
fati_OSL$session = as.numeric(fati_OSL$session)
#fati.cop = fati_OSL

# fix session factor levels
fati$session = factor(fati$VISIT)
fati$session = ifelse(fati$session == "Baseline",1,fati$session)
fati$session = ifelse(fati$session == 4,7,fati$session)
fati$session = ifelse(fati$session == 2,13,fati$session)
fati$session = ifelse(fati$session == 3,25,fati$session)
fati$eid = as.numeric(fati$patno)
fati$fatigue = fati %>% select(A,     B,     C,     D,     E,     F,     G,     H,     I) %>% rowMeans()
fati = fati %>% select(eid,session,fatigue)
fati = rbind(fati,fati_OSL)
long2 = merge(fati,long,by=c("eid","session"))
m = lmer(fatigue~nb_deviations+age+sex+(1|eid),long2)
summary(m)
effectsize::standardize_parameters(m)
m = lm(fatigue~nb_deviations+age+sex,long2%>%filter(session == 1))

# 2.1.4 number of deviations and age(ing)----
m = lm(nb_deviations~age+sex+TotalGrayVol,long%>%filter(session == 1))
m = lmer(nb_deviations~age+sex+TotalGrayVol+(1|eid),long)
summary(m)
effectsize::standardize_parameters(m)
#
#
# 2.2 Regional associations ----
# ---- Function to extract standardized coefficients manually
extract_std_coeffs <- function(data, predictor = "age") {
  regions <- data %>% select(ends_with("z_score")) %>% names()
  
  # Standardize predictor
  data <- data %>%
    mutate(std_predictor = scale(.data[[predictor]])[,1])
  
  std_coeffs <- sapply(regions, function(region) {
    # Standardize response
    data <- data %>%
      mutate(std_response = scale(.data[[region]])[,1])
    
    mod <- lmer(std_response ~ std_predictor + (1|eid), data = data)
    fixef(mod)["std_predictor"]
  })
  
  tibble(
    label = str_remove(regions, "_z_score"),
    std_coef = std_coeffs
  )
}

# ---- Function to plot
run_and_plot <- function(long, predictor) {
  coefs <- extract_std_coeffs(long, predictor)
  
  # Cortical
  test_cort <- z_cortical %>%
    filter(group == "MS") %>%
    left_join(coefs, by = "label")
  
  p1 <- ggplot(dk %>% as_tibble() %>%
                 left_join(test_cort %>% select(region, hemi, std_coef), by = c("region", "hemi"))) +
    geom_brain(atlas = dk, aes(fill = std_coef), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.35, 0.35)) +
    theme_void() +
    theme(legend.position = "none")
  
  # Subcortical
  test_sub <- z_subcortical %>%
    filter(group == "MS") %>%
    left_join(coefs, by = c("region.x" = "label"))
  
  p2 <- ggplot(test_sub) +
    geom_brain(atlas = aseg, side = "coronal", aes(fill = std_coef), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.35, 0.35)) +
    theme_void() +
    labs(fill = "Std.Coeff.")
  
  ggarrange(p1, p2,widths=c(2,.5))
}

# ---- Application of functions
age_plot     <- run_and_plot(long, "age")
edss_plot    <- run_and_plot(long, "edss")
pasat_plot   <- run_and_plot(long1, "PASAT")
fatigue_plot <- run_and_plot(long2, "fatigue")

large_plot = ggarrange(age_plot, edss_plot,
                       pasat_plot, fatigue_plot, 
                       ncol=1,labels=c("Age","EDSS","PASAT","Fatigue"),
                       hjust = c(0,0,0,0))
ggsave(paste(savepath,"Long_Effects_harmonised.pdf",sep=""),plot=large_plot, width = 10, height = 8)

# Function to also extract standardized coefficients for one predictor

# Helper: Standardize a variable
z <- function(x) {
  if (is.numeric(x)) {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  } else {
    return(x)
  }
}

# Function: Extract standardized coefficients by predictor
run_and_extract <- function(df, predictor) {
  regions <- df %>% select(ends_with("z_score")) %>% names()
  
  # Manually standardize all numeric variables
  df_std <- df %>%
    mutate(across(where(is.numeric), scale))  # z-score standardization
  
  extract_coef <- function(region) {
    tryCatch({
      formula_str <- paste0(region, " ~ ", predictor, " + (1|eid)")
      model <- lmer(as.formula(formula_str), data = df_std)
      coefs <- fixef(model)[predictor]
      tibble(region = gsub("_z_score", "", region),
             predictor = predictor,
             std_coef = as.numeric(coefs))
    }, error = function(e) {
      message("Skipping ", region, " due to error: ", e$message)
      return(NULL)
    })
  }
  
  results <- map_dfr(regions, extract_coef)
  return(results)
}
age_coefs     <- run_and_extract(long,  "age")
edss_coefs    <- run_and_extract(long,  "edss")
pasat_coefs   <- run_and_extract(long1, "PASAT")
fatigue_coefs <- run_and_extract(long2, "fatigue")

all_coefs <- bind_rows(age_coefs, edss_coefs, pasat_coefs, fatigue_coefs)

subcortical_wide <- all_coefs %>%
  pivot_wider(names_from = predictor, values_from = std_coef)
subcortical_wide[2:5] = round(subcortical_wide[2:5],2)
write.csv(x = subcortical_wide,paste(savepath,"long_associations_harmonised.csv",sep=""),row.names = FALSE)

run_and_extract_pvals <- function(df, predictor) {
  regions <- df %>% select(ends_with("z_score")) %>% names()
  
  # Manually standardize all numeric variables
  df_std <- df %>%
    mutate(across(where(is.numeric), scale))
  
  extract_pval <- function(region) {
    tryCatch({
      formula_str <- paste0(region, " ~ ", predictor, " + (1|eid)")
      model <- lmer(as.formula(formula_str), data = df_std)
      pval <- summary(model)$coefficients[predictor, "Pr(>|t|)"]
      tibble(region = gsub("_z_score", "", region),
             predictor = predictor,
             p_value = as.numeric(pval))
    }, error = function(e) {
      message("Skipping ", region, " due to error: ", e$message)
      return(NULL)
    })
  }
  
  results <- map_dfr(regions, extract_pval)
  return(results)
}
# apply functions and put it all together
age_pvals     <- run_and_extract_pvals(long,  "age")
edss_pvals    <- run_and_extract_pvals(long,  "edss")
pasat_pvals   <- run_and_extract_pvals(long1, "PASAT")
fatigue_pvals <- run_and_extract_pvals(long2, "fatigue")
all_pvals <- bind_rows(age_pvals, edss_pvals, pasat_pvals, fatigue_pvals)
#all_pvals$p_value = all_pvals$p_value * nrow(all_pvals) # can be used for Bongferroni correction

all_pvals_uncorrected = all_pvals
all_pvals$p_value = p.adjust(all_pvals$p_value,method = "fdr")
all_pvals_wide <- all_pvals %>%
  pivot_wider(names_from = predictor, values_from = p_value)
all_pvals_wide[2:5] = round(all_pvals_wide[2:5],3)
write.csv(x = all_pvals_wide,paste(savepath,"long_associations_pvalues_harmonised.csv",sep=""),row.names = FALSE)

all_pvals_wide1 <- all_pvals_uncorrected %>%
  pivot_wider(names_from = predictor, values_from = p_value)
all_pvals_wide1[2:5] = round(all_pvals_wide1[2:5],3)
write.csv(x = all_pvals_wide1,paste(savepath,"long_associations_pvalues_uncorrected_harmonised.csv",sep=""),row.names = FALSE)

# 1.4 Regional associations--------

# Helper: Standardize a variable
z <- function(x) {
  if (is.numeric(x)) {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  } else {
    return(x)
  }
}

# Extract standardized coefficients using lm()
extract_std_coeffs_cross <- function(data, predictor = "age") {
  regions <- data %>% select(ends_with("z_score")) %>% names()
  data <- data %>% mutate(std_predictor = scale(.data[[predictor]])[,1])
  
  std_coeffs <- sapply(regions, function(region) {
    data <- data %>% mutate(std_response = scale(.data[[region]])[,1])
    mod <- lm(std_response ~ std_predictor, data = data)
    coef(mod)["std_predictor"]
  })
  
  tibble(
    label = str_remove(regions, "_z_score"),
    std_coef = std_coeffs
  )
}

# Run and extract standardized coefficients (long format)
run_and_extract_cross <- function(df, predictor) {
  regions <- df %>% select(ends_with("z_score")) %>% names()
  df_std <- df %>% mutate(across(where(is.numeric), scale))
  
  map_dfr(regions, function(region) {
    tryCatch({
      model <- lm(as.formula(paste0(region, " ~ ", predictor)), data = df_std)
      tibble(region = gsub("_z_score", "", region),
             predictor = predictor,
             std_coef = coef(model)[predictor])
    }, error = function(e) {
      message("Skipping ", region, ": ", e$message)
      NULL
    })
  })
}

# Run and extract p-values (long format)
run_and_extract_pvals_cross <- function(df, predictor) {
  regions <- df %>% select(ends_with("z_score")) %>% names()
  df_std <- df %>% mutate(across(where(is.numeric), scale))
  
  map_dfr(regions, function(region) {
    tryCatch({
      model <- lm(as.formula(paste0(region, " ~ ", predictor)), data = df_std)
      pval <- summary(model)$coefficients[predictor, "Pr(>|t|)"]
      tibble(region = gsub("_z_score", "", region),
             predictor = predictor,
             p_value = pval)
    }, error = function(e) {
      message("Skipping ", region, ": ", e$message)
      NULL
    })
  })
}

# ---- Plotting Function
run_and_plot_cross <- function(df, predictor) {
  coefs <- extract_std_coeffs_cross(df, predictor)
  
  # Cortical
  test_cort <- z_cortical %>%
    filter(group == "MS") %>%
    left_join(coefs, by = "label")
  
  p1 <- ggplot(dk %>% as_tibble() %>%
                 left_join(test_cort %>% select(region, hemi, std_coef), by = c("region", "hemi"))) +
    geom_brain(atlas = dk, aes(fill = std_coef), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.35, 0.35)) +
    theme_void() +
    theme(legend.position = "none")
  
  # Subcortical
  test_sub <- z_subcortical %>%
    filter(group == "MS") %>%
    left_join(coefs, by = c("region.x" = "label"))
  
  p2 <- ggplot(test_sub) +
    geom_brain(atlas = aseg, side = "coronal", aes(fill = std_coef), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.35, 0.35)) +
    theme_void() +
    labs(fill = "Std.Coeff.")
  
  ggarrange(p1, p2, widths = c(2, 0.5))
}

# ---- Apply All Cross-sectional Analyses
age_plot     <- run_and_plot_cross(long %>% filter(session == 1), "age")
edss_plot    <- run_and_plot_cross(long %>% filter(session == 1), "edss")
pasat_plot   <- run_and_plot_cross(long1 %>% filter(session == 1), "PASAT")
fatigue_plot <- run_and_plot_cross(long2 %>% filter(session == 1), "fatigue")

# Combine plots
large_plot <- ggarrange(age_plot, edss_plot,pasat_plot, fatigue_plot, # did not compute: , pasat_plot, fatigue_plot
                        ncol = 1, labels = c("Age", "EDSS", "PASAT", "Fatigue"),
                        hjust = c(0, 0, 0, 0))
ggsave(paste0(savepath, "Cross_Effects_harmonised_harmonised.pdf"), plot = large_plot, width = 10, height = 8)

# ---- Extract Coefficients and P-values for Tables
age_coefs     <- run_and_extract_cross(long %>% filter(session == 1), "age")
edss_coefs    <- run_and_extract_cross(long %>% filter(session == 1), "edss")
pasat_coefs   <- run_and_extract_cross(long1 %>% filter(session == 1), "PASAT")
fatigue_coefs <- run_and_extract_cross(long2 %>% filter(session == 1), "fatigue")

all_coefs <- bind_rows(age_coefs, edss_coefs, pasat_coefs, fatigue_coefs) %>%
  pivot_wider(names_from = predictor, values_from = std_coef)
all_coefs[2:5] <- round(all_coefs[2:5], 2)
write.csv(all_coefs, paste0(savepath, "cross_associations_harmonised.csv"), row.names = FALSE)

age_pvals     <- run_and_extract_pvals_cross(long %>% filter(session == 1), "age")
edss_pvals    <- run_and_extract_pvals_cross(long %>% filter(session == 1), "edss")
pasat_pvals   <- run_and_extract_pvals_cross(long1 %>% filter(session == 1), "PASAT")
fatigue_pvals <- run_and_extract_pvals_cross(long2 %>% filter(session == 1), "fatigue")

all_pvals <- bind_rows(age_pvals, edss_pvals, pasat_pvals, fatigue_pvals)
all_pvals_uncorrected <- all_pvals

# FDR correction
all_pvals$p_value <- p.adjust(all_pvals$p_value, method = "fdr")

# Wide format tables
all_pvals_wide <- all_pvals %>%
  pivot_wider(names_from = predictor, values_from = p_value)
all_pvals_wide[2:5] <- round(all_pvals_wide[2:5], 3)
write.csv(all_pvals_wide, paste0(savepath, "cross_associations_pvalues_harmonised.csv"), row.names = FALSE)

all_pvals_uncorrected_wide <- all_pvals_uncorrected %>%
  pivot_wider(names_from = predictor, values_from = p_value)
all_pvals_uncorrected_wide[2:5] <- round(all_pvals_uncorrected_wide[2:5], 3)
write.csv(all_pvals_uncorrected_wide, paste0(savepath, "cross_associations_pvalues_uncorrected_harmonised.csv"), row.names = FALSE)

# ######### follow-up for CIs: Cross-sectional data
# ## Age associations
# x = (lm(scale(rh_caudalmiddlefrontal_volume)~scale(age),data = long %>% filter(session == 1)))
# summary(x)
# confint(x)
# ## EDSS assoctiations
# x = (lm(scale(Left.Thalamus)~scale(edss),data = long %>% filter(session == 1)))
# summary(x)
# confint(x)
# x = (lm(scale(Right.Thalamus)~scale(edss),data = long %>% filter(session == 1)))
# summary(x)
# confint(x)
# 
# x = (lm(scale(Right.Putamen)~scale(edss),data = long %>% filter(session == 1)))
# summary(x)
# confint(x)
# 
# # FSS
# x = (lm(scale(rh_superiorfrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
# summary(x)
# confint(x)
# x = (lm(scale(lh_superiorfrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
# summary(x)
# confint(x)
# x = (lm(scale(rh_caudalmiddlefrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
# summary(x)
# confint(x)
# x = (lm(scale(rh_paracentral_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
# summary(x)
# confint(x)
# x = (lm(scale(lh_supramarginal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
# summary(x)
# confint(x)

######### follow-up for CIs: Longitudinal data
## Age associations
#
#

# 3. Individual-level assessments----
#
# 3.1 Baseline ----
#
test = cross %>% filter(diagnosis == "MS") %>% select(ends_with("z_score")) %>%
  mutate_all(~ ifelse((.) <= -1.96, 1, 0)) %>% 
  colSums() / nrow(cross %>% filter(diagnosis == "MS"))
test[order(data.frame(test)$test)]

max(cross %>% filter(diagnosis == "HC") %>% select(ends_with("z_score")) %>%
      mutate_all(~ ifelse((.) <= -1.96, 1, 0)) %>% 
      colSums() / nrow(cross %>% filter(diagnosis == "MS")))


plot_indiv_c = merge(z_cortical %>% filter(group=="HC"),data.frame(perc = test, label = gsub("_z_score","",names(test)), by = "label"))
plot_indiv_c$perc=plot_indiv_c$perc*100
indi0 = ggplot(dk %>% as_tibble() %>% left_join(plot_indiv_c %>% select(region,hemi,perc) %>% as_tibble())) + 
  geom_brain(atlas = dk,aes(fill = perc),color="black")+
  scale_fill_gradient2(low = "white",mid = "green",high="black",limits = c(0,40)) +
  theme_void() + theme(legend.position="none")
plot_indiv = merge(z_subcortical,data.frame(perc = test, region.x = gsub("_z_score","",names(test)), by = "region.x"))
plot_indiv$perc=plot_indiv$perc*100
indi1 = ggplot(plot_indiv) + geom_brain(atlas = aseg, side = "coronal",aes(fill = perc),color="black")+
  scale_fill_gradient2(low = "white",mid = "green",high="black",limits = c(0,40)) + 
  #labs(title="Regional volume loss") + 
  theme_void() + labs(fill = "Deviation %")

baseline = ggarrange(indi0, indi1,widths=c(2,.5))
ggsave(paste(savepath,"regional_heterogeniety_harmonised_harmonised.pdf",sep=""),plot=baseline, width = 10, height = 2)
relative_numbers = data.frame(perc = test, label = gsub("_z_score","",names(test)))


# 3.2 Longitudinally ----
#
#
# 1. Identify z_score columns
z_cols <- long %>% select(ends_with("z_score")) %>% names()

# 2. Classify deviations: 1 if <= -1.96 else 0
long_dev <- long %>%
  select(eid, age, all_of(z_cols)) %>%
  mutate(across(all_of(z_cols), ~ ifelse(. <= -1.96, 1, 0)))

# 3. Pivot longer for region-wise analysis
long_dev_long <- long_dev %>%
  pivot_longer(cols = all_of(z_cols), names_to = "region", values_to = "deviation") %>%
  arrange(eid, region, age)  # important for transitions

# 4. Within-subject deviation proportion
within_subject <- long_dev_long %>%
  group_by(eid, region) %>%
  summarise(
    prop_deviation = mean(deviation, na.rm = TRUE),
    n_visits = n(),
    n_deviating = sum(deviation, na.rm = TRUE),
    .groups = "drop"
  )

# 5. Entropy of deviation state per region × subject
entropy_table <- long_dev_long %>%
  group_by(eid, region) %>%
  summarise(
    entropy = entropy::entropy(table(deviation), unit = "log2"),
    .groups = "drop"
  )

# 6. Transition counts and consistency (e.g., 0→1, 1→0)
transitions <- long_dev_long %>%
  group_by(eid, region) %>%
  summarise(
    transitions = sum(abs(diff(deviation))),  # number of state changes
    first = first(deviation),
    last = last(deviation),
    .groups = "drop"
  )

# 7. Merge all summaries into one table
summary_stats <- within_subject %>%
  left_join(entropy_table, by = c("eid", "region")) %>%
  left_join(transitions, by = c("eid", "region"))

# ----------- VISUALIZATIONS -----------

# A. Histogram of deviation consistency
ggplot(summary_stats, aes(x = prop_deviation)) +
  geom_histogram(binwidth = 0.1, fill = "steelblue") +
  facet_wrap(~region, scales = "free_y") +
  theme_minimal() +
  labs(title = "Proportion of visits with deviation", x = "Proportion", y = "Count")

# B. Entropy distribution
ggplot(summary_stats, aes(x = entropy)) +
  geom_histogram(binwidth = 0.2, fill = "darkorange") +
  facet_wrap(~region, scales = "free_y") +
  theme_minimal() +
  labs(title = "Entropy of deviation pattern", x = "Entropy (bits)", y = "Count")

# C. Transitions heatmap (optional: top N variable regions)
top_regions <- summary_stats %>%
  group_by(region) %>%
  summarise(avg_transitions = mean(transitions, na.rm = TRUE)) %>%
  top_n(9, avg_transitions) %>%
  pull(region)

ggplot(summary_stats %>% filter(region %in% top_regions),
       aes(x = reorder(eid, transitions), y = region, fill = transitions)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  labs(title = "Deviation transitions per subject", x = "Subject", y = "Region") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())







# --- Transition classification per subject × region ---
transition_breakdown <- long_dev_long %>%
  group_by(eid, region) %>%
  mutate(prev = lag(deviation)) %>%
  filter(!is.na(prev)) %>%
  summarise(
    worsening = sum(prev == 0 & deviation == 1, na.rm = TRUE),
    recovery  = sum(prev == 1 & deviation == 0, na.rm = TRUE),
    stable    = sum(prev == deviation, na.rm = TRUE),
    n_transitions = n(),
    .groups = "drop"
  )

# --- Summary tables (from previous steps) ---
within_subject <- long_dev_long %>%
  group_by(eid, region) %>%
  summarise(
    prop_deviation = mean(deviation, na.rm = TRUE),
    n_visits = n(),
    n_deviating = sum(deviation, na.rm = TRUE),
    .groups = "drop"
  )

entropy_table <- long_dev_long %>%
  group_by(eid, region) %>%
  summarise(
    entropy = entropy::entropy(table(deviation), unit = "log2"),
    .groups = "drop"
  )

# Combine everything
summary_stats <- within_subject %>%
  left_join(entropy_table, by = c("eid", "region")) %>%
  left_join(transition_breakdown, by = c("eid", "region"))

# --- Optional: visualize proportions of worsening/recovery ---
summary_stats_long <- summary_stats %>%
  pivot_longer(cols = c(worsening, recovery, stable),
               names_to = "transition_type", values_to = "count") %>%
  mutate(prop = count / n_transitions)

ggplot(summary_stats_long, aes(x = transition_type, y = prop, fill = transition_type)) +
  geom_boxplot(alpha = 0.8) +
  facet_wrap(~region, scales = "free_y") +
  theme_minimal() +
  labs(title = "Relative frequency of transition types",
       y = "Proportion", x = "Transition Type")


#########














# Threshold-based classification
binary_df <- long %>%
  select(eid, age, session, ends_with("z_score")) %>%
  mutate(across(ends_with("z_score"), ~ ifelse(. <= -1.96, 1, 0)))

# Helper to compute transitions
compute_transitions <- function(data) {
  regions <- names(data)[grepl("z_score", names(data))]
  
  transition_list <- lapply(regions, function(region) {
    temp <- data %>%
      select(eid, session, !!sym(region)) %>%
      arrange(eid, session) %>%
      group_by(eid) %>%
      mutate(
        prev = lag(!!sym(region)),
        transition = case_when(
          is.na(prev) ~ NA_character_,
          prev == 0 & !!sym(region) == 1 ~ "worsen",
          prev == 1 & !!sym(region) == 0 ~ "improve",
          prev == !!sym(region) ~ "stable",
          TRUE ~ NA_character_
        )
      )
    
    temp %>%
      filter(!is.na(transition)) %>%
      count(transition) %>%
      mutate(region = gsub("_z_score", "", region))
  })
  
  bind_rows(transition_list) %>%
    group_by(region, transition) %>%
    summarise(n = sum(n), .groups = "drop") %>%
    group_by(region) %>%
    mutate(prop = n / sum(n))
}

# Helper to compute entropy
compute_entropy <- function(data) {
  regions <- names(data)[grepl("z_score", names(data))]
  
  entropies <- sapply(regions, function(region) {
    tab <- table(data[[region]])
    probs <- tab / sum(tab)
    -sum(probs * log2(probs + 1e-9))  # Avoid log(0)
  })
  
  tibble(
    label = gsub("_z_score", "", names(entropies)),
    entropy = entropies
  )
}

# Helper to compute mean deviation proportions
compute_dev_prop <- function(data) {
  data %>%
    select(eid, session, ends_with("z_score")) %>%
    pivot_longer(cols = ends_with("z_score"), names_to = "region", values_to = "value") %>%
    filter(!is.na(value)) %>%
    group_by(region) %>%
    summarise(mean_prop = mean(value, na.rm = TRUE)) %>%
    mutate(label = gsub("_z_score", "", region))
}

# Compute all metrics
transitions <- compute_transitions(binary_df)
entropy_df <- compute_entropy(binary_df)
deviation_prop <- compute_dev_prop(binary_df)

# ---- Function to plot ggseg
plot_metric <- function(df, value_col,low,high, title = "") {
  # Cortical
  test_cort <- z_cortical %>%
    filter(group == "MS") %>%
    left_join(df, by = "label")
  
  p1 <- ggplot(dk %>% as_tibble() %>%
                 left_join(test_cort %>% select(region, hemi, all_of(value_col)), by = c("region", "hemi"))) +
    geom_brain(atlas = dk, aes_string(fill = value_col), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(low,high)) +
    theme_void() + ggtitle(title) + theme(legend.position = "none")
  
  # Subcortical
  df$region.x <- df$label
  test_sub <- z_subcortical %>%
    filter(group == "MS") %>%
    left_join(df, by = "region.x")
  
  p2 <- ggplot(test_sub) +
    geom_brain(atlas = aseg, side = "coronal", aes_string(fill = value_col), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",limits = c(low,high)) +
    theme_void() + 
    labs(fill = title)
  
  ggarrange(p1, p2, widths = c(2, 0.5))
}
# Plot each metric
entropy_plot <- plot_metric(entropy_df, "entropy", 0, 1, "Entropy")
deviation_prop = deviation_prop %>% select(label, mean_prop)
deviation_plot <- plot_metric(deviation_prop, "mean_prop",0,.25, "Mean Deviation Proportion")

# Plot each transition type
worsen_plot  <- plot_metric(transitions %>% filter(transition == "worsen") %>% rename(label = region), "prop",0,.075,  "Worsen Rate")
improve_plot <-  plot_metric(transitions %>% filter(transition == "improve") %>% rename(label = region), "prop",0,.075, "Improve Rate")
stable_plot  <- plot_metric(transitions %>% filter(transition == "stable") %>% rename(label = region),"prop",0.925,1,"Stable Rate")

# check stability closer
test = transitions %>% filter(transition == "stable")
test[order(test$prop),]
test[order(test$prop,decreasing = F),] # we are interested in the least stable regions
mean(test$prop)
sd(test$prop)
test = transitions %>% filter(transition == "worsen")
test[order(test$prop,decreasing = T),] # least stable AND decrasing is interesting
mean(test$prop)
sd(test$prop)
test = transitions %>% filter(transition == "improve")
test[order(test$prop,decreasing = T),] # least stable AND decrasing is interesting
mean(test$prop)
sd(test$prop)


