# Application of Brain Reference to cross sectional and longitudinal MS data
# Max Korbmacher, April 2025
#
# ----------------------------------- #
# --------------Structure------------ #
# ----------------------------------- #
# 0.1 Data wrangling----------------- #
# 0.2 Describing asymmetry in MS ---- #
# 1. Number of deviations------------ #
# 2. Longitudinal assessment--------- #
# 2.1 Development of deviations------ #
# 3. Cross-sectional associations---- #
# 4. Individual-level assessments---- #
# 4.1 At baseline-------------------- #
# 4.2 Longitudinally----------------- #
# ----------------------------------- #
# ----------------------------------- #
#
# wash your hands before eating
rm(list = ls(all.names = TRUE)) # clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.
# define the savepath
savepath = "/Users/max/Documents/Local/MS/NormativeModels/Asymmetry/results/"
#
# 0.1 Data wrangling----------------
# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,psych,effsize,ggseg,patchwork,rstatix,ggpubr,
               caret,lme4,lmerTest,haven,reshape2,stats,entropy,
               ggseg3d, ggsegExtra, mice, pwr,sjPlot)
# power
pwr.f2.test(u = 3, v = 43, f2 = NULL, sig.level = 0.05, power = .80)

# load data
cross = read.csv(paste(savepath,"Zscores_cross.csv",sep=""))
long = read.csv(paste(savepath,"Zscores_long.csv",sep=""))
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


long_surrogate = long[long$eid %in% long_no_na$eid,]
long_surrogate$eid <- as.integer(as.factor(long_surrogate$eid))  # ensures grouping is preserved

long <- long_surrogate %>%
  left_join(average_imputed_edss, by = c("eid", "session", "age")) %>%
  mutate(edss = if_else(is.na(edss), mean_edss, edss)) %>%
  select(-mean_edss)
long = na.omit(long)
long$eid = EID
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

# 0.2 Describing asymmetry in MS ----
asym = cross %>% select(eid, scanner, age, sex, diagnosis, data, EstimatedTotalIntraCranialVol,
                        contains("lh_") & !contains("z_score") & !contains("predicted") |                           
                          contains("Left") & !contains("z_score") & !contains("predicted"))
# estimate Cohen's d
plotit <- function(diagnosis_label) {
  cohen = cohen_low = cohen_high = c()
  
  # Extract only relevant data columns for given diagnosis
  dat = asym %>% 
    filter(diagnosis == diagnosis_label) %>%
    select(-c(eid, scanner, diagnosis, age, sex, data, EstimatedTotalIntraCranialVol))
  
  # Calculate Cohen's d for each variable
  for (i in seq_along(dat)) {
    d_val = cohen.d(dat[[i]], f = NA, mu = 0)
    cohen[i]      = as.numeric(d_val$estimate)
    cohen_low[i]  = as.numeric(d_val$conf.int[1])
    cohen_high[i] = as.numeric(d_val$conf.int[2])
  }
  
  # Assemble result table
  the_data = data.frame(
    label   = names(dat),
    d       = cohen,
    d_high  = cohen_high,
    d_low   = cohen_low,
    hemi    = "left"
  )
  
  # ──────────────────────────────────────────────
  # Cortical Plot
  # ──────────────────────────────────────────────
  plot_df = the_data %>% 
    filter(grepl("lh_", label)) %>%
    mutate(region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume", "", label)][1:34])
  
  cortical_HC = ggplot(
    dk %>% as_tibble() %>%
      left_join(plot_df, by = c("region", "hemi")) %>%
      filter(hemi == "left")
  ) + 
    geom_brain(atlas = dk, hemi = "left", aes(fill = d), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-2, 12)) +
    theme_void()
  
  # ──────────────────────────────────────────────
  # Subcortical Plot
  # ──────────────────────────────────────────────
  
  # Subcortical data: keep only "Left" labels
  plot_df1 = the_data %>%
    filter(grepl("^Left", label)) %>%
    arrange(d)
  
  # Fix naming: convert to match aseg atlas
  label_map_sub <- function(x) {
    x <- gsub("\\.", "-", x)
    if (x == "Left-Thalamus") x <- "Left-Thalamus-Proper"
    x
  }
  plot_df1$label <- vapply(plot_df1$label, label_map_sub, character(1))
  
  # Pull atlas regions that can be drawn
  coronal_brain_aseg <- as_tibble(aseg) %>%
    filter(side == "coronal") %>%
    select(label, geometry) %>%
    distinct()
  
  # Merge to keep only matched subcortical structures
  plot_df1 <- inner_join(plot_df1, coronal_brain_aseg, by = "label")
  
  # Optional warning if anything still failed
  missing_now <- setdiff(
    vapply(plot_df1$label, identity, character(1)),
    coronal_brain_aseg$label
  )
  if (length(missing_now) > 0) {
    warning("Still missing in atlas: ",
            paste(unique(missing_now), collapse = ", "))
  }
  
  # Subcortical brain plot
  subcortical_HC <- ggplot(plot_df1) +
    geom_brain(
      atlas = aseg,
      side  = "coronal",
      aes(fill = d),
      color = "black"
    ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         limits = c(-2, 12)) +
    theme_void()
  
  # Combine and return
  p_HC = ggpubr::ggarrange(
    cortical_HC, subcortical_HC,
    nrow = 1,
    widths = c(2, 0.5),
    common.legend = TRUE,
    legend = "right"
  )
  
  return(list(p_HC, plot_df, plot_df1))
}
leftwards_d = ggarrange(plotit("MS")[[1]],plotit("HC")[[1]],nrow=2, labels = c("MS","HC"))
ggsave(paste(savepath,"Cohens_d_leftwards_asym.pdf",sep=""),leftwards_d, width = 14, height = 6)
# inspect data
plotit("MS")[[2]][order(data.frame(plotit("MS")[[2]])$d),] # MS cortical
plotit("MS")[[3]][order(data.frame(plotit("MS")[[3]])$d),] # MS subcortical
plotit("HC")[[2]][order(data.frame(plotit("HC")[[2]])$d),] # HC cortical
plotit("HC")[[3]][order(data.frame(plotit("HC")[[3]])$d),] # HC subcortical


# 1.1 number of deviations---------
cross$nb_deviations = cross %>%
  select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(abs(.) > 1.96, 1, 0)) %>% # mutate_all(~ ifelse(abs(.) >= 1.96, 1, 0)) %>%
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
#diffs %>% summarise(mean_abs = mean(abs(difference)),sd_abs = sd(abs(difference)))
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
z_cortical <- z_long %>% filter(str_detect(region, "^lh_")) 
#%>%mutate(region = str_remove(region, "^lh_|^rh_"))
z_subcortical <- z_long %>% filter(!str_detect(region, "^lh_"))



# order the data frame
z_cortical = z_cortical[order(z_cortical$region),]
z_subcortical = z_subcortical[order(z_subcortical$region),] 

# label data correctly
z_cortical = (z_cortical %>% group_by(group)) %>% mutate(label = c(((brain_regions(dk)[brain_labels(dk) %in% gsub("_volume","",z_cortical$region)][1:34]))))
z_subcortical = (z_subcortical %>% group_by(group)) %>% mutate(label = (c(brain_labels(aseg)[grepl(c("halamus|allidum|mygdala|campus|utamen|audate"),brain_labels(aseg))],"Left-Cerebellum-Cortex","Right-Cerebellum-Cortex"))[order(c(brain_labels(aseg)[grepl(c("Left-Thalamus|Left-Pallidum|Left-Amygdala|Left-Hippocampus|Left-Putamen|Left-Caudate"),brain_labels(aseg))],"Left-Cerebellum-Cortex"))])
# add hemi
z_cortical$hemi = ifelse(grepl("lh_",z_cortical$region)==T,"left","right")
z_subcortical$hemi = ifelse(grepl("Left.",z_subcortical$region)==T,"left","right")

names(z_cortical) = c("label","group","Z","region","hemi")
names(z_subcortical) = c("region","group","Z","label","hemi")
#
# 2. Longitudinal assessment------- 
# 2.1 Nb of deviations, age, and clinical scores----
# 2.1.1 number of deviations and edss----
long$nb_deviations = long %>%
  select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(abs(.) >= 1.96, 1, 0)) %>%
  transmute(z_score_sum = rowSums(.)) %>%
  pull(z_score_sum)
m = lmer(edss~nb_deviations+age+sex+(1|eid),long)
summary(m)
m = lm(edss~nb_deviations+age+sex,long%>%filter(session == 1))
summary(m)
print("The number of deviations, measured by EDSS, does not tell us anything.")

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
pasat1$session = ifelse(pasat1$session == "BL_PASATcorrect",1,0)+ifelse(pasat1$session == "PASAT_24M",25,0)+ifelse(pasat1$session == "PASAT_OFAMS10",145,0)

long1 = merge(rbind(PASAT_OSL,pasat1),long,by=c("eid","session"))
m = lmer(PASAT~nb_deviations+age+sex+(1|eid),long1)
summary(m)
m = lm(PASAT~nb_deviations+age+sex,long1%>%filter(session == 1))
summary(m)
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
m = lm(fatigue~nb_deviations+age+sex,long2%>%filter(session == 1))
summary(m)

# 2.1.4 number of deviations and age(ing)----
m = lm(nb_deviations~age+sex,long%>%filter(session == 1))
m = lmer(nb_deviations~age+sex+(1|eid),long)
summary(m)
#
#
# 2.2 Regional associations ----
# Effects of Z-scores ----
# ---- Function to extract standardized coefficients manually using lmer()
plotit <- function(data, outcome = "age", regions) {
  message("Outcome (response) variable: ", outcome)
  data <- data %>%
    mutate(std_outcome = scale(.data[[outcome]])[, 1])
  
  # Initialize vectors for estimates, SEs, CIs, p-values
  estimates <- numeric(length(regions))
  ses       <- numeric(length(regions))
  ci_low    <- numeric(length(regions))
  ci_high   <- numeric(length(regions))
  p_values  <- numeric(length(regions))
  
  for (i in seq_along(regions)) {
    region <- regions[i]
    data <- data %>%
      mutate(std_predictor = scale(.data[[region]])[, 1])
    
    mod <- lmer(std_outcome ~ std_predictor + (1 | eid), data = data, REML = FALSE)
    
    est <- fixef(mod)["std_predictor"]
    se <- sqrt(diag(vcov(mod)))["std_predictor"]
    
    # 95% CI assuming normality
    ci_l <- est - 1.96 * se
    ci_h <- est + 1.96 * se
    
    # Wald z for approximate p-value
    z_val <- est / se
    p_val <- 2 * (1 - pnorm(abs(z_val)))
    
    estimates[i] <- est
    ses[i] <- se
    ci_low[i] <- ci_l
    ci_high[i] <- ci_h
    p_values[i] <- p_val
  }
  the_data <- tibble(
    label = str_remove(regions, "_z_score"),
    Std.Coeff. = estimates,
    SE = ses,
    CI_low = ci_low,
    CI_high = ci_high,
    p_value = p_values,
    hemi = "left"
  )
  # ──────────────────────────────────────────────
  # Cortical Plot
  # ──────────────────────────────────────────────
  cortical_df <- the_data %>% 
    filter(grepl("^lh_", label)) %>%
    mutate(region = brain_regions(dk)[brain_labels(dk) %in% gsub("_volume", "", label)][1:34])
  
  cortical_HC <- ggplot(
    dk %>% as_tibble() %>%
      left_join(cortical_df, by = c("region", "hemi")) %>%
      filter(hemi == "left")
  ) + 
    geom_brain(atlas = dk, hemi = "left", aes(fill = Std.Coeff.), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.25, 0.25), name = "Std.Coeff.") +
    theme_void()
  
  # ──────────────────────────────────────────────
  # Subcortical Plot
  # ──────────────────────────────────────────────
  subcortical_df <- the_data %>%
    filter(grepl("^Left", label)) %>%
    arrange(Std.Coeff.)
  
  label_map_sub <- function(x) {
    x <- gsub("\\.", "-", x)
    if (x == "Left-Thalamus") x <- "Left-Thalamus-Proper"
    x
  }
  subcortical_df$label <- vapply(subcortical_df$label, label_map_sub, character(1))
  
  coronal_brain_aseg <- as_tibble(aseg) %>%
    filter(side == "coronal") %>%
    select(label, geometry) %>%
    distinct()
  
  subcortical_df <- inner_join(subcortical_df, coronal_brain_aseg, by = "label")
  
  missing_now <- setdiff(subcortical_df$label, coronal_brain_aseg$label)
  if (length(missing_now) > 0) {
    warning("Still missing in atlas: ", paste(unique(missing_now), collapse = ", "))
  }
  
  subcortical_HC <- ggplot(subcortical_df) +
    geom_brain(atlas = aseg, side = "coronal", aes(fill = Std.Coeff.), color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits = c(-0.25, 0.25), name = "Std.Coeff.") +
    theme_void()
  
  combined_plot <- ggpubr::ggarrange(
    cortical_HC, subcortical_HC,
    nrow = 1,
    widths = c(2, 0.5),
    common.legend = TRUE,
    legend = "right"
  )
  
  return(list(combined_plot, cortical_df, subcortical_df))
}

# ---- Application of functions
## plotting
regions =  long %>% select(ends_with("z_score")) %>% names()
edss_plot = plotit(long, outcome = "edss", regions)[[1]]
age_plot = plotit(long, outcome = "age", regions)[[1]]
pasat_plot = plotit(long1, outcome = "PASAT", regions)[[1]]
fatigue_plot = plotit(long2, outcome = "fatigue", regions)[[1]]
large_plot = ggarrange(age_plot, edss_plot,
                       pasat_plot, fatigue_plot, 
                       ncol=1,labels=c("Age","EDSS","PASAT","Fatigue"),
                       hjust = c(0,0,0,0))
ggsave(paste(savepath,"Long_Z_Effects.pdf",sep=""),plot=large_plot, width = 10, height = 8)
#
## tables
edss_tab = rbind(plotit(long, outcome = "edss", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                 plotit(long, outcome = "edss", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
edss_tab$p.adjust = p.adjust(edss_tab$p_value,method = "fdr")
edss_tab %>% filter(p.adjust < 0.05)


age_tab = age_plot = rbind(plotit(long, outcome = "age", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                           plotit(long, outcome = "age", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
age_tab$p.adjust = p.adjust(age_tab$p_value,method = "fdr")
age_tab %>% filter(p.adjust < 0.05)

pasat_tab = rbind(plotit(long1, outcome = "PASAT", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                  plotit(long1, outcome = "PASAT", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))

pasat_tab$p.adjust = p.adjust(pasat_tab$p_value,method = "fdr")
pasat_tab %>% filter(p.adjust < 0.05)

fatigue_tab = rbind(plotit(long2, outcome = "fatigue", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                    plotit(long2, outcome = "fatigue", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
fatigue_tab$p.adjust = p.adjust(fatigue_tab$p_value,method = "fdr")
fatigue_tab %>% filter(p.adjust < 0.05)
sig = rbind(age_tab, edss_tab, pasat_tab, fatigue_tab) %>% filter(p.adjust < 0.05)
write.csv(sig, paste0(savepath, "Z_preds.csv"), row.names = FALSE)
#
# plot strongest effects of Zs ----
plot_lmer_with_lines <- function(dep_var, indep_var, random_effect, data,
                                 title = NULL, xlab_text = NULL, ylab_text = NULL,
                                 text_size = 14, flip = FALSE) {
  
  formula <- as.formula(paste(dep_var, "~", indep_var, "+ (1|", random_effect, ")"))
  model <- lmer(formula, data = data)
  data$predicted <- predict(model)
  
  p <- plot_model(model, type = "pred", terms = indep_var,title = NULL) +
    xlab(ifelse(is.null(xlab_text), indep_var, xlab_text)) +
    ylab(ifelse(is.null(ylab_text), dep_var, ylab_text)) +
    ggtitle(label = NULL) +
    #ggtitle(title %||% paste("Predicted", dep_var, "by", indep_var)) +
    theme_bw() +
    theme(
      legend.position = "none",
      text = element_text(size = text_size)
    ) +
    geom_line(
      data = data,
      aes_string(x = indep_var, y = dep_var, group = random_effect),
      color = "gray70", alpha = 0.5
    ) +
    geom_point(
      data = data,
      aes_string(x = indep_var, y = dep_var, group = random_effect),
      color = "black", size = 1, alpha = 0.6
    )
  
  if (flip) p <- p + coord_flip()
  
  return(p)
}
## Age
strong_a = plot_lmer_with_lines(
  dep_var = "age",
  indep_var = "lh_superiorparietal_volume_z_score",
  random_effect = "eid",
  data = long,
  #title = "Predicted Age",
  xlab_text = "Superior Parietal Asymmetry Z Score",
  ylab_text = "Age",
  text_size = 16,
  flip = T
)
ggsave(paste(savepath,"strong_a.pdf",sep=""),plot=strong_a, width = 8, height = 6)

strong_b = plot_lmer_with_lines(
  dep_var = "age",
  indep_var = "Left.Thalamus_z_score",
  random_effect = "eid",
  data = long,
  #title = "Predicted Age by Thalamus Asymmetry Z Score",
  xlab_text = "Thalamus Z Score",
  ylab_text = "Age",
  text_size = 16,
  flip = T
)
ggsave(paste(savepath,"strong_b.pdf",sep=""),plot=strong_b, width = 8, height = 6)

strong_c = plot_lmer_with_lines(
  dep_var = "age",
  indep_var = "lh_rostralanteriorcingulate_volume_z_score",
  random_effect = "eid",
  data = long,
  #title = "Predicted Age",
  xlab_text = "Rostral Anterior Cingulate Asymmetry Z Score",
  ylab_text = "Age",
  text_size = 16,
  flip = T
)
ggsave(paste(savepath,"strong_c.pdf",sep=""),plot=strong_c, width = 8, height = 6)

strong_d = plot_lmer_with_lines(
  dep_var = "edss",
  indep_var = "lh_rostralanteriorcingulate_volume_z_score",
  random_effect = "eid",
  data = long,
  #title = "EDSS",
  xlab_text = "Lateral Orbitofrontal Asymmetry Z Score",
  ylab_text = "EDSS",
  text_size = 16,
  flip = T
)
ggsave(paste(savepath,"strong_d.pdf",sep=""),plot=strong_d, width = 8, height = 6)

strong_e = plot_lmer_with_lines(
  dep_var = "PASAT",
  indep_var = "Left.Caudate_z_score",
  random_effect = "eid",
  data = long1 %>% select(eid,PASAT,Left.Caudate_z_score)%>%na.omit(),
  #title = "PASAT",
  xlab_text = "Caudate Asymmetry Z Score",
  ylab_text = "PASAT",
  text_size = 16,
  flip = T
)
ggsave(paste(savepath,"strong_e.pdf",sep=""),plot=strong_e, width = 8, height = 6)


#
# Effects of raw asymmetries ----
# let's try the same thing with the raw asymmetry 
regions = long %>% select(!contains("z_score")) %>% select(!contains("predicted")) %>%
  select(ends_with("_volume"), contains("halamus"),
         contains("allidum"), contains("mygdala"), 
         contains("campus"), contains("utamen"),
         contains("audate"), contains("bellum")) %>% names()
edss_plot = plotit(long, outcome = "edss", regions)[[1]]
age_plot = plotit(long, outcome = "age", regions)[[1]]
pasat_plot = plotit(long1, outcome = "PASAT", regions)[[1]]
fatigue_plot = plotit(long2, outcome = "fatigue", regions)[[1]]
large_plot = ggarrange(age_plot, edss_plot,
                       pasat_plot, fatigue_plot, 
                       ncol=1,labels=c("Age","EDSS","PASAT","Fatigue"),
                       hjust = c(0,0,0,0))
ggsave(paste(savepath,"Long_RAW_Effects.pdf",sep=""),plot=large_plot, width = 10, height = 8)
## tables
edss_tab = rbind(plotit(long, outcome = "edss", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                 plotit(long, outcome = "edss", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
edss_tab$p.adjust = p.adjust(edss_tab$p_value,method = "fdr")
edss_tab %>% filter(p.adjust < 0.05)

age_tab = age_plot = rbind(plotit(long, outcome = "age", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                           plotit(long, outcome = "age", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
age_tab$p.adjust = p.adjust(age_tab$p_value,method = "fdr")
age_tab %>% filter(p.adjust < 0.05)

pasat_tab = rbind(plotit(long1, outcome = "PASAT", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                  plotit(long1, outcome = "PASAT", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))

pasat_tab$p.adjust = p.adjust(pasat_tab$p_value,method = "fdr")
pasat_tab %>% filter(p.adjust < 0.05)

fatigue_tab = rbind(plotit(long2, outcome = "fatigue", regions)[[3]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value),
                    plotit(long2, outcome = "fatigue", regions)[[2]] %>% select(label, Std.Coeff., CI_low, CI_high, p_value))
fatigue_tab$p.adjust = p.adjust(fatigue_tab$p_value,method = "fdr")
fatigue_tab %>% filter(p.adjust < 0.05)
#
sig = rbind(age_tab, edss_tab, pasat_tab, fatigue_tab) %>% filter(p.adjust < 0.05)
write.csv(sig, paste0(savepath, "raw_preds.csv"), row.names = FALSE)

# plot strongest raw effects ----


# 3. Cross-sectional regional associations--------

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
    left_join(coefs, by = c("region" = "label"))
  
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
large_plot <- ggarrange(age_plot, edss_plot, pasat_plot, fatigue_plot,
                        ncol = 1, labels = c("Age", "EDSS", "PASAT", "Fatigue"),
                        hjust = c(0, 0, 0, 0))
ggsave(paste0(savepath, "Cross_Effects.pdf"), plot = large_plot, width = 10, height = 8)

# ---- Extract Coefficients and P-values for Tables
age_coefs     <- run_and_extract_cross(long %>% filter(session == 1), "age")
edss_coefs    <- run_and_extract_cross(long %>% filter(session == 1), "edss")
pasat_coefs   <- run_and_extract_cross(long1 %>% filter(session == 1), "PASAT")
fatigue_coefs <- run_and_extract_cross(long2 %>% filter(session == 1), "fatigue")

all_coefs <- bind_rows(age_coefs, edss_coefs, pasat_coefs, fatigue_coefs) %>%
  pivot_wider(names_from = predictor, values_from = std_coef)
all_coefs[2:5] <- round(all_coefs[2:5], 2)
write.csv(all_coefs, paste0(savepath, "cross_associations.csv"), row.names = FALSE)

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
write.csv(all_pvals_wide, paste0(savepath, "cross_associations_pvalues.csv"), row.names = FALSE)

all_pvals_uncorrected_wide <- all_pvals_uncorrected %>%
  pivot_wider(names_from = predictor, values_from = p_value)
all_pvals_uncorrected_wide[2:5] <- round(all_pvals_uncorrected_wide[2:5], 3)
write.csv(all_pvals_uncorrected_wide, paste0(savepath, "cross_associations_pvalues_uncorrected.csv"), row.names = FALSE)

######### follow-up for CIs: Cross-sectional data
## Age associations
x = (lm(scale(rh_caudalmiddlefrontal_volume)~scale(age),data = long %>% filter(session == 1)))
summary(x)
confint(x)
## EDSS assoctiations
x = (lm(scale(Left.Thalamus)~scale(edss),data = long %>% filter(session == 1)))
summary(x)
confint(x)
x = (lm(scale(Right.Thalamus)~scale(edss),data = long %>% filter(session == 1)))
summary(x)
confint(x)

x = (lm(scale(Right.Putamen)~scale(edss),data = long %>% filter(session == 1)))
summary(x)
confint(x)

# FSS
x = (lm(scale(rh_superiorfrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
summary(x)
confint(x)
x = (lm(scale(lh_superiorfrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
summary(x)
confint(x)
x = (lm(scale(rh_caudalmiddlefrontal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
summary(x)
confint(x)
x = (lm(scale(rh_paracentral_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
summary(x)
confint(x)
x = (lm(scale(lh_supramarginal_volume)~scale(fatigue),data = long2 %>% filter(session == 1)))
summary(x)
confint(x)

######### follow-up for CIs: Longitudinal data
## Age associations
#
#

# 4. Individual-level assessments----
#
# 4.1 Baseline ----
#
test = cross %>% filter(diagnosis == "MS") %>% select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(abs(.) > 1.96, 1, 0)) %>% 
  colSums() / nrow(cross %>% filter(diagnosis == "MS"))
test[order(data.frame(test)$test)]

max(cross %>% filter(diagnosis == "HC") %>% select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(abs(.) > 1.96, 1, 0)) %>% 
  colSums() / nrow(cross %>% filter(diagnosis == "HC")))
cross %>% filter(diagnosis == "HC") %>% select(ends_with("z_score")) %>%
  mutate_all(~ ifelse(abs(.) > 1.96, 1, 0)) %>% 
  colSums() / nrow(cross %>% filter(diagnosis == "HC"))


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
ggsave(paste(savepath,"regional_heterogeniety.pdf",sep=""),plot=baseline, width = 10, height = 2)
relative_numbers = data.frame(perc = test, label = gsub("_z_score","",names(test)))


# 4.2 Longitudinally ----
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


