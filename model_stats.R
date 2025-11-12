# Normative model performance
# 12 Nov 2025
# Max Korbmacher: max.korbmacher@gmail.com
#
# --------------------------- #
# ------- Contents ---------- #
# 0. Prep ------------------- #
# 1. Internal validation ---- #
# 2. External validation ---- #
# --------------------------- #

# 0. Prep --------------------
rm(list = ls(all.names = TRUE))
gc()

savepath <- "/Users/max/Documents/Local/MS/NormativeModels/Regular/results/"
model_path <- paste0(savepath, "models/")
zscore_file <- "/Users/max/Documents/Local/MS/NormativeModels/Regular/code/Zscores.csv"

if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, mfp, neuroCombat, MatchIt)

# Load data
df <- read.csv("/Users/max/Documents/Local/Data/Lifespan/cortical_subcortical.csv")
ms <- read.csv("/Users/max/Documents/Local/Data/Lifespan/MS_long.csv")

# Data wrangling
names(df)[names(df) == "Left.Thalamus.Proper"] <- "Left.Thalamus"
names(df)[names(df) == "Right.Thalamus.Proper"] <- "Right.Thalamus"

cross <- df %>%
  group_by(eid) %>%
  filter(session == min(session)) %>%
  ungroup() %>%
  filter(diagnosis == "HC")

cross$sex <- ifelse(cross$sex %in% c("F", "Female"), "female", "male")

ms$diagnosis <- "MS"
cross$diagnosis <- "HC"

ms_cross <- ms %>%
  group_by(eid) %>%
  filter(session == min(session)) %>%
  ungroup()

cross <- rbind(ms_cross %>% select(names(cross)), cross)
cross <- na.omit(cross)

# Matching
m.out <- matchit(factor(diagnosis) ~ age + sex + EstimatedTotalIntraCranialVol,
                 data = cross,
                 method = "nearest",
                 distance = "glm")
m.out <- match_data(m.out)

# Filter out test data
cross <- cross[!cross$eid %in% m.out$eid, ]
cross <- cross %>% filter(session != 2 & session != 3 & diagnosis == "HC")

# --------------------------- #
# Metric functions
compute_metrics <- function(y_true, preds) {
  y_true <- as.numeric(y_true)
  preds <- as.numeric(preds)
  rmse <- sqrt(mean((y_true - preds)^2))
  mae <- mean(abs(y_true - preds))
  r2 <- 1 - sum((y_true - preds)^2) / sum((y_true - mean(y_true))^2)
  corr <- cor(y_true, preds, method = "pearson")
  return(c(RMSE = rmse, MAE = mae, R2 = r2, Correlation = corr))
}

bootstrap_ci <- function(y_true, preds, n_boot = 1000) {
  y_true <- as.numeric(y_true)
  preds <- as.numeric(preds)
  metrics_list <- replicate(n_boot, {
    idx <- sample(seq_along(y_true), replace = TRUE)
    compute_metrics(y_true[idx], preds[idx])
  }, simplify = "array")
  
  ci <- apply(metrics_list, 1, function(x) quantile(x, probs = c(0.025, 0.975)))
  ci <- t(ci)
  rownames(ci) <- names(compute_metrics(y_true, preds))
  return(ci)
}

# --------------------------- #
# 1. Internal validation -----
model_files <- list.files(model_path, pattern = "\\.Rdata$", full.names = TRUE)

internal_stats <- data.frame(Model = character(), Region = character(), Sex = character(),
                             RMSE = numeric(), RMSE_Lower = numeric(), RMSE_Upper = numeric(),
                             MAE = numeric(), MAE_Lower = numeric(), MAE_Upper = numeric(),
                             R2 = numeric(), R2_Lower = numeric(), R2_Upper = numeric(),
                             Correlation = numeric(), Corr_Lower = numeric(), Corr_Upper = numeric(),
                             stringsAsFactors = FALSE)

for (file in model_files) {
  obj_names <- load(file)
  model <- get(obj_names[1])
  
  fname <- basename(file)
  fname_clean <- gsub("\\.Rdata$", "", fname)                # remove .Rdata
  region <- gsub("_(male|female)$", "", fname_clean)         # remove sex suffix
  sex0 <- ifelse(grepl("male", fname), "male", "female")
  
  if (!(region %in% names(cross))) {
    message(paste("Skipping model:", fname, "- region not found in data"))
    next
  }
  
  train_data <- cross %>% filter(sex == sex0)
  preds <- tryCatch(predict(model, newdata = train_data), error = function(e) NULL)
  
  if (is.null(preds)) {
    message(paste("Skipping model:", fname, "- prediction failed"))
    next
  }
  
  y_true <- as.numeric(train_data[[region]])
  preds <- as.numeric(preds)
  
  valid_idx <- which(!is.na(y_true) & !is.na(preds))
  y_true <- y_true[valid_idx]
  preds <- preds[valid_idx]
  
  if (length(y_true) > 1 && length(preds) == length(y_true)) {
    metrics <- compute_metrics(y_true, preds)
    ci <- bootstrap_ci(y_true, preds)
    
    internal_stats <- rbind(internal_stats,
                            data.frame(Model = fname, Region = region, Sex = sex0,
                                       RMSE = metrics["RMSE"], RMSE_Lower = ci["RMSE", 1], RMSE_Upper = ci["RMSE", 2],
                                       MAE = metrics["MAE"], MAE_Lower = ci["MAE", 1], MAE_Upper = ci["MAE", 2],
                                       R2 = metrics["R2"], R2_Lower = ci["R2", 1], R2_Upper = ci["R2", 2],
                                       Correlation = metrics["Correlation"], Corr_Lower = ci["Correlation", 1], Corr_Upper = ci["Correlation", 2]))
  } else {
    message(paste("Skipping model:", fname, "- incompatible lengths"))
  }
}

internal_stats = internal_stats %>%
  rowwise() %>%
  mutate(Mean_Original = mean(cross %>% filter(sex == Sex) %>% pull(Region), na.rm = TRUE)) %>%
  select(Model, Region, Sex, Mean_Original, everything()) %>% select(-Sex)

write.csv(internal_stats, file = paste0(savepath, "internal_validation_stats.csv"), row.names = FALSE)


internal_stats <- internal_stats %>%
  mutate(RMSE_pct = (RMSE / Mean_Original) * 100,
         MAE_pct = (MAE / Mean_Original) * 100)

# Summarize in a table
summary_long <- tibble(
  Metric = c("RMSE (%)", "MAE (%)", "R²", "Correlation"),
  Mean   = c(mean(internal_stats$RMSE_pct, na.rm = TRUE),
             mean(internal_stats$MAE_pct, na.rm = TRUE),
             mean(internal_stats$R2, na.rm = TRUE),
             mean(internal_stats$Correlation, na.rm = TRUE)),
  SD     = c(sd(internal_stats$RMSE_pct, na.rm = TRUE),
             sd(internal_stats$MAE_pct, na.rm = TRUE),
             sd(internal_stats$R2, na.rm = TRUE),
             sd(internal_stats$Correlation, na.rm = TRUE)),
  Min    = c(min(internal_stats$RMSE_pct, na.rm = TRUE),
             min(internal_stats$MAE_pct, na.rm = TRUE),
             min(internal_stats$R2, na.rm = TRUE),
             min(internal_stats$Correlation, na.rm = TRUE)),
  Max    = c(max(internal_stats$RMSE_pct, na.rm = TRUE),
             max(internal_stats$MAE_pct, na.rm = TRUE),
             max(internal_stats$R2, na.rm = TRUE),
             max(internal_stats$Correlation, na.rm = TRUE))
)

print(summary_long)


# Save summary table
write.csv(summary_long, paste0(savepath, "summary_metrics_internal.csv"), row.names = FALSE)

# # --------------------------- #
# # 2. External validation -----
# # External validation summary
# zscore_data <- read.csv(zscore_file)
# 
# # Identify region columns
# original_cols <- grep("_volume$", names(zscore_data), value = TRUE)
# predicted_cols <- grep("_volume_predicted$", names(zscore_data), value = TRUE)
# 
# # Ensure matching pairs
# regions <- gsub("_volume$", "", original_cols)
# 
# external_stats_list <- list()
# 
# for (region in regions) {
#   orig_col <- paste0(region, "_volume")
#   pred_col <- paste0(region, "_volume_predicted")
#   
#   y_true <- as.numeric(zscore_data[[orig_col]])
#   preds <- as.numeric(zscore_data[[pred_col]])
#   
#   valid_idx <- which(!is.na(y_true) & !is.na(preds))
#   y_true <- y_true[valid_idx]
#   preds <- preds[valid_idx]
#   
#   if (length(y_true) > 1 && length(preds) == length(y_true)) {
#     metrics <- compute_metrics(y_true, preds)
#     mean_orig <- mean(y_true, na.rm = TRUE)
#     
#     external_stats_list[[region]] <- tibble(
#       Region = region,
#       Mean_Original = mean_orig,
#       RMSE_pct = (metrics["RMSE"] / mean_orig) * 100,
#       MAE_pct = (metrics["MAE"] / mean_orig) * 100,
#       R2 = metrics["R2"],
#       Correlation = metrics["Correlation"]
#     )
#   }
# }
# 
# external_stats <- bind_rows(external_stats_list)
# 
# # Summarize across regions
# summary_external <- tibble(
#   Metric = c("RMSE (%)", "MAE (%)", "R²", "Correlation"),
#   Mean = c(mean(external_stats$RMSE_pct, na.rm = TRUE),
#            mean(external_stats$MAE_pct, na.rm = TRUE),
#            mean(external_stats$R2, na.rm = TRUE),
#            mean(external_stats$Correlation, na.rm = TRUE)),
#   SD = c(sd(external_stats$RMSE_pct, na.rm = TRUE),
#          sd(external_stats$MAE_pct, na.rm = TRUE),
#          sd(external_stats$R2, na.rm = TRUE),
#          sd(external_stats$Correlation, na.rm = TRUE)),
#   Min = c(min(external_stats$RMSE_pct, na.rm = TRUE),
#           min(external_stats$MAE_pct, na.rm = TRUE),
#           min(external_stats$R2, na.rm = TRUE),
#           min(external_stats$Correlation, na.rm = TRUE)),
#   Max = c(max(external_stats$RMSE_pct, na.rm = TRUE),
#           max(external_stats$MAE_pct, na.rm = TRUE),
#           max(external_stats$R2, na.rm = TRUE),
#           max(external_stats$Correlation, na.rm = TRUE))
# )
# 
# print(summary_external)
# 
# # Save outputs
# write.csv(external_stats, paste0(savepath, "external_regionwise_metrics.csv"), row.names = FALSE)
# write.csv(summary_external, paste0(savepath, "summary_metrics_external.csv"), row.names = FALSE)
# 
# 
# all_summary = rbind(summary_long,summary_external)
# all_summary$Data = c("training","training","training","training",
#                      "test", "test", "test", "test")
# write.csv(all_summary, paste0(savepath, "all_summary_metrics.csv"), row.names = FALSE)







# 2. External validation -----
# External validation summary
zscore_data <- read.csv(zscore_file)

# Identify original and predicted columns
original_cols <- grep("^(Left\\.|Right\\.)", names(zscore_data), value = TRUE)
predicted_cols <- grep("_predicted$", names(zscore_data), value = TRUE)

# Extract region names for predicted columns
regions_predicted <- gsub("_predicted$", "", predicted_cols)

external_stats_list <- list()

for (region in regions_predicted) {
  orig_col <- region
  pred_col <- paste0(region, "_predicted")
  
  y_true <- if (orig_col %in% names(zscore_data)) as.numeric(zscore_data[[orig_col]]) else NA
  preds <- if (pred_col %in% names(zscore_data)) as.numeric(zscore_data[[pred_col]]) else NA
  
  valid_idx <- which(!is.na(y_true) & !is.na(preds))
  y_true <- y_true[valid_idx]
  preds <- preds[valid_idx]
  
  if (length(y_true) > 1 && length(preds) == length(y_true)) {
    metrics <- compute_metrics(y_true, preds)
    mean_orig <- mean(y_true, na.rm = TRUE)
  } else {
    metrics <- c(RMSE = NA, MAE = NA, R2 = NA, Correlation = NA)
    mean_orig <- NA
  }
  
  external_stats_list[[region]] <- tibble(
    Region = region,
    Mean_Original = mean_orig,
    RMSE_pct = if (!is.na(metrics["RMSE"]) && !is.na(mean_orig)) (metrics["RMSE"] / mean_orig) * 100 else NA,
    MAE_pct = if (!is.na(metrics["MAE"]) && !is.na(mean_orig)) (metrics["MAE"] / mean_orig) * 100 else NA,
    R2 = metrics["R2"],
    Correlation = metrics["Correlation"]
  )
}

external_stats <- bind_rows(external_stats_list)

# Summarize across regions
summary_external <- tibble(
  Metric = c("RMSE (%)", "MAE (%)", "R²", "Correlation"),
  Mean = c(mean(external_stats$RMSE_pct, na.rm = TRUE),
           mean(external_stats$MAE_pct, na.rm = TRUE),
           mean(external_stats$R2, na.rm = TRUE),
           mean(external_stats$Correlation, na.rm = TRUE)),
  SD = c(sd(external_stats$RMSE_pct, na.rm = TRUE),
         sd(external_stats$MAE_pct, na.rm = TRUE),
         sd(external_stats$R2, na.rm = TRUE),
         sd(external_stats$Correlation, na.rm = TRUE)),
  Min = c(min(external_stats$RMSE_pct, na.rm = TRUE),
          min(external_stats$MAE_pct, na.rm = TRUE),
          min(external_stats$R2, na.rm = TRUE),
          min(external_stats$Correlation, na.rm = TRUE)),
  Max = c(max(external_stats$RMSE_pct, na.rm = TRUE),
          max(external_stats$MAE_pct, na.rm = TRUE),
          max(external_stats$R2, na.rm = TRUE),
          max(external_stats$Correlation, na.rm = TRUE))
)

print(summary_external)

# Save outputs
write.csv(external_stats, paste0(savepath, "external_regionwise_metrics.csv"), row.names = FALSE)
write.csv(summary_external, paste0(savepath, "summary_metrics_external.csv"), row.names = FALSE)

