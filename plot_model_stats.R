# Plot model stats
# 12 Nov 2025
# Max Korbmacher: max.korbmacher@gmail.com
#

library(tidyverse)

savepath <- "/Users/max/Documents/Local/MS/NormativeModels/Regular/results/"

# Load internal stats
internal_stats <- read.csv(paste0(savepath, "internal_validation_stats.csv"))
external_stats <- read.csv(paste0(savepath, "external_regionwise_metrics.csv"))
print(internal_stats)
# Compute percentage metrics if not already present
internal_stats <- internal_stats %>%
  mutate(RMSE_pct = (RMSE / Mean_Original) * 100,
         MAE_pct = (MAE / Mean_Original) * 100)

# Prepare long format for internal
internal_long <- internal_stats %>%
  select(Region, RMSE_pct, MAE_pct, R2, Correlation) %>%
  pivot_longer(cols = c(RMSE_pct, MAE_pct, R2, Correlation),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "Internal")

# Prepare long format for external
external_long <- external_stats %>%
  select(Region, RMSE_pct, MAE_pct, R2, Correlation) %>%
  pivot_longer(cols = c(RMSE_pct, MAE_pct, R2, Correlation),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "External")

# Combine
combined_long <- bind_rows(internal_long, external_long)

# Compute mean and SD for error bars
summary_data <- combined_long %>%
  group_by(Dataset, Region, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE), .groups = "drop")

# --------------------------- #
# Plot 1: RMSE% and MAE%
plot_rmse_mae <- summary_data %>%
  filter(Metric %in% c("RMSE_pct", "MAE_pct")) %>%
  ggplot(aes(x = Region, y = Mean, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~Dataset, ncol = 1) +
  labs(title = "RMSE% and MAE% by Region",
       y = "Percentage", x = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(savepath, "barplot_rmse_mae.png"), plot_rmse_mae, width = 12, height = 6)

# --------------------------- #
# Plot 2: R² and Correlation
plot_r_r2 <- summary_data %>%
  filter(Metric %in% c("R2", "Correlation")) %>%
  ggplot(aes(x = Region, y = Mean, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~Dataset, ncol = 1) +
  labs(title = "R² and Correlation by Region",
       y = "Value", x = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(savepath, "barplot_r_r2.png"), plot_r_r2, width = 12, height = 6)




















library(tidyverse)

savepath <- "/Users/max/Documents/Local/MS/NormativeModels/Regular/results/"

# Load internal stats
internal_stats <- read.csv(paste0(savepath, "internal_validation_stats.csv"))
external_stats <- read.csv(paste0(savepath, "external_regionwise_metrics.csv"))

# Compute percentage metrics for internal
internal_stats <- internal_stats %>%
  mutate(RMSE_pct = (RMSE / Mean_Original) * 100,
         MAE_pct = (MAE / Mean_Original) * 100)

# Prepare long format for internal
internal_long <- internal_stats %>%
  select(Region, RMSE_pct, MAE_pct, R2, Correlation) %>%
  pivot_longer(cols = c(RMSE_pct, MAE_pct, R2, Correlation),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "Internal")

# Prepare long format for external
external_long <- external_stats %>%
  select(Region, RMSE_pct, MAE_pct, R2, Correlation) %>%
  pivot_longer(cols = c(RMSE_pct, MAE_pct, R2, Correlation),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Dataset = "External")

# Combine
combined_long <- bind_rows(internal_long, external_long)

# Compute mean and SD for error bars
summary_data <- combined_long %>%
  group_by(Dataset, Region, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD = sd(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(SD = ifelse(is.na(SD) | SD == 0, 0.001, SD))  # Option 2 fix

# --------------------------- #
# Plot 1: RMSE% and MAE%
plot_rmse_mae <- summary_data %>%
  filter(Metric %in% c("RMSE_pct", "MAE_pct")) %>%
  ggplot(aes(x = Region, y = Mean, fill = Metric, group = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~Dataset, ncol = 1) +
  labs(title = "RMSE% and MAE% by Region",
       y = "Percentage", x = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(savepath, "barplot_rmse_mae.png"), plot_rmse_mae, width = 12, height = 6)

# --------------------------- #
# Plot 2: R² and Correlation
plot_r_r2 <- summary_data %>%
  filter(Metric %in% c("R2", "Correlation")) %>%
  ggplot(aes(x = Region, y = Mean, fill = Metric, group = Metric)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                position = position_dodge(width = 0.9), width = 0.3) +
  facet_wrap(~Dataset, ncol = 1) +
  labs(title = "R² and Correlation by Region",
       y = "Value", x = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(savepath, "barplot_r_r2.png"), plot_r_r2, width = 12, height = 6)
