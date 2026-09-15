# ==============================================================================
# Script Name: protease_peptide_assay_analysis.R
# Description: Processes kinetic microplate reader data for protease peptide 
#              probe assays (DTNB/TNB absorbance at 405nm), performs blank 
#              and baseline (T0) adjustments, fits OLS kinetic rate models, 
#              calculates specific enzymatic activity (pmol/min/µg), and 
#              exports dual-format (TIFF/SVG) publication graphics.
# Author: Dr Benjamin P Flynn
# Date: 2026-08-20
# ==============================================================================

# 1. Libraries -----------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(geomtextpath)
library(stringr)
library(svglite)

# 2. Configuration & Directory Setup -------------------------------------------
DATA_DIR   <- "data"
OUTPUT_DIR <- "output"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Physical & Enzymatic Reaction Constants
WELL_VOL_L     <- 0.1 / 1000  # Well reaction volume (100 µL converted to Liters)
EXT_COEF_M_CM  <- 13260       # Extinction coefficient of TNB at 405nm (M^-1 cm^-1)
PATH_LENGTH_CM <- 0.32        # Microplate optical path length (cm)
ENZYME_MASS_UG <- 0.05        # Total enzyme mass loaded per well (µg)
NO_OF_REPS     <- 3           # Number of technical replicates

# 3. Helper Functions ----------------------------------------------------------
save_ggplot_dual <- function(plot, filename_base, width = 8, height = 6) {
  # Save High-Resolution TIFF (Raster 300 DPI)
  tiff(
    filename = file.path(OUTPUT_DIR, paste0(filename_base, ".tiff")),
    units = "in", width = width, height = height, res = 300
  )
  print(plot)
  dev.off()
  
  # Save Scalable Vector Graphics (SVG)
  ggsave(
    filename = file.path(OUTPUT_DIR, paste0(filename_base, ".svg")),
    plot = plot, width = width, height = height, device = "svg"
  )
}

# Helper to round up to specific numeric increments without legacy 'plyr'
round_up_increment <- function(x, increment = 0.05) {
  ceiling(x / increment) * increment
}

round_down_increment <- function(x, increment = 0.05) {
  floor(x / increment) * increment
}

# 4. Data Ingestion & Preprocessing --------------------------------------------
plate_map_file <- file.path(DATA_DIR, "PlateMap.csv")
data_od_file   <- file.path(DATA_DIR, "Data.csv")

if (!file.exists(plate_map_file) || !file.exists(data_od_file)) {
  stop("Error: Required input files ('PlateMap.csv' and/or 'Data.csv') missing from data/ directory.")
}

data_ref <- read.csv(plate_map_file, header = TRUE)
data_od  <- read.csv(data_od_file, header = TRUE)

# Merge optical density measurements with plate metadata
df <- merge(data_od, data_ref, by = "Wells") %>%
  filter(!is.na(Tx) & Tx != "")

# 5. Baseline & Blank Adjustments ----------------------------------------------
rep_cols <- paste0("N", seq_len(NO_OF_REPS))

# Extract time-matched BLANK values
blank_vals <- df %>%
  filter(Tx == "BLANK") %>%
  select(Time, all_of(rep_cols)) %>%
  rename_with(~ paste0("BLANK.", .), starts_with("N"))

# Perform replicate-wise BLANK subtraction across time points
df <- df %>%
  left_join(blank_vals, by = "Time")

for (i in seq_len(NO_OF_REPS)) {
  df[[paste0("Adjusted.N", i)]] <- df[[paste0("N", i)]] - df[[paste0("BLANK.N", i)]]
}

# Perform Time-Zero (T0) baseline subtraction per treatment
zero_t_vals <- df %>%
  filter(Time == 0) %>%
  select(Tx, all_of(rep_cols)) %>%
  rename_with(~ paste0("ZeroT.", .), starts_with("N"))

df <- df %>%
  left_join(zero_t_vals, by = "Tx")

for (i in seq_len(NO_OF_REPS)) {
  df[[paste0("ZeroT.Adj.N", i)]] <- df[[paste0("N", i)]] - df[[paste0("ZeroT.N", i)]]
}

# Extract and subtract BLANK from T0-adjusted values
blank_t0_vals <- df %>%
  filter(Tx == "BLANK") %>%
  select(Time, paste0("ZeroT.Adj.N", seq_len(NO_OF_REPS)))
colnames(blank_t0_vals) <- c("Time", paste0("ZeroT.B.Adj.N", seq_len(NO_OF_REPS)))

df <- df %>%
  left_join(blank_t0_vals, by = "Time")

for (i in seq_len(NO_OF_REPS)) {
  df[[paste0("ZeroT.Adjusted.N", i)]] <- df[[paste0("ZeroT.Adj.N", i)]] - df[[paste0("ZeroT.B.Adj.N", i)]]
}

# Exclude BLANK wells from final analytical dataset
df <- df %>% filter(Tx != "BLANK")

# 6. Reshape into Long-Format Dataframe for Modeling & Graphics ----------------
df_ggplot_list <- list()

for (r in seq_len(NO_OF_REPS)) {
  temp_df <- df %>%
    mutate(
      REP = r,
      B.Adjusted.OD = .data[[paste0("Adjusted.N", r)]],
      T0.B.Adjusted.OD = .data[[paste0("ZeroT.Adjusted.N", r)]]
    ) %>%
    select(Tx, Time, Group, REP, B.Adjusted.OD, T0.B.Adjusted.OD)
  
  df_ggplot_list[[r]] <- temp_df
}

df_ggplot <- bind_rows(df_ggplot_list)

# 7. Kinetic Rate Modeling & Specific Activity Calculation ---------------------
treatments <- unique(df_ggplot$Tx)
lm_results <- list()

for (tx_id in treatments) {
  sub_df <- df_ggplot %>% filter(Tx == tx_id)
  
  model <- lm(B.Adjusted.OD ~ Time, data = sub_df)
  slope <- as.numeric(coef(model)["Time"])
  group_id <- unique(sub_df$Group)[1]
  
  lm_results[[tx_id]] <- data.frame(
    Tx = tx_id,
    Slope = slope,
    Group = group_id
  )
}

df_lm <- bind_rows(lm_results)

# Calculate molar enzymatic activity (pmol/min/µg)
# Formula: (Slope * Volume_L * 10^12) / (ExtCoef * PathLength * Mass_ug)
df_lm <- df_lm %>%
  mutate(
    Activity_pmol_min_ug = ceiling(
      (Slope * WELL_VOL_L * 1e12) / (EXT_COEF_M_CM * PATH_LENGTH_CM * ENZYME_MASS_UG)
    )
  )

write.csv(df_lm, file = file.path(OUTPUT_DIR, "Calculated_Protease_Activity.csv"), row.names = FALSE)

# 8. Downstream Kinetic & Cohort Visualizations --------------------------------
adj_types <- c("B.Adjusted.OD", "T0.B.Adjusted.OD")

for (adj_col in adj_types) {
  max_od   <- round_up_increment(max(df_ggplot[[adj_col]], na.rm = TRUE), 0.05)
  min_od   <- round_down_increment(min(df_ggplot[[adj_col]], na.rm = TRUE), 0.05)
  minor_od <- 0.05
  
  # Non-linear LOESS Curve Plot with Inline Labels
  plot_curve <- ggplot(df_ggplot, aes(x = Time, y = .data[[adj_col]])) +
    geom_point(aes(color = Tx), shape = 4, alpha = 0.8) +
    geom_textsmooth(
      aes(group = Tx, color = Tx, label = Tx),
      method = "loess", se = FALSE, size = 3.5, fontface = "bold", text_only = FALSE
    ) +
    scale_color_viridis_d(option = "D") +
    coord_cartesian(ylim = c(min_od, max_od)) +
    scale_y_continuous(breaks = seq(min_od, max_od, by = minor_od)) +
    labs(
      title = "Protease Kinetic Progress Curves (LOESS Fit)",
      x = "Time (min)",
      y = "Optical Density (405nm)",
      color = "Sample / Treatment"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  save_ggplot_dual(plot_curve, paste0(adj_col, "_CurveOfBestFit"))
  
  # OLS Linear Regression Kinetic Rate Plot
  plot_line <- ggplot(df_ggplot, aes(x = Time, y = .data[[adj_col]])) +
    geom_point(aes(color = Tx), shape = 4, alpha = 0.8) +
    geom_smooth(aes(group = Tx, color = Tx), method = "lm", se = FALSE) +
    scale_color_viridis_d(option = "D") +
    coord_cartesian(ylim = c(min_od, max_od)) +
    scale_y_continuous(breaks = seq(min_od, max_od, by = minor_od)) +
    labs(
      title = "Protease Reaction Velocities (OLS Linear Fit)",
      x = "Time (min)",
      y = "Optical Density (405nm)",
      color = "Sample / Treatment"
    ) +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
  
  save_ggplot_dual(plot_line, paste0(adj_col, "_LineOfBestFit"))
}

# Cohort Specific Activity Comparison Boxplot
plot_box <- ggplot(df_lm, aes(x = Group, y = Activity_pmol_min_ug, color = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_point(size = 2) +
  geom_text_repel(
    aes(label = Tx), size = 3, point.padding = 0.2, box.padding = 0.5,
    force = 2, segment.color = "grey50", show.legend = FALSE
  ) +
  scale_color_viridis_d(option = "D") +
  labs(
    title = "Protease Specific Activity Across Patient Cohorts",
    x = "Patient Cohort / Group",
    y = expression("Specific Activity (" * pmol/min/mu*g * ")"),
    color = "Cohort Group"
  ) +
  theme_light() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black")
  )

save_ggplot_dual(plot_box, "Protease_Specific_Activity_BoxPlot")