# ==============================================================================
# Complement Hemolytic Assay (CH50 / AP50) Analysis Pipeline
# ==============================================================================
# Description: Processes raw plate reader data to calculate 50% lysis (HAEM_50),
#              performs log-linear regression, and generates publication plots.
# Author: Dr Benjamin P Flynn
# Date: 2026-03-24
# ==============================================================================

library(tidyverse)

# ------------------------------------------------------------------------------
# 1. Project Configuration & Directory Setup
# ------------------------------------------------------------------------------
DATA_DIR    <- "data"
OUTPUT_DIR  <- "output"
INPUT_FILE  <- file.path(DATA_DIR, "CH50_assay_data.csv")

# Create output folder automatically if missing
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Setup & Data Loading
# ------------------------------------------------------------------------------
if (!file.exists(INPUT_FILE)) {
  stop(paste("Error: Input file not found at", INPUT_FILE, "- Please place raw CSV in the data/ folder."))
}

raw_data <- read.csv(INPUT_FILE, check.names = FALSE)

# Pivot dilution columns into long format
data_long <- raw_data %>%
  pivot_longer(
    cols = matches("^[0-9]"), 
    names_to = "Dilution", 
    values_to = "OD"
  ) %>%
  mutate(
    Dilution = as.numeric(Dilution),
    PlasmaFraction = 1 / Dilution
  ) %>%
  filter(!is.na(OD))

# Dynamic Assay Detection (CH50 typically extends down to ~128x dilution)
is_ch50 <- min(data_long$PlasmaFraction, na.rm = TRUE) < 0.02
assay_name <- if (is_ch50) "CH50" else "AP50"

# Calculate baseline negative (0% lysis) and positive (100% lysis) control means
neg_rows <- data_long %>% filter(Patient == "NegCon")
pos_rows <- data_long %>% filter(Patient == "PosCon")

neg_val <- mean(neg_rows$OD, na.rm = TRUE)
pos_val <- mean(pos_rows$OD, na.rm = TRUE)

# Calculate relative hemolysis fraction for test samples
test_samples <- data_long %>%
  filter(!Patient %in% c("NegCon", "PosCon")) %>%
  mutate(
    Haemolysis = (OD - neg_val) / (pos_val - neg_val)
  )

# ------------------------------------------------------------------------------
# 3. Model Fitting & HAEM_50 Calculation
# ------------------------------------------------------------------------------
results <- test_samples %>%
  group_by(Patient, Sample, ID) %>%
  reframe({
    model   <- lm(Haemolysis ~ log(PlasmaFraction))
    b_int   <- as.numeric(coef(model)[1]) 
    a_slope <- as.numeric(coef(model)[2])
    r_sq    <- summary(model)$r.squared
    
    # Solve for x where Haemolysis (y) = 0.5
    HAEM_50_fraction <- exp((0.5 - b_int) / a_slope) * 100
    
    data.frame(
      R_Squared = r_sq,
      Slope     = a_slope,
      Intercept = b_int,
      HAEM_50   = HAEM_50_fraction
    )
  })

# Handle calculation edge cases safely
results <- results %>%
  mutate(HAEM_50 = case_when(
    is.infinite(HAEM_50) ~ 100,
    is.nan(HAEM_50)      ~ 100,
    Slope <= 0           ~ 100,
    TRUE                 ~ HAEM_50
  ))

# Calculate reference cohort thresholds dynamically
normal_ref_val <- results %>% 
  filter(grepl("Normal", Patient, ignore.case = TRUE) | grepl("Ref.*2", Patient, ignore.case = TRUE)) %>% 
  pull(HAEM_50) %>% mean(na.rm = TRUE)

low_ref_val <- results %>% 
  filter(grepl("Low", Patient, ignore.case = TRUE) | grepl("Ref.*1", Patient, ignore.case = TRUE)) %>% 
  pull(HAEM_50) %>% mean(na.rm = TRUE)

# Safe fallback thresholds if reference strings are missing
if (is.nan(normal_ref_val) || is.na(normal_ref_val)) normal_ref_val <- 25
if (is.nan(low_ref_val) || is.na(low_ref_val)) low_ref_val <- 55

# Assign clinical activity status
results <- results %>%
  mutate(Status = case_when(
    HAEM_50 <= normal_ref_val ~ "Normal",
    HAEM_50 <= low_ref_val    ~ "Low/Borderline",
    TRUE                      ~ "Deficient"
  ))

# Export analytical results table
write.csv(results, file = file.path(OUTPUT_DIR, "Rresults.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Cohort Comparison Plot Generation (Inverted Y-Axis)
# ------------------------------------------------------------------------------
# Determine upper limit for drawing stems on inverted scale
max_y <- max(results$HAEM_50, na.rm = TRUE) * 1.05

comparison_plot <- ggplot(results, aes(x = Patient, y = HAEM_50)) +
  geom_segment(aes(x = Patient, xend = Patient, y = max_y, yend = HAEM_50), 
               linetype = "dotted", color = "gray") +
  geom_point(aes(color = Status), size = 5) +
  geom_label(aes(label = Sample), vjust = 1.6, size = 3, fill = "white") +
  scale_color_manual(values = c("Normal" = "#27ae60", "Deficient" = "#e74c3c", "Low/Borderline" = "#f39c12")) +
  # Invert Y-axis so lower values (higher activity) appear at the top
  scale_y_reverse(expand = expansion(mult = c(0.1, 0.15))) +
  labs(
    title = paste(assay_name, "Complement Activity (50% Haemolysis)"),
    subtitle = "Calculated using Plasma Fraction Regression",
    caption = "*Inverted axis: Lower values indicate increased complement activity",
    y = "Plasma Fraction for 50% Haemolysis (%) ↑ Higher Activity",
    x = "Patient ID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "top",
    plot.caption = element_text(face = "italic", color = "grey30")
  )

print(comparison_plot)

# Export Summary Plots (TIFF & SVG)
ggsave(file.path(OUTPUT_DIR, "HAEM_50Results.tiff"), plot = comparison_plot, width = 6, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "HAEM_50Results.svg"), plot = comparison_plot, width = 6, height = 5)

# ------------------------------------------------------------------------------
# 5. Individual Adaptive Dose-Response Curves
# ------------------------------------------------------------------------------
if (is_ch50) {
  x_limits <- c(-5.5, -1.0)
  x_breaks <- seq(-5.5, -1.0, 0.5)   
} else {
  x_limits <- c(-3.0, 0.0)
  x_breaks <- seq(-3.0, 0.0, 0.5)    
}

for (i in seq_len(nrow(results))) {
  curr_patient <- results$Patient[i]
  curr_sample  <- results$Sample[i]
  curr_id      <- results$ID[i]
  curr_slope   <- results$Slope[i]
  curr_int     <- results$Intercept[i]
  curr_HAEM_50 <- results$HAEM_50[i]
  
  curve_data <- test_samples %>% 
    filter(Patient %in% curr_patient, 
           Sample %in% curr_sample,
           ID %in% curr_id) %>%
    arrange(PlasmaFraction)
  
  if (nrow(curve_data) == 0) {
    message(paste("Skipping iteration", i, "| Patient:", curr_patient, "- No matching data found."))
    next
  }
  
  line_x <- seq(min(curve_data$PlasmaFraction), max(curve_data$PlasmaFraction), length.out = 100)
  line_y <- (curr_slope * log(line_x)) + curr_int
  reg_line_df <- data.frame(x = line_x, y = line_y)
  
  HAEM_50_coord <- curr_HAEM_50 / 100
  
  p <- ggplot(curve_data, aes(x = log(PlasmaFraction), y = Haemolysis)) +
    geom_line(color = "black", alpha = 0.2) +
    geom_point(size = 3, color = "#2c3e50") +
    geom_line(data = reg_line_df, aes(x = log(x), y = y), color = "#3498db", linewidth = 1.2)
  
  if (is.finite(HAEM_50_coord) && HAEM_50_coord > 0) {
    p <- p + 
      geom_segment(x = x_limits[1], xend = log(HAEM_50_coord), y = 0.5, yend = 0.5, 
                   linetype = "dashed", color = "#e74c3c") +
      geom_segment(x = log(HAEM_50_coord), xend = log(HAEM_50_coord), y = 0, yend = 0.5, 
                   linetype = "dashed", color = "#e74c3c") +
      annotate("label", x = log(HAEM_50_coord), y = 0.1, 
               label = paste0(assay_name, ": ", round(curr_HAEM_50, 2), "%"), 
               color = "white", fill = "#e74c3c", fontface = "bold", size = 3)
  }
  
  p <- p +
    scale_x_continuous(breaks = x_breaks) + 
    scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
    coord_cartesian(xlim = x_limits, ylim = c(-0.1, 1.2), expand = FALSE) +
    labs(
      title = paste(assay_name, "Haemolysis Curve:", curr_patient),
      subtitle = paste("Sample:", curr_sample, "| R²:", round(results$R_Squared[i], 3)),
      x = "Ln(Plasma fraction)", 
      y = "Haemolysis Fraction (0.0 - 1.0)"
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
  
  save_patient <- if(is.na(curr_patient)) "NA" else curr_patient
  save_sample  <- if(is.na(curr_sample)) "NA" else curr_sample
  
  out_filename_base <- file.path(OUTPUT_DIR, paste0(assay_name, "Curve_", save_patient, "_", save_sample))
  
  ggsave(paste0(out_filename_base, ".tiff"), plot = p, width = 5, height = 5, dpi = 300)
  ggsave(paste0(out_filename_base, ".svg"), plot = p, width = 5, height = 5)
}