# ==============================================================================
# Script Name: elisa_4pl_standard_curve_analysis.R
# Description: Processes raw ELISA plate reader data, performs blank adjustments,
#              fits a 4-parameter logistic (4PL) standard curve, interpolates 
#              sample concentrations, integrates multi-variable 3D organoid 
#              metadata, and exports publication-ready raster (TIFF) & vector
#              (SVG) figures.
# Author: Dr Benjamin P Flynn
# Date: 2026-05-13
# ==============================================================================

# 1. Libraries -----------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(drc)
library(forcats)
library(tidyr)
library(stringr)
library(ggrepel)
library(svglite)

# 2. Configuration & Directory Setup -------------------------------------------
# Create output directory if it doesn't exist
if (!dir.exists("output/")) {
  dir.create("output/", recursive = TRUE)
}

NO.OF.REPS <- 2
ELISA.Target  <- "Enzyme A"
DATE       <- "2026-05-13"

# 3. Helper Functions ----------------------------------------------------------
save_ggplot_tiff <- function(plot, filename, width = 8, height = 5) {
  # Save TIFF (Raster)
  tiff(
    filename = file.path("output/", paste0(filename, ".tiff")),
    units = "in", width = width, height = height, res = 300
  )
  print(plot)
  dev.off()
  
  # Save SVG (Vector)
  ggsave(
    filename = file.path("output/", paste0(filename, ".svg")),
    plot = plot, width = width, height = height,
    device = "svg"
  )
}

# 4. Load Data -----------------------------------------------------------------
df <- read.csv(file.path("data/", "Sample_ELISA_Results.csv"), header = TRUE)
REPS <- paste0("N", seq(1, NO.OF.REPS))

# 5. Data Preprocessing & Blank Adjustment -------------------------------------
# Blank adjustment applied BEFORE averaging replicates
BLANK <- as.numeric(rowMeans(
  df[df$Sample == "STD" & df$Conc == 0, c(REPS)], 
  na.rm = TRUE
))

for (i in REPS) {
  df[[i]] <- df[[i]] - BLANK
}

# Merge patient metadata
metadata <- read.csv(file.path("data/", "Sample_3D_Organoid_Phenotypes.csv"), header = TRUE)
df <- merge(df, metadata, by = "Sample", all.x = TRUE)

# Calculate average blank-adjusted OD
df$Average_OD <- rowMeans(df[, c(REPS)], na.rm = TRUE)

# Build replicate-level DataFrame for ggplot
COLNAMES <- c("Sample", "Dilution.Factor", "OD", "Spheroid.Profile", "Disease.Group")
df.ggplot <- data.frame(matrix(nrow = 0, ncol = length(COLNAMES)))
colnames(df.ggplot) <- COLNAMES

for (i in REPS) {
  df.inter <- df[, c("Sample", "Dilution.Factor", i, "Spheroid.Profile", "Disease.Group")]
  colnames(df.inter) <- COLNAMES
  df.ggplot <- rbind(df.ggplot, df.inter)
}
rm(df.inter)

df.ggplot <- df.ggplot[df.ggplot$Sample != "STD", ]
df.ggplot[is.na(df.ggplot[["Spheroid.Profile"]]), "Spheroid.Profile"] <- "Unknown"

# 6. Quality Control: OD Distribution Plot -------------------------------------
STD_OD <- df[df$Sample == "STD", "Average_OD"]

PLOT.OD <- df.ggplot %>%
  ggplot(aes(x = Sample, y = OD, color = Spheroid.Profile)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  geom_hline(yintercept = STD_OD, linetype = "dashed", color = "red") +
  theme_light() +
  scale_color_viridis_d(option = "H") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank()
  ) +
  labs(
    title = paste0("ELISA: ", ELISA.Target),
    x = "Sample",
    y = "Blank-adjusted OD450"
  )

print(PLOT.OD)
save_ggplot_tiff(PLOT.OD, "PLOT_OD_Values")

# 7. Standard Curve (4PL Model Fit) --------------------------------------------
STD <- df %>%
  filter(Sample == "STD") %>%
  dplyr::select(conc = Conc, OD = Average_OD) %>%
  transmute(
    conc = as.numeric(as.character(conc)),
    OD   = as.numeric(OD)
  ) 

fit <- drm(
  OD ~ conc, 
  data = STD,
  fct = LL.4(names = c("Slope", "Lower", "Upper", "EC50"))
)

# Save standard curve plot
tiff(
  filename = file.path("output/", "StandardCurve.tiff"),
  units = "in", width = 5, height = 5, res = 300
)
plot(
  fit, type = "all",
  xlab = "concentration (ng/mL)", 
  ylab = "OD450",
  main = paste0("ELISA: ", ELISA.Target, "\n", "4PL Standard Curve")
)
dev.off()

# 8. Estimate Concentrations ---------------------------------------------------
# Mean-level estimation
df$Est.Conc <- ED(fit, respLev = df$Average_OD, type = "absolute")[, 1] * df$Dilution.Factor
df$Est.Conc[is.na(df$Est.Conc) | df$Est.Conc < 0] <- 0

# Replicate-level estimation
df.ggplot$Est.Conc <- ED(fit, respLev = df.ggplot$OD, type = "absolute")[, 1] * df.ggplot$Dilution.Factor
df.ggplot$Est.Conc[is.na(df.ggplot$Est.Conc) | df.ggplot$Est.Conc < 0] <- 0

# Save analyzed dataframe
write.csv(df, file = file.path("output/", "RAnalysedResults.csv"), row.names = FALSE)

# 9. Downstream Visualization --------------------------------------------------
# Estimated concentrations plot
PLOT.EST.CONC <- df.ggplot %>%
  ggplot(aes(x = Sample, y = Est.Conc, color = Spheroid.Profile)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_light() +
  scale_color_viridis_d(option = "H") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank()
  ) +
  labs(
    title = paste0("ELISA: ", ELISA.Target),
    y = paste0("Estimated ", ELISA.Target, " concentration (ng/mL)"),
    x = "Sample"
  )

print(PLOT.EST.CONC)
save_ggplot_tiff(PLOT.EST.CONC, "PLOT_EstConc")

# Log2 concentration plot
PLOT.EST.CONC.LOG <- df.ggplot %>%
  ggplot(aes(x = Sample, y = log2(Est.Conc + 1), color = Spheroid.Profile)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point() +
  theme_light() +
  scale_color_viridis_d(option = "H") +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none",
    axis.line = element_line(color = "black"), 
    axis.ticks = element_line(color = "black"),
    strip.background = element_blank()
  ) +
  labs(
    title = paste0("ELISA: ", ELISA.Target),
    y = paste0("log2(Estimated ", ELISA.Target, " concentration + 1)"),
    x = "Sample"
  )

print(PLOT.EST.CONC.LOG)
save_ggplot_tiff(PLOT.EST.CONC.LOG, "PLOTConcLog2")