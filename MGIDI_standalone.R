################################################################################
#                                                                              #
#          Multi-trait Genotype-Ideotype Distance Index (MGIDI)                #
#                     Standalone R Script for Plant Breeding                   #
#                                                                              #
#  This script demonstrates MGIDI analysis for wheat breeding data             #
#  Reference: https://www.youtube.com/watch?v=8XBNPCXS5hQ                      #
#                                                                              #
################################################################################

# Clear R environment
rm(list = ls(all = TRUE))
graphics.off()

#===============================================================================
# SECTION 1: LOAD REQUIRED PACKAGES
#===============================================================================

# Install required packages if not already installed
required_packages <- c("metan", "readxl", "dplyr", "tidyr", "ggplot2", "lme4")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
}

# Install missing packages
invisible(lapply(required_packages, install_if_missing))

# Load required libraries
library(metan)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)

#===============================================================================
# SECTION 2: DATA INPUT AND PREPROCESSING
#===============================================================================

# Function to create sample wheat breeding data for demonstration
create_sample_data <- function(n_genotypes = 50, n_reps = 3, n_blocks = 2) {
  set.seed(123)
  
  # Generate genotype names
  genotypes <- paste0("G", sprintf("%02d", 1:n_genotypes))
  
  # Generate data frame
  data <- expand.grid(
    GEN = genotypes,
    REP = 1:n_reps,
    BLOCK = 1:n_blocks
  )
  
  # Add trait values with some genetic and environmental variation
  n_obs <- nrow(data)
  
  # Genetic effects (genotype-specific)
  gen_effects <- data.frame(
    GEN = genotypes,
    gen_yield = rnorm(n_genotypes, 0, 2),
    gen_height = rnorm(n_genotypes, 0, 5),
    gen_protein = rnorm(n_genotypes, 0, 1),
    gen_tkw = rnorm(n_genotypes, 0, 3),
    gen_days = rnorm(n_genotypes, 0, 2)
  )
  
  # Merge genetic effects with data
  data <- merge(data, gen_effects, by = "GEN")
  
  # Add environmental noise and calculate final trait values
  data$YIELD <- 45 + data$gen_yield + rnorm(n_obs, 0, 1.5)
  data$PLANT_HEIGHT <- 95 + data$gen_height + rnorm(n_obs, 0, 3)
  data$PROTEIN_CONTENT <- 12.5 + data$gen_protein + rnorm(n_obs, 0, 0.5)
  data$TKW <- 42 + data$gen_tkw + rnorm(n_obs, 0, 2)  # Thousand kernel weight
  data$DAYS_TO_MATURITY <- 125 + data$gen_days + rnorm(n_obs, 0, 1.5)
  data$SPIKE_LENGTH <- 10 + rnorm(n_obs, 0, 1)
  data$SPIKELETS_PER_SPIKE <- 18 + rnorm(n_obs, 0, 2)
  data$GRAINS_PER_SPIKE <- 45 + rnorm(n_obs, 0, 5)
  
  # Clean up - remove genetic effect columns
  data <- data[, !(names(data) %in% c("gen_yield", "gen_height", "gen_protein", 
                                       "gen_tkw", "gen_days"))]
  
  # Ensure positive values
  data$YIELD <- pmax(data$YIELD, 20)
  data$PLANT_HEIGHT <- pmax(data$PLANT_HEIGHT, 60)
  data$PROTEIN_CONTENT <- pmax(data$PROTEIN_CONTENT, 8)
  data$TKW <- pmax(data$TKW, 25)
  data$DAYS_TO_MATURITY <- pmax(data$DAYS_TO_MATURITY, 100)
  data$SPIKE_LENGTH <- pmax(data$SPIKE_LENGTH, 5)
  data$SPIKELETS_PER_SPIKE <- pmax(data$SPIKELETS_PER_SPIKE, 10)
  data$GRAINS_PER_SPIKE <- pmax(data$GRAINS_PER_SPIKE, 20)
  
  return(data)
}

# Load or create data
# For real analysis, replace this with:
# wheat_data <- read_excel("path/to/your/data.xlsx", sheet = "Sheet1")

cat("\n========================================\n")
cat("MGIDI Analysis - Wheat Breeding Data\n")
cat("========================================\n\n")

# Create sample data for demonstration
wheat_data <- create_sample_data()

cat("Data dimensions:", nrow(wheat_data), "observations x", ncol(wheat_data), "variables\n")
cat("Number of genotypes:", length(unique(wheat_data$GEN)), "\n")
cat("Number of replications:", length(unique(wheat_data$REP)), "\n")
cat("Number of blocks:", length(unique(wheat_data$BLOCK)), "\n\n")

# View data structure
str(wheat_data)
head(wheat_data)

#===============================================================================
# SECTION 3: DATA VALIDATION AND CLEANING
#===============================================================================

# Function to check and clean data
validate_data <- function(data) {
  cat("\n--- Data Validation ---\n")
  
  # Check for missing values
  na_counts <- colSums(is.na(data))
  if (any(na_counts > 0)) {
    cat("Missing values found:\n")
    print(na_counts[na_counts > 0])
  } else {
    cat("No missing values detected.\n")
  }
  
  # Check factor levels
  cat("\nFactor levels:\n")
  cat("GEN:", length(unique(data$GEN)), "levels\n")
  cat("REP:", length(unique(data$REP)), "levels\n")
  cat("BLOCK:", length(unique(data$BLOCK)), "levels\n")
  
  return(data)
}

# Validate data
wheat_data <- validate_data(wheat_data)

#===============================================================================
# SECTION 4: CALCULATE SUMMARY STATISTICS
#===============================================================================

# Function to calculate summary statistics for traits
calculate_summary_stats <- function(data, traits) {
  cat("\n--- Summary Statistics ---\n\n")
  
  summary_stats <- data.frame(
    Trait = character(),
    Mean = numeric(),
    SD = numeric(),
    Min = numeric(),
    Max = numeric(),
    CV = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (trait in traits) {
    if (trait %in% names(data)) {
      values <- data[[trait]]
      trait_mean <- mean(values, na.rm = TRUE)
      trait_sd <- sd(values, na.rm = TRUE)
      trait_min <- min(values, na.rm = TRUE)
      trait_max <- max(values, na.rm = TRUE)
      trait_cv <- (trait_sd / trait_mean) * 100
      
      summary_stats <- rbind(summary_stats, data.frame(
        Trait = trait,
        Mean = round(trait_mean, 2),
        SD = round(trait_sd, 2),
        Min = round(trait_min, 2),
        Max = round(trait_max, 2),
        CV = round(trait_cv, 2),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  print(summary_stats)
  return(summary_stats)
}

# Define trait columns
trait_columns <- c("YIELD", "PLANT_HEIGHT", "PROTEIN_CONTENT", "TKW", 
                   "DAYS_TO_MATURITY", "SPIKE_LENGTH", "SPIKELETS_PER_SPIKE", 
                   "GRAINS_PER_SPIKE")

# Calculate summary statistics
summary_stats <- calculate_summary_stats(wheat_data, trait_columns)

#===============================================================================
# SECTION 5: CONVERT VARIABLES TO APPROPRIATE TYPES
#===============================================================================

cat("\n--- Converting Variables to Factors ---\n")

# Convert categorical variables to factors using correct syntax
# Note: Using as.factor() without space between 'as.' and 'factor'
wheat_data$GEN <- as.factor(wheat_data$GEN)
wheat_data$REP <- as.factor(wheat_data$REP)
wheat_data$BLOCK <- as.factor(wheat_data$BLOCK)

cat("GEN converted to factor:", is.factor(wheat_data$GEN), "\n")
cat("REP converted to factor:", is.factor(wheat_data$REP), "\n")
cat("BLOCK converted to factor:", is.factor(wheat_data$BLOCK), "\n")

#===============================================================================
# SECTION 6: FIT MIXED MODELS AND EXTRACT BLUPS
#===============================================================================

cat("\n--- Fitting Mixed Models ---\n")

# Function to fit mixed model for a single trait
fit_trait_model <- function(data, trait_name) {
  formula_str <- paste(trait_name, "~ (1|GEN) + (1|REP) + (1|BLOCK)")
  model <- lmer(as.formula(formula_str), data = data)
  return(model)
}

# Function to extract BLUPs from mixed model
extract_blups <- function(model, trait_name) {
  blups <- ranef(model)$GEN
  blups_df <- data.frame(
    GEN = rownames(blups),
    BLUP = blups[, 1]
  )
  names(blups_df)[2] <- paste0(trait_name, "_BLUP")
  return(blups_df)
}

# Fit models and extract BLUPs for all traits
blup_list <- list()

for (trait in trait_columns) {
  cat("Fitting model for:", trait, "\n")
  model <- fit_trait_model(wheat_data, trait)
  blup_df <- extract_blups(model, trait)
  
  if (length(blup_list) == 0) {
    blup_list <- blup_df
  } else {
    blup_list <- merge(blup_list, blup_df, by = "GEN")
  }
}

cat("\nBLUP extraction complete.\n")
cat("BLUP data dimensions:", nrow(blup_list), "x", ncol(blup_list), "\n")

# Calculate adjusted means (BLUPs + overall mean)
for (trait in trait_columns) {
  blup_col <- paste0(trait, "_BLUP")
  adj_col <- paste0(trait, "_ADJ")
  overall_mean <- mean(wheat_data[[trait]], na.rm = TRUE)
  blup_list[[adj_col]] <- blup_list[[blup_col]] + overall_mean
}

head(blup_list)

#===============================================================================
# SECTION 7: PERFORM MGIDI ANALYSIS USING METAN PACKAGE
#===============================================================================

cat("\n========================================\n")
cat("MGIDI ANALYSIS\n")
cat("========================================\n\n")

# For MGIDI analysis, we need data in proper format
# Using the metan package's gamem() function

# Prepare data for gamem
cat("Preparing data for MGIDI analysis...\n")

# Fit mixed models using metan package
# gamem fits a mixed-effect model for multi-environment trials
# For single environment, we use it with genotype as random effect

tryCatch({
  # Fit the model using gamem (genotype as random effect)
  cat("Fitting gamem model...\n")
  
  model_gamem <- gamem(
    .data = wheat_data,
    gen = GEN,
    rep = REP,
    resp = c(YIELD, PLANT_HEIGHT, PROTEIN_CONTENT, TKW, 
             DAYS_TO_MATURITY, SPIKE_LENGTH, SPIKELETS_PER_SPIKE, 
             GRAINS_PER_SPIKE)
  )
  
  cat("Model fitting complete.\n\n")
  
  # Extract BLUPs using metan
  blups_metan <- gmd(model_gamem, "blupg")
  cat("BLUPs extracted successfully.\n")
  print(head(blups_metan, 10))
  
}, error = function(e) {
  cat("Note: gamem requires balanced data. Using alternative approach.\n")
  cat("Error:", conditionMessage(e), "\n")
})

#===============================================================================
# SECTION 8: COMPUTE MGIDI INDEX
#===============================================================================

cat("\n--- Computing MGIDI Index ---\n")

# Function to calculate MGIDI manually
# This is useful when the metan package approach doesn't work

calculate_mgidi <- function(data, traits, ideotype_direction = NULL) {
  # data: data frame with genotypes in rows and traits in columns
  # traits: vector of trait column names
  # ideotype_direction: named vector with "h" for higher is better, "l" for lower is better
  
  if (is.null(ideotype_direction)) {
    # Default: all traits higher is better
    ideotype_direction <- setNames(rep("h", length(traits)), traits)
  }
  
  # Extract trait data
  trait_data <- data[, traits, drop = FALSE]
  genotypes <- data$GEN
  
  # Rescale traits (0-100 scale)
  rescaled_data <- as.data.frame(lapply(names(trait_data), function(trait) {
    x <- trait_data[[trait]]
    x_min <- min(x, na.rm = TRUE)
    x_max <- max(x, na.rm = TRUE)
    
    if (ideotype_direction[trait] == "h") {
      # Higher is better: rescale so max = 100, min = 0
      rescaled <- ((x - x_min) / (x_max - x_min)) * 100
    } else {
      # Lower is better: rescale so min = 100, max = 0
      rescaled <- ((x_max - x) / (x_max - x_min)) * 100
    }
    return(rescaled)
  }))
  names(rescaled_data) <- traits
  
  # Define ideotype (100 for all traits after rescaling)
  ideotype <- rep(100, length(traits))
  
  # Calculate Euclidean distance from ideotype
  distances <- apply(rescaled_data, 1, function(row) {
    sqrt(sum((row - ideotype)^2, na.rm = TRUE))
  })
  
  # Calculate MGIDI (lower distance = higher MGIDI rank)
  # Normalize to 0-100 scale
  max_dist <- max(distances, na.rm = TRUE)
  min_dist <- min(distances, na.rm = TRUE)
  mgidi_values <- ((max_dist - distances) / (max_dist - min_dist)) * 100
  
  # Create result data frame
  result <- data.frame(
    GEN = genotypes,
    MGIDI = round(mgidi_values, 4),
    Distance = round(distances, 4),
    Rank = rank(-mgidi_values, ties.method = "first")
  )
  
  # Sort by rank
  result <- result[order(result$Rank), ]
  rownames(result) <- NULL
  
  return(result)
}

# Define ideotype direction for wheat traits
# h = higher is better, l = lower is better
ideotype_direction <- c(
  "YIELD_ADJ" = "h",           # Higher yield is better
  "PLANT_HEIGHT_ADJ" = "l",    # Shorter plants often preferred (lodging resistance)
  "PROTEIN_CONTENT_ADJ" = "h", # Higher protein is better
  "TKW_ADJ" = "h",             # Higher thousand kernel weight is better
  "DAYS_TO_MATURITY_ADJ" = "l",# Earlier maturity often preferred
  "SPIKE_LENGTH_ADJ" = "h",    # Longer spikes may indicate more grain
  "SPIKELETS_PER_SPIKE_ADJ" = "h", # More spikelets is better
  "GRAINS_PER_SPIKE_ADJ" = "h"     # More grains is better
)

# Calculate MGIDI using adjusted means (BLUPs + overall mean)
adj_traits <- paste0(trait_columns, "_ADJ")
mgidi_result <- calculate_mgidi(blup_list, adj_traits, ideotype_direction)

cat("\n--- MGIDI Results (Top 20 Genotypes) ---\n")
print(head(mgidi_result, 20))

#===============================================================================
# SECTION 9: SELECTION OF SUPERIOR GENOTYPES
#===============================================================================

cat("\n========================================\n")
cat("GENOTYPE SELECTION\n")
cat("========================================\n\n")

# Function to select top genotypes based on MGIDI
select_genotypes <- function(mgidi_result, selection_intensity = 0.20) {
  n_total <- nrow(mgidi_result)
  n_select <- ceiling(n_total * selection_intensity)
  
  selected <- mgidi_result[1:n_select, ]
  
  cat("Selection intensity:", selection_intensity * 100, "%\n")
  cat("Total genotypes:", n_total, "\n")
  cat("Selected genotypes:", n_select, "\n\n")
  
  return(selected)
}

# Select top 20% of genotypes
selected_genotypes <- select_genotypes(mgidi_result, selection_intensity = 0.20)

cat("--- Selected Genotypes ---\n")
print(selected_genotypes)

# Calculate selection differential
cat("\n--- Selection Differential ---\n")

selection_differential <- data.frame(
  Trait = character(),
  Population_Mean = numeric(),
  Selected_Mean = numeric(),
  Selection_Differential = numeric(),
  Selection_Differential_Pct = numeric(),
  stringsAsFactors = FALSE
)

for (trait in adj_traits) {
  pop_mean <- mean(blup_list[[trait]], na.rm = TRUE)
  sel_mean <- mean(blup_list[blup_list$GEN %in% selected_genotypes$GEN, trait], na.rm = TRUE)
  sel_diff <- sel_mean - pop_mean
  sel_diff_pct <- (sel_diff / pop_mean) * 100
  
  selection_differential <- rbind(selection_differential, data.frame(
    Trait = trait,
    Population_Mean = round(pop_mean, 2),
    Selected_Mean = round(sel_mean, 2),
    Selection_Differential = round(sel_diff, 2),
    Selection_Differential_Pct = round(sel_diff_pct, 2),
    stringsAsFactors = FALSE
  ))
}

print(selection_differential)

#===============================================================================
# SECTION 10: CONTRIBUTION OF TRAITS TO MGIDI
#===============================================================================

cat("\n--- Trait Contribution to Selection ---\n")

# Calculate trait contributions using factor analysis approach
calculate_trait_contributions <- function(data, traits) {
  # Standardize data
  std_data <- scale(data[, traits])
  
  # Calculate correlation matrix
  cor_matrix <- cor(std_data, use = "complete.obs")
  
  # Perform PCA to understand trait contributions
  pca_result <- prcomp(std_data, scale. = FALSE)
  
  # Calculate proportion of variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Calculate loadings contribution
  loadings <- abs(pca_result$rotation[, 1:min(3, ncol(pca_result$rotation))])
  
  # Weight by variance explained
  contributions <- rowSums(loadings * var_explained[1:ncol(loadings)])
  contributions <- contributions / sum(contributions) * 100
  
  result <- data.frame(
    Trait = traits,
    Contribution = round(contributions, 2)
  )
  result <- result[order(-result$Contribution), ]
  rownames(result) <- NULL
  
  return(result)
}

trait_contributions <- calculate_trait_contributions(blup_list, adj_traits)
print(trait_contributions)

#===============================================================================
# SECTION 11: VISUALIZATION
#===============================================================================

cat("\n--- Creating Visualizations ---\n")

# Plot 1: MGIDI values bar plot
plot_mgidi_bar <- function(mgidi_result, top_n = 20) {
  plot_data <- mgidi_result[1:top_n, ]
  plot_data$GEN <- factor(plot_data$GEN, levels = rev(plot_data$GEN))
  
  p <- ggplot(plot_data, aes(x = GEN, y = MGIDI, fill = MGIDI)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    labs(
      title = "MGIDI Values - Top 20 Genotypes",
      subtitle = "Multi-trait Genotype-Ideotype Distance Index",
      x = "Genotype",
      y = "MGIDI Value"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    )
  
  return(p)
}

# Plot 2: Trait contribution pie chart
plot_trait_contribution <- function(contributions) {
  p <- ggplot(contributions, aes(x = "", y = Contribution, fill = Trait)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y", start = 0) +
    labs(
      title = "Trait Contribution to MGIDI",
      fill = "Trait"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  return(p)
}

# Create plots
p1 <- plot_mgidi_bar(mgidi_result)
p2 <- plot_trait_contribution(trait_contributions)

# Display plots
print(p1)
print(p2)

# Save plots
tryCatch({
  ggsave("MGIDI_barplot.png", p1, width = 10, height = 8, dpi = 300)
  ggsave("MGIDI_contribution.png", p2, width = 8, height = 8, dpi = 300)
  cat("Plots saved successfully.\n")
}, error = function(e) {
  cat("Could not save plots:", conditionMessage(e), "\n")
})

#===============================================================================
# SECTION 12: EXPORT RESULTS
#===============================================================================

cat("\n--- Exporting Results ---\n")

# Create results summary
create_results_summary <- function(mgidi_result, selected_genotypes, 
                                   selection_differential, trait_contributions) {
  
  results <- list(
    mgidi_full = mgidi_result,
    selected = selected_genotypes,
    selection_diff = selection_differential,
    contributions = trait_contributions
  )
  
  return(results)
}

# Format output for export
format_for_export <- function(value, format_string = "%12.2f") {
  sprintf(format_string, value)
}

# Create formatted output table
output_table <- data.frame(
  Rank = mgidi_result$Rank,
  Genotype = mgidi_result$GEN,
  MGIDI = format_for_export(mgidi_result$MGIDI),
  Distance = format_for_export(mgidi_result$Distance),
  Selected = ifelse(mgidi_result$GEN %in% selected_genotypes$GEN, "Yes", "No")
)

cat("\n--- Formatted Output (Top 25) ---\n")
print(head(output_table, 25))

# Save results to CSV
results_summary <- create_results_summary(
  mgidi_result, 
  selected_genotypes, 
  selection_differential, 
  trait_contributions
)

tryCatch({
  write.csv(mgidi_result, "MGIDI_results.csv", row.names = FALSE)
  write.csv(selected_genotypes, "MGIDI_selected_genotypes.csv", row.names = FALSE)
  write.csv(selection_differential, "MGIDI_selection_differential.csv", row.names = FALSE)
  write.csv(trait_contributions, "MGIDI_trait_contributions.csv", row.names = FALSE)
  
  cat("\nResults exported successfully to CSV files:\n")
  cat("  - MGIDI_results.csv\n")
  cat("  - MGIDI_selected_genotypes.csv\n")
  cat("  - MGIDI_selection_differential.csv\n")
  cat("  - MGIDI_trait_contributions.csv\n")
}, error = function(e) {
  cat("Could not export results:", conditionMessage(e), "\n")
})

#===============================================================================
# SECTION 13: FINAL SUMMARY
#===============================================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n\n")

cat("Summary of MGIDI Analysis:\n")
cat("--------------------------\n")
cat("Total genotypes analyzed:", nrow(mgidi_result), "\n")
cat("Number of traits:", length(adj_traits), "\n")
cat("Selection intensity: 20%\n")
cat("Genotypes selected:", nrow(selected_genotypes), "\n\n")

cat("Top 5 genotypes:\n")
for (i in 1:5) {
  cat(sprintf("  %d. %s (MGIDI: %.2f)\n", 
              i, mgidi_result$GEN[i], mgidi_result$MGIDI[i]))
}

cat("\nAnalysis completed successfully.\n")
cat("Please check the output CSV files for detailed results.\n")

#===============================================================================
# END OF SCRIPT
#===============================================================================
