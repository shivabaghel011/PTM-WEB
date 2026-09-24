# ==============================================================================
# Automated Script to Generate a Single Normalized Alluvial Plot of PTMs
# ==============================================================================
# Input: alluvial_data.tsv, Master_Writer_Eraser_25Jun2026.tsv
# Output: alluvial_plots/all_modifications_alluvial.pdf (and .png)
# Preferred environment: RStudio
# ==============================================================================

# 1. Manage Package Dependencies
message("Checking and installing dependencies...")
cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "ggalluvial", "scales")
new_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])]
if (length(new_cran) > 0) {
  message("Installing missing CRAN packages: ", paste(new_cran, collapse = ", "))
  install.packages(new_cran, repos = "https://cloud.r-project.org")
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)
library(scales)

# Read the data
input_file <- "alluvial_data.tsv"
master_file <- "Master_Writer_Eraser_25Jun2026.tsv"

if (!file.exists(input_file)) {
  stop(paste("Input file not found:", input_file))
}
if (!file.exists(master_file)) {
  stop(paste("Master file not found:", master_file))
}

message("Reading data from ", input_file)
df <- read.csv(input_file, sep = "\t", stringsAsFactors = FALSE)

# ------------------------------------------------------------------------------
# 2. Extract Enzyme Counts from Master File for Species
# ------------------------------------------------------------------------------
message("Reading master file for enzyme count distribution...")
df_master <- read.delim(master_file, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Filter for target species (same as Perl script)
target_taxons <- c("9606", "10090", "10116", "3702", "559292", "285006", "83333", "83334", "331111", "199310", "585035", "574521", "585055")

df_master <- df_master %>%
  mutate(
    Taxon = case_when(
      grepl("TaxID:(\\d+)", `Organism ID`) ~ regmatches(`Organism ID`, regexpr("TaxID:(\\d+)", `Organism ID`)) %>% gsub("TaxID:", "", .),
      grepl("(\\d+)", `Organism ID`) ~ regmatches(`Organism ID`, regexpr("(\\d+)", `Organism ID`)),
      TRUE ~ ""
    )
  ) %>%
  filter(Taxon %in% target_taxons)

target_mods_list <- c("phosphorylation", "methylation", "acetylation", "ubiquitination", "glycosylation", "glycation")

df_enz_counts <- df_master %>%
  filter(`Gene Name` != "" & `Annotation` != "" & `Modification` != "") %>%
  mutate(`Gene Name` = toupper(trimws(`Gene Name`))) %>%
  rowwise() %>%
  mutate(
    Matched_PTMs = list(target_mods_list[sapply(target_mods_list, function(m) grepl(m, Modification, ignore.case = TRUE))])
  ) %>%
  unnest(Matched_PTMs) %>%
  mutate(
    Standard_PTM = case_when(
      Matched_PTMs == "phosphorylation" ~ "Phosphorylation",
      Matched_PTMs == "methylation"     ~ "Methylation",
      Matched_PTMs == "acetylation"     ~ "Acetylation",
      Matched_PTMs == "ubiquitination"  ~ "Ubiquitination",
      Matched_PTMs == "glycosylation"   ~ "Glycosylation",
      Matched_PTMs == "glycation"       ~ "Glycation"
    )
  )

writer_counts <- df_enz_counts %>%
  filter(grepl("writer", Annotation, ignore.case = TRUE)) %>%
  group_by(Standard_PTM) %>%
  summarise(W_Count = n_distinct(`Gene Name`), .groups = "drop")

eraser_counts <- df_enz_counts %>%
  filter(grepl("eraser", Annotation, ignore.case = TRUE)) %>%
  group_by(Standard_PTM) %>%
  summarise(E_Count = n_distinct(`Gene Name`), .groups = "drop")

# Ensure all 6 modifications are represented
all_mods_df <- data.frame(Standard_PTM = c("Phosphorylation", "Methylation", "Acetylation", "Ubiquitination", "Glycosylation", "Glycation"), stringsAsFactors=FALSE)
writer_counts <- all_mods_df %>% left_join(writer_counts, by="Standard_PTM") %>% replace_na(list(W_Count = 0))
eraser_counts <- all_mods_df %>% left_join(eraser_counts, by="Standard_PTM") %>% replace_na(list(E_Count = 0))

# ------------------------------------------------------------------------------
# 3. Data Processing & Axis Weight Calculations
# ------------------------------------------------------------------------------
# Unique substrates per modification
df_unique <- df %>%
  group_by(Modification, Interactor) %>%
  summarise(
    Type = first(Type),
    .groups = "drop"
  )

df_mod_totals <- df_unique %>%
  group_by(Modification) %>%
  summarise(Mod_Total = n(), .groups = "drop")

df_counts <- df_unique %>%
  group_by(Modification, Type) %>%
  summarise(Count = n(), .groups = "drop") %>%
  left_join(df_mod_totals, by = "Modification") %>%
  mutate(
    # Percentage within the PTM (for labels)
    Pct = (Count / Mod_Total) * 100,
    Type_Display = case_when(
      Type == "Writer_Only" ~ "Writer-Only",
      Type == "Eraser_Only" ~ "Eraser-Only",
      Type == "Shared"      ~ "Shared"
    ),
    # Only write the percentage, hide if less than 2.5% to prevent overlap
    Pct_Label = ifelse(Pct >= 2.5, sprintf("%.1f%%", Pct), "")
  )

# Add Writer and Eraser counts
df_counts <- df_counts %>%
  left_join(writer_counts, by = c("Modification" = "Standard_PTM")) %>%
  left_join(eraser_counts, by = c("Modification" = "Standard_PTM"))

# Grand Totals
total_writers <- sum(writer_counts$W_Count)
total_erasers <- sum(eraser_counts$E_Count)
total_interactors <- sum(df_counts$Count)

# Substrate count for normalization of within-PTM fractions on W/E sides
df_sums <- df_unique %>%
  group_by(Modification) %>%
  summarise(
    W_Shared_Substrates = sum(Type %in% c("Writer_Only", "Shared")),
    E_Shared_Substrates = sum(Type %in% c("Eraser_Only", "Shared")),
    .groups = "drop"
  )

df_counts <- df_counts %>%
  left_join(df_sums, by = "Modification") %>%
  mutate(
    # Writer side share
    Pct_W = (W_Count / total_writers) * 100,
    # Eraser side share
    Pct_E = (E_Count / total_erasers) * 100,
    # Interactors side share
    Pct_I = (Count / total_interactors) * 100
  )

# Print percentages and counts for the user
message("--- START STATISTICS TABLE ---")
cat(sprintf("%-18s | %-12s | %-12s | %-12s | %-12s\n", "Modification", "Total Substr", "Writer-Only", "Eraser-Only", "Shared (Pair)"))
for (m in unique(df_counts$Modification)) {
  sub_df <- df_counts %>% filter(Modification == m)
  tot <- max(sub_df$Mod_Total)
  w_only <- sub_df %>% filter(Type == "Writer_Only") %>% pull(Pct) %>% {ifelse(length(.), ., 0)}
  e_only <- sub_df %>% filter(Type == "Eraser_Only") %>% pull(Pct) %>% {ifelse(length(.), ., 0)}
  sh <- sub_df %>% filter(Type == "Shared") %>% pull(Pct) %>% {ifelse(length(.), ., 0)}
  
  w_only_c <- sub_df %>% filter(Type == "Writer_Only") %>% pull(Count) %>% {ifelse(length(.), ., 0)}
  e_only_c <- sub_df %>% filter(Type == "Eraser_Only") %>% pull(Count) %>% {ifelse(length(.), ., 0)}
  sh_c <- sub_df %>% filter(Type == "Shared") %>% pull(Count) %>% {ifelse(length(.), ., 0)}
  
  cat(sprintf("%-18s | %-12d | %5.1f%% (%-4d) | %5.1f%% (%-4d) | %5.1f%% (%-4d)\n", m, tot, w_only, w_only_c, e_only, e_only_c, sh, sh_c))
}
message("--- END STATISTICS TABLE ---")

# Calculate variable weights for the three axes (heights normalized to 100%)
df_counts <- df_counts %>%
  mutate(
    # Left axis weight (Writers)
    Weight_W = ifelse(Type == "Eraser_Only", 0, Pct_W * ifelse(W_Shared_Substrates == 0, 0, Count / W_Shared_Substrates)),
    # Middle axis weight (Interactors)
    Weight_I = Pct_I,
    # Right axis weight (Erasers)
    Weight_E = ifelse(Type == "Writer_Only", 0, Pct_E * ifelse(E_Shared_Substrates == 0, 0, Count / E_Shared_Substrates))
  )

# Determine vertical stacking order of modifications based on the shared percentage
df_shared_pcts <- df_counts %>%
  filter(Type == "Shared") %>%
  select(Modification, Pct) %>%
  complete(Modification = unique(df_counts$Modification), fill = list(Pct = 0)) %>%
  arrange(Pct) # Ascending order (lowest share first, highest at the top)

ordered_mods <- df_shared_pcts$Modification

# Build the wide format dataframe for ggalluvial
df_wide <- df_counts %>%
  mutate(
    Writer_Strata = ifelse(Type %in% c("Shared", "Writer_Only"), paste0(Modification, " (W)"), paste0("Invisible (W) for ", Modification)),
    Middle_Strata = paste0(Modification, " (", Type_Display, ")"),
    Eraser_Strata = ifelse(Type %in% c("Shared", "Eraser_Only"), paste0(Modification, " (E)"), paste0("Invisible (E) for ", Modification))
  )

# Convert to long format (Lodes format) for ggalluvial
df_long <- to_lodes_form(
  df_wide,
  key = "x",
  value = "stratum",
  id = "alluvium",
  axes = c("Writer_Strata", "Middle_Strata", "Eraser_Strata")
)

# Map columns back from the wide dataframe
df_long$Modification  <- df_wide$Modification[df_long$alluvium]
df_long$Type          <- df_wide$Type[df_long$alluvium]
df_long$Pct           <- df_wide$Pct[df_long$alluvium]
df_long$Pct_Label     <- df_wide$Pct_Label[df_long$alluvium]
df_long$Middle_Strata <- df_wide$Middle_Strata[df_long$alluvium]

# Assign the correct weight dynamically based on the axis
df_long$Weight <- case_when(
  df_long$x == "Writer_Strata" ~ df_wide$Weight_W[df_long$alluvium],
  df_long$x == "Middle_Strata" ~ df_wide$Weight_I[df_long$alluvium],
  df_long$x == "Eraser_Strata" ~ df_wide$Weight_E[df_long$alluvium]
)

# Map dynamic legend categories
# Legend label format: Phosphorylation_OnlyWriter, Phosphorylation_OnlyEraser, Phosphorylation_EnzymePair
df_long <- df_long %>%
  mutate(
    Legend_Category = paste0(
      Modification, "_",
      case_when(
        Type == "Writer_Only" ~ "OnlyWriter",
        Type == "Eraser_Only" ~ "OnlyEraser",
        Type == "Shared"      ~ "EnzymePair"
      )
    )
  )

# Define connection fills for left-to-middle and middle-to-right
df_long <- df_long %>%
  mutate(
    Fill_Left = ifelse(Type == "Eraser_Only", "Invisible", Legend_Category),
    Fill_Right = ifelse(Type == "Writer_Only", "Invisible", Legend_Category)
  )

# Define custom stratum factor levels to control the vertical stacking order without spacers
ordered_stratum_levels <- c()
for (m in ordered_mods) {
  ordered_stratum_levels <- c(
    ordered_stratum_levels,
    paste0("Invisible (W) for ", m),
    paste0(m, " (W)"),
    paste0("Invisible (E) for ", m),
    paste0(m, " (E)"),
    paste0(m, " (Writer-Only)"),
    paste0(m, " (Shared)"),
    paste0(m, " (Eraser-Only)")
  )
}
ordered_stratum_levels <- c(ordered_stratum_levels, "Invisible")

df_long$stratum <- factor(df_long$stratum, levels = ordered_stratum_levels)

# Construct legend categories in order of PTM stacking (highest sharing last/top)
legend_levels <- c()
for (m in ordered_mods) {
  legend_levels <- c(
    legend_levels,
    paste0(m, "_OnlyWriter"),
    paste0(m, "_OnlyEraser"),
    paste0(m, "_EnzymePair")
  )
}

# ------------------------------------------------------------------------------
# 4. Define the Color System (18 legend-mapped colors + transparent)
# ------------------------------------------------------------------------------
base_writer_colors <- c(
  "Phosphorylation" = "#1A6FC8",  # Blue
  "Acetylation"     = "#D96B07",  # Orange
  "Methylation"     = "#289128",  # Green
  "Glycosylation"   = "#C81D52",  # Red/Pink
  "Glycation"       = "#78249C",  # Purple
  "Ubiquitination"  = "#169C9C"   # Teal
)

base_eraser_colors <- c(
  "Phosphorylation" = "#85BDF2",  # Light Blue
  "Acetylation"     = "#F7BF8A",  # Light Orange
  "Methylation"     = "#87DE87",  # Light Green
  "Glycosylation"   = "#F28EAF",  # Light Pink
  "Glycation"       = "#D09BE3",  # Light Purple
  "Ubiquitination"  = "#7CEBEB"   # Light Teal
)

base_shared_colors <- c(
  "Phosphorylation" = "#4F95DD",  # Blend
  "Acetylation"     = "#E89548",  # Blend
  "Methylation"     = "#57B757",  # Blend
  "Glycosylation"   = "#DE5580",  # Blend
  "Glycation"       = "#A45FC0",  # Blend
  "Ubiquitination"  = "#3EC4C4"   # Blend
)

stratum_colors <- c("Invisible" = "#FFFFFF00")
stratum_border_colors <- c("Invisible" = "#FFFFFF00")

# Map colors to the 18 legend categories and stratum blocks
for (m in unique(df_counts$Modification)) {
  # Left Strata
  stratum_colors[paste0(m, " (W)")] <- base_writer_colors[m]
  stratum_border_colors[paste0(m, " (W)")] <- "#616161"
  
  stratum_colors[paste0("Invisible (W) for ", m)] <- "#FFFFFF00"
  stratum_border_colors[paste0("Invisible (W) for ", m)] <- "#FFFFFF00"
  
  # Right Strata
  stratum_colors[paste0(m, " (E)")] <- base_eraser_colors[m]
  stratum_border_colors[paste0(m, " (E)")] <- "#616161"
  
  stratum_colors[paste0("Invisible (E) for ", m)] <- "#FFFFFF00"
  stratum_border_colors[paste0("Invisible (E) for ", m)] <- "#FFFFFF00"
  
  # Middle Strata (Block categories)
  stratum_colors[paste0(m, " (Writer-Only)")] <- base_writer_colors[m]
  stratum_border_colors[paste0(m, " (Writer-Only)")] <- "#616161"
  
  stratum_colors[paste0(m, " (Shared)")] <- base_shared_colors[m]
  stratum_border_colors[paste0(m, " (Shared)")] <- "#616161"
  
  stratum_colors[paste0(m, " (Eraser-Only)")] <- base_eraser_colors[m]
  stratum_border_colors[paste0(m, " (Eraser-Only)")] <- "#616161"
  
  # Legend Categories
  stratum_colors[paste0(m, "_OnlyWriter")] <- base_writer_colors[m]
  stratum_colors[paste0(m, "_OnlyEraser")] <- base_eraser_colors[m]
  stratum_colors[paste0(m, "_EnzymePair")] <- base_shared_colors[m]
}

# ------------------------------------------------------------------------------
# 5. Generate the Alluvial Plot
# ------------------------------------------------------------------------------
message("Generating the single alluvial plot...")

# Reorder factor levels for the x-axis (from left to right)
df_long$x <- factor(df_long$x, levels = c("Writer_Strata", "Middle_Strata", "Eraser_Strata"))

p <- ggplot(df_long, aes(y = Weight, x = x, stratum = stratum, alluvium = alluvium)) +
  # First connection flow (Left -> Middle)
  geom_flow(
    data = filter(df_long, x %in% c("Writer_Strata", "Middle_Strata")),
    aes(fill = Fill_Left),
    width = 1/6,
    alpha = 0.85,
    color = "white",
    linewidth = 0.1,
    show.legend = TRUE # Let geom_flow handle the legend display
  ) +
  # Second connection flow (Middle -> Right)
  geom_flow(
    data = filter(df_long, x %in% c("Middle_Strata", "Eraser_Strata")),
    aes(fill = Fill_Right),
    width = 1/6,
    alpha = 0.85,
    color = "white",
    linewidth = 0.1,
    show.legend = FALSE # Duplicate, hide legend for second layer
  ) +
  # Strata blocks colored by their specific stratum names and mapped borders
  geom_stratum(aes(fill = stratum, color = stratum), width = 1/6, linewidth = 0.2, show.legend = FALSE) +
  # Color scale mappings
  scale_fill_manual(values = stratum_colors, name = "PTM Modification Categories", breaks = legend_levels) +
  scale_color_manual(values = stratum_border_colors) +
  # Set x-axis labels
  scale_x_discrete(
    limits = c("Writer_Strata", "Middle_Strata", "Eraser_Strata"),
    labels = c(
      "Writer_Strata" = "Writers\n(Left Group)",
      "Middle_Strata" = "Substrates\n(Middle)",
      "Eraser_Strata" = "Erasers\n(Right Group)"
    ),
    expand = c(.05, .05)
  ) +
  # Labels and themes
  labs(
    title = "PTM Substrate Sharing and Alluvial Networks",
    subtitle = "Varying weights: Left (Writers share), Right (Erasers share), Middle (Interactor share)",
    y = "Abundance & Distribution Percentage",
    x = ""
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10.5, face = "italic", hjust = 0.5),
    legend.position = "right", # Display the legend on the right side
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8.5),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(15, 15, 15, 15)
  )

# Create output directory and save the plots
dir.create("alluvial_plots", showWarnings = FALSE)
pdf_out <- "alluvial_plots/all_modifications_alluvial.pdf"
png_out <- "alluvial_plots/all_modifications_alluvial.png"

message("Saving PDF plot...")
ggsave(pdf_out, plot = p, width = 13, height = 13, device = "pdf")

message("Saving PNG plot...")
ggsave(png_out, plot = p, width = 13, height = 13, dpi = 300, device = "png")

message("Execution completed successfully!")
