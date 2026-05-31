# -------------------------------------------
# CLEAN SCRIPT: Process sharma2 Data Only
#
#
# -------------------------------------------
# This script does NOT load old .Rdata files.
# It builds the Phyloseq object from your ASV table and Mapping CSV.



# 1. Setup & Clean Environment
library('dplyr')
library('vegan')
library('phyloseq')
library('dada2')
library('ggplot2')

# Clear any variables that might exist from old sessions
# This ensures we don't accidentally use old variables
rm(list = ls()) 

# 2. Define Inputs
# Adjust these paths and names to match your actual files
asv_file <- 'sharma_dump/asv2.csv'   # Adjust to your actual CSV name
map_file <- "sharma_dump/sharmap2.csv"                   # Your mapping file with Group info
output_file <- "sharma_clean_physeq.Rdata"

# 3. Read Data
# Check if file exists before reading
if (!file.exists(asv_file)) {
  stop("ASV file not found! Please ensure your processed ASV table is in the same directory.")
}

# Read Mapping (Group info)
map_data <- read.csv2(map_file, stringsAsFactors = FALSE)

# Read ASV Table
# Ensure rownames are sample names (if not already in 'DNA_id' column)
asv_raw <- read.csv(asv_file, row.names = 1, check.names = FALSE)

# 4. Prepare Metadata (Match Sample Names)
# The first row of the mapping file is usually metadata, ensure it matches ASV rownames
sample_names <- asv_raw
metadata    <- map_data %>%
  # Rename to match ASV sample names if needed
  mutate(Sample_Name = row.names(asv_raw)) %>% 
  select(Sample_Name, Group, FASTQ.Run.Accession)       %>% 
  mutate(Sample_Name = row.names(asv_raw))

# Check for missing samples
missing_samples <- setdiff(row.names(asv_raw), unique(metadata$Sample_Name))
if (length(missing_samples) > 0) {
  warning("Missing samples in metadata: ", paste(missing_samples, collapse=", "))
}

# Ensure all columns have names (remove unnamed samples if any)
colnames(asv_raw)[is.na(colnames(asv_raw))] <- "unknown"

# 5. Construct Phyloseq Object
# Note: Ensure your ASV table has a column named 'Family' (or similar) for taxonomy
# If Family is a column in the ASV table (e.g. 'family'), use that.
# If not, you may need to run assignTaxonomy() beforehand.
# Assuming your ASV table has taxonomic columns:

sharmap2 <- read.csv2("~/sharma_dump/sharmap2.csv",
                      stringsAsFactors = FALSE)

# Use the run IDs directly as the sample IDs
sharmap2$DNA_id <- sharmap2$FASTQ.Run.Accession
rownames(sharmap2) <- sharmap2$DNA_id

# make a sorted map object matching the ASV table order
sharma.smap <- sharmap2[match(rownames(sharma.copytab_cut_filtered), sharmap2$DNA_id), , drop = FALSE]
rownames(sharma.smap) <- sharma.smap$DNA_id


# making a phyloseq object
sharma.OTU <- otu_table(sharma.copytab_cut_filtered, taxa_are_rows = FALSE)
sharma.SAMPLE <- sample_data(sharma.smap)

sharma.physeq1_cut_filter <- phyloseq(sharma.OTU, sharma.SAMPLE)


phy_data <- asv_raw
rownames(metadata) <- metadata$Sample_Name
metadata <- metadata[, colnames(phy_data)] # Keep only necessary columns

physeq <- phyloseq(
  otu_table(phy_data, taxa_are_rows = FALSE),
  sample_data(metadata)
)

load('sharma_dump.Rdata')
load('sharma_dump/pipeline_sharma.Rdata')
# Add Taxonomy (If Family column is in the ASV table)
# physeq <- add_taxon_table(physeq, "family_column_name") # Adjust if needed
# For now, we assume Family is in the phy_data dataframe directly

# 6. Normalize (Relative Abundance)
physeq <- transform_sample_counts(physeq, function(x) x / sum(x))

# 7. Process Families (Step 3 in your workflow)
phy_data <- psmelt(physeq)

# Step 3: Filter empty families
family_stats <- phy_data %>%
  group_by(Group, Family) %>%
  summarize(
    Mean = mean(Abundance, na.rm = TRUE),
    SD = sd(Abundance, na.rm = TRUE),
    N = n(),
    .groups = 'drop'
  )
# Step 4: Plotting (Dodged Bars)
plot <- family_stats %>%
  ggplot(aes(x = Family, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, alpha = 0.7) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.5, position = "dodge", color = "black") +
  scale_fill_manual(values = c(
    "Captive Chimps" = "steelblue",
    "Chimps" = "dodgerblue",
    "USA-Human" = "orange"
  )) +
  labs(
    x = "Family", 
    y = "Mean Relative Abundance", 
    title = "Average Relative Abundance of Families by Group",
    subtitle = "Dodged bars: 3 bars per family representing each group"
  ) +
  theme_minimal() +
  theme(
    # ADD face = "italic" HERE:
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), 
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    legend.position = "bottom"
  )


# Display the plot
print(plot)

# 8. Save the Clean Object
# Overwrite any old files to ensure you have the latest clean data
save(physeq, file = output_file)
cat("Clean object saved to:", output_file, "\n")

# -------------------------------------------
# End of Clean Script
# -------------------------------------------





## PCoA pipeline
library(ggplot2)
library(phyloseq)

# 1. Calculate the PCoA (Ordination)
pcoa_res <- ordinate(sharma.physeq_top_norm, 
                     method = "PCoA", 
                     distance = "bray")

# 2. Create a data frame of the coordinates
# This is the safest way to ensure the plot works perfectly
df_pcoa <- data.frame(
  PCoA1 = pcoa_res$vectors[,1],
  PCoA2 = pcoa_res$vectors[,2],
  Group = sample_data(sharma.physeq_top_norm)$Group
)

# 3. Plot using standard ggplot2
ggplot(df_pcoa, aes(x = PCoA1, y = PCoA2, color = Group, shape = Group)) +
  # Add the points
  geom_point(size = 4) + 
  
  # ADD THE ELLIPSES (This replaces the ggforce bubbles)
  # type = "t" creates a confidence ellipse
  stat_ellipse(level = 0.95, linewidth = 1) + 
  
  # Make it look professional
  theme_bw() +
  labs(title = "PCoA of Gut Microbiomes",
     #  subtitle = "Convergence of Captive Chimps toward Human profile",
       x = "PCoA1", 
       y = "PCoA2") +
  theme(panel.grid = element_blank(), 
        axis.line = element_line(color = "black"))

# Look for "Proportion of variance explained"
summary(pcoa_res) 

# 1. Extract the eigenvalues (the 'values' mentioned in your summary)
ev <- pcoa_res$values

# 2. Calculate the percentage of variance for each axis
# Variance % = (Axis Value / Sum of all Values) * 100
perc_var <- (ev / sum(ev)) * 100

# 3. Show the results for the first two axes
print(perc_var[1:2])



## The proximity plot!! 
# 1. Get the PCoA coordinates
pcoa_df <- data.frame(
  PCoA1 = pcoa_res$vectors[,1],
  PCoA2 = pcoa_res$vectors[,2],
  Group = sample_data(sharma.physeq_top_norm)$Group
)

# 2. Calculate Centroids (the "average" position of each group)
centroids <- pcoa_df %>% 
  group_by(Group) %>% 
  summarise(mean1 = mean(PCoA1), mean2 = mean(PCoA2))

wild_center <- centroids %>% filter(Group == "Wild Chimps")
human_center <- centroids %>% filter(Group == "USA Humans")

# 3. Calculate distances for Captive Chimps
captive_samples <- pcoa_df %>% filter(Group == "Captive Chimps")

captive_samples <- captive_samples %>% 
  rowwise() %>% 
  mutate(
    dist_to_wild = sqrt((PCoA1 - wild_center$mean1)^2 + (PCoA2 - wild_center$mean2)^2),
    dist_to_human = sqrt((PCoA1 - human_center$mean1)^2 + (PCoA2 - human_center$mean2)^2),
    # Positive value = closer to wild, Negative = closer to human
    shift_score = dist_to_human - dist_to_wild 
  )

# 4. Plot the "Shift Score"
ggplot(captive_samples, aes(x = Group, y = shift_score)) +
  geom_boxplot(fill = "#4DAF4A", alpha = 0.7) + 
  geom_jitter(width = 0.1) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + # The "Neutral" line
  theme_bw() +
  labs(title = "Proximity to Wild vs. Human Microbiomes",
       subtitle = "Negative values indicate a shift toward the Human profile",
       y = "Distance Difference (Human - Wild)", 
       x = "Captive Chimps")



## Wild signature decay plot
library(dplyr)
library(ggplot2)

# 1. Calculate the mean abundance for every Taxon in each group
# I'm using 'Genus' here as it's more precise, but you can change to 'Family'
signature_calc <- df_melted %>% 
  group_by(Genus) %>% 
  summarise(
    mean_wild = mean(Abundance[Group == "Wild Chimps"], na.rm = TRUE),
    mean_captive = mean(Abundance[Group == "Captive Chimps"], na.rm = TRUE),
    mean_human = mean(Abundance[Group == "USA Humans"], na.rm = TRUE)
  ) %>% 
  # Use a tiny number (0.0001) instead of 0 to avoid "division by zero" errors
  mutate(ratio_wild_human = mean_wild / (mean_human + 0.0001))

# 2. Pick the Top 10 taxa with the HIGHEST ratio (the most "Wild" taxa)
top_wild_taxa <- signature_calc %>% 
  slice_max(ratio_wild_human, n = 10) %>% 
  pull(Genus)

# --- VERIFICATION CHECK ---
print("The Wild Signature Taxa are:")
print(top_wild_taxa)  ## OUTPUT OF THIS STEP: --> [1] "The Wild Signature Taxa are:"

#[1] "Sarcina"                   "Candidatus Methanogranum"  "Coriobacteriaceae UCG-003"
#[4] "Clostridium"               "Succinivibrio"             "Methanobrevibacter"       
#[7] "Bacteroides"               "Collinsella"               "Flexilinea"               
#[10] "Ileibacterium"             "Treponema"                 NA         

# If this list is empty, your 'Genus' column might be named differently (e.g., 'genus')
# ---------------------------

# 3. Sum the abundance of ONLY these top 10 taxa for every sample
shift_df <- df_melted %>% 
  filter(Genus %in% top_wild_taxa) %>% 
  group_by(DNA_id, Group) %>% 
  summarise(Total_Wild_Signature = sum(Abundance), .groups = 'drop')

# 4. Plot as a boxplot to show the "Decay"
ggplot(shift_df, aes(x = Group, y = Total_Wild_Signature, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
  theme_bw() +
  labs(
    title = "Decay of the 'Wild Signature' Microbiome",
    subtitle = "Total abundance of the top 10 taxa most characteristic of wild chimps",
    x = "Population Group",
    y = "Total Relative Abundance (%)"
  ) +
  scale_fill_manual(values = c("Wild Chimps" = "#E41A1C", 
                               "Captive Chimps" = "#4DAF4A", 
                               "USA Humans" = "#377EB8")) +
  theme(legend.position = "none")




## Alpha diversity code
library(vegan)

# 2. Extract the OTU table as a matrix
# We need the data as a matrix for the vegan package to read it
otu_matrix <- as.matrix(otu_table(sharma.physeq_top25))

# 3. Calculate the Shannon Index directly
# This bypasses the Chao1 calculation entirely, so it won't crash
shannon_values <- diversity(otu_matrix, index = "shannon")

# 4. Combine these values into a data frame with your Group labels
alpha_df <- data.frame(
  Shannon = shannon_values,
  Group = sample_data(sharma.physeq_top25)$Group
)

# 5. Plot the results
ggplot(alpha_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) + 
  theme_bw() +
  labs(
    title = "Reduction in Microbial Diversity",
    x = "Population Group",
    y = "Shannon Diversity Index"
  ) +
  scale_fill_manual(values = c("Wild Chimps" = "#E41A1C", 
                               "Captive Chimps" = "#4DAF4A", 
                               "USA Humans" = "#377EB8")) +
  theme(legend.position = "none")

## make a heatmap 2nd try
library(pheatmap)
library(phyloseq)
library(dplyr)
library(tidyr)

# 1. Convert to data frame
df_melted <- psmelt(sharma.physeq_top_norm)

# 2. Filter for Treponema
# adding a check to make sure we only get samples that have treponema
target_taxa <- df_melted %>% 
  filter(Family == "Treponema" | Genus == "Treponema")

# 3. Create the abundance matrix
heatmap_data <- target_taxa %>% 
  group_by(DNA_id, Group) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  tidyr::pivot_wider(names_from = DNA_id, values_from = Abundance)

# Save group info and create matrix
group_info <- heatmap_data$Group
abundance_matrix <- as.matrix(heatmap_data[, -1])
rownames(abundance_matrix) <- group_info

# 4. THE "SAFE" Z-SCORE
# Transpose, scale, then transpose back
z_score_matrix <- t(scale(t(abundance_matrix)))

# --- THIS IS THE FIX ---
# Replace all NaNs (caused by dividing by zero) with 0
z_score_matrix[is.na(z_score_matrix)] <- 0 
# -----------------------

# 5. THE PLOT
pheatmap(z_score_matrix, 
         cluster_cols = FALSE, # Keep sample order
         cluster_rows = FALSE, # TURN THIS OFF - prevents the hclust error
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Z-Score Relative Abundance of Treponema succinifaciens",
         show_colnames = FALSE, 
         border_color = NA)





##TOTAL SHIFT group average heatmap without the NA label!! 
library(pheatmap)
library(dplyr)
library(tidyr)
library(tibble)

# 1. Convert and CLEAN both Samples and Taxa
df_melted_clean <- psmelt(sharma.physeq_top_norm) %>% 
  filter(!is.na(Group)) %>% # Removes samples without a group
  filter(!is.na(Genus))     # <--- THIS IS THE FIX: Removes any taxon without a Genus name

# 2. Calculate the MEAN abundance per group
group_means <- df_melted_clean %>% 
  group_by(Genus) %>% 
  summarise(
    Wild = mean(Abundance[Group == "Wild Chimps"], na.rm = TRUE),
    Captive = mean(Abundance[Group == "Captive Chimps"], na.rm = TRUE),
    Human = mean(Abundance[Group == "USA Humans"], na.rm = TRUE)
  )

# 3. Create the matrix
heatmap_matrix <- as.matrix(group_means[, -1]) 
rownames(heatmap_matrix) <- group_means$Genus

# 4. Z-score the data
z_score_matrix <- t(scale(t(heatmap_matrix)))
z_score_matrix[is.na(z_score_matrix)] <- 0 
ta
## fixing the italics issue for the plot underneath --------> still not fixed 


# 5. The Plot

newnames <- lapply(
  rownames(heatmap_matrix),
  function(x) bquote(italic(.(x))))

pheatmap(z_score_matrix,
         labels_row = (as.expression(newnames)),
         cluster_cols = FALSE, 
         cluster_rows = TRUE,   
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Taxonomic Shift: Wild vs. Captive vs. Human",
         show_colnames = TRUE, 
         border_color = NA)


## CONTRAST BATTLE BAR PLOT
# 1. Define your "Characters"
wild_taxa <- c("Coriobacteriaceae UCG-003", "Sarcina", "Candidatus Methanogranum") # Replace with your actual wild names
human_taxa <- c("Bacteroides", "Flexilinea", "Methanobrevibacter") # Replace with your actual human names

# 2. Filter the data for only these "Characters"
contrast_df <- df_melted %>% 
  filter(Genus %in% c(wild_taxa, human_taxa)) %>% 
  group_by(Genus, Group) %>% 
  summarise(MeanAbund = mean(Abundance, na.rm = TRUE))

# 3. Create a "Category" label for the legend
contrast_df <- contrast_df %>% 
  mutate(Signature = ifelse(Genus %in% wild_taxa, "Wild-Associated", "Human-Associated"))

# 4. The Plot
ggplot(contrast_df, aes(x = Group, y = MeanAbund, fill = Signature)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() +
  labs(title = "The Great Shift: Wild vs. Human Signatures",
       y = "Mean Relative Abundance (%)", 
       x = "") +
  scale_fill_manual(values = c("Wild-Associated" = "#E41A1C", "Human-Associated" = "#377EB8")) +
  theme(axis.text.x = element_text(face = "bold"))


##Summary table most abundant genera --> to actually 
# This creates a summary table of the most abundant Genera
top_taxa_summary <- df_melted %>% 
  group_by(Genus) %>% 
  summarise(
    mean_wild = mean(Abundance[Group == "Wild Chimps"], na.rm = TRUE),
    mean_captive = mean(Abundance[Group == "Captive Chimps"], na.rm = TRUE),
    mean_human = mean(Abundance[Group == "USA Humans"], na.rm = TRUE)
  ) %>% 
  # Sort by the ones that are most abundant in Wild Chimps
  arrange(desc(mean_wild)) %>% 
  head(20)

# View the table
print(top_taxa_summary)
