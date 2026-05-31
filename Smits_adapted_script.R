## Processing sharma and Hansen datasets
## Need to be in parallel until taxonomy because different variable regions

library(dada2); packageVersion("dada2")


# sharma first 
# where are the files can be anywhere!
sharma.path <- "/home/amisha/sharma_dump/"

# where shall I Put my outputs HAS TO BE IN YOUR HOME DIRECTORY
sharma.path.output <- "/home/amisha/sharma_dump/"

# define the file names
sharma.fnFs <- sort(list.files(sharma.path, pattern="1.fastq.gz", full.names = TRUE))
sharma.fnRs <- sort(list.files(sharma.path, pattern="2.fastq.gz", full.names = TRUE))

################## Code to subset only the filenames relevant to you ################ 
# Tell R where to find the mapping file, and use the mapping file to tell it which runs are wanted.
sharmap2 <- read.csv2("/home/amisha/sharma_dump/sharmap2.csv", stringsAsFactors = FALSE)
 
sharma.wanted_runs <- unique(sharmap2$FASTQ.Run.Accession)

#set the names right
sharma.runs_in_files <- sapply(strsplit(basename(sharma.fnFs), "_"), `[`, 1)

#Quick check if it was able to find all the files
sharma.missing <- setdiff(sharma.wanted_runs, sharma.runs_in_files)

if(length(sharma.missing) > 0){
  stop("Missing FASTQ files for: ", paste(sharma.missing, collapse=", "))
}
#tell it that you only want to keep the names present in "sharma.wanted_runs" 
sharma.keep <- sharma.runs_in_files %in% sharma.wanted_runs

sharma.fnFs <- sharma.fnFs[sharma.keep]
sharma.fnRs <- sharma.fnRs[sharma.keep]
sharma.sample.names <- sharma.runs_in_files[sharma.keep]
##################### End of subset samples code ##########################

## extract sample names (Already done in the sample selection block above)
sharma.sample.names <- c(sapply(strsplit(basename(sharma.fnFs), "_"), `[`, 1))

## define primers to cut (found em in the article) 

FWD <- "GTGCCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence (i.e. the point where you want to cut, we determined this earlier looking at alignments)
REV <- "GGACTACHVGGGTWTCTAAT"  ## CHANGE ME...


cutadapt <- "/home/james/miniconda2/envs/qiime2-2019.7/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine (server version)
#cutadapt <- "D:/script_test/cutadapt.exe" # James PC version

system2(cutadapt, args = "--version") # check that cutadapt is working (you should see a version number)

#DEFINING OUTPUT PATH

## we should all define different output paths in our own server
## This will give issues if you didn't change the path.output to YOUR OWN server, as you don't have permission in others' folders.
sharma.path.cut <- file.path(sharma.path.output, "cutadapt") # assumes that you have write access to path
if(!dir.exists(sharma.path.cut)) dir.create(sharma.path.cut)
sharma.fnFs.cut <- file.path(sharma.path.cut, basename(sharma.fnFs))
sharma.fnRs.cut <- file.path(sharma.path.cut, basename(sharma.fnRs))

# prepare flags for cutadapt
sharma.R1.flags <- paste("-g", FWD) 
sharma.R2.flags <- paste("-G", REV) 

# Run Cutadapt   IN LINE 53 removed "--trimmed-only" It negatively impacted the results.
start <- Sys.time()
## NB I added the "--trimmed-only" flag, this will remove all sequence (pairs) where the primer sequence is not found. This might break your code.
for(i in seq_along(sharma.fnFs)) {
  system2(cutadapt, args = c(sharma.R1.flags, sharma.R2.flags, "-n", 2, "-o", sharma.fnFs.cut[i], "-p", sharma.fnRs.cut[i],   sharma.fnFs[i], sharma.fnRs[i]))
}
end <- Sys.time()
end-start # 14 seconds

cbind(precut = sapply(sharma.fnFs, function(x) length(readLines(x))/4),postcut =sapply(sharma.fnFs.cut, function(x) length(readLines(x))/4))

## Now normal filter and dada pipeline

sharma.filtFs_cut <- file.path(sharma.path.output, "filtered_cut", paste0(sharma.sample.names, "_F_filtcut.fastq.gz"))
sharma.filtRs_cut <- file.path(sharma.path.output, "filtered_cut", paste0(sharma.sample.names, "_R_filtcut.fastq.gz"))
names(sharma.filtFs_cut) <- sharma.sample.names
names(sharma.filtRs_cut) <- sharma.sample.names

#You may need to adjust the truncLen= argument depending on which region you are looking at 
start <- Sys.time()
out <- filterAndTrim(sharma.fnFs.cut, sharma.filtFs_cut, sharma.fnRs.cut, sharma.filtRs_cut, truncLen=c(160,160),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=10)
end <- Sys.time()
end-start # 13 seconds for 6 samples (on single core, James' PC)

out # most reads retained??? YES! (for 160,160 and maxEE 1,1)

start <- Sys.time()
sharma.errF.1_cut <- learnErrors(sharma.filtFs_cut, multithread=10)
sharma.errR.1_cut <- learnErrors(sharma.filtRs_cut, multithread=10)
end <- Sys.time()
end-start # 1.4 minutes for 10 samples 18 threads 


start <- Sys.time()
sharma.dadaFs_cut.1 <- dada(sharma.filtFs_cut, err=sharma.errF.1_cut, multithread=10)
sharma.dadaRs_cut.1 <- dada(sharma.filtRs_cut, err=sharma.errR.1_cut, multithread=10)


end <- Sys.time()
end-start #~4 minutes for 6 samples

sharma.mergers_cut <- mergePairs(sharma.dadaFs_cut.1, sharma.filtFs_cut, sharma.dadaRs_cut.1, sharma.filtRs_cut, verbose=TRUE)

####################################### End of processing ##########################################################
## convert to sequence tables
sharma.seqtab_cut <- makeSequenceTable(sharma.mergers_cut) #warning for me: here sharma will run with the prefix sharma. for sharma.seqtab_cut for the first time.
dim(sharma.seqtab_cut)
table(nchar(getSequences(sharma.seqtab_cut))) ##credible length dist

## remove chimeras
start <- Sys.time()
sharma.seqtab_nochim <- removeBimeraDenovo(sharma.seqtab_cut, method="consensus", multithread=TRUE)
end <- Sys.time()
end-start # 10 seconds for 42 samples
dim(sharma.seqtab_nochim) # removed 60 ASVs

sharma.copytab_cut <- sharma.seqtab_nochim 

################## ASV ######################

# renaming ASVs for convenience
colnames(sharma.copytab_cut) <- paste0("ASV", 1:dim(sharma.copytab_cut)[2])

# how many reads per ASV?
plot(colSums(sharma.copytab_cut), log="y")

# what proportion of ASVs have less than 5 reads?
mean(colSums(sharma.copytab_cut) < 5)

# how many samples does each ASV appear in?
hist(as.numeric(colSums(sharma.copytab_cut > 0)),
     xlab = "Number of samples",
     ylab = "Number of ASVs",
     main = "Occurrences per ASV")

## and if we filter out all ASVs with read number < 10 --> replaced with 3

sharma.copytab_cut_filtered <- sharma.copytab_cut[, colSums(sharma.copytab_cut) > 3]
dim(sharma.copytab_cut_filtered) 
hist(as.numeric(colSums(sharma.copytab_cut_filtered > 0)),
     xlab = "Number of samples",
     ylab = "Number of ASVs",
     main = "Occurences per ASV (3+ reads)") 

## How many ASVs had fewer than 5 reads? maybe good for statistics since it can justify the filtering threshold of 10 reads.
mean(colSums(sharma.copytab_cut) < 3)

### Saving image to restart later
save.image("/home/amisha/sharma_dump.Rdata") #Path on the server "/home/[YOURNAME]/some_other_folder/goodfilename.Rdata")

View(sharma.copytab_cut_filtered)

write.csv(sharma.copytab_cut_filtered, 'sharma_dump/asv2.csv')

###################################### POST PROCESSING ################################
load("/home/amisha/Images/sharma_pipeline.Rdata")
library(vegan)
library(phyloseq)

# Read the mapping file 
sharmap2 <- read.csv2("~/sharma_dump/sharmap2.csv", 
                      stringsAsFactors = FALSE)

library(dplyr)

# assign the result back to sharmap2 and tell mutate WHICH data frame to use (using the %>% pipe)
sharmap2 <- sharmap2 %>% 
  mutate(Group = recode(Group, 
                        "Chimp" = "Wild Chimps", 
                        "Captive Chimps" = "Captive Chimps",   # Added this so it doesn't become NA
                        "USA-Human" = "USA Humans"))

sharmap2$Group <- factor(sharmap2$Group, levels = c("Wild Chimps", "Captive Chimps", "USA Humans"))

table(sharmap2$Group)
View(sharmap2)



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


# checking the range in reads per sample
range(rowSums(otu_table(sharma.physeq1_cut_filter)))

# Normalization (replaces rarefying) Normalizing to relative abundance (%)
sharma.physeq2_cut_filter <- transform_sample_counts(sharma.physeq1_cut_filter,  function(OTU) OTU / sum(OTU) * 100)

# checking the range in reads per sample
range(rowSums(otu_table(sharma.physeq2_cut_filter)))


norm_combined_physeq_top100 <- transform_sample_counts(sharma.physeq2_cut_filter, function(OTU) OTU/sum(OTU) * 100)

# checking the range in reads per sample
range(rowSums(otu_table(norm_combined_physeq_top100)))

plot_bar(norm_combined_physeq_top100, x = "DNA_id", fill = 'Group',) +
  facet_grid(~Group, scales="free_x", space = "free_x")+   #waarschijnlijk ~Population niet goed HOOOIIIII DEZE FF CHECKEN
  geom_bar(stat = "identity", color = NA) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

#aantekeningviir mezelf: plot iets dat logisch is: Waarschijnlijk is sharmap2 incompleet

#alternative code (James)
norm_combined_physeq_top100 <- transform_sample_counts(sharma.physeq1_cut_filter, function(OTU) OTU/sum(OTU) * 100)
norm_combined_physeq2_top100 <- transform_sample_counts(sharma.physeq2_cut_filter, function(OTU) OTU/sum(OTU) * 100)


plot_bar(norm_combined_physeq1_top100, x = "DNA_id", fill = "Phylum",) +
  facet_grid(~Population, scales="free_x", space = "free_x")+
  geom_bar(stat = "identity", color = NA) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# performing NMDS ordination 
sharma.mm <- metaMDS(vegdist(otu_table(sharma.physeq2_cut_filter)))

################## Plotting ##################
# plot NMDS
plot(sharma.mm, type = "n", display = "sites", main = "sharma")
sharma.mm

# color points by Group (grouping column)
sharma.sites <- scores(sharma.mm, display = "sites")
sharma.grp <- factor(sharma.smap$Group[match(rownames(sharma.sites), rownames(sharma.smap))])

sharma.cols <- c("red", "blue")
points(sharma.sites, pch = 16, col = sharma.cols[sharma.grp])

#the thing that creates the actual graph
ordihull(sharma.mm, sharma.grp, display="sites", draw="polygon", col=sharma.cols, label=TRUE)

# if you see all your points grouped according to studyID (or however else you have coded it in your mapping file) then it has worked!

  

##################### taxonomy #######################
library(dada2)
library(phyloseq)
library(ggplot2)

# choose the ASV table to annotate (nonchim is standard)
sharma.seqtab_use <- sharma.seqtab_nochim   # Check if the "seqtab_nochim" has "sharma." in front of it in the entire script using search function.

# pick top 100 most abundant ASVs across all samples
top25 <- order(colSums(sharma.seqtab_use)[1:25], decreasing = TRUE) #remove the 100 if you want everything, but that will take long.


# assign taxonomy - UPDATED SILVA FILE
#Silva reference file goes only up to genus level, not specie: "/home/james/tax/silva_nr99_v138.1_train_set.fa.gz"
sharmatop.taxa <- assignTaxonomy(
  sharma.seqtab_use[, top25, drop = FALSE],
  "/home/fakestudent/tax/silva_nr99_v138.2_toSpecies_trainset.fa.gz",
  multithread = 20
)


#########pipeline laatste stap
# make physeq subset of these ASVs
top.OTU2 <- otu_table(sharma.seqtab_use, taxa_are_rows = FALSE)

rownames(sharma.taxa) <- taxa_names(top.OTU2)
sharma.top <- phyloseq(top.OTU2, sharma.SAMPLE, tax_table(sharma.taxa))

library(ggplot2)

plot_bar(sharma.top, x = "DNA_id", fill = "Class") +
  geom_bar(stat = "identity", color = NA)


##################################################
library(dplyr)
library(phyloseq)
library(ggplot2)

# --- STEP 1: Fix the Mapping File (sharmap2) ---
sharmap2 <- as.data.frame(sharmap2) 
rownames(sharmap2) <- sharmap2$DNA_id

sharmap2 <- sharmap2 %>% 
  mutate(Group = recode(Group, 
                        "Chimp" = "Wild Chimps", 
                        "Captive" = "Captive Chimps", 
                        "USA-Human" = "USA Humans")) %>% 
  mutate(Group = factor(Group, levels = c("Wild Chimps", "Captive Chimps", "USA Humans")))
sharma.SAMPLE_NEW <- sample_data(sharmap2)

# --- STEP 2: Subset OTUs Safely (Fixes 'Subscript out of bounds') ---
asv_totals <- colSums(sharma.seqtab_use)
top <- order(asv_totals, decreasing = TRUE)[1:25] # Set to 25 or 100
sharmatop.OTU2 <- otu_table(sharma.seqtab_use[, top, drop = FALSE], taxa_are_rows = FALSE)

# --- STEP 3: Subset Taxonomy and Convert to Matrix (Fixes 'Taxa names do not match') ---
sharma.taxa_subset <- sharma.taxa[taxa_names(sharmatop.OTU2), ]
sharma.taxa_matrix <- as.matrix(sharma.taxa_subset)

# --- STEP 4: Build Phyloseq Object ---
sharma.physeq_top25 <- phyloseq(sharmatop.OTU2, 
                                sharma.SAMPLE_NEW, 
                                tax_table(sharma.taxa_matrix))

# --- STEP 5: Normalize and Plot ---
sharma.physeq_top_norm <- transform_sample_counts(
  sharma.physeq_top25,
  function(OTU) OTU / sum(OTU) * 100
)

# Using psmelt for the plot to avoid the 'Population/Group missing' error
df_melted <- psmelt(sharma.physeq_top_norm)

family_names = sort(unique(df_melted$Family))
newnames <- lapply(
  family_names,
  function(x) bquote(italic(.(x))))


ggplot(df_melted, aes(x = DNA_id, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", color = NA) +
  scale_fill_discrete("Family",
                       labels = as.expression(newnames))+ 
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  ylab("Relative abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.text = element_text(hjust = 0))


#iets met element_text(face = "italic")

############### end 
##oh no i want to make a table out of it!! WITH p-values!!!
library(dplyr)
library(tidyr)
library(writexl)

# 1. Calculate the Means per Group (as we did before)
means_table <- df_melted %>%
  group_by(Group, Family) %>%
  summarise(Mean_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = Group, values_from = Mean_Abundance) %>%
  mutate(across(-Family, ~ replace_na(.x, 0)))

# 2. Calculate P-Values for each Family
# We use Kruskal-Wallis because microbiome data is rarely normally distributed
p_values_table <- df_melted %>%
  group_by(Family) %>%
  summarise(
    p_value = kruskal.test(Abundance ~ Group)$p.value, 
    .groups = 'drop'
  ) %>%
  # Add a column for "Significance" stars to make the table look professional
  mutate(Significance = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE            ~ "ns" # ns = not significant
  ))

# 3. Merge the Means and the P-Values together
final_analysis_table <- left_join(means_table, p_values_table, by = "Family")

# 4. Reorder columns so p-value is at the end
final_analysis_table <- final_analysis_table %>% 
  select(Family, everything(), p_value, Significance)

# 5. Export to Excel
write_xlsx(final_analysis_table, "Family_Abundance_with_Pvalues.xlsx")

# View the result
print(final_analysis_table)

####################

# make phyloseq subset of these ASVs
sharmatop.OTU2 <- otu_table(sharma.seqtab_use[, top, drop = FALSE], taxa_are_rows = FALSE)

# make sure tax table rownames match OTU taxa names
# Added some code to remove empty spaces in the graph.
rownames(sharma.taxa) <- taxa_names(sharmatop.OTU2)

sharma.physeq_top25 <- phyloseq(sharmatop.OTU2, sharma.SAMPLE, tax_table(sharma.taxa))

sample_data(sharma.physeq_top25)$SampleIndex <- seq_len(nsamples(sharma.physeq_top25))

# normalize to relative abundance (%)
sharma.physeq_top_norm <- transform_sample_counts(
  sharma.physeq_top25,
  function(OTU) OTU / sum(OTU) * 100
)

p <- plot_bar(sharma.physeq_top_norm, x="DNA_id", fill="Family") +
  geom_bar(stat="identity", color=NA) +
  facet_grid(~Group, scales="free_x", space="free_x") +
  ylab("Relative abundance (%)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
print(p)

##screenshot but fancy
pdf("/home/amisha/sharma_dump/sharmatop25_species_barplot.pdf", width = 12, height = 6)
print(p)
dev.off()
# the extra stuff like "element_blank" is now included in the plot instead.

save.image("/home/amisha/sharma_dump/pipeline_sharma.Rdata") #Path on the server "/home/[YOURNAME]/some_other_folder/goodfilename.Rdata")


# done :)

############

library(phyloseq)
library(dplyr)
library(ggplot2)

# Load packages and your object
physeq <- sharma.physeq_top_norm  # Replace with your object name

# Check group names in your metadata
unique(sample_data(physeq)$Group)  # Replace if column name differs

# Step 1: Convert to relative abundance
physeq <- transform_sample_counts(physeq, function(x) x / sum(x))

# Step 2: Melt the phyloseq object
phy_data <- psmelt(physeq)

# Step 3: Keep only the Family column and your group column
# Adjust column name if needed (e.g. "Family", "Taxon", "fam")
family_stats <- phy_data %>%
  filter(Family != "") %>%
  group_by(Group, Family) %>%
  summarize(Mean_Rel_Abundance = mean(Abundance, na.rm = TRUE), .groups = 'drop')

family_stats <- phy_data %>%
  group_by(Group, Family) %>%
  summarize(
    Mean = mean(Abundance, na.rm = TRUE),
    SD = sd(Abundance, na.rm = TRUE),
    N = n(),
    .groups = 'drop'
  )

# Step 4: Create a dodged barplot (3 bars per family = 3 groups)
# This plots Captive Chimps, Chimps, and USA-Human side-by-side for each family
plot <- family_stats %>%
  ggplot(aes(x = Family, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.8, position = "dodge", color = "black") +
  scale_fill_manual(values = c("Captive Chimps" = "steelblue",
                               "Wild Chimps" = "dodgerblue",
                               "USA Humans" = "orange")) +
  labs(
    x = "Family", 
    y = "Mean Relative Abundance", 
    title = "Average Relative Abundance of Families by Group",
    subtitle = "Dodged bars: 3 bars per family representing each group"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10)
  )


# Display the plot
print(plot)



######
# Add standard deviation
family_stats <- phy_data %>%
  group_by(Group, Family) %>%
  summarize(
    Mean = mean(Abundance, na.rm = TRUE),
    SD = sd(Abundance, na.rm = TRUE),
    N = n(),
    .groups = 'drop'
  )

# Plot with error bars
plot <- family_stats %>%
  ggplot(aes(x = Family, y = Mean, fill = Group)) +
  scale_x_discrete(guide = guide_axis(angle = 45))+ 
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, position = "dodge", color = "black") +
  facet_wrap(~ Group, ncol = 1) +
  theme_minimal()

print(plot)
