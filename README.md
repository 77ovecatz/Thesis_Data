# Microbiome Analysis on Chimpanzees and Humans

The analyses are performed in RStudio.  

The project follows a standardized microbiome analysis workflow consisting of:
- Loading Data
- Pre-processing the Data
- Taxonomically assessing the processed Data 
- Performing Diversity analyses 
- Taxon investigation (Wen et al., 2023)

## Loading the Data
The data files are publicly available from the European Nucleotide Archive (ENA) browser.

## Pre-processing the Data
Initial processing is done with the following steps: 

### Quality Control
With the library `Cutadapt version 2.4` _(Martin, 2011)_
- Sequences have been identified and counted
- Adapter and primer sequences are removed
The sequences are now extracted and will be further cleaned and identified.

### ASV Identification 
With the libraries `DADA2` _(Callahan et al., 2016)_ and `phyloseq` _(McMurdie et al., 2013)_
- Sequences have been denoised
- Chimeric Sequences are removed
- Raw data is proceessed
- Amplicon Sequence Variants (ASV's) are identified
- The ASV data is integrated with asociated metadata into an Phyloseq object
The sequences have been processed, resulting in a ASV table, in an Phyloseq object which will be analyzed.

### Assessing the taxonomy
With data from the `SILVA` database _(Quastel et al., 2013)_ 
- The ASV table has been analysed
The output is a filter which exclusively selects the taxa within the relevant study groups. 

## Perfroming Diversity Analysis
With the library `vegan` _(Oksanen et al., 2022)_
The metrics that are calculated are:
- Alpha diversity
- Beta diversity
- Differential abundances
- Taxanomical analysis
The output of the analyses are used for my thesis.


## Data Visualisation
With the `ggplot2` as seen in _(Wickham, 2016)_ 
(placeholder images)

