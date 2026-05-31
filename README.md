# Microbiome Analysis on Chimpanzees and Humans

The analyses are performed in RStudio.  

The project follows a standardized microbiome analysis workflow consisting of:
- Data pre-processing
- Taxonomic assessment 
- Diversity analysis 
- Taxon investigation (Wen et al., 2023)


## Initial processing 
Initial processing includes: 
 ## quality testing 
 ## Cutadapt version 2.4 (Martin, 2011) for the removal of adapter and primer sequences.

Afterwards, the DADA2 package (Callahan et al., 2016) was used for 
denoising, 
chimera removal, 
and to process the raw sequencing data and 
identify amplicon sequence variants (ASVs). 
Afterwards, an ASV table was extracted from the sample set, and the SILVA database was used for the taxonomic classification (Quastel et al., 2013). 
A filter was added in the form of a mapping file (see appendix A). 
This reference database allowed for the assignment of bacterial taxonomies exclusively for the taxa within the relevant study groups. 

The phyloseq package (McMurdie et al., 2013) was used for data integration in generating an object containing the ASV table and associated metadata for statistical analyses. The vegan package (Oksanen et al., 2022) was used for the diversity metrics. 
For data visualization purposes, the ggplot2 package was used (Wickham, 2016). 
Sharma, mapping file
