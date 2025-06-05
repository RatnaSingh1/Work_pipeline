#for this analysis three files are required rooted-tree.qza, table.qza and taxonomy.qza
#convert the rooted-tree.qza, table.qza and taxonomy.qza into the tsv file using export command
#qiime tools export --input-path table.qza --output-path  table  # cd table
#biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
# qiime tools export --input-path rep-seqs-final_138.qza --output-path sequences   #sequences
#qiime tools export --input-path rooted-tree.qza --output-path rooted_phylo_tree     # for nwk format
#and converting biom to tsv using command
#biom convert -i feature-table.biom -o feature-table.tsv --to-tsv
#the three tsv files  feature-table.tsv, taxonomy.tsv and rooted.nwk with sampledata tsv will be use for the anaysis in R

#aparthis two more files required that are unweighted_unifrac_pcoa_results.qza, shannon.qza
lapply_vector <- c("DESeq2", "tidyverse", "pheatmap", "ggplot2", "RColorBrewer", "cowplot", "dplyr", "cluster", "phyloseq", "vegan", "ggplot2", "ANCOMBC", "dplyr", "ComplexHeatmap", "DECIPHER", "ape", "qiime2R", "microbiome")
lapply(lapply_vector, library, character.only = TRUE) 
library(openxlsx)
library(RColorBrewer)

# Define the getPalette function
colourCount <- 40  # Example number of colors you need
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))


otu <- read.table(file = "feature-table.tsv", sep = "\t", header = T)
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T)
metadata <- read.table(file = "metadata1.tsv", sep = "\t", header = T)

#importing shannon in qza format
shannon<-read_qza("shannon_vector.qza")

#the purpose of this code is to add a new column called "SampleID" to the data frame/tibble shannon$data. 
#The values in this column will be the row names of the original data frame/tibble, providing a unique identifier
#for each row/sample in the data.By executing this code, you will have a modified version of the shannon$data data 
#frame/tibble with the row names now appearing as a separate column called "SampleID".
shannon1<-shannon$data %>% rownames_to_column("SampleID")

# Convert "SampleID" column to integer in shannon data frame
shannon1 <- shannon$data %>%
  rownames_to_column("SampleID")

# Convert "SampleID" column to integer in metadata data frame joining metadata and shannon_entropy via SampleID
metadata1 <- metadata %>%
  left_join(shannon1, by = "SampleID")
# Analysing the data according to diversity


desired_order <-  c("Ctrl24", "EcN24h", "EPEC24h", "SK22D24h", "EcNEPEC24h", "SK22DEPEC24h", "EcNSupernat24h", "SK22DSupernat24h", "EcNSupernatEPEC24h", "SK22DSupernatEPEC24h", "EcNfed224h", "Sk22Dfed224h",
                                                         "Cntrl48", "EcN48h", "EPEC48h", "SK22D48h", "EcNEPEC48h", "SK22DEPEC48h", "EcNSupernat48h", "SK22DSupernat48h", "EcNSupernatEPEC48h", "SK22DSupernatEPEC48h")

metadata1$Sample_Name <- factor(metadata1$Sample_Name, levels = desired_order)
metadata1 %>%
  ggplot(aes(x = Sample_Name, y = shannon_entropy, color = Sample_Name)) +
  geom_boxplot() +
  geom_point(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.75)) +
  xlab("SampleID") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  facet_wrap(~ Days, scales = "free_x") +  # Facet by Days
  theme(legend.position = "center", 
        plot.title = element_text(hjust = 0.5, size = 12),  # Adjust title size
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11),   
        axis.text.y = element_text(size = 10),    
        axis.title.x = element_text(size = 10),  # Control x-axis title size
        axis.title.y = element_text(size = 10),  # Control y-axis title size
        legend.text = element_text(size = 10),   # Legend text size
        legend.title = element_text(size = 10),  # Legend title size
        strip.text = element_text(size = 10)) +  # Facet strip text size
  ggtitle("Alpha Diversity")


ggsave("./shannon/Shannon_by_treatment_box.png", height = 8, width = 10, device = "png")
ggsave("./shannon/Shannon_by_treatment_box.pdf", height = 5, width = 10, device = "pdf")

# Filter the data for specific groups in R
filtered_data <- metadata1[metadata1$Sample_Name %in% c("Ctrl24", 'EcN24h', "EPEC24h"), ]

filtered_data  %>%
  ggplot(aes(x = SampleID, y = shannon_entropy, color = SampleID)) +
  geom_boxplot() +
  geom_point(stat = "summary", fun.data = mean_se, position = position_dodge(width = 0.75)) +
  xlab("SampleID") +
  ylab("Shannon Diversity") +
  theme_classic()+
  facet_wrap(~ Days, scales = "free_x")+  # Facet by Days
  theme(legend.position = "center", 
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,  size = 12))+ # Rotate x-axis labels to 45 degrees
  # Increase the size for y-axis text
  ggtitle("Alpha Diversity within Samples") 



##PoCA plot

uwunifrac <- read_qza("unweighted_unifrac_pcoa_results.qza")



plot_data <- uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon1)

unique_treatments <- unique(plot_data$SampleID)
num_treatments <- length(unique_treatments)

# Adjust the shape values to match the number of unique treatments
shape_values <- 1:num_treatments

desired_order <-  c("Ctrl24", "EcN24h", "EPEC24h", "SK22D24h", "EcNEPEC24h", "SK22DEPEC24h", "EcNSupernat24h", "SK22DSupernat24h", "EcNSupernatEPEC24h", "SK22DSupernatEPEC24h", "EcNfed224h", "Sk22Dfed224h",
                    "Cntrl48", "EcN48h", "EPEC48h", "SK22D48h", "EcNEPEC48h", "SK22DEPEC48h", "EcNSupernat48h", "SK22DSupernat48h", "EcNSupernatEPEC48h", "SK22DSupernatEPEC48h")

plot_data$Sample_Name <- factor(plot_data$Sample_Name, levels = desired_order)

#clustering according to days
#filled with color and same shape,clustering according to days
ggplot(plot_data, aes(x = PC1, y = PC2, fill = Sample_Name, size = shannon_entropy)) +
  geom_point(shape = 21, alpha = 0.7, color = "black") + # shape 21 is a filled circle, color is for the outline
  scale_size_continuous(name = "Shannon Entropy") +  # Adjust the size range as desired
  xlab("PC1") +
  ylab("PC2")+
  theme_q2r() +
  theme(legend.position = "right",
        axis.text.x.bottom = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 10),              # Increase the size for y-axis text
        legend.text = element_text(size = 10),              # Increase the size for legend text
        strip.text = element_text(size = 10) )+
  facet_wrap(~ Days, scales = "free_x")
ggsave("./PoCA/PoCA_days.png", height = 8, width = 10, device = "png")
ggsave("./PoCA/PoCA_days.pdf", height = 8, width = 10, device = "pdf")


ggplot(plot_data, aes(x = PC1, y = PC2, fill = Sample_Name, size = shannon_entropy)) +
  geom_point(shape = 21, alpha = 0.7, color = "black") +
  scale_size_continuous(name = "Shannon Entropy") +
  xlab("PC1") +
  ylab("PC2") +
  theme_q2r() +
  theme(
    legend.position = "right",
    axis.text.x.bottom = element_text(angle = -90, size = 12),
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 10)
  ) +
  facet_wrap(~ Days, scales = "free_x") +
  scale_fill_manual(values = c("red", "blue", "green", "purple", "brown", "black", "magenta", "yellow", "cyan", "orange", "darkgreen",
                               "azure3", "aquamarine", "burlywood1", "chartreuse", "cyan4", "darkgoldenrod", "darkolivegreen2", "deepskyblue", "darkseagreen", "darkorchid4","brown3" ))


ggsave("./PoCA/PoCA_day_1.png", height = 8, width = 10, device = "png")
ggsave("./PoCA/PoCA_day_1.pdf", height = 5, width = 10, device = "pdf")






library(tidyverse)
library(qiime2R)
#Heatmap

metadata_r <- read_q2metadata("metadata1.tsv")
SVs <- read_qza("qiime/table.qza")$data
taxonomy <- read_qza("qiime/taxonomy.qza")$data

SVs <- apply(SVs, 2, function(x) x/sum(x)*100) # convert to percent

SVsToPlot <-  
  data.frame(MeanAbundance = rowMeans(SVs))%>% #find the average abundance of a SV
  rownames_to_column("Feature.ID") %>%
  arrange(desc(MeanAbundance)) %>%
  top_n(50, MeanAbundance) %>%
  pull(Feature.ID) # extract only the names from the table

joined <- SVs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")%>%
  gather(-Feature.ID, key = "SampleID", value = "Abundance")%>%
  mutate(Feature.ID = if_else(Feature.ID %in% SVsToPlot, Feature.ID, "Remainder")) %>% 
  group_by(SampleID, Feature.ID)%>%
  summarize(Abundance = sum(Abundance))%>%
  left_join(metadata_r) %>%
  mutate(NormAbundance = log10(Abundance + 0.001))%>% # do a log10 transformation after adding a 0.001% pseudocount. Could also add 1 read before transformation to percent
  left_join(taxonomy) %>%
  mutate(Feature = paste(Feature.ID, Taxon)) %>%
  mutate(Feature = gsub("[kpcofgs]__", "", Feature))

joined$SampleID <- factor(joined$SampleID, levels = c("Ctrl24", "EcN24h", "EPEC24h", "SK22D24h", "EcNEPEC24h", "SK22DEPEC24h", "EcNSupernat24h", "SK22DSupernat24h", "EcNSupernatEPEC24h", "SK22DSupernatEPEC24h", "EcNfed224h", "Sk22Dfed224h",
                                                      "Cntrl48", "EcN48h", "EPEC48h", "SK22D48h", "EcNEPEC48h", "SK22DEPEC48h", "EcNSupernat48h", "SK22DSupernat48h", "EcNSupernatEPEC48h", "SK22DSupernatEPEC48h"))

plot_heatmap <- ggplot(data = as.data.frame(joined ), aes(x = SampleID , y = Feature, fill = NormAbundance)) +
  geom_tile() +
  facet_grid(~`Days`, scales = "free_x") +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c(name = "log10(% Abundance)")

ggsave("./heatmap/heatmap1.png", plot_heatmap, height = 8, width = 20, device = "png")
ggsave("./heatmap/heatmap2_1.pdf", plot_heatmap, height = 8, width = 20, device = "pdf")




#Next phase analysis ,though some are common as above

otu <- read.table(file = "feature-table.tsv", sep = "\t", header = T, row.names = 1, 
                  skip = 1, comment.char = "")
taxonomy <- read.table(file = "taxonomy.tsv", sep = "\t", header = T ,row.names = 1)

# clean the taxonomy, Greengenes format
tax <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right")


tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "k__",""),
                        Phylum = str_replace(tax[,2], "p__",""),
                        Class = str_replace(tax[,3], "c__",""),
                        Order = str_replace(tax[,4], "o__",""),
                        Family = str_replace(tax[,5], "f__",""),
                        Genus = str_replace(tax[,6], "g__",""),
                        Species = str_replace(tax[,7], "s__",""),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""
tax.clean[tax.clean=="__"] <- ""

for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,7] != ""){
    tax.clean$Species[i] <- paste(tax.clean$Genus[i], tax.clean$Species[i], sep = " ")
  } else if (tax.clean[i,2] == ""){
    kingdom <- paste("Unclassified", tax.clean[i,1], sep = " ")
    tax.clean[i, 2:7] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Unclassified", tax.clean[i,2], sep = " ")
    tax.clean[i, 3:7] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Unclassified", tax.clean[i,3], sep = " ")
    tax.clean[i, 4:7] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Unclassified", tax.clean[i,4], sep = " ")
    tax.clean[i, 5:7] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Unclassified", tax.clean[i,5], sep = " ")
    tax.clean[i, 6:7] <- family
  } else if (tax.clean[i,7] == ""){
    tax.clean$Species[i] <- paste("Unclassified ",tax.clean$Genus[i], sep = " ")
  }
}

metadata <- read.table(file = "metadata1.tsv", sep = "\t", header = T, row.names = 1)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = tax_table(as.matrix(tax.clean))
SAMPLE <- sample_data(metadata)
TREE = read_tree("rooted_tree.nwk")
# merge the data
library(phytools)
ps1 <- phyloseq(SAMPLE, OTU, TAX, TREE)
#ps <- phyloseq(OTU, TAX, SAMPLE)
set.seed(111) # keep result reproductive
ps.rarefied = rarefy_even_depth(ps1, rngseed=1, sample.size=1103, replace=F)



#heatmap after ps1 normalization'


# Normalize by relative abundance
ps1_norm <- transform_sample_counts(ps1, function(x) x / sum(x))


# Get the mean abundance of each OTU across all samples
mean_abundance <- taxa_sums(ps1_norm) / nsamples(ps1_norm)

# Select the top 50 most abundant OTUs
top_50_otus <- names(sort(mean_abundance, decreasing = TRUE)[1:50])

# Prune the phyloseq object to only include the top 50 OTUs
ps1_top50 <- prune_taxa(top_50_otus, ps1_norm)
n <- plot_heatmap(ps1_top50)+
  geom_tile() +
  facet_grid(~`Days`, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 8),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)) 

ggsave("./heatmap/heatmap5.pdf", n, height = 12, width = 20, device = "pdf")








#heat map without normalization


# Get the total abundance of each OTU across all samples
total_abundance <- taxa_sums(ps1)

# Select the top 50 most abundant OTUs
top_50_otus <- names(sort(total_abundance, decreasing = TRUE)[1:50])

# Prune the phyloseq object to only include the top 50 OTUs
ps1_top50 <- prune_taxa(top_50_otus, ps1)


n <- plot_heatmap(ps1_top50)+
  geom_tile() +
  facet_grid(~`Days`, scales = "free_x") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)) 

ggsave("./heatmap/heatmap6_no_normalization.pdf", n, height = 12, width = 20, device = "pdf")


#a simple alpha diversity plot


# Reorder the Sample_Name variable in ps.rarefied according to your desired order
desired_order <- c("Ctrl24", "EcN24h", "EPEC24h", "SK22D24h", "EcNEPEC24h", "SK22DEPEC24h", "EcNSupernat24h", "SK22DSupernat24h", "EcNSupernatEPEC24h", "SK22DSupernatEPEC24h", "EcNfed224h", "Sk22Dfed224h",
                   "Cntrl48", "EcN48h", "EPEC48h", "SK22D48h", "EcNEPEC48h", "SK22DEPEC48h", "EcNSupernat48h", "SK22DSupernat48h", "EcNSupernatEPEC48h", "SK22DSupernatEPEC48h")
# Replace with your actual sample names
ps.rarefied@sam_data$Sample_Name <- factor(ps.rarefied@sam_data$Sample_Name, levels = desired_order)

# Now plot using the reordered Sample_Name
plot_richness(ps.rarefied, x = "Sample_Name", measures = c("Shannon"), color = "Sample_Name") +
  geom_boxplot() +  # Set fill to NA for transparency and color for outline
  facet_wrap(~ Days, scales = "free_x") +
  theme_q2r() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 20),
    axis.title = element_text(size = 14),  # Increase axis title size
    axis.text = element_text(size = 10),   # Increase axis text size
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    strip.text = element_text(size = 10)  # Increase facet label size
  )


 
ggsave("./alpha_diversity/alpha_diversity_shannon_days.png", height = 5, width = 9, device = "png")
ggsave("./alpha_diversity/alpha_diversity_shannon_days.pdf", height = 5, width = 8, device = "pdf")



#beta diversity though unifrac is more accurate than this one because this plotted after the rarefecation

dist = phyloseq::distance(ps.rarefied, method="bray")
ordination = ordinate(ps.rarefied, method="PCoA", distance=dist)

plot_ordination(ps.rarefied, ordination, color="Sample_Name") + 
  facet_wrap(~ Days, scales = "free_x" )
ggsave("./PoCA/PoCA_ps_rarefied.png", height = 5, width = 10, device = "png")


#specified color
dist <- phyloseq::distance(ps.rarefied, method = "bray")
ordination <- ordinate(ps.rarefied, method = "PCoA", distance = dist)

plot_ordination(ps.rarefied, ordination, color = "Sample_Name") +
  facet_wrap(~Days, scales = "free_x") +
  theme_classic()+
  scale_color_manual(values = c("red", "blue", "green", "purple", "brown", "black", "magenta", "yellow", "cyan", "orange", "darkgreen",
                                "azure3", "aquamarine", "burlywood1", "chartreuse", "cyan4", "darkgoldenrod", "darkolivegreen2", "deepskyblue", "darkseagreen", "darkorchid4", "brown3"))
ggsave("./PoCA/PoCA_ps_rarefied_colors.pdf", height = 5, width = 10, device = "pdf")


#abundance bar plot genus
# Normalize counts to relative abundance
ps.rel = transform_sample_counts(ps1, function(x) x/sum(x)*100)

# Agglomerate taxa at the Class level
glom <- tax_glom(ps.rel, taxrank = 'Class', NArm = FALSE)

# Melt the phyloseq object to long format
ps.melt <- psmelt(glom)

# Convert Class to character for easier manipulation
ps.melt$Class <- as.character(ps.melt$Class)

# Calculate median abundance per sample and class
ps.melt <- ps.melt %>%
  group_by(Sample_Name, Class) %>%
  mutate(median_abundance = median(Abundance))

# Identify classes with a median abundance greater than 1
keep <- unique(ps.melt$Class[ps.melt$median_abundance > 1])

# Assign "< 1%" to classes not in the 'keep' list
ps.melt$Class[!(ps.melt$Class %in% keep)] <- "< 1%"

# Summarize abundance by sample and class
ps.melt_sum <- ps.melt %>%
  group_by(Sample_Name,Days, Class) %>%
  summarise(Abundance = sum(Abundance))

# Plot the summarized data
ggplot(ps.melt_sum, aes(x = Sample_Name, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y = "%") +
  facet_wrap(~ Days, scales = "free_x" )+
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 12),              
        strip.text = element_text(size = 14))
ggsave("./abundance_plot/abundance_Class_days.pdf", height = 14, width = 25, device = "pdf")

filtered_data_ps <- ps.melt_sum[ps.melt_sum$Sample_Name %in% c("Ctrl24", "EcN24h", "EPEC24h", "SK22D24h", "Cntrl48", 'EcN48h', "EPEC48h", "SK22D48h"), ]
ggplot(filtered_data_ps, aes(x = Sample_Name, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y = "%") +
  facet_wrap(~ Days, scales = "free_x" )+
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 12),              
        strip.text = element_text(size = 14))

ggsave("./abundance_plot/abundance_Class_days_single_strain.pdf", height = 8, width = 12, device = "pdf")


filtered_data_ps1 <- ps.melt_sum[ps.melt_sum$Sample_Name %in% c("Ctrl24", "Cntrl48", "EcNEPEC24h", "SK22DEPEC24h", "Ctrl48","EcNEPEC48h", "SK22DEPEC48h"), ]
ggplot(filtered_data_ps1, aes(x = Sample_Name, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y = "%") +
  facet_wrap(~ Days, scales = "free_x" )+
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 12),              
        strip.text = element_text(size = 14))

ggsave("./abundance_plot/abundance_Class_days_double_strain.pdf", height = 8, width = 12, device = "pdf")



filtered_data_ps2 <- ps.melt_sum[ps.melt_sum$Sample_Name %in% c("Ctrl24", "Cntrl48", "EcNSupernatEPEC24h", "SK22DSupernatEPEC24h","EcNSupernatEPEC48h", "SK22DSupernatEPEC48h"), ]
ggplot(filtered_data_ps2, aes(x = Sample_Name, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y = "%") +
  facet_wrap(~ Days, scales = "free_x" )+
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 12),              
        strip.text = element_text(size = 14))

ggsave("./abundance_plot/abundance_Class_days_supernatent double_strain.pdf", height = 8, width = 12, device = "pdf")


filtered_data_ps3 <- ps.melt_sum[ps.melt_sum$Sample_Name %in% c("Ctrl24", "Cntrl48", "EcNSupernat24h", "SK22DSupernat24h","EcNSupernat48h", "SK22DSupernat48h"), ]
ggplot(filtered_data_ps3, aes(x = Sample_Name, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity", aes(fill = Class)) + 
  labs(x = "", y = "%") +
  facet_wrap(~ Days, scales = "free_x" )+
  theme_classic() + 
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = -90, size = 12),
        axis.text.y = element_text(size = 12),              
        legend.text = element_text(size = 12),              
        strip.text = element_text(size = 14))

ggsave("./abundance_plot/abundance_Class_days_supernatent_single_strain.pdf", height = 8, width = 12, device = "pdf")

#sample_data(ps1)$treatment <- as.factor(sample_data(ps1)$treatment) # factorize for DESeq2
ps.taxa <- tax_glom(ps1, taxrank = 'Species', NArm = FALSE)

sample_data <- sample_data(ps.taxa)
sample_data$treatment <- as.factor(sample_data$Sample_Name)
sample_data$Days <- as.numeric(sample_data$Days)

#Instead of Degs the log2 ratio calculation
# 1 "Ctrl24", "EcN24h


subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "EcN24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)



# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$EcN24h <- otu_table$EcN24h + 1




# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN24h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# store in top_otus_ECN_Ctrl_24 

top_otus_ECN_Ctrl_24 <- merged_table

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECN_Ctrl_24 <- top_otus_ECN_Ctrl_24[!is.na(top_otus_ECN_Ctrl_24$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECN_Ctrl_24  <- top_otus_clean_ECN_Ctrl_24 [abs(top_otus_clean_ECN_Ctrl_24 $Log2_Ratio) >= 5, ]




# Remove OTUs where both Sample1 and Sample2 are less than 10
#filtered_otus_top_otus_clean_ECN_Ctrl_24 <- filtered_otus_top_otus_clean_ECN_Ctrl_24[!(filtered_otus_top_otus_clean_ECN_Ctrl_24$`Ctrl24` < 100 & filtered_otus_top_otus_clean_ECN_Ctrl_24$`EcN24h` < 100), ]


# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_ECN_Ctrl_24 <- filtered_otus_top_otus_clean_ECN_Ctrl_24[!(filtered_otus_top_otus_clean_ECN_Ctrl_24$`Total_Abundance` < 100) , ]


# Define the getPalette function
colourCount <- 50  # Example number of colors you need
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))

# Plotting the log2 ratios for the top OTUs by Class
ggplot(filtered_otus_top_otus_clean_ECN_Ctrl_24, aes(x = Family, y = Log2_Ratio, fill = Family)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcN24h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between EcN24h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )
write.csv(filtered_otus_top_otus_clean_ECN_Ctrl_24, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_ECN_Ctrl_24.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/ECN24_Ctrl24_family.pdf", height = 8, width = 15, device = "pdf")





#2 Ctrl24 vs Epec24h

subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "EPEC24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$EPEC24h <- otu_table$EPEC24h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EPEC24h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)



# Order by total abundance and select top N most abundant OTUs

top_otus_EPEC_Ctrl_24 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EPEC_Ctrl_24 <- top_otus_EPEC_Ctrl_24[!is.na(top_otus_EPEC_Ctrl_24$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EPEC_Ctrl_24  <- top_otus_clean_EPEC_Ctrl_24 [abs(top_otus_clean_EPEC_Ctrl_24 $Log2_Ratio) >= 5, ]

#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_EPEC_Ctrl_24 <- filtered_otus_top_otus_clean_EPEC_Ctrl_24[!(filtered_otus_top_otus_clean_EPEC_Ctrl_24$`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EPEC_Ctrl_24, aes(x = Family, y = Log2_Ratio, fill = Family)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EPEC24h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between EPEC24h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
      limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
      breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
      expand = c(0, 0)           # Remove extra space at the axis edges
  )
write.csv(filtered_otus_top_otus_clean_EPEC_Ctrl_24, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EPEC_Ctrl_24.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/EPEC24h_Ctrl24_family.pdf", height = 8, width = 15, device = "pdf")






# 3 Ctrl24 vs EcNEPEC24h


subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "EcNEPEC24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$EcNEPEC24h <- otu_table$EcNEPEC24h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcNEPEC24h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EcNEPEC_ctrl_24h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]


# Remove rows with NA in Log2_Ratio
top_otus_clean_EcNEPEC_ctrl_24h <- top_otus_EcNEPEC_ctrl_24h[!is.na(top_otus_EcNEPEC_ctrl_24h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h  <- top_otus_clean_EcNEPEC_ctrl_24h [abs(top_otus_clean_EcNEPEC_ctrl_24h $Log2_Ratio) >= 5, ]

#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h <- filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h[!(filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h$`Total_Abundance` < 100) , ]


# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcNEPEC24h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between EcNEPEC24h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
        ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EcNEPEC_ctrl_24h.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/EcNEPEC24h_Ctrl24.pdf", height = 6, width = 10, device = "pdf")


#




# 4 Ctrl24 vs SK22DEPEC24h



subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "SK22DEPEC24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$SK22DEPEC24h <- otu_table$SK22DEPEC24h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22DEPEC24h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22DEPEC_ctrl_24h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22DEPEC_ctrl_24h <- top_otus_SK22DEPEC_ctrl_24h[!is.na(top_otus_SK22DEPEC_ctrl_24h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h  <- top_otus_clean_SK22DEPEC_ctrl_24h [abs(top_otus_clean_SK22DEPEC_ctrl_24h $Log2_Ratio) >= 5, ]

#apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h <- filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h[!(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h$`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22DEPEC24h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22DEPEC24h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/SK22DEPEC24h_Ctrl24.pdf", height = 6, width = 10, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_SK22DEPEC_ctrl_24h.csv", row.names = FALSE, quote = FALSE)


# 5 Ctrl24 vs SK22D24h


subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "SK22D24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$SK22D24h <- otu_table$SK22D24h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22D24h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22D_ctrl_24h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22D_ctrl_24h <- top_otus_SK22D_ctrl_24h[!is.na(top_otus_SK22D_ctrl_24h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22D_ctrl_24h  <- top_otus_clean_SK22D_ctrl_24h [abs(top_otus_clean_SK22D_ctrl_24h $Log2_Ratio) >= 5, ]

#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_SK22D_ctrl_24h <- filtered_otus_top_otus_clean_SK22D_ctrl_24h[!(filtered_otus_top_otus_clean_SK22D_ctrl_24h$`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22D_ctrl_24h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22D24h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22D24h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/SK22D24h_Ctrl24.pdf", height = 6, width = 10, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_SK22D_ctrl_24h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_SK22D_ctrl_24h.csv", row.names = FALSE, quote = FALSE)




# EcNfed224h vs 24h
subset_indices <- sample_data$Sample_Name %in% c("Ctrl24", "EcNfed224h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1
otu_table$EcNfed224h <- otu_table$EcNfed224h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcNfed224h"  # Replace with the actual sample name
sample2 <- "Ctrl24"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EcNfed224h_ctrl_24h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EcNfed224h_ctrl_24h <- top_otus_EcNfed224h_ctrl_24h[!is.na(top_otus_EcNfed224h_ctrl_24h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h  <- top_otus_clean_EcNfed224h_ctrl_24h [abs(top_otus_clean_EcNfed224h_ctrl_24h $Log2_Ratio) >= 5, ]

#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h <- filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h[!(filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h$`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h, aes(x = Family, y = Log2_Ratio, fill = Family)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcNfed224h / Ctrl24)") +
  ggtitle("Top OTUs Log2 Ratio Between EcNfed224h and Ctrl24 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/EcNfed224h_Ctrl24_family.pdf", height = 8, width = 12, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EcNfed224h_ctrl_24h.csv", row.names = FALSE, quote = FALSE)






####48hrs



# 1 "Ctrl48", "EcN48h


subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "EcN48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 90% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$EcN48h <- otu_table$EcN48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN48h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_ECN_Ctrl_48 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECN_Ctrl_48 <- top_otus_ECN_Ctrl_48[!is.na(top_otus_ECN_Ctrl_48$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECN_Ctrl_48  <- top_otus_clean_ECN_Ctrl_48 [abs(top_otus_clean_ECN_Ctrl_48 $Log2_Ratio) >= 5, ]
#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_ECN_Ctrl_48 <- filtered_otus_top_otus_clean_ECN_Ctrl_48[!(filtered_otus_top_otus_clean_ECN_Ctrl_48$`Total_Abundance` < 100) , ]


# Define the getPalette function
colourCount <- 50  # Example number of colors you need
getPalette <- colorRampPalette(brewer.pal(12, "Set3"))

# Plotting the log2 ratios for the top OTUs by Class
ggplot(filtered_otus_top_otus_clean_ECN_Ctrl_48, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcN48h / Ctrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between EcN48h and Ctrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )
write.csv(filtered_otus_top_otus_clean_ECN_Ctrl_48, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_ECN_Ctrl_48.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/ECN48_Ctrl48.pdf", height = 7, width = 11, device = "pdf")







#2 Ctrl48 vs Epec48h
subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "EPEC48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)



# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

#pseudocounts
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$EPEC48h <- otu_table$EPEC48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EPEC48h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)



# Order by total abundance and select top N most abundant OTUs

top_otus_EPEC_Ctrl_48 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EPEC_Ctrl_48 <- top_otus_EPEC_Ctrl_48[!is.na(top_otus_EPEC_Ctrl_48$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EPEC_Ctrl_48  <- top_otus_clean_EPEC_Ctrl_48 [abs(top_otus_clean_EPEC_Ctrl_48 $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_EPEC_Ctrl_48 <- filtered_otus_top_otus_clean_EPEC_Ctrl_48[!(filtered_otus_top_otus_clean_EPEC_Ctrl_48$`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EPEC_Ctrl_48, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EPEC48h / Ctrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between EPEC48h and Ctrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),              
        legend.text = element_text(size = 8),              
        strip.text = element_text(size = 8)) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)           # Remove extra space at the axis edges
  )
write.csv(filtered_otus_top_otus_clean_EPEC_Ctrl_48, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EPEC_Ctrl_48.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/EPEC48h_Ctrl48.pdf", height = 6, width = 10, device = "pdf")






# 3 Ctrl48 vs EcNEPEC48h


subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "EcNEPEC48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

#pseudocount
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$EcNEPEC48h <- otu_table$EcNEPEC48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcNEPEC48h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EcNEPEC_ctrl_48h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]


# Remove rows with NA in Log2_Ratio
top_otus_clean_EcNEPEC_ctrl_48h <- top_otus_EcNEPEC_ctrl_48h[!is.na(top_otus_EcNEPEC_ctrl_48h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h  <- top_otus_clean_EcNEPEC_ctrl_48h [abs(top_otus_clean_EcNEPEC_ctrl_48h $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h <- filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h[!(filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h$`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcNEPEC48h / Ctrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between EcNEPEC48h and Ctrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EcNEPEC_ctrl_48h.csv", row.names = FALSE, quote = FALSE)
ggsave("./abundance_logratio/new_logratio_pseudo/EcNEPEC48h_Ctrl48.pdf", height = 6, width = 10, device = "pdf")



#




# 4 Cntrl48 vs SK22DEPEC48h



subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "SK22DEPEC48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

#pseudocount
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$SK22DEPEC48h <- otu_table$SK22DEPEC48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22DEPEC48h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22DEPEC_ctrl_48h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22DEPEC_ctrl_48h <- top_otus_SK22DEPEC_ctrl_48h[!is.na(top_otus_SK22DEPEC_ctrl_48h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h  <- top_otus_clean_SK22DEPEC_ctrl_48h [abs(top_otus_clean_SK22DEPEC_ctrl_48h $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h <- filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h[!(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h$`Total_Abundance` < 100) , ]


# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22DEPEC48h / Ctrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22DEPEC48h and Ctrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),              
        legend.text = element_text(size = 8),              
        strip.text = element_text(size = 8)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/SK22DEPEC48h_Ctrl48.pdf", height = 6, width = 10, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_SK22DEPEC_ctrl_48h.csv", row.names = FALSE, quote = FALSE)



# 5 Ctrl48 vs SK22D48h


subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "SK22D48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$SK22D48h <- otu_table$SK22D48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22D48h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22D_ctrl_48h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22D_ctrl_48h <- top_otus_SK22D_ctrl_48h[!is.na(top_otus_SK22D_ctrl_48h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22D_ctrl_48h  <- top_otus_clean_SK22D_ctrl_48h [abs(top_otus_clean_SK22D_ctrl_48h $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_SK22D_ctrl_48h <- filtered_otus_top_otus_clean_SK22D_ctrl_48h[!(filtered_otus_top_otus_clean_SK22D_ctrl_48h$`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22D_ctrl_48h, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22D48h / Ctrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22D48h and Ctrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size =10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/SK22D48h_Ctrl48.pdf", height = 6, width = 10, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_SK22D_ctrl_48h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_SK22D_ctrl_48h.csv", row.names = FALSE, quote = FALSE)



# EcNfed224h vs 48h
subset_indices <- (sample_data$Sample_Name %in% c("Cntrl48") & sample_data$Days == 48) |
(sample_data$Sample_Name %in% c("EcNfed224h") & sample_data$Days == 24)

ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$EcNfed224h <- otu_table$EcNfed224h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcNfed224h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EcNfed224h_ctrl_48h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EcNfed224h_ctrl_48h <- top_otus_EcNfed224h_ctrl_48h[!is.na(top_otus_EcNfed224h_ctrl_48h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h  <- top_otus_clean_EcNfed224h_ctrl_48h [abs(top_otus_clean_EcNfed224h_ctrl_48h $Log2_Ratio) >= 5, ]

#Now apply the filter where sample1 and sample2 counts are >= 10
filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h <- filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h[!(filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h$`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h, aes(x = Family, y = Log2_Ratio, fill = Family)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcNfed224h / Cntrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between EcNfed224h and Cntrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/EcNfed224h_Ctrl48_family.pdf", height = 8, width = 12, device = "pdf")

# Save the data frame to a CSV file
write.csv(filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_EcNfed224h_ctrl_48h.csv", row.names = FALSE, quote = FALSE)







# 24hrs vs 48hrs


# 1 Ctrl24 vs Ctrl48


subset_indices <- sample_data$Sample_Name %in% c("Cntrl48", "Ctrl24")
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$Ctrl24 <- otu_table$Ctrl24 + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "Ctrl24"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name


# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_control <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_control <- top_otus_control[!is.na(top_otus_control$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_control  <- top_otus_clean_control [abs(top_otus_clean_control $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_control <- filtered_otus_top_otus_clean_control[!(filtered_otus_top_otus_clean_control$`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_control, aes(x = reorder(Class, -Log2_Ratio), y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (Ctrl24 / Cntrl48)") +
  ggtitle("Top OTUs Log2 Ratio Between Ctrl24 and Cntrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/24_48/ctrl24_Ctrl48.pdf", height = 6, width = 10, device = "pdf")





# 2 ECN24 vs ECN48


subset_indices <- sample_data$Sample_Name %in% c("EcN24h", "EcN48h")
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

#Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EcN24h <- otu_table$EcN24h + 1
otu_table$EcN48h <- otu_table$EcN48h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN24h"  # Replace with the actual sample name
sample2 <- "EcN48h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EcN <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EcN <- top_otus_EcN[!is.na(top_otus_EcN$Log2_Ratio), ]


# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EcN <- top_otus_clean_EcN[abs(top_otus_clean_EcN$Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_EcN  <- filtered_otus_top_otus_clean_EcN [!(filtered_otus_top_otus_clean_EcN $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EcN, aes(x = Genus, y = Log2_Ratio, fill = Genus)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcN24h / EcN48h)") +
  ggtitle("Top OTUs Log2 Ratio Between EcN24h and EcN48h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-50, 100),         # Set y-axis limits (adjust as needed)
    breaks = seq(-50, 100, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/24_48/EcN24_EcN48_factor_5_genus.pdf", height = 8, width = 10, device = "pdf")



# 3 EPEC24h vs EPEC48h



subset_indices <- sample_data$Sample_Name %in% c("EPEC24h", "EPEC48h")
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)


# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EPEC24h <- otu_table$EPEC24h + 1
otu_table$EPEC48h <- otu_table$EPEC48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EPEC24h"  # Replace with the actual sample name
sample2 <- "EPEC48h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_EPEC <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_EPEC <- top_otus_EPEC[!is.na(top_otus_EPEC$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_EPEC   <- top_otus_clean_EPEC  [abs(top_otus_clean_EPEC $Log2_Ratio) >= 5, ]


# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_EPEC <- filtered_otus_top_otus_clean_EPEC [!(filtered_otus_top_otus_clean_EPEC $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_EPEC , aes(x = reorder(Genus, -Log2_Ratio), y = Log2_Ratio, fill = Genus)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EPEC24h / EPEC48h)") +
  ggtitle("Top OTUs Log2 Ratio Between EPEC24h and EPEC48h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 200, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/24_48/EPEC24_EPEC48_Genus.pdf", height = 9, width = 10, device = "pdf")

#4  SK22D24h vs SK22D48h



subset_indices <- sample_data$Sample_Name %in% c("SK22D24h", "SK22D48h")
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$SK22D24h <- otu_table$SK22D24h + 1
otu_table$SK22D48h <- otu_table$SK22D48h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22D24h"  # Replace with the actual sample name
sample2 <- "SK22D48h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22D <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22D <- top_otus_SK22D[!is.na(top_otus_SK22D$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22D <- top_otus_clean_SK22D [abs(top_otus_clean_SK22D $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_SK22D <- filtered_otus_top_otus_clean_SK22D [!(filtered_otus_top_otus_clean_SK22D $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22D, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22D24h / SK22D48h)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22D24h and SK22D48h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/24_48/SK22D24_SK22D48_class.pdf", height = 6, width = 10, device = "pdf")


#5 SK22DEPEC24h vs SK22DEPEC48h


subset_indices <- sample_data$Sample_Name %in% c("SK22DEPEC24h", "SK22DEPEC48h")
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$SK22DEPEC48h <- otu_table$SK22DEPEC48h + 1
otu_table$SK22DEPEC24h <- otu_table$SK22DEPEC24h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22DEPEC24h"  # Replace with the actual sample name
sample2 <- "SK22DEPEC48h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22DEPEC <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22DEPEC <- top_otus_SK22DEPEC[!is.na(top_otus_SK22DEPEC$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22DEPEC  <- top_otus_clean_SK22DEPEC [abs(top_otus_clean_SK22DEPEC $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_SK22DEPEC <- filtered_otus_top_otus_clean_SK22DEPEC [!(filtered_otus_top_otus_clean_SK22DEPEC $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22DEPEC, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (SK22DEPEC24h / SK22DEPEC48h)") +
  ggtitle("Top OTUs Log2 Ratio Between SK22DEPEC24h and SK22DEPEC48h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/24_48/SK22DEPEC24_SK22DEPEC48_class.pdf", height = 6, width = 10, device = "pdf")




# 3 EcN24h vs EPEC24h



subset_indices <- sample_data$Sample_Name %in% c("EcN24h", "EPEC24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EcN24h <- otu_table$EcN24h + 1
otu_table$EPEC24h <- otu_table$EPEC24h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN24h"  # Replace with the actual sample name
sample2 <- "EPEC24h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_ECN_EPEC24 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECN_EPEC24 <- top_otus_ECN_EPEC24[!is.na(top_otus_ECN_EPEC24$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECN_EPEC24  <- top_otus_clean_ECN_EPEC24 [abs(top_otus_clean_ECN_EPEC24 $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_ECN_EPEC24<- filtered_otus_top_otus_clean_ECN_EPEC24 [!(filtered_otus_top_otus_clean_ECN_EPEC24 $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_ECN_EPEC24, aes(x =Genus, y = Log2_Ratio, fill = Genus)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio ( ECN24h/ EPEC24h )") +
  ggtitle("Top OTUs Log2 Ratio Between ECN24h and EPEC24h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/ECN24_EPEC24_genus.pdf", height = 6, width = 10, device = "pdf")




# 3 EPEC24h vs EcNEPEC24h



subset_indices <- sample_data$Sample_Name %in% c("EPEC24h", "EcNEPEC24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EPEC24h <- otu_table$EPEC24h + 1
otu_table$EcNEPEC24h <- otu_table$EcNEPEC24h + 1

# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcNEPEC24h"  # Replace with the actual sample name
sample2 <- "EPEC24h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_ECNEPEC24_EPEC24 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECNEPEC24_EPEC24<- top_otus_ECNEPEC24_EPEC24[!is.na(top_otus_ECNEPEC24_EPEC24$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECNEPEC24_EPEC24  <- top_otus_clean_ECNEPEC24_EPEC24 [abs(top_otus_clean_ECNEPEC24_EPEC24 $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_ECNEPEC24_EPEC24 <- filtered_otus_top_otus_clean_ECNEPEC24_EPEC24 [!(filtered_otus_top_otus_clean_ECNEPEC24_EPEC24 $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_ECNEPEC24_EPEC24 , aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (EcNEPEC24h / EPEC24h)") +
  ggtitle("Top OTUs Log2 Ratio Between EcNEPEC24h and EPEC24h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/ECNEPEC24_EPEC24_.pdf", height = 6, width = 10, device = "pdf")

# ECN vs  SKD22 24hrs


subset_indices <- sample_data$Sample_Name %in% c("EcN24h", "SK22D24h") & sample_data$Days == 24
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EcN24h <- otu_table$EcN24h + 1
otu_table$SK22D24h <- otu_table$SK22D24h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN24h"  # Replace with the actual sample name
sample2 <- "SK22D24h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_ECN_SK22D24h <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECN_SK22D24h <- top_otus_ECN_SK22D24h[!is.na(top_otus_ECN_SK22D24h$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECN_SK22D24h  <- top_otus_clean_ECN_SK22D24h [abs(top_otus_clean_ECN_SK22D24h$Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_ECN_SK22D24h<- filtered_otus_top_otus_clean_ECN_SK22D24h [!(filtered_otus_top_otus_clean_ECN_SK22D24h $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_ECN_SK22D24h, aes(x = Genus, y = Log2_Ratio, fill = Genus)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio ( ECN24h/ SK22D24h )") +
  ggtitle("Top OTUs Log2 Ratio Between ECN24h and SK22D24h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/ECN24_SK22D24h.pdf", height = 6, width = 10, device = "pdf")
write.csv(filtered_otus_top_otus_clean_ECN_SK22D24h, file = "./abundance_logratio/new_logratio_pseudo/filtered_otus_top_otus_clean_ECN_SK22D24h.csv", row.names = FALSE, quote = FALSE)


# Sk224h vs ctrl 48 hr


subset_indices <- sample_data$Sample_Name %in% c("SK22D24h") & sample_data$Days == 24  |
  sample_data$Sample_Name %in% c("Cntrl48") & sample_data$Days == 48

ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$Cntrl48 <- otu_table$Cntrl48 + 1
otu_table$SK22D24h <- otu_table$SK22D24h + 1



# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "SK22D24h"  # Replace with the actual sample name
sample2 <- "Cntrl48"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_SK22D24h_ctrl48 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_SK22D24h_ctrl48 <- top_otus_SK22D24h_ctrl48[!is.na(top_otus_SK22D24h_ctrl48$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_SK22D24h_ctrl48  <- top_otus_clean_SK22D24h_ctrl48 [abs(top_otus_clean_SK22D24h_ctrl48$Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_SK22D24h_ctrl48<- filtered_otus_top_otus_clean_SK22D24h_ctrl48 [!(filtered_otus_top_otus_clean_SK22D24h_ctrl48 $`Total_Abundance` < 100) , ]

# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_SK22D24h_ctrl48, aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (  SK22D24h/Ctrl48 )") +
  ggtitle("Top OTUs Log2 Ratio Between SK22D24h/Cntrl48 ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-100, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-100, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )

ggsave("./abundance_logratio/new_logratio_pseudo/SK22D24h_ctrl48hrs.pdf", height = 6, width = 10, device = "pdf")




# 3 EcN48h vs EPEC48h



subset_indices <- sample_data$Sample_Name %in% c("EcN48h", "EPEC48h") & sample_data$Days == 48
ps.taxa.sub <- prune_samples(subset_indices, ps.taxa)

# filter sparse features, with > 80% zeros
ps.taxa.pse.sub <- prune_taxa(rowSums(otu_table(ps.taxa.sub) == 0) < ncol(otu_table(ps.taxa.sub)) * 0.9, ps.taxa.sub)

# Extract the OTU table from the phyloseq object
otu_table <- as.data.frame(otu_table(ps.taxa.pse.sub))

# Add a pseudocount of 1 to all OTU counts in sample1 and sample2 to avoid division by zero
otu_table <- otu_table
otu_table$EcN48h <- otu_table$EcN48h + 1
otu_table$EPEC48h <- otu_table$EPEC48h + 1


# Ensure the OTU IDs are preserved
otu_table$OTU <- rownames(otu_table)

# Assuming 'otu_table' is already in the correct format
# Define sample names
sample1 <- "EcN48h"  # Replace with the actual sample name
sample2 <- "EPEC48h"  # Replace with the actual sample name

# Calculate the ratio (Sample1 / Sample2)
otu_table$Ratio <- otu_table[, sample1] / otu_table[, sample2]

# Log-transform the ratio
otu_table$Log2_Ratio <- log2(otu_table$Ratio)

# Handle potential division by zero by replacing Inf and NaN with NA
otu_table$Log2_Ratio[is.infinite(otu_table$Log2_Ratio) | is.nan(otu_table$Log2_Ratio)] <- NA

# Calculate total abundance across both samples for each OTU
otu_table$Total_Abundance <- rowSums(otu_table[, c(sample1, sample2)], na.rm = TRUE)

otu_table$OTU <- rownames(otu_table)

#  'otu_table' and 'tax_table' are already in the correct format
# Ensure the OTU IDs are preserved in tax_table as rownames

tax_table <- as.data.frame(tax_table(ps.taxa.pse.sub))
tax_table$OTU <- rownames(tax_table)

# Merge the otu_table with tax_table based on OTU ID
merged_table <- merge(otu_table, tax_table, by.x = "OTU", by.y = "row.names", all.x = TRUE)


# Order by total abundance and select top N most abundant OTUs

top_otus_ECN_EPEC48 <- merged_table[order(merged_table$Total_Abundance, decreasing = TRUE), ]

# Remove rows with NA in Log2_Ratio
top_otus_clean_ECN_EPEC48 <- top_otus_ECN_EPEC48[!is.na(top_otus_ECN_EPEC48$Log2_Ratio), ]

# Filter OTUs where the absolute Log2_Ratio is greater than 1 (i.e., difference is more than a factor of 2)
filtered_otus_top_otus_clean_ECN_EPEC48   <- top_otus_clean_ECN_EPEC48  [abs(top_otus_clean_ECN_EPEC48  $Log2_Ratio) >= 5, ]

# Remove OTUs where OTU in total_abundance less than 100

filtered_otus_top_otus_clean_ECN_EPEC48 <- filtered_otus_top_otus_clean_ECN_EPEC48 [!(filtered_otus_top_otus_clean_ECN_EPEC48 $`Total_Abundance` < 100) , ]



# Plotting the log2 ratios for the top OTUs 
ggplot(filtered_otus_top_otus_clean_ECN_EPEC48  , aes(x = Class, y = Log2_Ratio, fill = Class)) +
  geom_bar(stat = "identity") +
  xlab("OTU") +
  ylab("Log2 Ratio (ECN48h / EPEC48h  )") +
  ggtitle("Top OTUs Log2 Ratio Between ECN48h and EPEC48h ") +
  theme_classic() +
  theme(legend.position = "right", 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle = 70, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),              
        legend.text = element_text(size = 10),              
        strip.text = element_text(size = 10)
  ) +
  scale_fill_brewer(palette = "Set3")+
  scale_fill_manual(values = getPalette(colourCount))+
  scale_y_continuous(
    limits = c(-130, 20),         # Set y-axis limits (adjust as needed)
    breaks = seq(-130, 20, by = 10),  # Set breaks for y-axis ticks
    expand = c(0, 0)            # Remove extra space at the axis edges
  )
ggsave("./abundance_logratio/new_logratio_pseudo/ECN48_EPEC48.pdf", height = 6, width = 10, device = "pdf")


