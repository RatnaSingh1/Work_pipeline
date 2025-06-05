library("data.table")   # Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("tidyverse")
library("DECIPHER")
library("phangorn")
library("ggplot2")
#loading folder
picrust2 = "picrust2_out_pipeline2"
list.files(picrust2, recursive = TRUE)

#building output file

p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")


mapfile = "picrust2_out_pipeline1/description_mapfiles"
list.files(mapfile, recursive = TRUE)

mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

#Load map files
mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")

#Load picrust2 output files
p2EC = as.data.frame(fread(p2_EC))
# Example first row containing column names
column_names <- c("function", "X361", "X362", "X363", "X364", "X365", "X366", "X367", "X368", "X369", "X370", "X371", "X372", "X373", "X374", "X375", "X376", "X377", "X378", "X379", "X380", "X381")
# Assign the first row as column names
colnames(p2EC) <- column_names
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[-1,-1])
p2EC = round(p2EC)


p2KO = as.data.frame(fread(p2_KO))
# Example first row containing column names
column_names <- c("function", "X361", "X362", "X363", "X364", "X365", "X366", "X367", "X368", "X369", "X370", "X371", "X372", "X373", "X374", "X375", "X376", "X377", "X378", "X379", "X380", "X381")
# Assign the first row as column names
colnames(p2KO) <- column_names
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[-1,-1])
p2KO = round(p2KO)


p2PW = as.data.frame(fread(p2_PW))
column_names <- c("pathway", "X361", "X362", "X363", "X364", "X365", "X366", "X367", "X368", "X369", "X370", "X371", "X372", "X373", "X374", "X375", "X376", "X377", "X378", "X379", "X380", "X381")
# Assign the first row as column names
colnames(p2PW) <- column_names

rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[-1,-1])
p2PW = round(p2PW)



p2EC1 <- (p2EC[, c("X367", "X368", "X369","X364", "X365", "X366")])
p2KO1 <- (p2KO[, c("X367", "X368", "X369","X364", "X365", "X366")])
p2PW1 <- (p2PW[, c("X367", "X368", "X369","X364", "X365", "X366")])


#pca analysis
# Assuming p2EC1 is the subset of the original dataset as mentioned
library(ggfortify)
install.packages("factoextra")
install.packages("ggplot2")
# Perform PCA analysis
pca_result <- prcomp(p2EC1, scale = TRUE)

# Access the results
summary(pca_result) # Provides summary of the PCA results
pca_result$rotation # Principal component loadings (eigenvectors)
pca_result$x # Transformed data (scores)

autoplot(pca_result)
ggsave("images/metabolite_pca_EC1_2.png", width = 5, height = 4)
# Define vectors for treated and untreated samples


# Assuming pca_result is a list containing the rotation matrix

# Create a data frame with the samples X364 to X369 and their loadings on PC1 and PC2
samples <- c("X364", "X365", "X366", "X367", "X368", "X369")
sample_indices <- match(samples, rownames(pca_result$rotation))

samples_data <- data.frame(
  SampleID = samples,
  PC1 = pca_result$rotation[sample_indices, "PC1"],
  PC2 = pca_result$rotation[sample_indices, "PC2"]
)

# Create a new column in samples_data to differentiate treated and untreated samples
samples_data$Treatment <- ifelse(samples_data$SampleID %in% c("X367", "X368", "X369"), "Treated", "Untreated")

# Plot the samples X364 to X369 from PC1 and PC2 based on their loadings
ggplot(samples_data, aes(x = PC1, y = PC2, color = Treatment)) +
  geom_point(size = 3) +
  labs(x = "PC1", y = "PC2") +
  ggtitle("Samples X364 to X369 in PC1 and PC2") +
  theme_minimal() +
  geom_text(aes(label = SampleID), vjust = -0.5) +
  scale_color_manual(values = c("green", "blue"), labels = c("Treated", "Untreated")) +
  guides(color = guide_legend(title = "Sample Group"))



ggsave("images/metabolite_pca_EC1.png", width = 5, height = 4)




library(ggplot2)

# Assuming pca_result is already calculated

# Convert PCA results to a data frame
pca_data <- as.data.frame(pca_result$x)
colnames(pca_data) <- c("PC1", "PC2")  # Rename the columns for clarity

# Assuming p2EC contains the original dataset with sample metadata
# Add sample labels and color based on "Treatment" column (Treated vs. Untreated)
pca_data$SampleID <- rownames(p2EC1)
pca_data$Treatment <- ifelse(pca_data$SampleID %in% c("X367", "X368", "X369"), "Treated", "Untreated")

# Plot PCA results using ggplot2
ggplot() +
  geom_point(data = pca_data[pca_data$Treatment == "Treated", ], aes(x = PC1, y = PC2, color = Treatment), size = 3) +
  geom_point(data = pca_data[pca_data$Treatment == "Untreated", ], aes(x = PC1, y = PC2, color = Treatment), size = 3) +
  labs(x = "PC1", y = "PC2") +
  ggtitle("Enzymes") +
  scale_color_manual(values = c("green", "blue"))

ggsave("images/metabolite.png", width = 5, height = 4)



