#install.packages("webshot")
#if ( !require(RSelenium) ) {
  install.packages("kaleido", repos = "https://cloud.r-project.org/")
#}

#p %>%
#  export(file = "filename.svg",
#         selenium = RSelenium::rsDriver(browser = "chrome"))
library(RSelenium)
library(webshot)
  library(data.table)
  library(plotly)
  library(ggplot)
  library(dplyr)
#install.packages("kaleido")
# Step 1: Data Preparation
# Assuming you have a data frame named 'abundance_data' with columns: 'MetaCyc_ID', 'Treated_Replicate1', 'Treated_Replicate2', 'Untreated_Replicate1', 'Untreated_Replicate2'
# Make sure the replicate columns contain abundance values

# Step 1: Data Preparation
abundance_data <- read.table('path_abun_unstrat.tsv', sep = '\t', header = TRUE)


mapfile = "description_mapfiles"
list.files(mapfile, recursive = TRUE)

mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")


#Load map files
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")
mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapKO1 <- mapKO %>%
  rename(function1 = "function")





  # Make sure abundance_data is a data frame
  abundance_data_df <- as.data.frame(abundance_data)
  
  
  # #Merge with map file data, Perform the inner join and arrange the data
  abundance_data <- abundance_data_df %>%
    tibble::rownames_to_column(var = "Pathway") %>%
    inner_join(mapPW, by = c("pathway" = "pathway")) %>%
    arrange(Pathway)
 # Print the first few rows of abundance_data1 to check the content
  print(head(abundance_data))
  

#Merge with map file data



# Step 2: Data Analysis (Normalization using Z-score)
# Aggregate the replicates
abundance_data$Treated <- rowMeans(abundance_data[, c("X367", "X368", "X369")], na.rm = TRUE)
abundance_data$Untreated <- rowMeans(abundance_data[, c("X364", "X365", "X366")], na.rm = TRUE)

# Calculate Z-scores for Treated and Untreated
abundance_data$Treated_zscore <- scale(abundance_data$Treated)
abundance_data$Untreated_zscore <- scale(abundance_data$Untreated)

# Step 3: Filter Rows Based on Condition
# Discard rows with treated z-score values below 100
#abundance_data_filtered <- abundance_data[abundance_data$Treated_zscore >= 100, ]
abundance_data_filtered <- abundance_data
treated_untreated_matrix <- abundance_data_filtered[, c("pathway", "Treated_zscore", "description", "Untreated_zscore")]

#Step 4: Visualization
# Load the required packages
#with plotly

write.csv(treated_untreated_matrix, "pathway/treated_untreated_matrix_pathway.csv", row.names = TRUE)

top20_data <- head(treated_untreated_matrix[order(abundance_data_filtered$Treated_zscore, decreasing = TRUE), ], 30)
write.csv(top20_data, "pathway/top20_data_zcore.csv")

plot_obj <- plot_ly(data = top20_data, x = ~pathway) %>%
  add_trace(y = ~Treated_zscore, name = "Treated", type = "bar", marker = list(color = "blue")) %>%
  add_trace(y = ~Untreated_zscore, name = "Untreated", type = "bar", marker = list(color = "red")) %>%
  layout(xaxis = list(title = "Pathway"), yaxis = list(title = "Normalized Abundance (Z-score)"), title = "Normalized Abundance in Treated Samples (Filtered)")

plotly::export(plot_obj, #the graph to export
               file = "pathway/treated_untreated_matrix_pathway.png") #the name and type of file (can be .png, .jpeg, etc.)




# for EC

abundance_data_EC <- read.table('pred_metagenome_unstrat.tsv', sep = '\t', header = TRUE)

abundance_data_EC1 <- abundance_data_EC %>%
  rename("function" = function1)

# Make sure abundance_data is a data frame
abundance_data_EC_df <- as.data.frame(abundance_data_EC1)


# #Merge with map file data, Perform the inner join and arrange the data
abundance_data_EC_df <- abundance_data_EC_df %>%
  tibble::rownames_to_column(var = "Enzyme") %>%
  inner_join(mapEC, by = c("function" = "function")) %>%
  arrange(Enzyme)
# Print the first few rows of abundance_data1 to check the content
print(head(abundance_data_EC_df))

# Step 2: Data Analysis (Normalization using Z-score)
# Aggregate the replicates
abundance_data_EC_df$Treated <- rowMeans(abundance_data_EC_df[, c("X367", "X368", "X369")], na.rm = TRUE)
abundance_data_EC_df$Untreated <- rowMeans(abundance_data_EC_df[, c("X364", "X365", "X366")], na.rm = TRUE)

# Calculate Z-scores for Treated and Untreated
abundance_data_EC_df$Treated_zscore <- scale(abundance_data_EC_df$Treated)
abundance_data_EC_df$Untreated_zscore <- scale(abundance_data_EC_df$Untreated)

# Step 3: Filter Rows Based on Condition
# Discard rows with treated z-score values below 100
#abundance_data_filtered <- abundance_data[abundance_data$Treated_zscore >= 100, ]
abundance_data_EC_df_filtered <- abundance_data_EC_df
treated_untreated_matrix_EC <- abundance_data_EC_df[, c("function", "Treated_zscore", "description", "Untreated_zscore")]

#Step 4: Visualization
# Load the required packages
#with plotly

write.csv(treated_untreated_matrix_EC, "enzyme/treated_untreated_matrix_pathway.csv", row.names = TRUE)

top20_data_EC <- head(treated_untreated_matrix_EC[order(treated_untreated_matrix_EC$Treated_zscore, decreasing = TRUE), ], 30)

write.csv(top20_data_EC, "enzyme/top20_data_zcore_enzyme.csv")

plot_obj_EC <- plot_ly(data = top20_data_EC, x = ~description) %>%
  add_trace(y = ~Treated_zscore, name = "Treated", type = "bar", marker = list(color = "blue")) %>%
  add_trace(y = ~Untreated_zscore, name = "Untreated", type = "bar", marker = list(color = "red")) %>%
  layout(xaxis = list(title = "Enzyme"), yaxis = list(title = "Normalized Abundance (Z-score)"), title = "Normalized Abundance in Treated Samples (Filtered)")

plotly::export(plot_obj_EC, #the graph to export
               file = "enzyme/treated_untreated_matrix_enzyme.png") #the name and type of file (can be .png, .jpeg, etc.)




#for KO


abundance_data_KO <- read.table('pred_metagenome_unstrat_KO.tsv', sep = '\t', header = TRUE)

abundance_data_KO1 <- abundance_data_KO %>%
  rename("function1" = "function.")

# Make sure abundance_data is a data frame
abundance_data_KO_df <- as.data.frame(abundance_data_KO1)


# #Merge with map file data, Perform the inner join and arrange the data
abundance_data_KO_df <- abundance_data_KO_df %>%
  tibble::rownames_to_column(var = "KO") %>%
  inner_join(mapKO1, by = c("function1" = "function1")) %>%
  arrange(KO)
# Print the first few rows of abundance_data1 to check the content
print(head(abundance_data_KO_df))

# Step 2: Data Analysis (Normalization using Z-score)
# Aggregate the replicates
abundance_data_KO_df$Treated <- rowMeans(abundance_data_KO_df[, c("X367", "X368", "X369")], na.rm = TRUE)
abundance_data_KO_df$Untreated <- rowMeans(abundance_data_KO_df[, c("X364", "X365", "X366")], na.rm = TRUE)

# Calculate Z-scores for Treated and Untreated
abundance_data_KO_df$Treated_zscore <- scale(abundance_data_KO_df$Treated)
abundance_data_KO_df$Untreated_zscore <- scale(abundance_data_KO_df$Untreated)

# Step 3: Filter Rows Based on Condition
# Discard rows with treated z-score values below 100
#abundance_data_filtered <- abundance_data[abundance_data$Treated_zscore >= 100, ]
abundance_data_KO_df_filtered <- abundance_data_KO_df
treated_untreated_matrix_KO <- abundance_data_KO_df[, c("function1", "description", "Treated_zscore",  "Untreated_zscore")]

#Step 4: Visualization
# Load the required packages
#with plotly

write.csv(treated_untreated_matrix_KO, "KO/treated_untreated_matrix_KO.csv", row.names = TRUE)

top20_data_KO <- head(treated_untreated_matrix_KO[order(treated_untreated_matrix_KO$Treated_zscore, decreasing = TRUE), ], 30)

write.csv(top20_data_KO, "KO/top20_data_zcore_KO.csv")

plot_obj_KO <- plot_ly(data = top20_data_KO, x = ~description) %>%
  add_trace(y = ~Treated_zscore, name = "Treated", type = "bar", marker = list(color = "blue")) %>%
  add_trace(y = ~Untreated_zscore, name = "Untreated", type = "bar", marker = list(color = "red")) %>%
  layout(xaxis = list(title = "Enzyme"), yaxis = list(title = "Normalized Abundance (Z-score)"), title = "Normalized Abundance in Treated Samples (Filtered)")

plotly::export(plot_obj_KO, #the graph to export
               file = "KO/treated_untreated_matrix_KO.png") #the name and type of file (can be .png, .jpeg, etc.)



#for whole KO


#for pathway


ggplot(treated_untreated_matrix, aes(x = description)) +
  geom_bar(aes(y = Treated_zscore, fill = "Treated"), stat = "identity", width = 0.5) +
  geom_bar(aes(y = Untreated_zscore, fill = "Untreated"), stat = "identity", width = 0.5) +
  labs(x = "Pathway", y = "Normalized Abundance (log2)", title = "Normalized Abundance in Treated and Untreated Samples") +
  scale_fill_manual(values = c("Treated" = "blue", "Untreated" = "red"), labels = c("Treated", "Untreated")) +
  theme_bw()
ggsave("pathway/pathway_whole.png", height = 3, width = 8, device = "png")

#for whole KO

ggplot(treated_untreated_matrix_KO, aes(x = description)) +
  geom_bar(aes(y = Treated_zscore, fill = "Treated"), stat = "identity", width = 0.5) +
  geom_bar(aes(y = Untreated_zscore, fill = "Untreated"), stat = "identity", width = 0.5) +
  labs(x = "KO ID", y = "Normalized Abundance (log2)", title = "Normalized Abundance in Treated and Untreated Samples") +
  scale_fill_manual(values = c("Treated" = "blue", "Untreated" = "red"), labels = c("Treated", "Untreated")) +
  theme_bw()
ggsave("KO/KO_whole.png", height = 3, width = 8, device = "png")



#for whole enzyme
treated_untreated_matrix_EC
ggplot(treated_untreated_matrix_EC, aes(x = description)) +
  geom_bar(aes(y = Treated_zscore, fill = "Treated"), stat = "identity", width = 0.5) +
  geom_bar(aes(y = Untreated_zscore, fill = "Untreated"), stat = "identity", width = 0.5) +
  labs(x = "Enzymes", y = "Normalized Abundance (log2)", title = "Normalized Abundance in Treated and Untreated Samples") +
  scale_fill_manual(values = c("Treated" = "blue", "Untreated" = "red"), labels = c("Treated", "Untreated")) +
  theme_bw()
ggsave("enzyme/enzyme_whole.png", height = 3, width = 8, device = "png")























# with ggplot2

#  top20_data, Treated_zscore, and Untreated_zscore

# Combine the data into a single data frame
# Merge top20_data with the original data frame to get the description column back

# Create a data frame with pathway, Treated_zscore, and Untreated_zscore columns
plot_data <- data.frame(
  pathway = top20_data$pathway,
  description = top20_data$description,
  Treated_zscore = top20_data$Treated_zscore,
  Untreated_zscore = top20_data$Untreated_zscore
)

# Melt the data for plotting in ggplot
plot_data_melted <- reshape2::melt(plot_data, id.vars = c("pathway", "description"), variable.name = "Group", value.name = "Z-score")


# Plot using ggplot2
ggplot(plot_data_melted, aes(x = description, y = `Z-score`, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pathway", y = "Normalized Abundance (Z-score)",
       title = "Normalized Abundance in Treated and Untreated Samples (Filtered)") +
  scale_fill_manual(values = c("blue", "red"), labels = c("Treated", "Untreated")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank())  # Remove grid lines

# Save the ggplot2 object using ggsave
ggsave("pathway/treated_untreated_matrix_pathway2.png", width = 10, height = 6)























# Step 4: Visualization for a Specific Gene
# Select a specific gene (change 'your_gene_name' to the desired gene)
your_gene_name <- "ASPASN-PWY"

# Subset the data for the specific gene (replace "your_gene_name" with the actual gene name)
gene_data <- subset(plot_data_melted, pathway == your_gene_name)

# Subset the data for the specific gene (replace "your_gene_name" with the actual gene name)
gene_data <- subset(plot_data_melted, pathway == your_gene_name)

# Create a vertical bar plot for the specific gene with all replicate values
ggplot(gene_data, aes(x = description, y = `Z-score`, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Pathway", y = "Normalized Abundance", title = your_gene_name) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Set x-axis labels to be horizontal and centered
        plot.title = element_text(angle = 0))  # Keep title as it is (no rotation)
ggsave("pathway/gene.png", width = 10, height = 6)

