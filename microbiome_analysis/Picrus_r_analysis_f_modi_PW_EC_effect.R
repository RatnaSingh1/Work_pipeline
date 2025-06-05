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


mapfile = "picrust2_out_pipeline2/description_mapfiles"
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
column_names <- c("function", "X364", "X365", "X366", "X367", "X368", "X369")
# Assign the first row as column names
colnames(p2EC) <- column_names
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[-1,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
# Example first row containing column names
column_names <- c("function", "X364", "X365", "X366", "X367", "X368", "X369")
# Assign the first row as column names
colnames(p2KO) <- column_names
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[-1,-1])
p2KO = round(p2KO)


p2PW = as.data.frame(fread(p2_PW))
column_names <- c("pathway", "X364", "X365", "X366", "X367", "X368", "X369")
# Assign the first row as column names
colnames(p2PW) <- column_names
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[-1,-1])
p2PW = round(p2PW)



p2EC1 <- (p2EC[, c("X367", "X368", "X369","X364", "X365", "X366")])
p2KO1 <- (p2KO[, c("X367", "X368", "X369","X364", "X365", "X366")])
p2PW1 <- (p2PW[, c("X367", "X368", "X369","X364", "X365", "X366")])

treated <- (p2EC1[, c("X367", "X368", "X369")])
untreated <- (p2EC1[, c("X364", "X365", "X366")])             
# Create a factor vector representing the conditions
conditions <- c(rep("treated", 3), rep("untreated", 3))

# Perform the ALDEx2 analysis for EC
set.seed(12345)
system.time({
  aldex2_EC1 = aldex(p2EC1, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", 
                     conditions = conditions, verbose = TRUE)
})



# pathway
p2PW1 <- (p2PW[, c("X367", "X368", "X369","X364", "X365", "X366")])

treated <- (p2PW1[, c("X367", "X368", "X369")])
untreated <- (p2PW1[, c("X364", "X365", "X366")])
conditions <- c(rep("treated", 3), rep("untreated", 3))
set.seed(12345)
system.time({
  aldex2_PW1 = aldex(p2PW1, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", 
                     conditions = conditions, verbose = TRUE)
})
head(aldex2_PW1, 10)

#ALDEx2 authors suggest that an effect size of 1 or greater can be used as significance cutoff



png("images/ALDEx2_picrust2_effect_1_treated_vs_untreated2.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,1))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())

# Assuming filtered_data_PW1 is your filtered data frame containing pathways with effect size > 1

# Load the required library
library(ggplot2)



# Filter the data to retain only rows with effect size greater than 1
filtered_data_PW1 <- aldex2_PW1[aldex2_PW1$effect > 1, c( "effect", "we.ep")]

# Write the filtered data to a CSV file
write.csv(filtered_data_PW1, "images/data/filtered_data_PW1.csv", row.names = TRUE)

# Filter the data to retain only rows with effect size greater than 1
filtered_data_EC1 <- aldex2_EC1[aldex2_EC1$effect > 1, c("effect", "we.ep")]
# Write the filtered data to a CSV file
write.csv(filtered_data_EC1, "images/data/filtered_data_EC1.csv", row.names = TRUE)





#MW and MA plots

png("images/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,3))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Dispersion", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Dispersion", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.6, rare.cex = 0.6, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Dispersion", ylab = "Difference",col = "blue", rare.col = "red", called.col = "green")
title(main = "(PW) MW Plot")

m <- aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(PW) MA Plot")
invisible(dev.off())



#Relationship between effect, difference, and P values

png("images/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(2,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")


plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())



# Subset-1
df_EC1 = aldex2_EC1 %>% tibble::rownames_to_column(var = "EC") %>% 
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)


df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)


write.table(df_EC1, file = "images/data/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

write.table(df_PW1, file = "images/data/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)



df_effect_filter <- df_EC1[, c( "EC", "description", "effect", "we.ep")]
write.table(df_EC1, file = "images/data/ALDEx2_picrust2_EC_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)


df_effect_filter_PW <- df_PW1[, c( "Pathway", "description", "effect", "we.ep")]
write.table(df_PW1, file = "images/data/ALDEx2_picrust2_PW_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)


#Pathway effect size
top300_data_PW <- df_effect_filter_PW[order(df_effect_filter_PW$effect, decreasing = TRUE), ]


top30_data_PW2 <- top30_data_PW %>%
  filter(effect > 1 & we.ep <= 0.05) %>%
  arrange(desc(effect))


ggplot(top30_data_PW2, aes(y = effect, x = description)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(y = "Effect Size", x = "Pathway Names") +
  ggtitle("Pathways") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 25),  # Adjust size here
    axis.text.y = element_text(size = 16),  # Adjust size here
    axis.title = element_text(size = 18),  # Adjust size here
    plot.title = element_text(size = 18, hjust = 0.5),  # Adjust size and center title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

ggsave("images/pathway_effect1.png", width = 20, height = 20)
write.table(top30_data_PW1, file = "outfiles/top30_data_PW2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#Enzymes
df_effect_filter_EC <- df_EC1[, c( "EC", "description", "effect", "we.ep")]
write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

top30_data_effect_EC <- head(df_effect_filter_EC[order(df_effect_filter_EC$effect, decreasing = TRUE), ], 300)

top30_data_effect_EC1 <- top30_data_effect_EC %>%
  filter(effect > 5 & we.ep <= 0.01) %>%
  arrange(desc(effect))



ggplot(top30_data_effect_EC1, aes(y = effect, x = description)) +
  geom_bar(stat = "identity", fill = "blue")  +
  labs(y = "Effect Size", x = "Enzyme Names") +
  ggtitle("Enzymes ") +
  theme_minimal() +


  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 25),  # Adjust size here
    axis.text.y = element_text(size = 12),  # Adjust size here
    axis.title = element_text(size = 14),  # Adjust size here
    plot.title = element_text(size = 20, hjust = 0.5),  # Adjust size and center title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

ggsave("images/EC_effect2.png", , width = 18, height = 14)
write.table(top30_data_effect_EC1, file = "images/data/top30_data_EC.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)



#diff.win



#plotted based on diff.btw

#Pathway
df_diff.win.treated_filter_PW <- df_PW1[, c( "Pathway", "description", "diff.win", "we.ep")]
top30_data_diff.win.treated_PW <- df_diff.win.treated_filter_PW[order(df_diff.win.treated_filter_PW$diff.win, decreasing = TRUE), ]
top30_data_diff.win.treated_PW1 <- df_diff.win.treated_filter_PW %>%
  filter(diff.win > 0.5 & we.ep <= 0.5)



ggplot(top30_data_diff.win.treated_PW1, aes(y = diff.win, x = description)) +
  geom_bar(stat = "identity", fill = "blue")  +
  labs(y = "diff.win", x = "Pathway") +
  ggtitle("Pathway") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16),  # Adjust size here
    axis.text.y = element_text(size = 12),  # Adjust size here
    axis.title = element_text(size = 14),  # Adjust size here
    plot.title = element_text(size = 16, hjust = 0.5),  # Adjust size and center title
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()  # Remove minor grid lines
  )

ggsave("images/pathway_diff.win.treated2.png",  width = 15, height = 14)







library(ggfortify)

# Perform PCA analysis
pca_result <- prcomp(p2EC1, scale = TRUE)

# Access the results
summary(pca_result) # Provides summary of the PCA results

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


