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


# Perform the ALDEx2 analysis for KO
p2KO = round(p2KO)   # importnat to make  integer round
p2KO1 <- (p2KO[, c("X367", "X368", "X369","X364", "X365", "X366")])

treated <- (p2KO1[, c("X367", "X368", "X369")])
untreated <- (p2KO1[, c("X364", "X365", "X366")])
conditions <- c(rep("treated", 3), rep("untreated", 3))

set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO1, mc.samples = 500, test = "t", 
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
head(aldex2_EC1, 10)

#ALDEx2 authors suggest that an effect size of 1 or greater can be used as significance cutoff



png("images/ALDEx2_picrust2_effect_1_treated_vs_untreated.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())

# Assuming filtered_data_PW1 is your filtered data frame containing pathways with effect size > 1

# Load the required library
library(ggplot2)


# Filter the data to retain only rows with effect size greater than 1
filtered_data_PW1 <- aldex2_PW1[aldex2_PW1$effect > 1, c("rab.all", "effect")]

# Write the filtered data to a CSV file
write.csv(filtered_data_PW1, "images/filtered_data_PW1.csv", row.names = TRUE)

# Filter the data to retain only rows with effect size greater than 1
filtered_data_EC1 <- aldex2_EC1[aldex2_EC1$effect > 1, c("rab.all", "effect")]
# Write the filtered data to a CSV file
write.csv(filtered_data_EC1, "images/filtered_data_EC1.csv", row.names = TRUE)

# Filter the data to retain only rows with effect size greater than 1
filtered_data_KO1 <- aldex2_KO1[aldex2_KO1$effect > 1, c("rab.all", "effect")]

# Write the filtered data to a CSV file
write.csv(filtered_data_KO1, "images/filtered_data_KO.csv", row.names = TRUE)



#MW and MA plots

png("images/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,3))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.6, rare.cex = 0.6, 
           called.cex = 0.6, cutoff = 0.5, xlab = "Dispersion", ylab = "Difference",col = "blue", rare.col = "red", called.col = "green")
title(main = "(PW) MW Plot")

m <- aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference", col = "blue", rare.col = "red", called.col = "green")
title(main = "(PW) MA Plot")
invisible(dev.off())



#Relationship between effect, difference, and P values

png("images/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
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

df_KO1 = aldex2_KO1 %>% tibble::rownames_to_column(var = "KO") %>% 
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)


write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_KO1, file = "outfiles/ALDEx2_picrust2_KO_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_PW1, file = "outfiles/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)



df_effect_filter <- df_EC1[, c( "EC", "description", "effect", "we.ep")]
write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

df_effect_filter_KO <- df_KO1[, c( "KO", "description", "effect", "we.ep")]
write.table(df_KO1, file = "outfiles/ALDEx2_picrust2_KO_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

df_effect_filter_PW <- df_PW1[, c( "Pathway", "description", "effect", "we.ep")]
write.table(df_PW1, file = "outfiles/ALDEx2_picrust2_PW_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)


#Pathway
top30_data_PW <- head(df_effect_filter_PW[order(df_effect_filter_PW$effect, decreasing = TRUE), ], 40)

top30_data_PW2 <- top30_data_PW %>%
  filter(we.ep <= 0.06) %>%
  head(40)

ggplot(top30_data_PW2, aes(y = effect, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Effect Size", x = "Pathway Names") +
  ggtitle("Pathways") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/pathway_effect2.png", width = 10, height = 8)
write.table(top30_data_PW2, file = "outfiles/top30_data_PW2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#Enzymes
df_effect_filter_EC <- df_EC1[, c( "EC", "description", "effect", "we.ep")]
write.table(df_EC1, file = "outfiles/ALDEx2_picrust2_EC_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
top30_data_EC1 <- df_effect_filter_EC %>%
  filter(effect > 1) %>%
  arrange(desc(effect))%>%
 head(40)

top30_data_EC2 <- top30_data_EC1 %>%
  filter(we.ep <= 0.06) %>%
  head(40)

ggplot(top30_data_EC2, aes(y = effect, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Effect Size", x = "Enzyme Names") +
  ggtitle("Enzymes ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/EC_effect1.png", , width = 10, height = 8)
write.table(top30_data_EC2, file = "outfiles/top30_data_EC2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
#KO

df_effect_filter_KO <- df_KO1[, c( "KO", "description", "effect", "we.ep")]
write.table(df_KO1, file = "outfiles/ALDEx2_picrust2_KO_results_filter.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
top30_data_KO1 <- df_effect_filter_KO %>%
  filter(effect > 1) %>%
  arrange(desc(effect))%>%
  head(40)

top30_data_KO2 <- top30_data_KO1 %>%
  filter(we.ep <= 0.06) %>%
  head(40)

ggplot(top30_data_KO2, aes(y = effect, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Effect Size", x = "Enzyme Names") +
  ggtitle("KEGG-KO ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/KO_effect1.png",, width = 10, height = 8)
write.table(top30_data_KO2, file = "outfiles/top30_data_KO.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)




#plotted based on rab.win

#Pathway
df_rab.win.treated_filter_PW <- df_PW1[, c( "Pathway", "description", "rab.win.treated", "we.ep")]
top30_data_rab.win.treated <- df_rab.win.treated_filter_PW %>%
  filter(rab.win.treated > 1) %>%
  arrange(desc(rab.win.treated))%>%
  head(40)

top30_data_rab.win.treated2 <- top30_data_rab.win.treated

ggplot(top30_data_rab.win.treated2, aes(y = rab.win.treated, x = Pathway)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "Pathway Names") +
  ggtitle("Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/pathway_rab.win.treated1.png",  width = 8, height = 4)


#Enzymes
df_rab.win.treated_filter_EC <- df_EC1[, c( "EC", "description", "rab.win.treated", "we.ep")]
top30_data_rab.win.treated <- df_rab.win.treated_filter_EC %>%
  filter(rab.win.treated > 1) %>%
  arrange(desc(rab.win.treated))%>%
  head(40)

top30_data_rab.win.treated2 <- top30_data_rab.win.treated
  
ggplot(top30_data_rab.win.treated2, aes(y = rab.win.treated, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "Enzyme Names") +
  ggtitle("Enzymes ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/EC_rab.win.treated.png",  width = 9, height = 6)

#KO

df_rab.win.treated_filter_KO <- df_KO1[, c( "KO", "description", "rab.win.treated", "we.ep")]
top30_data_rab.win.treated <- df_rab.win.treated_filter_KO %>%
  filter(rab.win.treated > 1) %>%
  arrange(desc(rab.win.treated))%>%
  head(40)

top30_data_rab.win.treated2 <- top30_data_rab.win.treated

ggplot(top30_data_rab.win.treated2, aes(y = rab.win.treated, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "KEGG KO Names") +
  ggtitle("Enzymes ") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/KO_rab.win.treated.png", width = 8, height = 4)






#plotted based on diff.win

#Pathway
df_diff.win.treated_filter_PW <- df_PW1[, c( "Pathway", "description", "diff.win", "we.ep")]
top30_data_diff.win.treated <- df_diff.win.treated_filter_PW %>%
  filter(diff.win > 1) %>%
  arrange(desc(diff.win))%>%
  head(40)

top30_data_diff.win.treated2 <- top30_data_diff.win.treated

ggplot(top30_data_diff.win.treated2, aes(y = diff.win, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "Pathway") +
  ggtitle("Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/pathway_diff.win.treated.png",  width = 8, height = 4)


#Enzymes

df_diff.win.treated_filter_EC <- df_EC1[, c( "EC", "description", "diff.win", "we.ep")]
top30_data_diff.win.treated <- df_diff.win.treated_filter_EC %>%
  filter(diff.win > 1) %>%
  arrange(desc(diff.win))%>%
  head(40)

top30_data_diff.win.treated2 <- top30_data_diff.win.treated

ggplot(top30_data_diff.win.treated2, aes(y = diff.win, x = EC)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "Pathway Names") +
  ggtitle("Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/EC_diff.win.treated1_ec1.png",  width = 8, height = 4)
write.table(top30_data_diff.win.treated2, file = "images/top30_data_EC2.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
#KO

df_diff.win.treated_filter_KO <- df_KO1[, c( "KO", "description", "diff.win", "we.ep")]
top30_data_diff.win.treated <- df_diff.win.treated_filter_KO %>%
  filter(diff.win > 1) %>%
  arrange(desc(diff.win))%>%
  head(40)

top30_data_diff.win.treated2 <- top30_data_diff.win.treated

ggplot(top30_data_diff.win.treated2, aes(y = diff.win, x = description)) +
  geom_point(color = "blue") +
  labs(y = "Raw Abundance", x = "Pathway Names") +
  ggtitle("Pathway") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),  # Rotate x-axis labels vertically
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank())  # Remove minor grid lines

ggsave("images/KO_diff.win.treated1.png", width = 8, height = 4)







