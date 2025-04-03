
#Complete script for methylation analysis and Plots generation


#Load Libraries
library(NanoMethViz)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rtracklayer)
library(impute)
library(edgeR)
library(reshape2)
library(pheatmap)
library(umap)
library(data.table)
library(patchwork)
library(purrr)

# Set Bam files and path
ModSamples <- c("Barcode01.bam","Barcode02.bam",
                "Barcode03.bam","Barcode04.bam", 
                "Barcode05.bam", "Barcode06.bam",
                "Barcode07.bam", "Barcode08.bam", 
                "Barcode09.bam","Barcode10.bam",
                "Barcode11.bam","Barcode12.bam")

#Set paths to .bam files
ModPaths <- c("D:/MtMet26Feb/Indexed/Barcode01.bam",
              "D:/MtMet26Feb/Indexed/Barcode02.bam",
              "D:/MtMet26Feb/Indexed/Barcode03.bam",
              "D:/MtMet26Feb/Indexed/Barcode04.bam",
              "D:/MtMet26Feb/Indexed/Barcode05.bam",
              "D:/MtMet26Feb/Indexed/Barcode06.bam",
              "D:/MtMet26Feb/Indexed/Barcode07.bam",
              "D:/MtMet26Feb/Indexed/Barcode08.bam",
              "D:/MtMet26Feb/Indexed/Barcode09.bam",
              "D:/MtMet26Feb/Indexed/Barcode10.bam",
              "D:/MtMet26Feb/Indexed/Barcode11.bam",
              "D:/MtMet26Feb/Indexed/Barcode12.bam")


#Set Sample Conditions
DataSample <- data.frame(sample = ModSamples,
                         group = c("LoDec", "LoDec", "LoDec", "HiDec", "HiDec", 
                                   "HiDec", "LoDMSO", "LoDMSO", "LoDMSO", 
                                   "HiDMSO", "HiDMSO", "HiDMSO"), 
                         stringsAsFactors = F)

DataSample2 <- data.frame(sample= ModSamples,
                          group= c("Dec", "Dec", "Dec", "Dec", "Dec", "Dec", 
                                   "DMSO", "DMSO", "DMSO", "DMSO", "DMSO", 
                                   "DMSO"), stringsAsFactors = F)

DataSample3 <- data.frame(sample = ModSamples, 
                          group = c( "Lo", "Lo", "Lo", "Hi", "Hi", "Hi", "Lo", 
                                     "Lo", "Lo", "Hi", "Hi", "Hi"), 
                          stringsAsFactors=F)

DataSample4 <- data.frame (sample= ModSamples[c(4:6, 10:12)],
                           group= c("HiDec", "HiDec", "HiDec", "HiDMSO", 
                                    "HiDMSO", "HiDMSO"), stringsAsFactors=F)

DataSample5 <- data.frame (sample= ModSamples[c(1:3, 7:9)],
                           group= c("LoDec", "LoDec", "LoDec", "LoDMSO", 
                                    "LoDMSO", "LoDMSO"), stringsAsFactors=F)


# Convert Bam files to ModBamFiles object
ModBams <- ModBamFiles(samples = ModSamples, paths = ModPaths)

ModBams4 <- ModBamFiles(samples=ModSamples[c(4:6, 10:12)], 
                        paths = ModPaths[c(4:6, 10:12)])

ModBams5 <- ModBamFiles(samples = ModSamples[c(1:3, 7:9)], 
                        paths = ModPaths[c(1:3, 7:9)])


# Exons annotations
ExonsHoSap <- get_exons_hg38()


#ModBamResult object
Result <- ModBamResult(ModBams, DataSample, ExonsHoSap, mod_code = "m")

Result2 <- ModBamResult(ModBams, DataSample2, ExonsHoSap, mod_code = "m")

Result3 <- ModBamResult(ModBams, DataSample3, ExonsHoSap, mod_code = "m")

Result4 <- ModBamResult(ModBams4, DataSample4, ExonsHoSap, mod_code = "m")

Result5 <- ModBamResult(ModBams5, DataSample5, ExonsHoSap, mod_code = "m")

# Write Tabix File for further analysis
modbam_to_tabix(Result, "D:/MtMet26Feb/Tabix.tsv.bgz", mod_code = "m")
modbam_to_tabix(Result2, "D:/MtMet26Feb/Tabix2.tsv.bgz", mod_code = "m")
modbam_to_tabix(Result3, "D:/MtMet26Feb/Tabix3.tsv.bgz", mod_code = "m")
modbam_to_tabix(Result4, "D:/MtMet26Feb/Tabix4.tsv.bgz", mod_code = "m")
modbam_to_tabix(Result5, "D:/MtMet26Feb/Tabix5.tsv.bgz", mod_code = "m")



#Creation of BSSEQ object
MetSeq <- methy_to_bsseq("D:/MtMet26Feb/Tabix.tsv.bgz", 
                         out_folder = tempdir(), verbose = TRUE)

MetSeq2 <- methy_to_bsseq("D:/MtMet26Feb/Tabix2.tsv.bgz", 
                         out_folder = tempdir(), verbose = TRUE)

MetSeq3 <- methy_to_bsseq("D:/MtMet26Feb/Tabix3.tsv.bgz", 
                         out_folder = tempdir(), verbose = TRUE)

MetSeq4 <- methy_to_bsseq("D:/MtMet26Feb/Tabix4.tsv.bgz", 
                         out_folder = tempdir(), verbose = TRUE)

MetSeq5 <- methy_to_bsseq("D:/MtMet26Feb/Tabix5.tsv.bgz", 
                         out_folder = tempdir(), verbose = TRUE)


#Tables of chromosomes Annotations
ChromInfo <- getChromInfoFromUCSC("hg38")
ChromInfo <- subset(ChromInfo, chrom %in% c(paste0("chr", 1:22), "chrX","chrM"))
ChromInfo$start <- 0
ChromInfo$end <- ChromInfo$size
ChromInfo <- ChromInfo[, c(1,5,6)]
rownames(ChromInfo) <- ChromInfo$chrom 
ChromMit <- subset(ChromInfo, chrom == "chrM")



#Creation of a Log-Methylation Ratio matrix
EdgerMet <- bsseq_to_log_methy_ratio(MetSeq, ChromInfo)
EdgerMet <- apply(EdgerMet, 2, function(x) {
  x[is.na(x)] <- 0  # Replace NA with 0
  1 / (1 + exp(-x)) # Apply sigmoid conversion
})
EdgerMetT <- data.frame(t(EdgerMet))
colnames(EdgerMetT) <- c(paste0("chr", 1:22), "chrX", "chrM")
EdgerMetT$Cond <- DataSample$group[match(rownames(EdgerMetT), 
                                         DataSample$sample)]

EdgerMet2 <- bsseq_to_log_methy_ratio(MetSeq2, ChromInfo)
EdgerMet2 <- apply(EdgerMet2, 2, function(x) {
  x[is.na(x)] <- 0  # Replace NA with 0
  1 / (1 + exp(-x)) # Apply sigmoid conversion
})
EdgerMetT2 <- data.frame(t(EdgerMet2))
colnames(EdgerMetT2) <- c(paste0("chr", 1:22), "chrX", "chrM")
EdgerMetT2$Cond <- DataSample2$group[match(rownames(EdgerMetT2), 
                                           DataSample2$sample)]

EdgerMet3 <- bsseq_to_log_methy_ratio(MetSeq3, ChromInfo)
EdgerMet3 <- apply(EdgerMet3, 2, function(x) {
  x[is.na(x)] <- 0  # Replace NA with 0
  1 / (1 + exp(-x)) # Apply sigmoid conversion
})
EdgerMetT3 <- data.frame(t(EdgerMet3))
colnames(EdgerMetT3) <- c(paste0("chr", 1:22), "chrX", "chrM")
EdgerMetT3$Cond <- DataSample3$group[match(rownames(EdgerMetT3),
                                           DataSample3$sample)]

EdgerMet4 <- bsseq_to_log_methy_ratio(MetSeq4, ChromInfo)
EdgerMet4 <- apply(EdgerMet4, 2, function(x) {
  x[is.na(x)] <- 0  # Replace NA with 0
  1 / (1 + exp(-x)) # Apply sigmoid conversion
})
EdgerMetT4 <- data.frame(t(EdgerMet4))
colnames(EdgerMetT4) <- c(paste0("chr", 1:22), "chrX", "chrM")
EdgerMetT4$Cond <- DataSample4$group[match(rownames(EdgerMetT4), 
                                           DataSample4$sample)]
EdgerMetT4 <- na.omit(EdgerMetT4)

EdgerMet5 <- bsseq_to_log_methy_ratio(MetSeq5, ChromInfo)
EdgerMet5 <- apply(EdgerMet5, 2, function(x) {
  x[is.na(x)] <- 0  # Replace NA with 0
  1 / (1 + exp(-x)) # Apply sigmoid conversion
})
EdgerMetT5 <- data.frame(t(EdgerMet5))
colnames(EdgerMetT5) <- c(paste0("chr", 1:22), "chrX", "chrM")
EdgerMetT5$Cond <- DataSample5$group[match(rownames(EdgerMetT5), 
                                           DataSample5$sample)]
EdgerMetT5 <- na.omit(EdgerMetT5)




#Long Format  DataFrame
Metlong <- reshape2::melt(EdgerMetT, id.vars = "Cond")

# Boxplot

#Sort Conditions
Metlong$Cond <- factor(Metlong$Cond, levels = c("HiDMSO", "HiDec", "LoDMSO", 
                                                "LoDec"))

# Define inset region(chrM)
Zoom <- "chrM"


#List for means comparisons
Comparisons <- list(
  c("LoDec", "HiDec"),
  c("LoDec", "LoDMSO"),
  c("LoDec", "HiDMSO"),
  c("HiDec", "LoDMSO"),
  c("HiDec", "HiDMSO"),
  c("LoDMSO", "HiDMSO")
)

Comparisons2 <- list(
  c("DMSO", "Dec"))

Comparisons3 <- list(
  c("Hi", "Lo"))

Comparisons4 <- list(
  c("HiDec", "HiDMSO"))

Comparisons5 <- list(
  c("LoDMSO", "LoDec"))

# Main Plot
P1 <- ggplot(Metlong, aes(x = variable, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Chromosome", y = "Log-Methyl Ratio",
       title = "Methylation Distribution by Chromosome") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  stat_compare_means( method = "kruskal.test", 
                     label = "p.signif")+
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",  
    "HiDec" = "#EB58F8",  
    "LoDMSO" = "#98BDf4",  
    "HiDMSO" = "#5757F0")) +
  # Conection Lines
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.275), 
               linetype = "solid", color = "black", size = 1) + 
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.076), 
               linetype = "solid", color = "black", size = 1)

# Inset (chrM)
P1_2 <- ggplot(subset(Metlong, variable == Zoom), aes(x = Cond, y = value, 
                                                      fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(x = "", y = "", title = "ChrM") + 
  stat_compare_means(comparisons = Comparisons, label = "p.signif", 
                     method= "kruskal.test")+
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",  
    "HiDec" = "#EB58F8",  
    "LoDMSO" = "#98BDf4",  
    "HiDMSO" = "#5757F0")) +
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))  

# Match Plots
  P1 + inset_element(P1_2, left = 0.6, bottom = 0.06, right = 0.94, top = 0.4)


# Delete Condition column

MetMatrix <- as.matrix(EdgerMetT[, -ncol(EdgerMetT)])
MetMatrix <- MetMatrix[, colnames(MetMatrix) != "chrM"]

# Heatmap
pheatmap(MetMatrix, scale = "row", 
         annotation_row = data.frame(Cond = EdgerMetT$Cond,
                                     row.names = rownames(EdgerMetT)), 
         cluster_rows = TRUE, cluster_cols = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Methylation Heatmap by chromosome",
         show_rownames = F,
         fontsize_row = 6,
         fontsize_col = 8)


# Data Scale
MetMatrixScale <- scale(MetMatrix)

#Check Rows
NSamples <- nrow(MetMatrixScale)

# Adjust Neighbors
NNeighbors <- min(12, NSamples - 1)  # Samples -1

# Umap
MetUmap <- umap(MetMatrixScale, n_neighbors = NNeighbors)

# Convert to dataframe
UmapDf <- as.data.frame(MetUmap$layout)
colnames(UmapDf) <- c("UMAP1", "UMAP2")

# Add condition
UmapDf$Cond <- EdgerMetT$Cond

# Umap Plot with ggplot2
ggplot(UmapDf, aes(x = UMAP1, y = UMAP2, color = Cond)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMap by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green3", "purple")) 



#Long Format  DataFrame
Metlong2 <- reshape2::melt(EdgerMetT2, id.vars = "Cond")

#Sort Conditions
Metlong2$Cond <- factor(Metlong2$Cond, levels = c("DMSO", "Dec"))

# Boxplot


# Define inset region(chrM)
Zoom <- "chrM"

# Main Plot
P2 <- ggplot(Metlong2, aes(x = variable, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Chromosome", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Chromosome") + 
  stat_compare_means( method ="wilcox.test", label="p.signif")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(
    "DMSO" = "#626256",  
    "Dec" = "#EB58F8")) +
  # Conection Lines
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.271), 
               linetype = "solid", color = "black", size = 1) + 
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.076), 
               linetype = "solid", color = "black", size = 1)

# Inset (chrM)
P2_2 <- ggplot(subset(Metlong2, variable == Zoom), aes(x = Cond, y = value, 
                                                      fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(x = "", y = "", title = "ChrM") +  
  scale_fill_manual(values = c(
    "DMSO" = "#626256",  
    "Dec" = "#EB58F8")) +
  stat_compare_means(comparisons=Comparisons2, label="p.signif", 
                     method = "wilcox.test",label.y=0.0019)+
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))  

# Match Plots
P2 + inset_element(P2_2, left = 0.6, bottom = 0.06, right = 0.94, top = 0.4)


# Delete Condition column
MetMatrix2 <- as.matrix(EdgerMetT2[, -ncol(EdgerMetT2)])
MetMatrix2 <- MetMatrix2[, colnames(MetMatrix2) != "chrM"]
# Heatmap
pheatmap(MetMatrix2, scale = "row", 
         annotation_row = data.frame(Cond = EdgerMetT2$Cond, 
                                     row.names = rownames(EdgerMetT2)), 
         cluster_rows = TRUE, cluster_cols = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Methylation Heatmap by chromosome",
         show_rownames = F)

# Data Scale
MetMatrixScale2 <- scale(MetMatrix2)

#Check Rows
NSamples2 <- nrow(MetMatrixScale2)

# Adjust Neighbors
NNeighbors2 <- min(12, NSamples2 - 1)  # Samples -1

# Umap
MetUmap2 <- umap(MetMatrixScale2, n_neighbors = NNeighbors2)

# Convert to dataframe
UmapDf2 <- as.data.frame(MetUmap2$layout)
colnames(UmapDf2) <- c("UMAP1", "UMAP2")

# Add condition
UmapDf2$Cond <- EdgerMetT2$Cond

# Umap Plot with ggplot2
ggplot(UmapDf2, aes(x = UMAP1, y = UMAP2, color = Cond)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMap by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green3", "purple"))


#Long Format  DataFrame
Metlong3 <- reshape2::melt(EdgerMetT3, id.vars = "Cond")

#Sort Conditions
Metlong3$Cond <- factor(Metlong3$Cond, levels = c("Hi", "Lo"))

# Boxplot


# Define inset region(chrM)
Zoom <- "chrM"

# Main Plot
P3 <- ggplot(Metlong3, aes(x = variable, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Chromosome", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Chromosome") + 
  stat_compare_means( method ="wilcox.test", label="p.signif")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(
    "Hi" = "#FF0000",  
    "Lo" = "#FFC201")) +
  # Conection Lines
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.271), 
               linetype = "solid", color = "black", size = 1) + 
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.076), 
               linetype = "solid", color = "black", size = 1)

# Inset (chrM)
P3_2 <- ggplot(subset(Metlong3, variable == Zoom), aes(x = Cond, y = value, 
                                                       fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(x = "", y = "", title = "ChrM") +  
  scale_fill_manual(values = c(
    "Hi" = "#FF0000",  
    "Lo" = "#FFC201")) +
  stat_compare_means(comparisons=Comparisons3, label="p.signif", 
                     method = "wilcox.test",label.y=0.0019)+
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))  

# Match Plots
P3 + inset_element(P3_2, left = 0.6, bottom = 0.06, right = 0.94, top = 0.4)

# Delete Condition column
MetMatrix3 <- as.matrix(EdgerMetT3[, -ncol(EdgerMetT3)])
MetMatrix3 <- MetMatrix3[, colnames(MetMatrix3) != "chrM"]
# Heatmap
pheatmap(MetMatrix3, scale = "row", 
         annotation_row = data.frame(Cond = EdgerMetT3$Cond, row.names = rownames(EdgerMetT3)), 
         cluster_rows = TRUE, cluster_cols = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Methylation Heatmap by chromosome",
         show_rownames = F)

# Data Scale
MetMatrixScale3 <- scale(MetMatrix3)

#Check Rows
NSamples3 <- nrow(MetMatrixScale3)

# Adjust Neighbors
NNeighbors3 <- min(12, NSamples3 - 1)  # Samples -1

# Umap
MetUmap3 <- umap(MetMatrixScale3, n_neighbors = NNeighbors3)

# Convert to dataframe
UmapDf3 <- as.data.frame(MetUmap3$layout)
colnames(UmapDf3) <- c("UMAP1", "UMAP2")

# Add condition
UmapDf3$Cond <- EdgerMetT3$Cond

# Umap Plot with ggplot2
ggplot(UmapDf3, aes(x = UMAP1, y = UMAP2, color = Cond)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMap by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green3", "purple"))




#Long Format  DataFrame
Metlong4 <- reshape2::melt(EdgerMetT4, id.vars = "Cond")
#Sort Conditions
Metlong4$Cond <- factor(Metlong4$Cond, levels = c("HiDMSO", "HiDec"))

# Boxplot


# Define inset region(chrM)
Zoom <- "chrM"

# Main Plot
P4 <- ggplot(Metlong4, aes(x = variable, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Chromosome", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Chromosome") + 
  stat_compare_means( method ="wilcox.test", label="p.signif")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(
    "HiDec" = "#EB58F8",  
    "HiDMSO" = "#5757F0")) +
  # Conection Lines
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.271), 
               linetype = "solid", color = "black", size = 1) + 
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.076), 
               linetype = "solid", color = "black", size = 1)

# Inset (chrM)
P4_2 <- ggplot(subset(Metlong4, variable == Zoom), aes(x = Cond, y = value, 
                                                       fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(x = "", y = "", title = "ChrM") +  
  scale_fill_manual(values = c(
    "HiDec" = "#EB58F8",  
    "HiDMSO" = "#5757F0")) +
  stat_compare_means(comparisons=Comparisons4, label="p.signif", 
                     method = "wilcox.test", label.y=0.00186)+
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))  

# Match Plots
P4 + inset_element(P4_2, left = 0.6, bottom = 0.06, right = 0.94, top = 0.4)



# Delete Condition column
MetMatrix4 <- as.matrix(EdgerMetT4[, -ncol(EdgerMetT4)])
MetMatrix4 <- MetMatrix4[, colnames(MetMatrix4) != "chrM"]
# Heatmap
pheatmap(MetMatrix4, scale = "row", 
         annotation_row = data.frame(Cond = EdgerMetT4$Cond, row.names = rownames(EdgerMetT4)), 
         cluster_rows = TRUE, cluster_cols = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Methylation Heatmap by chromosome",
         show_rownames = F)


# Data Scale
MetMatrixScale4 <- scale(MetMatrix4)

#Check Rows
NSamples4 <- nrow(MetMatrixScale4)

# Adjust Neighbors
NNeighbors4 <- min(6, NSamples4 - 1)  # Samples -1

# Umap
MetUmap4 <- umap(MetMatrixScale4, n_neighbors = NNeighbors4)

# Convert to dataframe
UmapDf4 <- as.data.frame(MetUmap4$layout)
colnames(UmapDf4) <- c("UMAP1", "UMAP2")

# Add condition
UmapDf4$Cond <- EdgerMetT4$Cond

# Umap Plot with ggplot2
ggplot(UmapDf4, aes(x = UMAP1, y = UMAP2, color = Cond)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMap by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue"))






#Long Format  DataFrame
Metlong5 <- reshape2::melt(EdgerMetT5, id.vars = "Cond")

#Sort Conditions
Metlong5$Cond <- factor(Metlong5$Cond, levels = c("LoDMSO", "LoDec"))

# Boxplot


# Define inset region(chrM)
Zoom <- "chrM"

# Main Plot
P5 <- ggplot(Metlong5, aes(x = variable, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) + 
  theme_minimal() +
  labs(x = "Chromosome", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Chromosome") + 
  stat_compare_means( method ="wilcox.test", label="p.signif")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",  
    "LoDMSO" = "#98BDF4")) +
  # Conection Lines
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.271), 
               linetype = "solid", color = "black", size = 1) + 
  geom_segment(aes(x = Zoom, xend = "chrX", 
                   y = 0.002, yend = 0.076), 
               linetype = "solid", color = "black", size = 1)

# Inset (chrM)
P5_2 <- ggplot(subset(Metlong5, variable == Zoom), aes(x = Cond, y = value, 
                                                       fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(x = "", y = "", title = "ChrM") +  
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",  
    "LoDMSO" = "#98BDF4")) +
  stat_compare_means(comparisons=Comparisons5, label="p.signif", 
                     method = "wilcox.test", label.y=0.00186)+
  theme(legend.position = "none", 
        panel.border = element_rect(color = "black", fill = NA, size = 1.5))  

# Match Plots
P5 + inset_element(P5_2, left = 0.6, bottom = 0.06, right = 0.94, top = 0.4)


# Delete Condition column
MetMatrix5 <- as.matrix(EdgerMetT5[, -ncol(EdgerMetT5)])
MetMatrix5 <- MetMatrix5[, colnames(MetMatrix5) != "chrM"]
# Heatmap
pheatmap(MetMatrix5, scale = "row", 
         annotation_row = data.frame(Cond = EdgerMetT5$Cond, row.names = rownames(EdgerMetT5)), 
         cluster_rows = TRUE, cluster_cols = F, 
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         main = "Methylation Heatmap by chromosome",
         show_rownames = F)


# Data Scale
MetMatrixScale5 <- scale(MetMatrix5)

#Check Rows
NSamples5 <- nrow(MetMatrixScale5)

# Adjust Neighbors
NNeighbors5 <- min(6, NSamples5 - 1)  # Samples -1

# Umap
MetUmap5 <- umap(MetMatrixScale5, n_neighbors = NNeighbors5)

# Convert to dataframe
UmapDf5 <- as.data.frame(MetUmap5$layout)
colnames(UmapDf5) <- c("UMAP1", "UMAP2")

# Add condition
UmapDf5$Cond <- EdgerMetT5$Cond

# Umap Plot with ggplot2
ggplot(UmapDf5, aes(x = UMAP1, y = UMAP2, color = Cond)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMap by condition", x = "UMAP 1", y = "UMAP 2") +
  theme_minimal() +
  scale_color_manual(values = c("red", "blue", "green3", "purple"))


# Combined boxplot-violinplot by condition 



ggplot(Metlong, aes(x = Cond, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,1.5)+
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) +
  theme_minimal() +
  labs(x = "Condition", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons= Comparisons, label="p.signif", 
                     method = "wilcox.test", label.x= 1.5)+
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",  
    "HiDec" = "#EB58F8",  
    "LoDMSO" = "#98BDf4",  
    "HiDMSO" = "#5757F0"   
  ))


ggplot(Metlong2, aes(x = Cond, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) +
  ylim(0,1)+
  theme_minimal() +
  labs(x = "Condition", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons=Comparisons2, label="p.signif", 
                     method = "wilcox.test", label.x= 1.5)+
  scale_fill_manual(values = c(
    "DMSO" = "#626256",  
    "Dec" = "#EB58F8"))





ggplot(Metlong3, aes(x = Cond, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) +
  ylim(0,1)+
  theme_minimal() +
  labs(x = "Condition", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons= Comparisons3, label="p.signif", 
                     method = "wilcox.test", label.x= 1.5)+
  scale_fill_manual(values = c(
    "Hi" = "#FF0000",  
    "Lo" = "#FFC201"))




ggplot(Metlong4, aes(x = Cond, y = value, fill = Cond)) +
  geom_boxplot(outlier.shape = NA) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) +
  ylim(0,1)+
  theme_minimal() +
  labs(x = "Condition", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons= Comparisons4, label="p.signif", 
                     method = "wilcox.test", label.x= 1.5)+
  scale_fill_manual(values = c(
    "HiDec" = "#EB58F8",  
    "HiDMSO" = "#5757F0"   
  ))




ggplot(Metlong5, aes(x = Cond, y = value, fill = Cond)) +
  geom_violin(alpha = 0.5, position = position_dodge(0.9)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(0,1)+
  theme_minimal() +
  labs(x = "Condition", y = "Log-Methyl Ratio", 
       title = "Methylation Distribution by Condition") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  stat_compare_means(comparisons= Comparisons5, label="p.signif", 
                     method = "wilcox.test", label.x= 1.5)+
  scale_fill_manual(values = c(
    "LoDec" = "#EDADED",
    "LoDMSO" = "#98BDf4"
  ))
  


#Wilcoxon Test and p-values by chromosome


# Obtain chromosomes
Chromosomes <- setdiff(colnames(EdgerMetT), "Cond")

# Generate all comparisons
Compar <- expand.grid(Cond1 = unique(EdgerMetT$Cond),
                           Cond2 = unique(EdgerMetT$Cond),
                           stringsAsFactors = FALSE) %>%
  filter(Cond1 != Cond2)  # Eliminar comparaciones redundantes

# Wilcoxon iterative function
WilByChr <- function(chr, cond1, cond2, data) {
  x <- data %>% filter(Cond == cond1) %>% pull(!!sym(chr))
  y <- data %>% filter(Cond == cond2) %>% pull(!!sym(chr))
  
  if (length(x) > 1 & length(y) > 1) {
    test <- wilcox.test(x, y)
    return(data.frame(Chromosome = chr, 
                      Group1 = cond1, 
                      Group2 = cond2, 
                      W = test$statistic, 
                      p_value = test$p.value))
  } else {
    return(NULL)  # Omit comparisons with not enough data
  }
}

# Apply function to all chromosomes and conditions
Results <- map_dfr(Chromosomes, function(chr) {
  map_dfr(seq_len(nrow(Compar)), function(i) {
    WilByChr(chr, Compar$Cond1[i], Compar$Cond2[i], EdgerMetT)
  })
})


# Obtain chromosomes
Chromosomes2 <- setdiff(colnames(EdgerMetT2), "Cond")

# Generate all comparisons
Compar2 <- expand.grid(Cond1 = unique(EdgerMetT2$Cond),
                      Cond2 = unique(EdgerMetT2$Cond),
                      stringsAsFactors = FALSE) %>%
  filter(Cond1 != Cond2)  # Eliminar comparaciones redundantes

# Wilcoxon iterative function
WilByChr2 <- function(chr, cond1, cond2, data) {
  x <- data %>% filter(Cond == cond1) %>% pull(!!sym(chr))
  y <- data %>% filter(Cond == cond2) %>% pull(!!sym(chr))
  
  if (length(x) > 1 & length(y) > 1) {
    test <- wilcox.test(x, y)
    return(data.frame(Chromosome = chr, 
                      Group1 = cond1, 
                      Group2 = cond2, 
                      W = test$statistic, 
                      p_value = test$p.value))
  } else {
    return(NULL)  # Omit comparisons with not enough data
  }
}

# Apply function to all chromosomes and conditions
Results2 <- map_dfr(Chromosomes2, function(chr) {
  map_dfr(seq_len(nrow(Compar2)), function(i) {
    WilByChr2(chr, Compar2$Cond1[i], Compar2$Cond2[i], EdgerMetT2)
  })
})

# Obtain chromosomes
Chromosomes3 <- setdiff(colnames(EdgerMetT3), "Cond")

# Generate all comparisons
Compar3 <- expand.grid(Cond1 = unique(EdgerMetT3$Cond),
                      Cond2 = unique(EdgerMetT3$Cond),
                      stringsAsFactors = FALSE) %>%
  filter(Cond1 != Cond2)  # Eliminar comparaciones redundantes

# Wilcoxon iterative function
WilByChr3 <- function(chr, cond1, cond2, data) {
  x <- data %>% filter(Cond == cond1) %>% pull(!!sym(chr))
  y <- data %>% filter(Cond == cond2) %>% pull(!!sym(chr))
  
  if (length(x) > 1 & length(y) > 1) {
    test <- wilcox.test(x, y)
    return(data.frame(Chromosome = chr, 
                      Group1 = cond1, 
                      Group2 = cond2, 
                      W = test$statistic, 
                      p_value = test$p.value))
  } else {
    return(NULL)  # Omit comparisons with not enough data
  }
}

# Apply function to all chromosomes and conditions
Results3 <- map_dfr(Chromosomes3, function(chr) {
  map_dfr(seq_len(nrow(Compar3)), function(i) {
    WilByChr3(chr, Compar3$Cond1[i], Compar3$Cond2[i], EdgerMetT3)
  })
})



# Obtain chromosomes
Chromosomes4 <- setdiff(colnames(EdgerMetT4), "Cond")

# Generate all comparisons
Compar4 <- expand.grid(Cond1 = unique(EdgerMetT4$Cond),
                      Cond2 = unique(EdgerMetT4$Cond),
                      stringsAsFactors = FALSE) %>%
  filter(Cond1 != Cond2)  # Eliminar comparaciones redundantes

# Wilcoxon iterative function
WilByChr4 <- function(chr, cond1, cond2, data) {
  x <- data %>% filter(Cond == cond1) %>% pull(!!sym(chr))
  y <- data %>% filter(Cond == cond2) %>% pull(!!sym(chr))
  
  if (length(x) > 1 & length(y) > 1) {
    test <- wilcox.test(x, y)
    return(data.frame(Chromosome = chr, 
                      Group1 = cond1, 
                      Group2 = cond2, 
                      W = test$statistic, 
                      p_value = test$p.value))
  } else {
    return(NULL)  # Omit comparisons with not enough data
  }
}

# Apply function to all chromosomes and conditions
Results4 <- map_dfr(Chromosomes4, function(chr) {
  map_dfr(seq_len(nrow(Compar4)), function(i) {
    WilByChr(chr, Compar4$Cond1[i], Compar4$Cond2[i], EdgerMetT4)
  })
})



# Obtain chromosomes
Chromosomes5 <- setdiff(colnames(EdgerMetT5), "Cond")

# Generate all comparisons
Compar5 <- expand.grid(Cond1 = unique(EdgerMetT5$Cond),
                      Cond2 = unique(EdgerMetT5$Cond),
                      stringsAsFactors = FALSE) %>%
  filter(Cond1 != Cond2)  # Eliminar comparaciones redundantes

# Wilcoxon iterative function
WilByChr5 <- function(chr, cond1, cond2, data) {
  x <- data %>% filter(Cond == cond1) %>% pull(!!sym(chr))
  y <- data %>% filter(Cond == cond2) %>% pull(!!sym(chr))
  
  if (length(x) > 1 & length(y) > 1) {
    test <- wilcox.test(x, y)
    return(data.frame(Chromosome = chr, 
                      Group1 = cond1, 
                      Group2 = cond2, 
                      W = test$statistic, 
                      p_value = test$p.value))
  } else {
    return(NULL)  # Omit comparisons with not enough data
  }
}

# Apply function to all chromosomes and conditions
Results5 <- map_dfr(Chromosomes5, function(chr) {
  map_dfr(seq_len(nrow(Compar5)), function(i) {
    WilByChr5(chr, Compar5$Cond1[i], Compar5$Cond2[i], EdgerMetT5)
  })
})



#Wilcoxon test for global metehylation by condition

WilResults <- compare_means(value ~ Cond, data = Metlong, 
                                method = "wilcox.test")

WilResults2 <- compare_means(value ~ Cond, data = Metlong2, 
                            method = "wilcox.test")

WilResults3 <- compare_means(value ~ Cond, data = Metlong3, 
                            method = "wilcox.test")

WilResults4 <- compare_means(value ~ Cond, data = Metlong4, 
                            method = "wilcox.test")

WilResults5 <- compare_means(value ~ Cond, data = Metlong5, 
                                                                               method = "wilcox.test")

Results <- subset(Results, Group1 == "LoDMSO" | Group1 == "HiDMSO")

Results2 <- subset(Results2, Group1 == "DMSO")

Results3 <- subset(Results3, Group1 == "Lo")

Results4 <- subset(Results4, Group1 == "HiDMSO")

Results5 <- subset(Results5, Group1 == "LoDMSO")


#Generate log-methyl Ratio (sigmoid converted) table
MetTable <- EdgerMetT

MetTable$Sample <- rownames(MetTable)

#Store in Excel file
writexl::write_xlsx(list("LogMethylRatios(Sigmoid)" = MetTable, 
                         "GlobalMethylation" = WilResults, 
                         "GlobalMethDmsoVsDec" = WilResults2,
                         "GlobalMethHiVsLo" = WilResults3, 
                         "GlobalMethHiDmsoVsLoDec" = WilResults4,
                         "GlobalMethLoDmsoVsLoDec" = WilResults5,
                         "MethByChr" = Results,
                         "MethByChrDmsoVsDec" = Results2,
                         "MethByChrHiVsLo" = Results3,
                         "MethByChrHiDmsoVsHiDec" = Results4,
                         "MethByChrLoDmsoVsLoDec" = Results5), 
           path = "D:/MtMet26Feb/PostAnalisis/FinalResults.xlsx")



Tabix <- "D:/MtMet26Feb/Tabix.tsv.bgz"
Tabix2 <- "D:/MtMet26Feb/Tabix2.tsv.bgz"
Tabix3 <- "D:/MtMet26Feb/Tabix3.tsv.bgz"
Tabix4 <- "D:/MtMet26Feb/Tabix4.tsv.bgz"
Tabix5 <- "D:/MtMet26Feb/Tabix5.tsv.bgz"




Dmr <- NanoMethResult(Tabix, DataSample)
Dmr2 <- NanoMethResult(Tabix2, DataSample2)
Dmr3 <- NanoMethResult(Tabix3, DataSample3)
Dmr4 <- NanoMethResult(Tabix4, DataSample4)
Dmr5 <- NanoMethResult(Tabix5, DataSample5)



plot_region(Dmr,  "chrM", 0, 16569,heatmap=F) + 
  ggplot2::scale_color_manual(values = c("#F700FF","#001AFF", "#FCB0FF", 
                                         "#82C5FF"))+
  ggplot2::coord_cartesian(ylim = c(0, 0.1))
  

plot_region_heatmap(Dmr,  "chrM", 0, 16569) + 
  ggplot2::scale_color_gradient(low="#00CD00", high = "#FFFF00", na.value = NA)
plot_region(Dmr2, "chrM", 0, 16569,heatmap=T)
plot_region(Dmr3, "chrM", 0, 16569, heatmap=T)
plot_region(Dmr4, "chrM", 0, 16569, heatmap=T)
plot_region(Dmr5, "chrM", 0, 16569, heatmap =T)







#