
#This script is for mutation analysis of 13513 position mitochondrial 
#chromosome in .bam files obtained from Nanopore sequencing 

#Load libraries
library(Rsamtools)
library(VariantAnnotation)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)

# .bam files list
BamFiles <- list.files(path = "route/to/files.bam", pattern = "*.bam$", 
                        full.names = TRUE)

# Positon of interest in mitochondrial genome
MitoChr <- "chrM"  # Or "MT", depending of BAM index previously created
Pos <- 13513



# Pileups parameters
PileupParams <- PileupParam(max_depth = 10000, min_base_quality = 1)

#  Extracting bases in 13513 position
BaseCounts <- function(bam) {
  Result <- pileup(bam, scanBamParam = ScanBamParam(which = GRanges(MitoChr, 
                  IRanges(Pos, Pos))),
                   pileupParam = PileupParams)
  
# If there are no readings in the position, return empty row with sample name
  if (nrow(Result) == 0) {
    return(data.frame(seqnames = MitoChr, pos = Pos, nucleotide = NA, 
                      count = 0))
  }
  
  return(Result)
}

# Apply function to all .bam files
BCounts <- lapply(BamFiles, BaseCounts)

# Add sample name to all results
BCounts <- Map(function(df, sample) {
  df$Sample <- sample
  return(df)
}, BCounts, basename(BamFiles))

# Define columns
ReqCols <- c("seqnames", "pos", "strand", "nucleotide", "count", "which_label",
             "Sample")

# Homogenize dataframes
BCounts <- lapply(BCounts, function(df) {
  # Add missing columns
  MissCols <- setdiff(ReqCols, colnames(df))
  for (col in MissCols) df[[col]] <- NA  # Complete with NAs
  
  # Convertir nusleotide column to factor
  if (!is.factor(df$nucleotide)) {
    df$nucleotide <- factor(df$nucleotide, levels = c("A", "C", "G", "T", "N", 
                                                      "X", "Y", "Z")) 
  }
  
  return(df[, ReqCols, drop = FALSE])  # Make sure correct order
})

# Bind in an Data Frame
BaseResults <- do.call(rbind, BCounts)
BaseResults <- BaseResults %>%
  filter(strand == "+") %>%  # Filter by positive strand
  select(1,2,4,5,7) %>%  # Select specific columns
  filter(!nucleotide %in% c("-", NA))  # Eliminate rows with "-" or NAs 



# Counts base per sample
DfSummary <- BaseResults %>%
  group_by(Sample, nucleotide) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  tidyr::spread(nucleotide, count, fill = 0)

#Vector of condition for every sample
DfSummary$Cond <- c("LoDe", "LoDe", "LoDe", "HiDe", "HiDe", 
                    "HiDe", "LoDMSO", "LoDMSO", "LoDMSO", "HiDMSO", "HiDMSO", 
                    "HiDMSO")

#Convert dataframe to long format
DfSummaryLong <- DfSummary %>%
  pivot_longer(cols = c(A, G), names_to = "Nucleotide", values_to = "Count")

#stacked bar plot
ggplot(DfSummaryLong, aes(x = Cond, y = Count, fill = Nucleotide)) +
  geom_bar(stat = "identity")+
  labs(title = "Mitochondrial Chromosome (13513 position)",
       x = "Condition", y = "Mutation counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#boxplot
ggplot(DfSummaryLong, aes(x = Cond, y = Count, fill = Nucleotide)) +
  geom_boxplot()  +
  labs(title = "Mitochondrial Chromosome (13513 position)",
       x = "Sample", y = "Mutation Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#bar plot
ggplot(DfSummaryLong, aes(x = Cond, y = Count, fill = Nucleotide)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mitochondrial Chromosome (13513 position)",
       x = "Sample", y = "Mutation Counts")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#write results in excel file
writexl::write_xlsx(DfSummaryLong, "d:/MtMet26Feb/MitoMutationCounts.xlsx")
