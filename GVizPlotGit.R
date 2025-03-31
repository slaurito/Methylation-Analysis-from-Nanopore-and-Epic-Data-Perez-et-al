
#This script is intended to plot the Chromosome location of ESR1 gene, 
#exon annotation, and Illumina EPIC methylation array probes localization in 
#gene coordinates, along with the average of beta values of each probe, 
#for every condition depicted (i.e. "Low" vs "High")
#This script is valid for another examples


#Load needed libraries
library(Gviz)
library(GenomicRanges)
library(biomaRt)
library(dplyr)
library(tidyr)
library(rtracklayer)
library(ggplot2)
library(ggbio)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicFeatures)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(patchwork)
library(cowplot)
library(ggplotify)


#The beta values table (probesID in rows, and the columns must be "ProbeID",
#"chromosome", "start", "end", "Low.1","Low.2", " High.1", "High.2", etc) 
#could be a excel file or a .txt (for a .txt,
#replace readxl::read_xlsx for read.table() or read.delim())


#Load table
Data <- readxl::read_xlsx("path/to/your/file.xlsx")

#Name columns
colnames(Data) <- c("ProbeId", "chromosome", "start", "end", "13L.1", "13L.2", "
                    13L.3", "8L.1", "8L.2", "8L.3", "13H.1", "13H.2", 
                    "13H.3", "8H.1", "8H.2", "8H.3" )


# Convert to GRanges Object
DataMet <- GRanges(
  seqnames = paste0("chr", Data$chromosome),  # Add prefix "chr"to chromosome 
  ranges = IRanges(start = Data$start, end = Data$end),
  IlmnID = Data$ProbeId  # save ProbeID
)

# Add columns with Beta-values
mcols(DataMet) <- Data[, -c(1:4)]  # Exclude first 5 columns
                                   #(coordinates and strand)


#Download Genomic data from Esnsembl using BioMart
Ensembl <- useEnsembl(biomart = "ENSEMBL_MART_ENSEMBL", 
                      dataset = "hsapiens_gene_ensembl")



GeneInfo <- getBM(attributes = c("chromosome_name", "exon_chrom_start", 
                                  "exon_chrom_end", "strand", 
                                  "ensembl_gene_id",
                                 "ensembl_exon_id",
                                  "ensembl_transcript_id",
                                  "external_gene_name"),
                   filters = "external_gene_name", 
                   values = "ESR1", 
                   mart = Ensembl)




colnames(GeneInfo) <- c("chromosome", "start", "end", "strand", 
                        "gene", "exon", "transcript", "symbol" )

GeneInfo$strand <- ifelse(GeneInfo$strand ==1, "+", "-")

GeneInfo$chromosome <- paste0("chr", GeneInfo$chromosome)

#Select transcripts corresponding to ESR1 gene
GeneInfo <- subset(GeneInfo, transcript == "ENST00000440973" | transcript == 
                     "ENST000000443427" | transcript == "ENST00000206249" |
                     transcript == "ENST000004275312") 
#Create Chromosome Track
IdeoTrack <- IdeogramTrack(chromosome = 6, 
                           genome = "hg38", 
                           name = "Chromosome 6")

#Create Axix Track (Coordinates Annotation)
AxisTrack <- GenomeAxisTrack()

#Create Gene Track (ESR1 transcripts and exon annotation)
GeneTrack <- GeneRegionTrack(GeneInfo,
                             name= "ESR1",
                             rotation.title=T,
                             genome = "hg38", 
                             chromosome = 6, 
                             transcript = GeneInfo$transcript,
                             transcriptAnnotation = "transcript",
                             exon=GeneInfo$exon,
                             collapseTranscripts=F,
                             shape="arrow",
                             fill = "green4")


#Create details for Data Track
GroupColors <- c("High" = "#FF0000", "Low" = "#FFC201")
# Assign Group color
GroupVector <- c("Low", "Low", "Low", "Low", "Low", "Low", 
                  "High", "High", "High", "High", "High", "High")

#create Data Track with beta values corresponding to  EPIC probes BetaValues
#included in ESR1 coordinatescfor High andl Low samples
DTrack <- DataTrack(DataMet, 
                    name = "Beta Values",
                    groups= GroupVector,
                    legend = T,
                    type = c("p", "smooth"),
                    aggregateGroups = T,
                    col = GroupColors,
                    lwd=3,
                    ) 



#Plot Completed
plotTracks(list(IdeoTrack,AxisTrack,GeneTrack, DTrack), 
           background.panel = "#FFFEDB",
           background.title = "darkblue",
           cex.title=1.5)












