#### cgMLST comparison ####

# clean environment
rm(list=ls())

# set working directory
setwd("/PathTo/RstudioWorkingDirectory/cgMLSTcomparison")

# install regular packages
install.packages("ggplot2")
install.packages("plyr")
install.packages("reshape2")
install.packages("tidyr")
install.packages("AER")
install.packages("spgs")
install.packages("usethis")

# install devtool packages
install_github("vqv/ggbiplot")

# call library
library(ggplot2)
library(plyr)
library(reshape2)
library(tidyr)
library(AER)
library(spgs)
library(usethis)
library(devtools)
library(ggbiplot)

# read dataframe
data = read.table("Listeria-cgMLST-Additional-file-3.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")

# check dimension
dim(data)
# => [1] 2856   63

# check 10 first lines
head(data, 10)

# check nature of variables (integer or factor)
str(data)

# check and reorganize levels variables
levels(data$reference_strain)
levels(data$sample_origin)
data$sample_origin = factor(data$sample_origin, levels=c("original_DNA", "fifth_culture", "tenth_culture"))
levels(data$sample_origin)
levels(data$DNA_extraction_replicate)
data$DNA_extraction_replicate = factor(data$DNA_extraction_replicate, levels=c("original_DNA", "extraction_A", "extraction_B", "extraction_C"))
levels(data$DNA_extraction_replicate)
levels(data$sequencing_replicate)
levels(data$targeted_depth)
data$targeted_depth = factor(data$targeted_depth, levels=c("Dr100-Dk75", "Dr90-Dk68", "Dr80-Dk60", "Dr70-Dk53", "Dr60-Dk45", "Dr50-Dk38", "Dr40-Dk31", "Dr30-Dk23", "Dr20-Dk16", "Dr10-Dk8"))
levels(data$targeted_depth)
levels(data$workflow)
data$workflow = factor(data$workflow, levels=c("BIGSdb", "INNUENDO", "GENPAT", "SeqSphere", "BioNumericsAB", "BioNumericsAF", "MentaLiST"))
levels(data$workflow)

#### BioNumericsAB versus BioNumericsAF ####

# subset BioNumerics
data_BioNumerics = subset(data,data$workflow %in% c("BioNumericsAB","BioNumericsAF"))
dim(data_BioNumerics)
# => [1] 840  63

# Wilcoxon tests
## switch from long to short dataframe
str(data_BioNumerics)
data_BioNumerics_short <- dcast(data_BioNumerics, formula = sample+reference_strain+targeted_depth~workflow, value.var = "identical_alleles_against_reference")
str(data_BioNumerics_short)
dim(data_BioNumerics_short)
# => [1] 420   5
## Wilcoxon tests for all read depth coverage
wilcox.test(data_BioNumerics_short$BioNumericsAB, data_BioNumerics_short$BioNumericsAF, paired = TRUE)
# => p-value < 2.2e-16

## subset each targeted_depth
levels(data_BioNumerics_short$targeted_depth)
data_BioNumerics_short_Dr100 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr100-Dk75"))
dim(data_BioNumerics_short_Dr100)
# => [1] 42  5
data_BioNumerics_short_Dr90 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr90-Dk68"))
dim(data_BioNumerics_short_Dr90)
# => [1] 42  5
data_BioNumerics_short_Dr80 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr80-Dk60"))
dim(data_BioNumerics_short_Dr80)
# => [1] 42  5
data_BioNumerics_short_Dr70 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr70-Dk53"))
dim(data_BioNumerics_short_Dr70)
# => [1] 42  5
data_BioNumerics_short_Dr60 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr60-Dk45"))
dim(data_BioNumerics_short_Dr60)
# => [1] 42  5
data_BioNumerics_short_Dr50 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr50-Dk38"))
dim(data_BioNumerics_short_Dr50)
# => [1] 42  5
data_BioNumerics_short_Dr40 = subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr40-Dk31"))
dim(data_BioNumerics_short_Dr40)
# => [1] 42  5
data_BioNumerics_short_Dr30= subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr30-Dk23"))
dim(data_BioNumerics_short_Dr30)
# => [1] 42  5
data_BioNumerics_short_Dr20= subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr20-Dk16"))
dim(data_BioNumerics_short_Dr20)
# => [1] 42  5
data_BioNumerics_short_Dr10= subset(data_BioNumerics_short,data_BioNumerics_short$targeted_depth %in% c("Dr10-Dk8"))
dim(data_BioNumerics_short_Dr10)
# => [1] 42  5

## Wilcoxon signed rank tests for each targeted_depth (in case of dependent sampling, which is  the case)
wilcox.test(data_BioNumerics_short_Dr100$BioNumericsAB, data_BioNumerics_short_Dr100$BioNumericsAF, paired = TRUE)
# => p-value = 4.112e-07
wilcox.test(data_BioNumerics_short_Dr90$BioNumericsAB, data_BioNumerics_short_Dr90$BioNumericsAF, paired = TRUE)
# => p-value = 4.184e-07
wilcox.test(data_BioNumerics_short_Dr80$BioNumericsAB, data_BioNumerics_short_Dr80$BioNumericsAF, paired = TRUE)
# => p-value = 1.936e-06
wilcox.test(data_BioNumerics_short_Dr70$BioNumericsAB, data_BioNumerics_short_Dr70$BioNumericsAF, paired = TRUE)
# => p-value = 7.891e-07
wilcox.test(data_BioNumerics_short_Dr60$BioNumericsAB, data_BioNumerics_short_Dr60$BioNumericsAF, paired = TRUE)
# => p-value = 4.012e-06
wilcox.test(data_BioNumerics_short_Dr50$BioNumericsAB, data_BioNumerics_short_Dr50$BioNumericsAF, paired = TRUE)
# => p-value = 2.189e-07
wilcox.test(data_BioNumerics_short_Dr40$BioNumericsAB, data_BioNumerics_short_Dr40$BioNumericsAF, paired = TRUE)
# => p-value = 9.359e-07
wilcox.test(data_BioNumerics_short_Dr30$BioNumericsAB, data_BioNumerics_short_Dr30$BioNumericsAF, paired = TRUE)
# => p-value = 2.649e-07
wilcox.test(data_BioNumerics_short_Dr20$BioNumericsAB, data_BioNumerics_short_Dr20$BioNumericsAF, paired = TRUE)
# => p-value = 9.564e-08
wilcox.test(data_BioNumerics_short_Dr10$BioNumericsAB, data_BioNumerics_short_Dr10$BioNumericsAF, paired = TRUE)
# => p-value = 9.379e-08

# plot identical_alleles_against_reference for BioNumerics
p = ggplot(data = data_BioNumerics, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of BioNumerics cgMLST workflows, workflow-based (AB: n=420) \n or combining workflow-based and -free (AF: n=420) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021).") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("BioNumerics-A-reference.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("BioNumerics-A-reference.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# zoom plot identical_alleles_against_reference for BioNumerics all references merged
p = ggplot(data = data_BioNumerics, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of BioNumerics cgMLST workflows, workflow-based (AB: n=420) \n or combining workflow-based and -free (AF: n=420) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021).") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B") +  
  facet_wrap(~ workflow)
p
plot(p)
ggsave("BioNumerics-B.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("BioNumerics-B.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# zoom plot identical_alleles_against_reference for BioNumerics
p = ggplot(data = data_BioNumerics, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1710, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of BioNumerics cgMLST workflows, workflow-based (AB: n=420) \n or combining workflow-based and -free (AF: n=420) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021).") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("C") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("BioNumerics-C-reference-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("BioNumerics-C-reference-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# zoom plot identical_alleles_against_reference for BioNumerics all references merged
p = ggplot(data = data_BioNumerics, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 8, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1730, 1750), breaks = c(1730,1735,1740,1745,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of BioNumerics cgMLST workflows, workflow-based (AB: n=420) \n or combining workflow-based and -free (AF: n=420) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021).") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("D") +
  facet_wrap(~ workflow)
p
plot(p)
ggsave("BioNumerics-D-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("BioNumerics-D-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# keep BioNumericsAF output for comparison with the other cgMLST workflow
## rename levels of a factor
levels(data$workflow)
data$workflow = revalue(data$workflow, c(
  "BIGSdb" = "BIGSdb", 
  "INNUENDO" = "INNUENDO", 
  "GENPAT" = "GENPAT", 
  "SeqSphere" = "SeqSphere", 
  "BioNumericsAB" = "BioNumericsAB", 
  "BioNumericsAF" = "BioNumerics",
  "MentaLiST" = "MentaLiST"))
levels(data$workflow)
## subset without BioNumericsAB
data_cgMLST = subset(data,data$workflow %in% c("BIGSdb","INNUENDO","GENPAT","SeqSphere","BioNumerics","MentaLiST"))
dim(data)
# => [1] 2856   63
dim(data_cgMLST)
# => [1] 2436   63
levels(data_cgMLST$workflow)
str(data_cgMLST)

#### Principal Component Analysis (PCA) ####

# preparation of the dataframe

## simplify names of 'numerical columns (num and int) with acronyms
colnames(data_cgMLST)
dataACRONYMS  <- data_cgMLST
names(dataACRONYMS)[names(dataACRONYMS) == "BBMap_read_depth"] <- "DEPTH"
names(dataACRONYMS)[names(dataACRONYMS) == "BBMap_read_breadth"] <- "BREADTH"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_0bp"] <- "C0"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_1000bp"] <- "C1000"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_5000bp"] <- "C5000"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_10000bp"] <- "C10000"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_25000bp"] <- "C25000"
names(dataACRONYMS)[names(dataACRONYMS) == "contigs_50000bp"] <- "C50000"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_0bp"] <- "TL0"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_1000bp"] <- "TL1000"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_5000bp"] <- "TL5000"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_10000bp"] <- "TL10000"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_25000bp"] <- "TL25000"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length_50000bp"] <- "TL50000"
names(dataACRONYMS)[names(dataACRONYMS) == "largest_contig"] <- "LC"
names(dataACRONYMS)[names(dataACRONYMS) == "total_length"] <- "TL"
names(dataACRONYMS)[names(dataACRONYMS) == "GC."] <- "GC"
names(dataACRONYMS)[names(dataACRONYMS) == "N50"] <- "N50"
names(dataACRONYMS)[names(dataACRONYMS) == "NG50"] <- "NG50"
names(dataACRONYMS)[names(dataACRONYMS) == "N75"] <- "N75"
names(dataACRONYMS)[names(dataACRONYMS) == "NG75"] <- "NG75"
names(dataACRONYMS)[names(dataACRONYMS) == "L50"] <- "L50"
names(dataACRONYMS)[names(dataACRONYMS) == "LG50"] <- "LG50"
names(dataACRONYMS)[names(dataACRONYMS) == "L75"] <- "L75"
names(dataACRONYMS)[names(dataACRONYMS) == "LG75"] <- "LG75"
names(dataACRONYMS)[names(dataACRONYMS) == "misassemblies"] <- "MA"
names(dataACRONYMS)[names(dataACRONYMS) == "misassembled_contigs"] <- "MAC"
names(dataACRONYMS)[names(dataACRONYMS) == "misassembled_contigs_length"] <- "MACL"
names(dataACRONYMS)[names(dataACRONYMS) == "local_misassemblies"] <- "LMA"
names(dataACRONYMS)[names(dataACRONYMS) == "scaffold_gap_ext_mis"] <- "SQEM"
names(dataACRONYMS)[names(dataACRONYMS) == "scaffold_gap_loc_mis"] <- "SQLM"
names(dataACRONYMS)[names(dataACRONYMS) == "unaligned_mis_contigs"] <- "UAMC"
names(dataACRONYMS)[names(dataACRONYMS) == "unaligned_contigs"] <- "UAC"
names(dataACRONYMS)[names(dataACRONYMS) == "unaligned_contigs_partial"] <- "UACP"
names(dataACRONYMS)[names(dataACRONYMS) == "unaligned_length"] <- "UAL"
names(dataACRONYMS)[names(dataACRONYMS) == "genome_fraction_."] <- "GF"
names(dataACRONYMS)[names(dataACRONYMS) == "Duplication_ratio"] <- "DR"
names(dataACRONYMS)[names(dataACRONYMS) == "Ns_per_100_kbp"] <- "N100"
names(dataACRONYMS)[names(dataACRONYMS) == "mismatches_per_100_kbp"] <- "MM100"
names(dataACRONYMS)[names(dataACRONYMS) == "indels_per_100_kbp"] <- "ID100"
names(dataACRONYMS)[names(dataACRONYMS) == "largest_alignment"] <- "LA"
names(dataACRONYMS)[names(dataACRONYMS) == "total_aligned_length"] <- "TAL"
names(dataACRONYMS)[names(dataACRONYMS) == "NA50"] <- "NA50"
names(dataACRONYMS)[names(dataACRONYMS) == "NGA50"] <- "NGA50"
names(dataACRONYMS)[names(dataACRONYMS) == "NA75"] <- "NA75"
names(dataACRONYMS)[names(dataACRONYMS) == "NGA75"] <- "NGA75"
names(dataACRONYMS)[names(dataACRONYMS) == "LA50"] <- "LA50"
names(dataACRONYMS)[names(dataACRONYMS) == "LGA50"] <- "LGA50"
names(dataACRONYMS)[names(dataACRONYMS) == "LA75"] <- "LA75"
names(dataACRONYMS)[names(dataACRONYMS) == "LGA75"] <- "LGA75"
names(dataACRONYMS)[names(dataACRONYMS) == "unidentified_alleles_against_schema"] <- "UIAAS"
names(dataACRONYMS)[names(dataACRONYMS) == "misidentified_alleles_against_reference"] <- "MIAAR"
names(dataACRONYMS)[names(dataACRONYMS) == "identified_alleles_against_schema"] <- "IAAS"
names(dataACRONYMS)[names(dataACRONYMS) == "identical_alleles_against_reference"] <- "IAAR"
colnames(dataACRONYMS)

## convert the values in a column into row names in a copy of dataframe
datarowname <- dataACRONYMS
row.names(datarowname) <- datarowname$index
datarowname[1] <- NULL
## check absence of index as variable
str(dataACRONYMS)
str(datarowname)
## check row names
rownames(dataACRONYMS)
rownames(datarowname)
## check variables
dimnames(dataACRONYMS)
dimnames(datarowname)

## identify a constant/zero column to not include it in PCA
which(apply(dataACRONYMS, 2, var)==0)
which(apply(datarowname, 2, var)==0)

## remove columns with NAs from the dataframe (i.e. all cgMLST workflows without assembly columns)
datarownameNAcol <- datarowname[ , colSums(is.na(datarowname)) == 0]
dim(datarowname)
# => [1] 2436   62
dim(datarownameNAcol)
# => [1] 2436   14
str(datarownameNAcol)

## remove rows with NAs from the dataframe (i.e. MentaLIST rows with NA in assembly columns)
datarownameNArow = datarowname %>% drop_na()
dim(datarowname)
# => [1] 2436   62
dim(datarownameNArow)
# => [1] 2016   62
str(datarownameNArow)

# PCAs from datarownameNAcol (i.e. all cgMLST workflows without assembly columns)

## get variable numbers
dimnames(datarownameNAcol)
comment <- scan(what="character")
[1] "sample"                  
[2] "reference_strain"        
[3] "sample_origin"           
[4] "DNA_extraction_replicate"
[5] "sequencing_replicate"    
[6] "targeted_depth"          
[7] "DEPTH"                      
[8] "BREADTH"                      
[9] "INNUca_read_depth"       
[10] "workflow"                
[11] "UIAAS"                   
[12] "MIAAR"                   
[13] "IAAS"                    
[14] "IAAR" 

rm(comment)

## compute the principal components
datarownameNAcol.pca <- prcomp(datarownameNAcol[,c(7:8,13:14)], center = TRUE,scale. = TRUE)

## check PCA summary
summary(datarownameNAcol.pca)
comment <- scan(what="character")
Importance of components:
  PC1    PC2    PC3     PC4
Standard deviation     1.3954 1.0102 0.8765 0.51365
Proportion of Variance 0.4868 0.2551 0.1921 0.06596
Cumulative Proportion  0.4868 0.7420 0.9340 1.00000

rm(comment)

## call PCA objects
str(datarownameNAcol.pca)
comment <- scan(what="character")
List of 5
$ sdev    : num [1:4] 1.395 1.01 0.877 0.514
$ rotation: num [1:4, 1:4] 0.627 0.637 0.142 0.426 -0.148 ...
..- attr(*, "dimnames")=List of 2
.. ..$ : chr [1:4] "RD" "BBMap_read_breadth" "IAAS" "IAAR"
.. ..$ : chr [1:4] "PC1" "PC2" "PC3" "PC4"
$ center  : Named num [1:4] 57.7 99.3 1743.8 1683.9
..- attr(*, "names")= chr [1:4] "RD" "BBMap_read_breadth" "IAAS" "IAAR"
$ scale   : Named num [1:4] 28.4108 0.0759 8.6182 272.1378
..- attr(*, "names")= chr [1:4] "RD" "BBMap_read_breadth" "IAAS" "IAAR"
$ x       : num [1:2436, 1:4] -2.326 -1.505 -1.332 0.434 -0.482 ...
..- attr(*, "dimnames")=List of 2
.. ..$ : chr [1:2436] "1" "2" "3" "4" ...
.. ..$ : chr [1:4] "PC1" "PC2" "PC3" "PC4"
- attr(*, "class")= chr "prcomp"

rm(comment)

## plot PCA
ggbiplot(datarownameNAcol.pca)

### add labels (not necessary because too many data)
ggbiplot(datarownameNAcol.pca, labels=rownames(datarownameNAcol))

### add ellipse from dataframe
datarownameNAcol.refence <- datarownameNAcol$reference_strain
ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.refence)

### add circle (not necessary because confusing)
ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  circle=TRUE, groups=datarownameNAcol.refence)

### change of principal components (not necessary because few information)
ggbiplot(datarownameNAcol.pca,ellipse=TRUE, choices=c(3,4), groups=datarownameNAcol.refence)

### scale the samples (obs.scale) and the variables (var.scale) (not necessary because already well defined)
ggbiplot(datarownameNAcol.pca,ellipse=TRUE, obs.scale = 2, var.scale = 1, groups=datarownameNAcol.refence)

### all categorical parameters of interest
#### reference_strain
datarownameNAcol.refence <- datarownameNAcol$reference_strain
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.refence) +
  scale_color_discrete(name = 'reference genomes') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A")
p
plot(p)
ggsave("PCA-no-assembly-A-reference.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-A-reference.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### sample_origin
datarownameNAcol.plating <- datarownameNAcol$sample_origin
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.plating) +
  scale_color_discrete(name = 'successive plating') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B")
p
plot(p)
ggsave("PCA-no-assembly-B-plating.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-B-plating.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### DNA_extraction_replicate
datarownameNAcol.DNA <- datarownameNAcol$DNA_extraction_replicate
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.DNA) +
  scale_color_discrete(name = 'DNA extraction replicate') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("C")
p
plot(p)
ggsave("PCA-no-assembly-C-DNA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-C-DNA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### sequencing_replicate
datarownameNAcol.NGS <- datarownameNAcol$sequencing_replicate
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.NGS) +
  scale_color_discrete(name = 'sequencing replicate') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("D")
p
plot(p)
ggsave("PCA-no-assembly-D-NGS.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-D-NGS.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### targeted_depth
datarownameNAcol.DrDk <- datarownameNAcol$targeted_depth
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.DrDk) +
  scale_color_discrete(name = 'targeted depth') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("E")
p
plot(p)
ggsave("PCA-no-assembly-E-DrDk.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-E-DrDk.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow
datarownameNAcol.workflow <- datarownameNAcol$workflow
p=ggbiplot(datarownameNAcol.pca,ellipse=TRUE,  groups=datarownameNAcol.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("F")
p
plot(p)
ggsave("PCA-no-assembly-F-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-F-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow focusing on "INNUENDO","GENPAT","SeqSphere"
dim(datarownameNAcol)
dataSubset = subset(datarownameNAcol,datarownameNAcol$workflow %in% c("INNUENDO","GENPAT","SeqSphere"))
dim(data_noSeqSphere)
dataSubset.pca <- prcomp(dataSubset[,c(7:8,13:14)], center = TRUE,scale. = TRUE)
dataSubset.workflow <- dataSubset$workflow
p=ggbiplot(dataSubset.pca,ellipse=TRUE,  groups=dataSubset.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("G")
p
plot(p)
ggsave("PCA-no-assembly-G-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-G-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow focusing on "BIGSdb","INNUENDO","GENPAT"
dim(datarownameNAcol)
dataSubset = subset(datarownameNAcol,datarownameNAcol$workflow %in% c("BIGSdb","INNUENDO","GENPAT"))
dim(data_noSeqSphere)
dataSubset.pca <- prcomp(dataSubset[,c(7:8,13:14)], center = TRUE,scale. = TRUE)
dataSubset.workflow <- dataSubset$workflow
p=ggbiplot(dataSubset.pca,ellipse=TRUE,  groups=dataSubset.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("H")
p
plot(p)
ggsave("PCA-no-assembly-H-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-assembly-H-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# PCAs from datarownameNArow (i.e. all cgMLST workflows excepted MentaLIST which does not present assembly output)

## get variable numbers
dimnames(datarownameNArow)
comment <- scan(what="character")
[1] "sample"                  
[2] "reference_strain"        
[3] "sample_origin"           
[4] "DNA_extraction_replicate"
[5] "sequencing_replicate"    
[6] "targeted_depth"          
[7] "DEPTH"                      
[8] "BREADTH"                      
[9] "INNUca_read_depth"       
[10] "workflow"                
[11] "C0"                      
[12] "C1000"                   
[13] "C5000"                   
[14] "C10000"                  
[15] "C25000"                  
[16] "C50000"                  
[17] "TL0"                     
[18] "TL1000"                  
[19] "TL5000"                  
[20] "TL10000"                 
[21] "TL25000"                 
[22] "TL50000"                 
[23] "LC"                      
[24] "TL"                      
[25] "GC"                      
[26] "N50"                     
[27] "NG50"                    
[28] "N75"                     
[29] "NG75"                    
[30] "L50"                     
[31] "LG50"                    
[32] "L75"                     
[33] "LG75"                    
[34] "MA"                      
[35] "MAC"                     
[36] "MACL"                    
[37] "LMA"                     
[38] "SQEM"                    
[39] "SQLM"                    
[40] "UAMC"                    
[41] "UAC"                     
[42] "UACP"                    
[43] "UAL"                     
[44] "GF"                      
[45] "DR"                      
[46] "N100"                    
[47] "MM100"                   
[48] "ID100"                   
[49] "LA"                      
[50] "TAL"                     
[51] "NA50"                    
[52] "NGA50"                   
[53] "NA75"                    
[54] "NGA75"                   
[55] "LA50"                    
[56] "LGA50"                   
[57] "LA75"                    
[58] "LGA75"                   
[59] "UIAAS"                   
[60] "MIAAR"                   
[61] "IAAS"                    
[62] "IAAR"

rm(comment)

## compute the principal components
datarownameNArow.pca <- prcomp(datarownameNArow[,c(7:8,11:58,61:62)], center = TRUE,scale. = TRUE)

## check PCA summary
summary(datarownameNArow.pca)

## call PCA objects
str(datarownameNArow.pca)

## plot PCA
p=ggbiplot(datarownameNArow.pca) +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("All numerical parameters")
p
plot(p)
ggsave("PCA-no-MentaLIST-all-parameters.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-all-parameters.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

## grouped manually numerical parameters

### group 1
comment <- scan(what="character")
[11] "C0"                      
[12] "C1000"                   
[13] "C5000"                   
[14] "C10000"                  
[15] "C25000"                  
[16] "C50000" 

rm(comment)

### group 2
comment <- scan(what="character")
[25] "GC" 
[17] "TL0"                     
[18] "TL1000"                  
[19] "TL5000"                  
[20] "TL10000"                 
[21] "TL25000"                 
[22] "TL50000"  
[24] "TL"
[50] "TAL" 
[36] "MACL"  

rm(comment)

### group 3
comment <- scan(what="character")
[26] "N50"                     
[27] "NG50"                    
[28] "N75"                     
[29] "NG75" 
[38] "SQEM" 
[51] "NA50" 
[52] "NGA50"                   
[53] "NA75"                    
[54] "NGA75"  
[49] "LA"  

rm(comment)

### group 4
comment <- scan(what="character")
[30] "L50"                     
[31] "LG50"                    
[32] "L75"                     
[33] "LG75"

rm(comment)

### group 5
comment <- scan(what="character")
[55] "LA50"                    
[56] "LGA50"                   
[57] "LA75"   
[58] "LGA75" 

rm(comment)

### group 6
comment <- scan(what="character")
[7] "DEPTH" 
[44] "GF"  

rm(comment)

### group 7
comment <- scan(what="character")
[37] "LMA"
[43] "UAL"
[47] "MM100"  
[39] "SQLM" 

rm(comment)

### group 8
comment <- scan(what="character")
[45] "DR"  
[46] "N100"  
[41] "UAC" 

rm(comment)

### group 9
comment <- scan(what="character")
[34] "MA"  
[35] "MAC"

rm(comment)

### ungrouped
comment <- scan(what="character")
[8] "BREADTH"       
[23] "LC"                      
[40] "UAMC"                    
[42] "UACP"                    
[48] "ID100"                   
[61] "IAAS"                    
[62] "IAAR"

rm(comment)

## compute the principal components
datarownameNArow.pca.subset <- prcomp(datarownameNArow[,c(7,8,12,14,18,20,25,26,30,34,37,38,39,40,42,45,48,49,51,55,61,62)], center = TRUE,scale. = TRUE)

## plot subset PCA
p=ggbiplot(datarownameNArow.pca.subset) +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("Subset numerical parameters")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

### all categorical parameters of interest
#### reference_strain
datarownameNArow.refence <- datarownameNArow$reference_strain
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.refence) +
  scale_color_discrete(name = 'reference genomes') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-A-reference.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-A-reference.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### sample_origin
datarownameNArow.plating <- datarownameNArow$sample_origin
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.plating) +
  scale_color_discrete(name = 'successive plating') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-B-plating.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-B-plating.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### DNA_extraction_replicate
datarownameNArow.DNA <- datarownameNArow$DNA_extraction_replicate
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.DNA) +
  scale_color_discrete(name = 'DNA extraction replicate') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("C")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-C-DNA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-C-DNA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### sequencing_replicate
datarownameNArow.NGS <- datarownameNArow$sequencing_replicate
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.NGS) +
  scale_color_discrete(name = 'sequencing replicate') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("D")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-D-NGS.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-D-NGS.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### targeted_depth
datarownameNArow.DrDk <- datarownameNArow$targeted_depth
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.DrDk) +
  scale_color_discrete(name = 'targeted depth') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("E")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-E-DrDk.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-E-DrDk.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow
datarownameNArow.workflow <- datarownameNArow$workflow
p=ggbiplot(datarownameNArow.pca.subset,ellipse=TRUE,  groups=datarownameNArow.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("F")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-F-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-F-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow focusing on "INNUENDO","GENPAT","SeqSphere"
dim(datarownameNArow)
dataSubset = subset(datarownameNArow,datarownameNArow$workflow %in% c("INNUENDO","GENPAT","SeqSphere"))
dim(data_noSeqSphere)
dim(dataSubset)
which(apply(dataSubset, 2, var)==0)
dataSubset.pca.subset <- prcomp(dataSubset[,c(7,8,12,14,18,20,25,26,30,34,37,38,39,42,45,48,49,51,55,61,62)], center = TRUE,scale. = TRUE)
dataSubset.workflow <- dataSubset$workflow
p=ggbiplot(dataSubset.pca.subset,ellipse=TRUE,  groups=dataSubset.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("G")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-G-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-G-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()
#### workflow focusing on "BIGSdb","INNUENDO","GENPAT"
dim(datarownameNArow)
dataSubset = subset(datarownameNArow,datarownameNArow$workflow %in% c("BIGSdb","INNUENDO","GENPAT"))
dim(data_noSeqSphere)
dim(dataSubset)
which(apply(dataSubset, 2, var)==0)
dataSubset.pca.subset <- prcomp(dataSubset[,c(7,8,12,14,18,20,25,26,30,34,37,38,39,42,45,48,49,51,55,61,62)], center = TRUE,scale. = TRUE)
dataSubset.workflow <- dataSubset$workflow
p=ggbiplot(dataSubset.pca.subset,ellipse=TRUE,  groups=dataSubset.workflow) +
  scale_color_discrete(name = 'cgMLST workflows') +
  theme(legend.direction = 'vertical', legend.position = 'right') +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("H")
p
plot(p)
ggsave("PCA-no-MentaLIST-subset-parameters-H-cgMLST.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("PCA-no-MentaLIST-subset-parameters-H-cgMLST.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()         

# GLMs from datarownameNAcol (i.e. all cgMLST workflows without assembly columns)
## test distribution of the variable to explain
### Gaussian
shapiro.test(datarownameNAcol$IAAR)
# => W = 0.2398, p-value < 2.2e-16
# => different than gaussian
### Poisson
poisson.test(x = sum(datarownameNAcol$IAAR), T = length(datarownameNAcol$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 4102009, time base = 2436, p-value < 2.2e-16
# => different than two side poisson distribution
poisson.test(x = sum(datarownameNAcol$IAAR), T = length(datarownameNAcol$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 4102009, time base = 2436, p-value = 1
# => one side poisson distribution with lower hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≥ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 < 0)
poisson.test(x = sum(datarownameNAcol$IAAR), T = length(datarownameNAcol$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 4102009, time base = 2436, p-value < 2.2e-16
# => different than one side poisson distribution with upper hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≤ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 > 0)
### uniform
chisq.test(datarownameNAcol$IAAR)
# => X-squared = 107092, df = 2435, p-value < 2.2e-16
# => different than uniform
## check distribution visualy
hist(datarownameNAcol$IAAR)

## rename variables for GLM
datarownameNAcolGLM <- datarownameNAcol
colnames(datarownameNAcolGLM)
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "reference_strain"] <- "REFERENCE"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "sample_origin"] <- "PLATING"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "DNA_extraction_replicate"] <- "DNA"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "sequencing_replicate"] <- "SEQUENCING"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "targeted_depth"] <- "DrDk"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "workflow"] <- "WORKFLOW"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "DEPTH"] <- "DEPTH"
names(datarownameNAcolGLM)[names(datarownameNAcolGLM) == "BREADTH"] <- "BREADTH"
colnames(datarownameNAcolGLM)
str(datarownameNAcolGLM)

## select the variable to explain and the variables of interest (do not add DrDk because already linear with DEPTH)
formulaAddition <- IAAR~REFERENCE+PLATING+DNA+SEQUENCING+DEPTH+BREADTH+WORKFLOW+IAAS
formulaMultiplication <- IAAR~REFERENCE*PLATING*DNA*SEQUENCING*DEPTH*BREADTH*WORKFLOW*IAAS

## run GLM
GLM1 <- glm(formulaAddition, data=datarownameNAcolGLM, family = poisson)
GLM2 <- glm(formulaAddition, data=datarownameNAcolGLM, family = poisson(link = "log"))
GLM3 <- glm(formulaAddition, data=datarownameNAcolGLM, family = quasipoisson)
GLM4 <- glm(formulaMultiplication, data=datarownameNAcolGLM, family = poisson)
# => glm.fit: algorithm did not converge
GLM5 <- glm(formulaMultiplication, data=datarownameNAcolGLM, family = poisson(link = "log"))
# => glm.fit: algorithm did not converge
GLM6 <- glm(formulaMultiplication, data=datarownameNAcolGLM, family = quasipoisson)
# => glm.fit: algorithm did not converge

### test over dispersion  (i.e. alpha > 0 with p<5%)
dispersiontest(GLM1,trafo=1)
# => z = 11.069, p-value < 2.2e-16
# => sample estimates: alpha 39.17442
dispersiontest(GLM2,trafo=1)
# => z = 11.069, p-value < 2.2e-16
# => sample estimates: alpha 39.17442
dispersiontest(GLM3,trafo=1)
# => only Poisson GLMs can be tested
# => presence of  GLM overdispersion => keep quasipoisson (GLM3)

## check GLM results

summary(GLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAddition, family = poisson, data = datarownameNAcolGLM)
Deviance Residuals: 
  Min       1Q   Median       3Q      Max  
-51.924   -1.574    0.446    2.844    9.400  
Coefficients: (1 not defined because of singularities)
Estimate Std. Error  z value Pr(>|z|)    
(Intercept)          -4.939e+01  1.270e+00  -38.874  < 2e-16 ***
REFERENCEATCC19115   -3.015e-02  1.590e-03  -18.959  < 2e-16 ***
REFERENCEATCCBAA679  -3.815e-02  1.675e-03  -22.772  < 2e-16 ***
PLATINGfifth_culture -1.767e-02  1.726e-03  -10.237  < 2e-16 ***
PLATINGtenth_culture -2.394e-02  1.760e-03  -13.604  < 2e-16 ***
DNAextraction_A      -7.447e-03  1.318e-03   -5.652 1.59e-08 ***
DNAextraction_B      -4.528e-03  1.311e-03   -3.454 0.000552 ***
DNAextraction_C              NA         NA       NA       NA    
SEQUENCINGNextSeq_B  -2.256e-02  1.100e-03  -20.508  < 2e-16 ***
DEPTH                 8.124e-04  2.942e-05   27.613  < 2e-16 ***
BREADTH               5.613e-01  1.273e-02   44.082  < 2e-16 ***
WORKFLOWINNUENDO     -2.005e-02  1.768e-03  -11.342  < 2e-16 ***
WORKFLOWGENPAT        1.680e-04  1.660e-03    0.101 0.919392    
WORKFLOWSeqSphere     4.500e-05  1.652e-03    0.027 0.978260    
WORKFLOWBioNumerics  -2.872e-02  1.884e-03  -15.243  < 2e-16 ***
WORKFLOWMentaLiST    -1.851e-01  1.733e-03 -106.763  < 2e-16 ***
IAAS                  6.301e-04  8.098e-05    7.781 7.21e-15 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 185666  on 2435  degrees of freedom
Residual deviance: 154801  on 2420  degrees of freedom
AIC: 177133
Number of Fisher Scoring iterations: 4

rm(comment)

summary(GLM2)
comment <- scan(what="character")
Call: glm(formula = formulaAddition, family = poisson(link = "log"), data = datarownameNAcolGLM)
Deviance Residuals: 
  Min       1Q   Median       3Q      Max  
-51.924   -1.574    0.446    2.844    9.400  
Coefficients: (1 not defined because of singularities)
Estimate Std. Error  z value Pr(>|z|)    
(Intercept)          -4.939e+01  1.270e+00  -38.874  < 2e-16 ***
REFERENCEATCC19115   -3.015e-02  1.590e-03  -18.959  < 2e-16 ***
REFERENCEATCCBAA679  -3.815e-02  1.675e-03  -22.772  < 2e-16 ***
PLATINGfifth_culture -1.767e-02  1.726e-03  -10.237  < 2e-16 ***
PLATINGtenth_culture -2.394e-02  1.760e-03  -13.604  < 2e-16 ***
DNAextraction_A      -7.447e-03  1.318e-03   -5.652 1.59e-08 ***
DNAextraction_B      -4.528e-03  1.311e-03   -3.454 0.000552 ***
DNAextraction_C              NA         NA       NA       NA    
SEQUENCINGNextSeq_B  -2.256e-02  1.100e-03  -20.508  < 2e-16 ***
DEPTH                 8.124e-04  2.942e-05   27.613  < 2e-16 ***
BREADTH               5.613e-01  1.273e-02   44.082  < 2e-16 ***
WORKFLOWINNUENDO     -2.005e-02  1.768e-03  -11.342  < 2e-16 ***
WORKFLOWGENPAT        1.680e-04  1.660e-03    0.101 0.919392    
WORKFLOWSeqSphere     4.500e-05  1.652e-03    0.027 0.978260    
WORKFLOWBioNumerics  -2.872e-02  1.884e-03  -15.243  < 2e-16 ***
WORKFLOWMentaLiST    -1.851e-01  1.733e-03 -106.763  < 2e-16 ***
IAAS                  6.301e-04  8.098e-05    7.781 7.21e-15 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 185666  on 2435  degrees of freedom
Residual deviance: 154801  on 2420  degrees of freedom
AIC: 177133
Number of Fisher Scoring iterations: 4

rm(comment)

summary(GLM3)
comment <- scan(what="character")
Call: glm(formula = formulaAddition, family = quasipoisson, data = datarownameNAcolGLM)
Deviance Residuals: 
  Min       1Q   Median       3Q      Max  
-51.924   -1.574    0.446    2.844    9.400  
Coefficients: (1 not defined because of singularities)
Estimate Std. Error t value Pr(>|t|)    
(Intercept)          -4.939e+01  8.079e+00  -6.113 1.14e-09 ***
REFERENCEATCC19115   -3.015e-02  1.011e-02  -2.981 0.002899 ** 
REFERENCEATCCBAA679  -3.815e-02  1.065e-02  -3.581 0.000349 ***
PLATINGfifth_culture -1.767e-02  1.097e-02  -1.610 0.107568    
PLATINGtenth_culture -2.394e-02  1.119e-02  -2.139 0.032522 *  
DNAextraction_A      -7.447e-03  8.379e-03  -0.889 0.374225    
DNAextraction_B      -4.528e-03  8.335e-03  -0.543 0.587039    
DNAextraction_C              NA         NA      NA       NA    
SEQUENCINGNextSeq_B  -2.256e-02  6.995e-03  -3.225 0.001277 ** 
DEPTH                 8.124e-04  1.871e-04   4.342 1.47e-05 ***
BREADTH               5.613e-01  8.097e-02   6.932 5.31e-12 ***
WORKFLOWINNUENDO     -2.005e-02  1.124e-02  -1.784 0.074616 .  
WORKFLOWGENPAT        1.680e-04  1.056e-02   0.016 0.987305    
WORKFLOWSeqSphere     4.500e-05  1.050e-02   0.004 0.996581    
WORKFLOWBioNumerics  -2.872e-02  1.198e-02  -2.397 0.016605 *  
WORKFLOWMentaLiST    -1.851e-01  1.102e-02 -16.789  < 2e-16 ***
IAAS                  6.301e-04  5.150e-04   1.224 0.221255    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for quasipoisson family taken to be 40.44085)
Null deviance: 185666  on 2435  degrees of freedom
Residual deviance: 154801  on 2420  degrees of freedom
AIC: NA
Number of Fisher Scoring iterations: 4

rm(comment)

## save GLM results
write.csv(summary(GLM1)['coefficients'], file="GLM1.csv")
write.csv(summary(GLM2)['coefficients'], file="GLM2.csv")
write.csv(summary(GLM3)['coefficients'], file="GLM3.csv")

# GLMs from datarownameNArow (i.e. all cgMLST workflows excepted MentaLIST which does not present assembly output)
## test distribution of the variable to explain
### Gaussian
shapiro.test(datarownameNArow$IAAR)
# => W = 0.18834, p-value < 2.2e-16
# => different than gaussian
### Poisson
poisson.test(x = sum(datarownameNArow$IAAR), T = length(datarownameNArow$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 3492649, time base = 2016, p-value < 2.2e-16
# => different than two side poisson distribution
poisson.test(x = sum(datarownameNArow$IAAR), T = length(datarownameNArow$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 3492649, time base = 2016, p-value = 1
# => one side poisson distribution with lower hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≥ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 < 0)
poisson.test(x = sum(datarownameNArow$IAAR), T = length(datarownameNArow$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 3492649, time base = 2016, p-value < 2.2e-16
# => different than one side poisson distribution with upper hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≤ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 > 0)
### uniform
chisq.test(datarownameNArow$IAAR)
# => X-squared = 5913.9, df = 2015, p-value < 2.2e-16
# => different than uniform
## check distribution visualy
hist(datarownameNArow$IAAR)

## rename variables for GLM
datarownameNArowGLM <- datarownameNArow
colnames(datarownameNArowGLM)
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "reference_strain"] <- "REFERENCE"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "sample_origin"] <- "PLATING"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "DNA_extraction_replicate"] <- "DNA"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "sequencing_replicate"] <- "SEQUENCING"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "targeted_depth"] <- "DrDk"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "workflow"] <- "WORKFLOW"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "RD"] <- "DEPTH"
names(datarownameNArowGLM)[names(datarownameNArowGLM) == "RB"] <- "BREADTH"
colnames(datarownameNArowGLM)
str(datarownameNArowGLM)

## select the variable to explain and the variables of interest (do not add DrDk because already linear with others)
formulaAdditionAssembly <- IAAR~REFERENCE+PLATING+DNA+SEQUENCING+DEPTH+BREADTH+WORKFLOW+IAAS+C0+C1000+C5000+C10000+C25000+C50000+TL0+TL1000+TL5000+TL10000+TL25000+TL50000+LC+TL+GC+N50+NG50+N75+NG75+L50+LG50+L75+LG75+MA+MAC+MACL+LMA+SQEM+SQLM+UAMC+UAC+UACP+UAL+GF+DR+N100+MM100+ID100+LA+TAL+NA50+NGA50+NA75+NGA75+LA50+LGA50+LA75+LGA75
formulaMultiplicationAssembly <- IAAR~REFERENCE*PLATING*DNA*SEQUENCING*DEPTH*BREADTH*WORKFLOW*IAAS**C0*C1000*C5000*C10000*C25000*C50000*TL0*TL1000*TL5000*TL10000*TL25000*TL50000*LC*TL*GC*N50*NG50*N75*NG75*L50*LG50*L75*LG75*MA*MAC*MACL*LMA*SQEM*SQLM*UAMC*UAC*UACP*UAL*GF*DR*N100*MM100*ID100*LA*TAL*NA50*NGA50*NA75*NGA75*LA50*LGA50*LA75*LGA75

## run GLM
AssemblyGLM1 <- glm(formulaAdditionAssembly, data=datarownameNArowGLM, family = poisson)
AssemblyGLM2 <- glm(formulaAdditionAssembly, data=datarownameNArowGLM, family = poisson(link = "log"))
AssemblyGLM3 <- glm(formulaAdditionAssembly, data=datarownameNArowGLM, family = quasipoisson)
AssemblyGLM4 <- glm(formulaMultiplicationAssembly, data=datarownameNArowGLM, family = poisson)
# => glm.fit: algorithm did not converge
AssemblyGLM5 <- glm(formulaMultiplicationAssembly, data=datarownameNArowGLM, family = poisson(link = "log"))
# => glm.fit: algorithm did not converge
AssemblyGLM6 <- glm(formulaMultiplicationAssembly, data=datarownameNArowGLM, family = quasipoisson)
# => glm.fit: algorithm did not converge

### test over dispersion (i.e. alpha > 0 with p<5%)
dispersiontest(AssemblyGLM1,trafo=1)
# => z = -454.29, p-value = 1
# => sample estimates: alpha -0.9875226 
dispersiontest(AssemblyGLM2,trafo=1)
# => z = -454.29, p-value = 1
# => sample estimates: alpha -0.9875226 
dispersiontest(AssemblyGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of GLM overdispersion => keep poisson (AssemblyGLM1)

## check GLM results

summary(AssemblyGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssembly, family = poisson, data = datarownameNArowGLM)
Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-1.37988  -0.03322   0.00203   0.03558   1.51349  
Coefficients: (2 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)    
(Intercept)           5.343e+00  5.993e+00   0.892  0.37262    
REFERENCEATCC19115    3.125e-02  4.366e-02   0.716  0.47414    
REFERENCEATCCBAA679   3.784e-02  5.153e-02   0.734  0.46283    
PLATINGfifth_culture  6.164e-05  1.921e-03   0.032  0.97440    
PLATINGtenth_culture  9.518e-05  1.973e-03   0.048  0.96153    
DNAextraction_A       1.432e-04  1.468e-03   0.098  0.92231    
DNAextraction_B       1.176e-05  1.442e-03   0.008  0.99349    
DNAextraction_C              NA         NA      NA       NA    
SEQUENCINGNextSeq_B  -2.146e-04  1.245e-03  -0.172  0.86317    
DEPTH                -4.393e-06  3.533e-05  -0.124  0.90104    
BREADTH               4.290e-03  1.486e-02   0.289  0.77284    
WORKFLOWINNUENDO     -2.036e-03  3.454e-03  -0.589  0.55560    
WORKFLOWGENPAT       -2.665e-03  2.717e-03  -0.981  0.32656    
WORKFLOWSeqSphere    -3.396e-03  6.904e-03  -0.492  0.62285    
WORKFLOWBioNumerics  -5.518e-05  3.714e-03  -0.015  0.98815    
IAAS                  5.810e-04  1.256e-04   4.627 3.71e-06 ***
C0                   -2.787e-04  2.585e-04  -1.078  0.28102    
C1000                 1.129e-03  1.438e-03   0.785  0.43238    
C5000                -9.297e-04  5.988e-03  -0.155  0.87662    
C10000               -7.646e-04  7.817e-03  -0.098  0.92209    
C25000               -6.099e-04  8.036e-03  -0.076  0.93950    
C50000                1.248e-03  4.911e-03   0.254  0.79946    
TL0                   6.088e-06  2.275e-06   2.676  0.00746 ** 
TL1000               -1.202e-06  1.267e-06  -0.949  0.34284    
TL5000                1.116e-07  1.158e-06   0.096  0.92328    
TL10000              -3.424e-08  9.753e-07  -0.035  0.97200    
TL25000               6.223e-09  3.656e-07   0.017  0.98642    
TL50000              -4.691e-08  1.313e-07  -0.357  0.72094    
LC                   -1.443e-08  5.600e-08  -0.258  0.79666    
TL                           NA         NA      NA       NA    
GC                   -7.733e-02  1.319e-01  -0.586  0.55756    
N50                   3.594e-08  7.767e-08   0.463  0.64354    
NG50                 -5.392e-08  9.261e-08  -0.582  0.56042    
N75                  -8.383e-08  1.900e-07  -0.441  0.65907    
NG75                  7.817e-08  1.916e-07   0.408  0.68337    
L50                   5.598e-03  1.249e-02   0.448  0.65388    
LG50                 -6.861e-03  1.270e-02  -0.540  0.58911    
L75                  -2.874e-03  1.014e-02  -0.283  0.77698    
LG75                  2.585e-03  9.505e-03   0.272  0.78568    
MA                   -5.408e-03  9.188e-03  -0.589  0.55615    
MAC                   6.786e-03  1.049e-02   0.647  0.51785    
MACL                 -1.304e-09  9.470e-09  -0.138  0.89047    
LMA                  -2.237e-05  7.345e-04  -0.030  0.97571    
SQEM                 -6.377e-05  3.750e-03  -0.017  0.98643    
SQLM                 -1.817e-05  7.703e-04  -0.024  0.98119    
UAMC                 -2.647e-03  1.730e-02  -0.153  0.87836    
UAC                   6.976e-04  6.603e-04   1.056  0.29075    
UACP                  4.574e-03  4.425e-03   1.034  0.30133    
UAL                  -6.501e-06  2.425e-06  -2.681  0.00734 ** 
GF                   -8.626e-03  3.350e-02  -0.258  0.79678    
DR                   -3.048e-01  3.068e+00  -0.099  0.92087    
N100                 -3.284e-04  1.827e-05 -17.970  < 2e-16 ***
MM100                -2.434e-04  1.190e-03  -0.204  0.83801    
ID100                 3.193e-03  6.246e-03   0.511  0.60918    
LA                    1.630e-08  5.580e-08   0.292  0.77022    
TAL                  -3.298e-06  1.806e-06  -1.826  0.06786 .  
NA50                 -6.894e-09  1.033e-07  -0.067  0.94678    
NGA50                 2.060e-08  1.139e-07   0.181  0.85647    
NA75                  1.486e-08  1.088e-07   0.137  0.89142    
NGA75                -1.093e-08  1.019e-07  -0.107  0.91462    
LA50                 -3.245e-03  1.308e-02  -0.248  0.80414    
LGA50                 3.945e-03  1.294e-02   0.305  0.76052    
LA75                  1.266e-03  8.574e-03   0.148  0.88265    
LGA75                -1.412e-03  7.644e-03  -0.185  0.85342    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 6542.130  on 2015  degrees of freedom
Residual deviance:   25.142  on 1954  degrees of freedom
AIC: 18886
Number of Fisher Scoring iterations: 3

rm(comment)

summary(AssemblyGLM2)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssembly, family = poisson(link = "log"),data = datarownameNArowGLM)
Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-1.37988  -0.03322   0.00203   0.03558   1.51349  
Coefficients: (2 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)    
(Intercept)           5.343e+00  5.993e+00   0.892  0.37262    
REFERENCEATCC19115    3.125e-02  4.366e-02   0.716  0.47414    
REFERENCEATCCBAA679   3.784e-02  5.153e-02   0.734  0.46283    
PLATINGfifth_culture  6.164e-05  1.921e-03   0.032  0.97440    
PLATINGtenth_culture  9.518e-05  1.973e-03   0.048  0.96153    
DNAextraction_A       1.432e-04  1.468e-03   0.098  0.92231    
DNAextraction_B       1.176e-05  1.442e-03   0.008  0.99349    
DNAextraction_C              NA         NA      NA       NA    
SEQUENCINGNextSeq_B  -2.146e-04  1.245e-03  -0.172  0.86317    
DEPTH                -4.393e-06  3.533e-05  -0.124  0.90104    
BREADTH               4.290e-03  1.486e-02   0.289  0.77284    
WORKFLOWINNUENDO     -2.036e-03  3.454e-03  -0.589  0.55560    
WORKFLOWGENPAT       -2.665e-03  2.717e-03  -0.981  0.32656    
WORKFLOWSeqSphere    -3.396e-03  6.904e-03  -0.492  0.62285    
WORKFLOWBioNumerics  -5.518e-05  3.714e-03  -0.015  0.98815    
IAAS                  5.810e-04  1.256e-04   4.627 3.71e-06 ***
C0                   -2.787e-04  2.585e-04  -1.078  0.28102    
C1000                 1.129e-03  1.438e-03   0.785  0.43238    
C5000                -9.297e-04  5.988e-03  -0.155  0.87662    
C10000               -7.646e-04  7.817e-03  -0.098  0.92209    
C25000               -6.099e-04  8.036e-03  -0.076  0.93950    
C50000                1.248e-03  4.911e-03   0.254  0.79946    
TL0                   6.088e-06  2.275e-06   2.676  0.00746 ** 
TL1000               -1.202e-06  1.267e-06  -0.949  0.34284    
TL5000                1.116e-07  1.158e-06   0.096  0.92328    
TL10000              -3.424e-08  9.753e-07  -0.035  0.97200    
TL25000               6.223e-09  3.656e-07   0.017  0.98642    
TL50000              -4.691e-08  1.313e-07  -0.357  0.72094    
LC                   -1.443e-08  5.600e-08  -0.258  0.79666    
TL                           NA         NA      NA       NA    
GC                   -7.733e-02  1.319e-01  -0.586  0.55756    
N50                   3.594e-08  7.767e-08   0.463  0.64354    
NG50                 -5.392e-08  9.261e-08  -0.582  0.56042    
N75                  -8.383e-08  1.900e-07  -0.441  0.65907    
NG75                  7.817e-08  1.916e-07   0.408  0.68337    
L50                   5.598e-03  1.249e-02   0.448  0.65388    
LG50                 -6.861e-03  1.270e-02  -0.540  0.58911    
L75                  -2.874e-03  1.014e-02  -0.283  0.77698    
LG75                  2.585e-03  9.505e-03   0.272  0.78568    
MA                   -5.408e-03  9.188e-03  -0.589  0.55615    
MAC                   6.786e-03  1.049e-02   0.647  0.51785    
MACL                 -1.304e-09  9.470e-09  -0.138  0.89047    
LMA                  -2.237e-05  7.345e-04  -0.030  0.97571    
SQEM                 -6.377e-05  3.750e-03  -0.017  0.98643    
SQLM                 -1.817e-05  7.703e-04  -0.024  0.98119    
UAMC                 -2.647e-03  1.730e-02  -0.153  0.87836    
UAC                   6.976e-04  6.603e-04   1.056  0.29075    
UACP                  4.574e-03  4.425e-03   1.034  0.30133    
UAL                  -6.501e-06  2.425e-06  -2.681  0.00734 ** 
GF                   -8.626e-03  3.350e-02  -0.258  0.79678    
DR                   -3.048e-01  3.068e+00  -0.099  0.92087    
N100                 -3.284e-04  1.827e-05 -17.970  < 2e-16 ***
MM100                -2.434e-04  1.190e-03  -0.204  0.83801    
ID100                 3.193e-03  6.246e-03   0.511  0.60918    
LA                    1.630e-08  5.580e-08   0.292  0.77022    
TAL                  -3.298e-06  1.806e-06  -1.826  0.06786 .  
NA50                 -6.894e-09  1.033e-07  -0.067  0.94678    
NGA50                 2.060e-08  1.139e-07   0.181  0.85647    
NA75                  1.486e-08  1.088e-07   0.137  0.89142    
NGA75                -1.093e-08  1.019e-07  -0.107  0.91462    
LA50                 -3.245e-03  1.308e-02  -0.248  0.80414    
LGA50                 3.945e-03  1.294e-02   0.305  0.76052    
LA75                  1.266e-03  8.574e-03   0.148  0.88265    
LGA75                -1.412e-03  7.644e-03  -0.185  0.85342    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 6542.130  on 2015  degrees of freedom
Residual deviance:   25.142  on 1954  degrees of freedom
AIC: 18886
Number of Fisher Scoring iterations: 3

rm(comment)

summary(AssemblyGLM3)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssembly, family = quasipoisson,data = datarownameNArowGLM)
Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-1.37988  -0.03322   0.00203   0.03558   1.51349  
Coefficients: (2 not defined because of singularities)
Estimate Std. Error  t value Pr(>|t|)    
(Intercept)           5.343e+00  6.799e-01    7.858 6.40e-15 ***
REFERENCEATCC19115    3.125e-02  4.954e-03    6.308 3.48e-10 ***
REFERENCEATCCBAA679   3.784e-02  5.847e-03    6.471 1.23e-10 ***
PLATINGfifth_culture  6.164e-05  2.179e-04    0.283 0.777307    
PLATINGtenth_culture  9.518e-05  2.239e-04    0.425 0.670800    
DNAextraction_A       1.432e-04  1.666e-04    0.860 0.390146    
DNAextraction_B       1.176e-05  1.636e-04    0.072 0.942675    
DNAextraction_C              NA         NA       NA       NA    
SEQUENCINGNextSeq_B  -2.146e-04  1.413e-04   -1.519 0.128944    
DEPTH                -4.393e-06  4.008e-06   -1.096 0.273244    
BREADTH               4.290e-03  1.686e-03    2.544 0.011034 *  
WORKFLOWINNUENDO     -2.036e-03  3.919e-04   -5.195 2.26e-07 ***
WORKFLOWGENPAT       -2.665e-03  3.082e-04   -8.647  < 2e-16 ***
WORKFLOWSeqSphere    -3.396e-03  7.834e-04   -4.335 1.53e-05 ***
WORKFLOWBioNumerics  -5.518e-05  4.214e-04   -0.131 0.895829    
IAAS                  5.810e-04  1.425e-05   40.781  < 2e-16 ***
C0                   -2.787e-04  2.933e-05   -9.501  < 2e-16 ***
C1000                 1.129e-03  1.631e-04    6.920 6.11e-12 ***
C5000                -9.297e-04  6.794e-04   -1.368 0.171367    
C10000               -7.646e-04  8.869e-04   -0.862 0.388784    
C25000               -6.099e-04  9.118e-04   -0.669 0.503656    
C50000                1.248e-03  5.573e-04    2.239 0.025268 *  
TL0                   6.088e-06  2.582e-07   23.581  < 2e-16 ***
TL1000               -1.202e-06  1.438e-07   -8.360  < 2e-16 ***
TL5000                1.116e-07  1.314e-07    0.849 0.396110    
TL10000              -3.424e-08  1.107e-07   -0.309 0.757059    
TL25000               6.223e-09  4.149e-08    0.150 0.880779    
TL50000              -4.691e-08  1.490e-08   -3.148 0.001668 ** 
LC                   -1.443e-08  6.354e-09   -2.271 0.023259 *  
TL                           NA         NA       NA       NA    
GC                   -7.733e-02  1.496e-02   -5.169 2.60e-07 ***
N50                   3.594e-08  8.813e-09    4.078 4.72e-05 ***
NG50                 -5.392e-08  1.051e-08   -5.131 3.16e-07 ***
N75                  -8.383e-08  2.156e-08   -3.889 0.000104 ***
NG75                  7.817e-08  2.174e-08    3.595 0.000333 ***
L50                   5.598e-03  1.417e-03    3.952 8.03e-05 ***
LG50                 -6.861e-03  1.441e-03   -4.760 2.07e-06 ***
L75                  -2.874e-03  1.151e-03   -2.497 0.012624 *  
LG75                  2.585e-03  1.078e-03    2.397 0.016641 *  
MA                   -5.408e-03  1.043e-03   -5.187 2.35e-07 ***
MAC                   6.786e-03  1.191e-03    5.699 1.39e-08 ***
MACL                 -1.304e-09  1.074e-09   -1.214 0.225021    
LMA                  -2.237e-05  8.334e-05   -0.268 0.788435    
SQEM                 -6.377e-05  4.255e-04   -0.150 0.880883    
SQLM                 -1.817e-05  8.740e-05   -0.208 0.835378    
UAMC                 -2.647e-03  1.963e-03   -1.349 0.177541    
UAC                   6.976e-04  7.492e-05    9.311  < 2e-16 ***
UACP                  4.574e-03  5.021e-04    9.110  < 2e-16 ***
UAL                  -6.501e-06  2.751e-07  -23.630  < 2e-16 ***
GF                   -8.626e-03  3.800e-03   -2.270 0.023338 *  
DR                   -3.048e-01  3.481e-01   -0.875 0.381415    
N100                 -3.284e-04  2.073e-06 -158.378  < 2e-16 ***
MM100                -2.434e-04  1.351e-04   -1.802 0.071736 .  
ID100                 3.193e-03  7.086e-04    4.506 7.00e-06 ***
LA                    1.630e-08  6.331e-09    2.574 0.010117 *  
TAL                  -3.298e-06  2.049e-07  -16.093  < 2e-16 ***
NA50                 -6.894e-09  1.172e-08   -0.588 0.556419    
NGA50                 2.060e-08  1.292e-08    1.594 0.111084    
NA75                  1.486e-08  1.235e-08    1.203 0.229077    
NGA75                -1.093e-08  1.156e-08   -0.945 0.344834    
LA50                 -3.245e-03  1.485e-03   -2.186 0.028957 *  
LGA50                 3.945e-03  1.469e-03    2.686 0.007286 ** 
LA75                  1.266e-03  9.728e-04    1.301 0.193442    
LGA75                -1.412e-03  8.673e-04   -1.628 0.103624    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for quasipoisson family taken to be 0.01287379)
Null deviance: 6542.130  on 2015  degrees of freedom
Residual deviance:   25.142  on 1954  degrees of freedom
AIC: NA
Number of Fisher Scoring iterations: 3

rm(comment)

## save GLM results
write.csv(summary(AssemblyGLM1)['coefficients'], file="AssemblyGLM1.csv")
write.csv(summary(AssemblyGLM2)['coefficients'], file="AssemblyGLM2.csv")
write.csv(summary(AssemblyGLM3)['coefficients'], file="AssemblyGLM3.csv")

# GLMs from datarownameNArow for each cgMLST workflow (i.e. excepted MentaLIST which does not present assembly output)
## subset
datarownameNArowGLM.BIGSdb=subset(datarownameNArowGLM,datarownameNArowGLM$WORKFLOW %in% c("BIGSdb"))
datarownameNArowGLM.INNUENDO=subset(datarownameNArowGLM,datarownameNArowGLM$WORKFLOW %in% c("INNUENDO"))
datarownameNArowGLM.GENPAT=subset(datarownameNArowGLM,datarownameNArowGLM$WORKFLOW %in% c("GENPAT"))
datarownameNArowGLM.SeqSphere=subset(datarownameNArowGLM,datarownameNArowGLM$WORKFLOW %in% c("SeqSphere"))
datarownameNArowGLM.BioNumerics=subset(datarownameNArowGLM,datarownameNArowGLM$WORKFLOW %in% c("BioNumerics"))

## test distribution of the variable to explain
### Gaussian
shapiro.test(datarownameNArowGLM.BIGSdb$IAAR)
# => W = 0.71443, p-value < 2.2e-16
shapiro.test(datarownameNArowGLM.INNUENDO$IAAR)
# => W = 0.65059, p-value < 2.2e-16
shapiro.test(datarownameNArowGLM.GENPAT$IAAR)
# => W = 0.69449, p-value < 2.2e-16
shapiro.test(datarownameNArowGLM.SeqSphere$IAAR)
# => W = 0.71305, p-value < 2.2e-16
shapiro.test(datarownameNArowGLM.BioNumerics$IAAR)
# => W = 0.44714, p-value < 2.2e-16
# => different than gaussian
### Poisson
poisson.test(x = sum(datarownameNArowGLM.BIGSdb$IAAR), T = length(datarownameNArowGLM.BIGSdb$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 733243, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.INNUENDO$IAAR), T = length(datarownameNArowGLM.INNUENDO$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 586080, time base = 336, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.GENPAT$IAAR), T = length(datarownameNArowGLM.GENPAT$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 732442, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.SeqSphere$IAAR), T = length(datarownameNArowGLM.SeqSphere$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 733276, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.BioNumerics$IAAR), T = length(datarownameNArowGLM.BioNumerics$IAAR), r = 1,
             alternative = "two.sided",
             conf.level = 0.95)
# => number of events = 707608, time base = 420, p-value < 2.2e-16
# => different than two side poisson distribution
poisson.test(x = sum(datarownameNArowGLM.BIGSdb$IAAR), T = length(datarownameNArowGLM.BIGSdb$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 733243, time base = 420, p-value = 1
poisson.test(x = sum(datarownameNArowGLM.INNUENDO$IAAR), T = length(datarownameNArowGLM.INNUENDO$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 586080, time base = 336, p-value = 1
poisson.test(x = sum(datarownameNArowGLM.GENPAT$IAAR), T = length(datarownameNArowGLM.GENPAT$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 732442, time base = 420, p-value = 1
poisson.test(x = sum(datarownameNArowGLM.SeqSphere$IAAR), T = length(datarownameNArowGLM.SeqSphere$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 733276, time base = 420, p-value = 1
poisson.test(x = sum(datarownameNArowGLM.BioNumerics$IAAR), T = length(datarownameNArowGLM.BioNumerics$IAAR), r = 1,
             alternative = "less",
             conf.level = 0.95)
# => number of events = 707608, time base = 420, p-value = 1
# => one side poisson distribution with lower hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≥ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 < 0)
poisson.test(x = sum(datarownameNArowGLM.BIGSdb$IAAR), T = length(datarownameNArowGLM.BIGSdb$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 733243, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.INNUENDO$IAAR), T = length(datarownameNArowGLM.INNUENDO$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 586080, time base = 336, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.GENPAT$IAAR), T = length(datarownameNArowGLM.GENPAT$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 732442, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.SeqSphere$IAAR), T = length(datarownameNArowGLM.SeqSphere$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 733276, time base = 420, p-value < 2.2e-16
poisson.test(x = sum(datarownameNArowGLM.BioNumerics$IAAR), T = length(datarownameNArowGLM.BioNumerics$IAAR), r = 1,
             alternative = "greater",
             conf.level = 0.95)
# => number of events = 707608, time base = 420, p-value < 2.2e-16
# => different than one side poisson distribution with upper hypothesis (𝐻𝐻0: 𝜆𝜆2 − 𝜆𝜆1 ≤ 0 vs. 𝐻𝐻𝐴𝐴: 𝜆𝜆2 − 𝜆𝜆1 > 0)
### uniform
chisq.test(datarownameNArowGLM.BIGSdb$IAAR)
# => X-squared = 1.9744, df = 419, p-value = 1
# => uniform
chisq.test(datarownameNArowGLM.INNUENDO$IAAR)
# => X-squared = 0.18608, df = 335, p-value = 1
# => uniform
chisq.test(datarownameNArowGLM.GENPAT$IAAR)
# => X-squared = 0.60701, df = 419, p-value = 1
# => uniform
chisq.test(datarownameNArowGLM.SeqSphere$IAAR)
# => X-squared = 1.7409, df = 419, p-value = 1
# => uniform
chisq.test(datarownameNArowGLM.BioNumerics$IAAR)
# => X-squared = 5359.9, df = 419, p-value < 2.2e-16
# => different than uniform
### dooble check uniform
datarownameNArowGLM.BIGSdb$IAAR.rate <- datarownameNArowGLM.BIGSdb$IAAR / 1748
chisq.unif.test(datarownameNArowGLM.BIGSdb$IAAR.rate, min.bin.size=10)
# => Error in chisq.test(counts, ...) : 'x' must at least have 2 elements
datarownameNArowGLM.INNUENDO$IAAR.rate <- datarownameNArowGLM.INNUENDO$IAAR / 1748
chisq.unif.test(datarownameNArowGLM.INNUENDO$IAAR.rate, min.bin.size=10)
# => Error in chisq.test(counts, ...) : 'x' must at least have 2 elements
datarownameNArowGLM.GENPAT$IAAR.rate <- datarownameNArowGLM.GENPAT$IAAR / 1748
chisq.unif.test(datarownameNArowGLM.GENPAT$IAAR.rate, min.bin.size=10)
# => Error in chisq.test(counts, ...) : 'x' must at least have 2 elements
datarownameNArowGLM.SeqSphere$IAAR.rate <- datarownameNArowGLM.SeqSphere$IAAR / 1748
chisq.unif.test(datarownameNArowGLM.SeqSphere$IAAR.rate, min.bin.size=10)
# => Error in chisq.test(counts, ...) : 'x' must at least have 2 elements
datarownameNArowGLM.BioNumerics$IAAR.rate <- datarownameNArowGLM.BioNumerics$IAAR / 1748
chisq.unif.test(datarownameNArowGLM.BioNumerics$IAAR.rate, min.bin.size=10)
# => X-squared = 607.81, df = 2, a = 0, b = 1, p-value < 2.2e-16
# => different than uniform

## check distribution visually
hist(datarownameNArowGLM.BIGSdb$IAAR)
hist(datarownameNArowGLM.INNUENDO$IAAR)
hist(datarownameNArowGLM.GENPAT$IAAR)
hist(datarownameNArowGLM.SeqSphere$IAAR)
hist(datarownameNArowGLM.BioNumerics$IAAR)

## select the variable to explain and the variables of interest (do not add DrDk because already linear with others)
formulaAdditionAssemblyAlone <- IAAR~REFERENCE+PLATING+DNA+SEQUENCING+DEPTH+BREADTH+IAAS+C0+C1000+C5000+C10000+C25000+C50000+TL0+TL1000+TL5000+TL10000+TL25000+TL50000+LC+TL+GC+N50+NG50+N75+NG75+L50+LG50+L75+LG75+MA+MAC+MACL+LMA+SQEM+SQLM+UAMC+UAC+UACP+UAL+GF+DR+N100+MM100+ID100+LA+TAL+NA50+NGA50+NA75+NGA75+LA50+LGA50+LA75+LGA75
formulaMultiplicationAssemblyAlone <- IAAR~REFERENCE*PLATING*DNA*SEQUENCING*DEPTH*BREADTH*IAAS*C0*C1000*C5000*C10000*C25000*C50000*TL0*TL1000*TL5000*TL10000*TL25000*TL50000*LC*TL*GC*N50*NG50*N75*NG75*L50*LG50*L75*LG75*MA*MAC*MACL*LMA*SQEM*SQLM*UAMC*UAC*UACP*UAL*GF*DR*N100*MM100*ID100*LA*TAL*NA50*NGA50*NA75*NGA75*LA50*LGA50*LA75*LGA75

## run GLM
BIGSdbGLM1 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BIGSdb, family = poisson)
BIGSdbGLM2 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BIGSdb, family = poisson(link = "log"))
BIGSdbGLM3 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BIGSdb, family = quasipoisson)
INNUENDOGLM1 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.INNUENDO, family = poisson)
INNUENDOGLM2 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.INNUENDO, family = poisson(link = "log"))
INNUENDOGLM3 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.INNUENDO, family = quasipoisson)
GENPATGLM1 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.GENPAT, family = poisson)
GENPATGLM2 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.GENPAT, family = poisson(link = "log"))
GENPATGLM3 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.GENPAT, family = quasipoisson)
SeqSphereGLM1 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.SeqSphere, family = poisson)
SeqSphereGLM2 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.SeqSphere, family = poisson(link = "log"))
SeqSphereGLM3 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.SeqSphere, family = quasipoisson)
BioNumericsGLM1 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BioNumerics, family = poisson)
BioNumericsGLM2 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BioNumerics, family = poisson(link = "log"))
BioNumericsGLM3 <- glm(formulaAdditionAssemblyAlone, data=datarownameNArowGLM.BioNumerics, family = quasipoisson)

### test over dispersion  (i.e. alpha > 0 with p<5%)
dispersiontest(BIGSdbGLM1,trafo=1)
# => z = -22447, p-value = 1
# => sample estimates: alpha -0.9768956 
dispersiontest(BIGSdbGLM2,trafo=1)
# => z = -22447, p-value = 1
# => sample estimates: alpha -0.9997595 
dispersiontest(BIGSdbGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (BIGSdbGLM1)
dispersiontest(INNUENDOGLM1,trafo=1)
# => z = -121031, p-value = 1
# => sample estimates: alpha -0.9999817 
dispersiontest(INNUENDOGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (INNUENDOGLM1)
dispersiontest(GENPATGLM1,trafo=1)
# => z = -36069, p-value = 1
# => sample estimates: alpha -0.9998695  
dispersiontest(GENPATGLM2,trafo=1)
# => z = -36069, p-value = 1
# => sample estimates: alpha -0.9998695  
dispersiontest(GENPATGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (GENPATGLM1)
dispersiontest(SeqSphereGLM1,trafo=1)
# => z = -23298, p-value = 1
# => sample estimates: alpha -0.9997593
dispersiontest(SeqSphereGLM2,trafo=1)
# => z = -23298, p-value = 1
# => sample estimates: alpha -0.9997593
dispersiontest(SeqSphereGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (SeqSphereGLM1)
dispersiontest(BioNumericsGLM1,trafo=1)
# => z = -215.63, p-value = 1
# => sample estimates: alpha -0.9769445 
dispersiontest(BioNumericsGLM2,trafo=1)
# => z = -215.63, p-value = 1
# => sample estimates: alpha -0.9769445  
dispersiontest(BioNumericsGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (BioNumericsGLM1)

## check significant parameters (<0.01)
summary(BIGSdbGLM1) # no significant parameters
summary(BIGSdbGLM2) # no significant parameters
summary(BIGSdbGLM3) # 18 significant parameters
summary(INNUENDOGLM1) # no significant parameters
summary(INNUENDOGLM2) # no significant parameters
summary(INNUENDOGLM3) # 8 significant parameters
summary(GENPATGLM1) # no significant parameters
summary(GENPATGLM2) # no significant parameters
summary(GENPATGLM3) # 18 significant parameters
summary(SeqSphereGLM1) # no significant parameters
summary(SeqSphereGLM2) # no significant parameters
summary(SeqSphereGLM3) # 20 significant parameters
summary(BioNumericsGLM1) # 1 significant parameters: N100
summary(BioNumericsGLM2) # 1 significant parameters: N100
summary(BioNumericsGLM3) # 17 significant parameters

## keep GLM1
summary(BIGSdbGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssemblyAlone, family = poisson, data = datarownameNArowGLM.BIGSdb)
Deviance Residuals: 
  Min         1Q     Median         3Q        Max  
-0.075249  -0.006135   0.000345   0.006807   0.084397  
Coefficients: (4 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)
(Intercept)           7.896e+00  1.751e+01   0.451    0.652
REFERENCEATCC19115    1.225e-02  1.355e-01   0.090    0.928
REFERENCEATCCBAA679   1.285e-02  1.639e-01   0.078    0.938
PLATINGfifth_culture  6.360e-06  4.441e-03   0.001    0.999
PLATINGtenth_culture  9.992e-05  4.520e-03   0.022    0.982
DNAextraction_A      -5.232e-05  3.387e-03  -0.015    0.988
DNAextraction_B      -1.244e-05  3.346e-03  -0.004    0.997
DNAextraction_C              NA         NA      NA       NA
SEQUENCINGNextSeq_B  -2.148e-05  2.842e-03  -0.008    0.994
DEPTH                 1.174e-06  8.097e-05   0.015    0.988
BREADTH               2.061e-04  3.467e-02   0.006    0.995
IAAS                         NA         NA      NA       NA
C0                    8.394e-05  4.111e-03   0.020    0.984
C1000                -2.246e-04  5.777e-03  -0.039    0.969
C5000                -2.914e-04  1.747e-02  -0.017    0.987
C10000               -3.289e-03  3.630e-02  -0.091    0.928
C25000                3.985e-03  3.672e-02   0.109    0.914
C50000               -5.310e-04  1.578e-02  -0.034    0.973
TL0                   1.482e-07  1.065e-05   0.014    0.989
TL1000                9.151e-08  7.966e-06   0.011    0.991
TL5000                8.282e-09  3.426e-06   0.002    0.998
TL10000               1.535e-07  3.413e-06   0.045    0.964
TL25000              -2.132e-07  1.974e-06  -0.108    0.914
TL50000               8.800e-09  4.087e-07   0.022    0.983
LC                    1.542e-09  1.282e-07   0.012    0.990
TL                           NA         NA      NA       NA
GC                    7.757e-03  3.567e-01   0.022    0.983
N50                  -9.614e-09  6.272e-07  -0.015    0.988
NG50                  5.865e-09  6.147e-07   0.010    0.992
N75                  -5.721e-10  3.636e-07  -0.002    0.999
NG75                  1.376e-09  3.843e-07   0.004    0.997
L50                  -8.734e-04  5.458e-02  -0.016    0.987
LG50                  6.651e-04  5.134e-02   0.013    0.990
L75                  -7.303e-04  3.370e-02  -0.022    0.983
LG75                  1.072e-03  2.398e-02   0.045    0.964
MA                   -1.227e-03  1.597e-02  -0.077    0.939
MAC                   2.995e-04  2.298e-02   0.013    0.990
MACL                  8.335e-10  2.212e-08   0.038    0.970
LMA                  -1.447e-04  1.666e-03  -0.087    0.931
SQEM                 -7.250e-05  6.236e-03  -0.012    0.991
SQLM                  1.429e-05  2.292e-03   0.006    0.995
UAMC                         NA         NA      NA       NA
UAC                   1.195e-04  1.228e-02   0.010    0.992
UACP                  1.557e-04  2.016e-02   0.008    0.994
UAL                  -2.330e-07  1.284e-05  -0.018    0.986
GF                   -7.258e-03  1.239e-01  -0.059    0.953
DR                   -7.264e-01  1.208e+01  -0.060    0.952
N100                 -1.025e-05  6.023e-04  -0.017    0.986
MM100                 5.257e-05  2.624e-03   0.020    0.984
ID100                -7.324e-04  1.575e-02  -0.047    0.963
LA                   -3.570e-10  1.218e-07  -0.003    0.998
TAL                   3.874e-08  3.765e-06   0.010    0.992
NA50                  7.652e-09  3.974e-07   0.019    0.985
NGA50                -5.206e-09  3.717e-07  -0.014    0.989
NA75                  8.314e-09  4.666e-07   0.018    0.986
NGA75                -6.751e-09  4.212e-07  -0.016    0.987
LA50                  4.765e-04  3.654e-02   0.013    0.990
LGA50                -1.908e-04  3.140e-02  -0.006    0.995
LA75                  1.461e-03  3.241e-02   0.045    0.964
LGA75                -1.571e-03  1.898e-02  -0.083    0.934
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 1.97544  on 419  degrees of freedom
Residual deviance: 0.10102  on 364  degrees of freedom
AIC: 4019.3
Number of Fisher Scoring iterations: 3

rm(comment)

summary(INNUENDOGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssemblyAlone, family = poisson, data = datarownameNArowGLM.INNUENDO)
Deviance Residuals: 
  Min          1Q      Median          3Q         Max  
-0.0208000  -0.0009894   0.0005833   0.0022093   0.0086684  
Coefficients: (12 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)
(Intercept)           7.301e+00  1.399e+01   0.522    0.602
REFERENCEATCC19115    1.690e-03  4.076e-01   0.004    0.997
REFERENCEATCCBAA679   1.695e-03  4.982e-01   0.003    0.997
PLATINGfifth_culture  5.527e-05  5.150e-03   0.011    0.991
PLATINGtenth_culture  9.074e-05  5.195e-03   0.017    0.986
DNAextraction_A       6.807e-06  4.045e-03   0.002    0.999
DNAextraction_B      -3.643e-05  3.895e-03  -0.009    0.993
DNAextraction_C              NA         NA      NA       NA
SEQUENCINGNextSeq_B   1.230e-05  3.342e-03   0.004    0.997
DEPTH                -5.146e-07  9.506e-05  -0.005    0.996
BREADTH              -3.956e-05  3.629e-02  -0.001    0.999
IAAS                         NA         NA      NA       NA
C0                   -2.789e-05  2.738e-03  -0.010    0.992
C1000                 3.049e-05  4.581e-03   0.007    0.995
C5000                 3.779e-05  1.628e-02   0.002    0.998
C10000               -4.928e-05  3.457e-02  -0.001    0.999
C25000                1.344e-04  3.665e-02   0.004    0.997
C50000               -1.182e-04  2.116e-02  -0.006    0.996
TL0                  -4.280e-08  4.495e-05  -0.001    0.999
TL1000               -5.878e-08  5.950e-06  -0.010    0.992
TL5000                2.284e-09  3.434e-06   0.001    0.999
TL10000               5.307e-09  3.314e-06   0.002    0.999
TL25000              -1.160e-09  2.132e-06  -0.001    1.000
TL50000               2.507e-09  5.926e-07   0.004    0.997
LC                   -1.858e-07  9.248e-05  -0.002    0.998
TL                           NA         NA      NA       NA
GC                    2.015e-03  3.668e-01   0.005    0.996
N50                  -1.768e-07  1.613e-03   0.000    1.000
NG50                  1.672e-07  1.613e-03   0.000    1.000
N75                   1.131e-06  6.585e-03   0.000    1.000
NG75                  7.210e-07  6.563e-03   0.000    1.000
L50                  -4.766e-01  1.423e+02  -0.003    0.997
LG50                 -2.872e-05  4.864e-02  -0.001    1.000
L75                   2.207e-01  6.592e+01   0.003    0.997
LG75                 -6.569e-05  3.224e-02  -0.002    0.998
MA                    1.864e-04  6.964e-02   0.003    0.998
MAC                          NA         NA      NA       NA
MACL                 -4.741e-07  1.426e-04  -0.003    0.997
LMA                   1.918e-05  3.373e-03   0.006    0.995
SQEM                         NA         NA      NA       NA
SQLM                         NA         NA      NA       NA
UAMC                         NA         NA      NA       NA
UAC                   5.786e-05  2.289e-02   0.003    0.998
UACP                         NA         NA      NA       NA
UAL                  -1.160e-08  4.952e-05   0.000    1.000
GF                    7.233e-04  4.317e-01   0.002    0.999
DR                           NA         NA      NA       NA
N100                         NA         NA      NA       NA
MM100                -1.150e-04  5.558e-03  -0.021    0.983
ID100                 3.549e-04  1.821e-02   0.019    0.984
LA                    1.858e-07  9.248e-05   0.002    0.998
TAL                   9.861e-08  4.671e-05   0.002    0.998
NA50                  1.772e-07  1.612e-03   0.000    1.000
NGA50                -1.678e-07  1.612e-03   0.000    1.000
NA75                 -1.131e-06  6.585e-03   0.000    1.000
NGA75                -7.212e-07  6.563e-03   0.000    1.000
LA50                  4.766e-01  1.423e+02   0.003    0.997
LGA50                        NA         NA      NA       NA
LA75                 -2.206e-01  6.592e+01  -0.003    0.997
LGA75                        NA         NA      NA       NA
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 0.186103  on 335  degrees of freedom
Residual deviance: 0.006147  on 288  degrees of freedom
AIC: 3221.5
Number of Fisher Scoring iterations: 2

rm(comment)

summary(GENPATGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssemblyAlone, family = poisson, data = datarownameNArowGLM.GENPAT)
Deviance Residuals: 
  Min         1Q     Median         3Q        Max  
-0.064526  -0.003717   0.000834   0.004754   0.083201  
Coefficients: (4 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)
(Intercept)           6.783e+00  3.001e+01   0.226    0.821
REFERENCEATCC19115    9.840e-03  2.181e-01   0.045    0.964
REFERENCEATCCBAA679   1.164e-02  2.661e-01   0.044    0.965
PLATINGfifth_culture  4.429e-05  4.613e-03   0.010    0.992
PLATINGtenth_culture  1.309e-04  4.753e-03   0.028    0.978
DNAextraction_A      -1.601e-05  3.520e-03  -0.005    0.996
DNAextraction_B      -8.998e-05  3.351e-03  -0.027    0.979
DNAextraction_C              NA         NA      NA       NA
SEQUENCINGNextSeq_B  -2.235e-05  3.145e-03  -0.007    0.994
DEPTH                -6.306e-07  8.220e-05  -0.008    0.994
BREADTH               3.938e-04  3.491e-02   0.011    0.991
IAAS                         NA         NA      NA       NA
C0                    1.341e-05  1.843e-03   0.007    0.994
C1000                -2.087e-04  4.953e-03  -0.042    0.966
C5000                -6.371e-05  2.532e-02  -0.003    0.998
C10000               -4.919e-04  3.077e-02  -0.016    0.987
C25000               -1.124e-04  3.059e-02  -0.004    0.997
C50000                2.031e-04  2.082e-02   0.010    0.992
TL0                   3.096e-07  1.130e-05   0.027    0.978
TL1000                9.020e-08  4.577e-06   0.020    0.984
TL5000                1.652e-08  4.832e-06   0.003    0.997
TL10000              -4.766e-08  4.022e-06  -0.012    0.991
TL25000              -1.596e-08  1.476e-06  -0.011    0.991
TL50000              -1.642e-08  5.410e-07  -0.030    0.976
LC                   -1.847e-08  5.032e-07  -0.037    0.971
TL                           NA         NA      NA       NA
GC                   -4.230e-04  3.549e-01  -0.001    0.999
N50                   3.486e-08  1.572e-06   0.022    0.982
NG50                 -4.054e-08  1.550e-06  -0.026    0.979
N75                   9.434e-08  6.840e-06   0.014    0.989
NG75                 -1.090e-07  6.824e-06  -0.016    0.987
L50                   3.116e-04  1.167e-01   0.003    0.998
LG50                 -2.565e-03  9.615e-02  -0.027    0.979
L75                  -6.914e-04  1.051e-01  -0.007    0.995
LG75                 -3.275e-04  7.958e-02  -0.004    0.997
MA                    5.123e-04  3.038e-02   0.017    0.987
MAC                  -1.003e-03  3.260e-02  -0.031    0.975
MACL                 -4.557e-12  4.233e-08   0.000    1.000
LMA                  -1.702e-05  2.054e-03  -0.008    0.993
SQEM                  3.058e-05  1.314e-02   0.002    0.998
SQLM                 -3.357e-05  3.400e-03  -0.010    0.992
UAMC                         NA         NA      NA       NA
UAC                  -8.489e-05  4.930e-03  -0.017    0.986
UACP                 -1.035e-05  1.884e-02  -0.001    1.000
UAL                  -1.777e-07  1.330e-05  -0.013    0.989
GF                   -5.516e-03  2.201e-01  -0.025    0.980
DR                    2.873e-01  2.777e+01   0.010    0.992
N100                 -9.251e-06  8.831e-04  -0.010    0.992
MM100                -7.275e-05  3.076e-03  -0.024    0.981
ID100                -4.763e-04  1.619e-02  -0.029    0.977
LA                    1.879e-08  5.031e-07   0.037    0.970
TAL                  -2.330e-08  7.268e-06  -0.003    0.997
NA50                 -3.401e-08  7.954e-07  -0.043    0.966
NGA50                 3.791e-08  7.660e-07   0.049    0.961
NA75                 -5.804e-08  6.695e-06  -0.009    0.993
NGA75                 6.901e-08  6.671e-06   0.010    0.992
LA50                 -5.179e-04  5.892e-02  -0.009    0.993
LGA50                 2.983e-03  4.136e-02   0.072    0.942
LA75                  1.618e-03  1.019e-01   0.016    0.987
LGA75                -5.452e-04  7.188e-02  -0.008    0.994
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 0.607339  on 419  degrees of freedom
Residual deviance: 0.054826  on 364  degrees of freedom
AIC: 4018.8
Number of Fisher Scoring iterations: 3

rm(comment)

summary(SeqSphereGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssemblyAlone, family = poisson, data = datarownameNArowGLM.SeqSphere)
Deviance Residuals: 
  Min         1Q     Median         3Q        Max  
-0.080994  -0.004923   0.000630   0.006767   0.065263  
Coefficients: (7 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)
(Intercept)           6.691e+00  1.623e+01   0.412    0.680
REFERENCEATCC19115    1.212e-02  2.106e-01   0.058    0.954
REFERENCEATCCBAA679   1.251e-02  2.569e-01   0.049    0.961
PLATINGfifth_culture  1.178e-04  4.435e-03   0.027    0.979
PLATINGtenth_culture  1.123e-04  4.594e-03   0.024    0.980
DNAextraction_A      -2.905e-05  3.461e-03  -0.008    0.993
DNAextraction_B      -6.532e-06  3.312e-03  -0.002    0.998
DNAextraction_C              NA         NA      NA       NA
SEQUENCINGNextSeq_B  -8.202e-05  2.926e-03  -0.028    0.978
DEPTH                 5.583e-07  8.585e-05   0.007    0.995
BREADTH               7.943e-04  3.456e-02   0.023    0.982
IAAS                         NA         NA      NA       NA
C0                    1.141e-06  6.963e-04   0.002    0.999
C1000                -1.653e-04  3.523e-03  -0.047    0.963
C5000                 7.474e-04  1.230e-02   0.061    0.952
C10000               -5.917e-04  1.460e-02  -0.041    0.968
C25000               -7.701e-05  1.680e-02  -0.005    0.996
C50000               -6.763e-05  1.126e-02  -0.006    0.995
TL0                   1.843e-07  1.054e-05   0.017    0.986
TL1000                3.435e-08  3.391e-06   0.010    0.992
TL5000               -1.396e-07  2.323e-06  -0.060    0.952
TL10000               8.708e-08  1.883e-06   0.046    0.963
TL25000               8.174e-09  7.035e-07   0.012    0.991
TL50000              -1.802e-09  2.822e-07  -0.006    0.995
LC                   -5.069e-09  2.651e-07  -0.019    0.985
TL                           NA         NA      NA       NA
GC                    1.356e-03  3.184e-01   0.004    0.997
N50                  -8.591e-09  4.157e-07  -0.021    0.984
NG50                  1.775e-08  4.975e-07   0.036    0.972
N75                   2.634e-08  5.345e-07   0.049    0.961
NG75                 -2.608e-08  5.253e-07  -0.050    0.960
L50                  -1.463e-03  5.470e-02  -0.027    0.979
LG50                  1.828e-03  5.938e-02   0.031    0.975
L75                   1.651e-03  2.138e-02   0.077    0.938
LG75                 -1.385e-03  2.204e-02  -0.063    0.950
MA                   -4.431e-04  1.044e-02  -0.042    0.966
MAC                          NA         NA      NA       NA
MACL                  1.038e-09  4.460e-08   0.023    0.981
LMA                  -8.148e-05  1.970e-03  -0.041    0.967
SQEM                         NA         NA      NA       NA
SQLM                         NA         NA      NA       NA
UAMC                         NA         NA      NA       NA
UAC                  -3.166e-05  1.111e-03  -0.028    0.977
UACP                 -2.383e-04  1.621e-02  -0.015    0.988
UAL                  -1.134e-07  1.043e-05  -0.011    0.991
GF                   -1.315e-03  2.185e-01  -0.006    0.995
DR                   -5.209e-02  4.214e+00  -0.012    0.990
N100                 -1.628e-03  1.782e-01  -0.009    0.993
MM100                 1.742e-05  3.416e-03   0.005    0.996
ID100                -3.297e-04  1.638e-02  -0.020    0.984
LA                    5.964e-09  2.709e-07   0.022    0.982
TAL                   1.074e-07  9.643e-06   0.011    0.991
NA50                  7.097e-09  4.381e-07   0.016    0.987
NGA50                -1.694e-08  5.132e-07  -0.033    0.974
NA75                 -1.133e-08  1.851e-07  -0.061    0.951
NGA75                 1.007e-08  1.751e-07   0.058    0.954
LA50                  1.275e-03  5.626e-02   0.023    0.982
LGA50                -1.320e-03  6.076e-02  -0.022    0.983
LA75                 -1.403e-03  1.641e-02  -0.085    0.932
LGA75                 9.068e-04  1.590e-02   0.057    0.955
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 1.74158  on 419  degrees of freedom
Residual deviance: 0.10112  on 367  degrees of freedom
AIC: 4013.4
Number of Fisher Scoring iterations: 3

rm(comment)

summary(BioNumericsGLM1)
comment <- scan(what="character")
Call: glm(formula = formulaAdditionAssemblyAlone, family = poisson, data = datarownameNArowGLM.BioNumerics)
Deviance Residuals: 
  Min        1Q    Median        3Q       Max  
-0.93420  -0.05068   0.00220   0.05196   0.91224  
Coefficients: (3 not defined because of singularities)
Estimate Std. Error z value Pr(>|z|)    
(Intercept)           9.897e-01  1.756e+01   0.056    0.955    
REFERENCEATCC19115    7.269e-02  2.041e-01   0.356    0.722    
REFERENCEATCCBAA679   8.133e-02  2.438e-01   0.334    0.739    
PLATINGfifth_culture  4.000e-04  4.461e-03   0.090    0.929    
PLATINGtenth_culture  5.177e-04  4.670e-03   0.111    0.912    
DNAextraction_A       3.699e-05  3.449e-03   0.011    0.991    
DNAextraction_B       5.515e-05  3.519e-03   0.016    0.987    
DNAextraction_C              NA         NA      NA       NA    
SEQUENCINGNextSeq_B  -8.419e-04  2.948e-03  -0.286    0.775    
DEPTH                -9.969e-06  9.293e-05  -0.107    0.915    
BREADTH               1.219e-02  3.520e-02   0.346    0.729    
IAAS                         NA         NA      NA       NA    
C0                   -3.867e-04  4.134e-03  -0.094    0.925    
C1000                 2.550e-03  6.321e-03   0.403    0.687    
C5000                -2.001e-04  2.102e-02  -0.010    0.992    
C10000               -9.046e-03  2.764e-02  -0.327    0.743    
C25000                5.862e-03  2.563e-02   0.229    0.819    
C50000                3.391e-03  2.115e-02   0.160    0.873    
TL0                   4.637e-06  1.053e-05   0.441    0.660    
TL1000               -1.257e-06  8.673e-06  -0.145    0.885    
TL5000                1.647e-08  3.550e-06   0.005    0.996    
TL10000               6.984e-07  3.409e-06   0.205    0.838    
TL25000              -4.332e-07  1.179e-06  -0.367    0.713    
TL50000              -6.539e-08  5.910e-07  -0.111    0.912    
LC                   -2.624e-07  2.674e-07  -0.981    0.326    
TL                           NA         NA      NA       NA    
GC                    1.482e-02  3.742e-01   0.040    0.968    
N50                   1.273e-07  1.703e-07   0.747    0.455    
NG50                 -2.131e-07  2.282e-07  -0.934    0.351    
N75                  -6.402e-08  1.226e-06  -0.052    0.958    
NG75                  1.000e-07  1.251e-06   0.080    0.936    
L50                   2.327e-02  2.992e-02   0.778    0.437    
LG50                 -2.139e-02  3.025e-02  -0.707    0.479    
L75                  -8.710e-03  3.399e-02  -0.256    0.798    
LG75                  5.437e-03  3.444e-02   0.158    0.875    
MA                   -4.188e-03  3.664e-02  -0.114    0.909    
MAC                   2.545e-03  3.894e-02   0.065    0.948    
MACL                  3.081e-08  3.785e-08   0.814    0.416    
LMA                  -2.297e-04  1.654e-03  -0.139    0.890    
SQEM                  1.543e-02  3.497e-02   0.441    0.659    
SQLM                  5.528e-04  1.625e-03   0.340    0.734    
UAMC                 -2.897e-03  1.866e-02  -0.155    0.877    
UAC                   8.678e-04  4.357e-03   0.199    0.842    
UACP                 -1.666e-04  2.063e-02  -0.008    0.994    
UAL                  -2.795e-06  1.202e-05  -0.233    0.816    
GF                   -2.163e-02  2.185e-01  -0.099    0.921    
DR                    1.324e+00  8.974e+00   0.148    0.883    
N100                 -3.406e-04  3.367e-05 -10.116   <2e-16 ***
  MM100                -8.987e-05  2.703e-03  -0.033    0.973    
ID100                 2.433e-03  1.668e-02   0.146    0.884    
LA                    2.687e-07  2.749e-07   0.977    0.328    
TAL                  -1.735e-06  7.909e-06  -0.219    0.826    
NA50                  1.116e-08  3.431e-07   0.033    0.974    
NGA50                 8.267e-08  3.903e-07   0.212    0.832    
NA75                  1.270e-08  2.726e-07   0.047    0.963    
NGA75                -4.426e-08  2.421e-07  -0.183    0.855    
LA50                 -1.209e-02  3.940e-02  -0.307    0.759    
LGA50                 6.819e-03  3.991e-02   0.171    0.864    
LA75                  4.740e-03  1.786e-02   0.265    0.791    
LGA75                -4.077e-03  1.702e-02  -0.240    0.811    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
(Dispersion parameter for poisson family taken to be 1)
Null deviance: 5835.822  on 419  degrees of freedom
Residual deviance:    9.688  on 363  degrees of freedom
AIC: 4014.1
Number of Fisher Scoring iterations: 3

rm(comment)

## save GLM results
write.csv(summary(BIGSdbGLM1)['coefficients'], file="BIGSdbGLM1.csv")
write.csv(summary(BIGSdbGLM2)['coefficients'], file="BIGSdbGLM2.csv")
write.csv(summary(BIGSdbGLM3)['coefficients'], file="BIGSdbGLM3.csv")
write.csv(summary(INNUENDOGLM1)['coefficients'], file="INNUENDOGLM1.csv")
write.csv(summary(INNUENDOGLM2)['coefficients'], file="INNUENDOGLM2.csv")
write.csv(summary(INNUENDOGLM3)['coefficients'], file="INNUENDOGLM3.csv")
write.csv(summary(GENPATGLM1)['coefficients'], file="GENPATGLM1.csv")
write.csv(summary(GENPATGLM2)['coefficients'], file="GENPATGLM2.csv")
write.csv(summary(GENPATGLM3)['coefficients'], file="GENPATGLM3.csv")
write.csv(summary(SeqSphereGLM1)['coefficients'], file="SeqSphereGLM1.csv")
write.csv(summary(SeqSphereGLM2)['coefficients'], file="SeqSphereGLM2.csv")
write.csv(summary(SeqSphereGLM3)['coefficients'], file="SeqSphereGLM3.csv")
write.csv(summary(BioNumericsGLM1)['coefficients'], file="BioNumericsGLM1.csv")
write.csv(summary(BioNumericsGLM2)['coefficients'], file="BioNumericsGLM2.csv")
write.csv(summary(BioNumericsGLM3)['coefficients'], file="BioNumericsGLM3.csv")


## Additional GLMs with parameters restricted to assembly
formulaAdditionAssemblyRestricted <- IAAR~IAAS+C0+C1000+C5000+C10000+C25000+C50000+TL0+TL1000+TL5000+TL10000+TL25000+TL50000+LC+TL+GC+N50+NG50+N75+NG75+L50+LG50+L75+LG75+MA+MAC+MACL+LMA+SQEM+SQLM+UAMC+UAC+UACP+UAL+GF+DR+N100+MM100+ID100+LA+TAL+NA50+NGA50+NA75+NGA75+LA50+LGA50+LA75+LGA75
formulaAdditionAssemblyAloneRestricted <- IAAR~IAAS+C0+C1000+C5000+C10000+C25000+C50000+TL0+TL1000+TL5000+TL10000+TL25000+TL50000+LC+TL+GC+N50+NG50+N75+NG75+L50+LG50+L75+LG75+MA+MAC+MACL+LMA+SQEM+SQLM+UAMC+UAC+UACP+UAL+GF+DR+N100+MM100+ID100+LA+TAL+NA50+NGA50+NA75+NGA75+LA50+LGA50+LA75+LGA75

## run GLM
RestrictedAssemblyGLM1 <- glm(formulaAdditionAssemblyRestricted, data=datarownameNArowGLM, family = poisson)
RestrictedAssemblyGLM2 <- glm(formulaAdditionAssemblyRestricted, data=datarownameNArowGLM, family = poisson(link = "log"))
RestrictedAssemblyGLM3 <- glm(formulaAdditionAssemblyRestricted, data=datarownameNArowGLM, family = quasipoisson)
RestrictedBIGSdbGLM1 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BIGSdb, family = poisson)
RestrictedBIGSdbGLM2 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BIGSdb, family = poisson(link = "log"))
RestrictedBIGSdbGLM3 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BIGSdb, family = quasipoisson)
RestrictedINNUENDOGLM1 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.INNUENDO, family = poisson)
RestrictedINNUENDOGLM2 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.INNUENDO, family = poisson(link = "log"))
RestrictedINNUENDOGLM3 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.INNUENDO, family = quasipoisson)
RestrictedGENPATGLM1 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.GENPAT, family = poisson)
RestrictedGENPATGLM2 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.GENPAT, family = poisson(link = "log"))
RestrictedGENPATGLM3 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.GENPAT, family = quasipoisson)
RestrictedSeqSphereGLM1 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.SeqSphere, family = poisson)
RestrictedSeqSphereGLM2 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.SeqSphere, family = poisson(link = "log"))
RestrictedSeqSphereGLM3 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.SeqSphere, family = quasipoisson)
RestrictedBioNumericsGLM1 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BioNumerics, family = poisson)
RestrictedBioNumericsGLM2 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BioNumerics, family = poisson(link = "log"))
RestrictedBioNumericsGLM3 <- glm(formulaAdditionAssemblyAloneRestricted, data=datarownameNArowGLM.BioNumerics, family = quasipoisson)

### test over dispersion  (i.e. alpha > 0 with p<5%)
dispersiontest(RestrictedAssemblyGLM1,trafo=1)
# => z = -443.57, p-value = 1
# => sample estimates: alpha -0.9863838
dispersiontest(RestrictedAssemblyGLM2,trafo=1)
# => z = -443.57, p-value = 1
# => sample estimates: alpha -0.9863838
dispersiontest(RestrictedAssemblyGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (BIGSdbGLM1)
dispersiontest(RestrictedBIGSdbGLM1,trafo=1)
# => z = -21547, p-value = 1
# => sample estimates: alpha -0.9997529
dispersiontest(RestrictedBIGSdbGLM2,trafo=1)
# => z = -21547, p-value = 1
# => sample estimates: alpha -0.9997529 
dispersiontest(RestrictedBIGSdbGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (BIGSdbGLM1)
dispersiontest(RestrictedINNUENDOGLM1,trafo=1)
# => z = -109620, p-value = 1
# => sample estimates: alpha -0.99998
dispersiontest(RestrictedINNUENDOGLM2,trafo=1)
# => z = -109620, p-value = 1
# => sample estimates: alpha -0.99998
dispersiontest(RestrictedINNUENDOGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (INNUENDOGLM1)
dispersiontest(RestrictedGENPATGLM1,trafo=1)
# => z = -34046, p-value = 1
# => sample estimates: alpha -9998608   
dispersiontest(RestrictedGENPATGLM2,trafo=1)
# => z = -34046, p-value = 1
# => sample estimates: alpha -9998608   
dispersiontest(RestrictedGENPATGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (GENPATGLM1)
dispersiontest(RestrictedSeqSphereGLM1,trafo=1)
# => z = -22235, p-value = 1
# => sample estimates: alpha -0.999751 
dispersiontest(RestrictedSeqSphereGLM2,trafo=1)
# => z = -22235, p-value = 1
# => sample estimates: alpha -0.999751 
dispersiontest(RestrictedSeqSphereGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (SeqSphereGLM1)
dispersiontest(RestrictedBioNumericsGLM1,trafo=1)
# => z = -205.91, p-value = 1
# => sample estimates: alpha -0.9762896
dispersiontest(RestrictedBioNumericsGLM2,trafo=1)
# => z = -205.91, p-value = 1
# => sample estimates: alpha -0.9762896
dispersiontest(RestrictedBioNumericsGLM3,trafo=1)
# => only Poisson GLMs can be tested
# => absence of  GLM overdispersion => keep poisson (BioNumericsGLM1)

## check significant parameters (<0.01)
summary(RestrictedAssemblyGLM1) # 5 significant parameters (IAAS, TL0, UAL, N100, TAL)
summary(RestrictedBIGSdbGLM1) # 0 significant parameters
summary(RestrictedINNUENDOGLM1) # 0 significant parameters
summary(RestrictedGENPATGLM1) # 0 significant parameters
summary(RestrictedSeqSphereGLM1) # 0 significant parameters
summary(RestrictedBioNumericsGLM1) # 1 significant parameters: N100

# => no differences of GLM conclusions with or without initial parameters in the GLM formula (REFERENCE+PLATING+DNA+SEQUENCING+DEPTH+BREADTH)

## save GLM results
write.csv(summary(RestrictedAssemblyGLM1)['coefficients'], file="RestrictedAssemblyGLM1.csv")
write.csv(summary(RestrictedAssemblyGLM2)['coefficients'], file="RestrictedAssemblyGLM2.csv")
write.csv(summary(RestrictedAssemblyGLM3)['coefficients'], file="RestrictedAssemblyGLM3.csv")
write.csv(summary(RestrictedBIGSdbGLM1)['coefficients'], file="RestrictedBIGSdbGLM1.csv")
write.csv(summary(RestrictedBIGSdbGLM2)['coefficients'], file="RestrictedBIGSdbGLM2.csv")
write.csv(summary(RestrictedBIGSdbGLM3)['coefficients'], file="RestrictedBIGSdbGLM3.csv")
write.csv(summary(RestrictedINNUENDOGLM1)['coefficients'], file="RestrictedINNUENDOGLM1.csv")
write.csv(summary(RestrictedINNUENDOGLM2)['coefficients'], file="RestrictedINNUENDOGLM2.csv")
write.csv(summary(RestrictedINNUENDOGLM3)['coefficients'], file="RestrictedINNUENDOGLM3.csv")
write.csv(summary(RestrictedGENPATGLM1)['coefficients'], file="RestrictedGENPATGLM1.csv")
write.csv(summary(RestrictedGENPATGLM2)['coefficients'], file="RestrictedGENPATGLM2.csv")
write.csv(summary(RestrictedGENPATGLM3)['coefficients'], file="RestrictedGENPATGLM3.csv")
write.csv(summary(RestrictedSeqSphereGLM1)['coefficients'], file="RestrictedSeqSphereGLM1.csv")
write.csv(summary(RestrictedSeqSphereGLM2)['coefficients'], file="RestrictedSeqSphereGLM2.csv")
write.csv(summary(RestrictedSeqSphereGLM3)['coefficients'], file="RestrictedSeqSphereGLM3.csv")
write.csv(summary(RestrictedBioNumericsGLM1)['coefficients'], file="RestrictedBioNumericsGLM1.csv")
write.csv(summary(RestrictedBioNumericsGLM2)['coefficients'], file="RestrictedBioNumericsGLM2.csv")
write.csv(summary(RestrictedBioNumericsGLM3)['coefficients'], file="RestrictedBioNumericsGLM3.csv")

#### Graphical confirmation: BIGSdb versus INNUENDO versus GENPAT versus SeqSphere versus BioNumerics versus MentaLiST ####

# plot identified_alleles_against_schema for reference
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identified alleles against schema (restricted scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-A-schema-ATCC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-A-schema-ATCC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identified_alleles_against_schema for plating
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identified alleles against schema (restricted scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of successive plating and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B") +
  facet_grid(sample_origin ~ workflow)
p
plot(p)
ggsave("cgMLST-B-schema-plating.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-B-schema-plating.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identified_alleles_against_schema for extraction
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identified alleles against schema (restricted scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of DNA extraction replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("C") +
  facet_grid(DNA_extraction_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-C-schema-extraction.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-C-schema-extraction.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identified_alleles_against_schema for sequencing
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identified alleles against schema (restricted scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("D") +
  facet_grid(sequencing_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-D-schema-sequencing.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-D-schema-sequencing.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identical_alleles_against_reference for reference
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("E") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-E-reference-ATCC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-E-reference-ATCC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identical_alleles_against_reference for plating
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of successive plating and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("F") +
  facet_grid(sample_origin ~ workflow)
p
plot(p)
ggsave("cgMLST-F-schema-plating.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-F-schema-plating.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identical_alleles_against_reference for extraction
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of DNA extraction replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("G") +
  facet_grid(DNA_extraction_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-G-reference-extraction.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-G-reference-extraction.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identical_alleles_against_reference for sequencing
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("H") +
  facet_grid(sequencing_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-H-reference-sequencing.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-H-reference-sequencing.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom identical_alleles_against_reference for reference
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("I") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-I-reference-ATCC-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-I-reference-ATCC-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom identical_alleles_against_reference for plating
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1730, 1750), breaks = c(1730,1735,1740,1745,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of successive plating and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("J") +
  facet_grid(sample_origin ~ workflow)
p
plot(p)
ggsave("cgMLST-J-schema-plating.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-J-schema-plating.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom identical_alleles_against_reference for extraction
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1730, 1750), breaks = c(1730,1735,1740,1745,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of DNA extraction replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("K") +
  facet_grid(DNA_extraction_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-K-reference-extraction-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-K-reference-extraction-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom identical_alleles_against_reference for sequencing
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1730, 1750), breaks = c(1730,1735,1740,1745,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("L") +
  facet_grid(sequencing_replicate ~ workflow)
p
plot(p)
ggsave("cgMLST-L-reference-sequencing-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-L-reference-sequencing-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot identical_alleles_against_reference for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-A-reference.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-A-reference.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom identical_alleles_against_reference for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = identical_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Identical alleles against reference (restricted scale)", limits = c(1700, 1750), breaks = c(1700,1710,1720,1730,1740,1750)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-B-reference-zoom.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-B-reference-zoom.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot misidentified_alleles_against_reference for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = misidentified_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Misidentified alleles against reference (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-A-reference-misidentified.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-A-reference-misidentified.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom misidentified_alleles_against_reference for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = misidentified_alleles_against_reference)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Misidentified alleles against reference (restricted scale)", limits = c(0, 30), breaks = c(0,10,20,30)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-B-reference-zoom-misidentified.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-B-reference-zoom-misidentified.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot unidentified_alleles_against_schema for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = unidentified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Unidentified alleles against schema (extended scale)", limits = c(-1, 1750), breaks = c(0,500,1000,1500)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Impact of sequencing replicate and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("A") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-A-reference-unidentified.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-A-reference-unidentified.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# plot zoom unidentified_alleles_against_schema for reference_strain in figure
p = ggplot(data = data_cgMLST, aes(x = targeted_depth, y = unidentified_alleles_against_schema)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.2), size = 2, color = "#000000", alpha = 0.7, shape = ".") +
  theme(axis.text.x = element_text (color = "#000000", size = 6, angle = 90, vjust = 0.5)) +
  scale_y_continuous(name = "Unidentified alleles against schema (restricted scale)", limits = c(0, 50), breaks = c(0,10,20,30,40,50)) +
  scale_x_discrete(name = "Targeted read (Dr) and kmer (Dk) depth (X)") +
  theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
        strip.text.x = element_text(size=8, face = "bold"),
        strip.text.y = element_text(size=8, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Zoom in impact of reference genome and downsampling of paired-end reads (i.e. 2x150bp) \n of Listeria monocytogene (ATCC19114, ATCC19115 and ATCCBAA679) \n on outcomes of cgMLST workflows (n=420 each) performed with BIGSdb Pasteur schema \n (1748 alleles downloaded on March 08, 2021). \n INNUca assemblies (n=336) from INNUENDO cannot be performed for Dr=10X and Dr=20X") +
  theme(plot.title = element_text(size = 30, face = "bold")) +
  ggtitle("B") +
  facet_grid(reference_strain ~ workflow)
p
plot(p)
ggsave("cgMLST-merged-B-reference-zoom-unidentified.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("cgMLST-merged-B-reference-zoom-unidentified.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
dev.off()

# Numerical confirmation: Precision (i.e. IAAR/1748*100)
## add precision (%)
data_cgMLST$precision <- data_cgMLST$identical_alleles_against_reference / 1748 * 100
## pass from long to short dataframe
data_cgMLST_short <- dcast(data_cgMLST, formula = sample+reference_strain+targeted_depth~workflow, value.var = "precision")
## check variables
str(data_cgMLST_short)
## check dimensions
dim(data_cgMLST_short)
# => [1] 420   9
## per ATCC
data_cgMLST_short_ATCC19114=subset(data_cgMLST_short,reference_strain == "ATCC19114")
table_ATCC19114=ddply(data_cgMLST_short_ATCC19114, .(targeted_depth), summarize, 
      BIGSdb_mean=mean(BIGSdb), 
      BIGSdb_sd=sd(BIGSdb), 
      INNUENDO_mean=mean(INNUENDO), 
      INNUENDO_sd=sd(INNUENDO),
      GENPAT_mean=mean(GENPAT), 
      GENPAT_sd=sd(GENPAT),
      SeqSphere_mean=mean(SeqSphere), 
      SeqSphere_sd=sd(SeqSphere),
      BioNumerics_mean=mean(BioNumerics), 
      BioNumerics_sd=sd(BioNumerics),
      MentaLiST_mean=mean(MentaLiST), 
      MentaLiST_sd=sd(MentaLiST))
write.table(table_ATCC19114, file = "Precision_ATCC19114.csv", sep = ",", col.names = NA)
data_cgMLST_short_ATCC19115=subset(data_cgMLST_short,reference_strain == "ATCC19115")
table_ATCC19115=ddply(data_cgMLST_short_ATCC19115, .(targeted_depth), summarize, 
       BIGSdb_mean=mean(BIGSdb), 
       BIGSdb_sd=sd(BIGSdb), 
       INNUENDO_mean=mean(INNUENDO), 
       INNUENDO_sd=sd(INNUENDO),
       GENPAT_mean=mean(GENPAT), 
       GENPAT_sd=sd(GENPAT),
       SeqSphere_mean=mean(SeqSphere), 
       SeqSphere_sd=sd(SeqSphere),
       BioNumerics_mean=mean(BioNumerics), 
       BioNumerics_sd=sd(BioNumerics),
       MentaLiST_mean=mean(MentaLiST), 
       MentaLiST_sd=sd(MentaLiST))
write.table(table_ATCC19115, file = "Precision_ATCC19115.csv", sep = ",", col.names = NA)
data_cgMLST_short_ATCCBAA679=subset(data_cgMLST_short,reference_strain == "ATCCBAA679")
table_ATCCBAA679=ddply(data_cgMLST_short_ATCCBAA679, .(targeted_depth), summarize, 
       BIGSdb_mean=mean(BIGSdb), 
       BIGSdb_sd=sd(BIGSdb), 
       INNUENDO_mean=mean(INNUENDO), 
       INNUENDO_sd=sd(INNUENDO),
       GENPAT_mean=mean(GENPAT), 
       GENPAT_sd=sd(GENPAT),
       SeqSphere_mean=mean(SeqSphere), 
       SeqSphere_sd=sd(SeqSphere),
       BioNumerics_mean=mean(BioNumerics), 
       BioNumerics_sd=sd(BioNumerics),
       MentaLiST_mean=mean(MentaLiST), 
       MentaLiST_sd=sd(MentaLiST))
write.table(table_ATCCBAA679, file = "Precision_ATCCBAA679.csv", sep = ",", col.names = NA)

## all ATCC included
table_ATCC=ddply(data_cgMLST_short, .(targeted_depth), summarize, 
       BIGSdb_mean=mean(BIGSdb), 
       BIGSdb_sd=sd(BIGSdb), 
       INNUENDO_mean=mean(INNUENDO), 
       INNUENDO_sd=sd(INNUENDO),
       GENPAT_mean=mean(GENPAT), 
       GENPAT_sd=sd(GENPAT),
       SeqSphere_mean=mean(SeqSphere), 
       SeqSphere_sd=sd(SeqSphere),
       BioNumerics_mean=mean(BioNumerics), 
       BioNumerics_sd=sd(BioNumerics),
       MentaLiST_mean=mean(MentaLiST), 
       MentaLiST_sd=sd(MentaLiST))
write.table(table_ATCC, file = "Precision_ATCC.csv", sep = ",", col.names = NA)

## all DrDk included
summary(data_cgMLST_short$BIGSdb)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.26   99.71  100.00   99.88  100.00  100.00 
summary(data_cgMLST_short$INNUENDO)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# => 99.66   99.71   99.83   99.79   99.83   99.83      84
summary(data_cgMLST_short$GENPAT)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.37   99.71   99.83   99.77   99.83   99.83
summary(data_cgMLST_short$SeqSphere)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.31   99.71  100.00   99.88  100.00  100.00 
summary(data_cgMLST_short$BioNumerics)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 66.19   97.77   99.89   96.38  100.00  100.00 
summary(data_cgMLST_short$MentaLiST)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 0.00   93.82   99.71   83.00   99.94  100.00

## DrDk >= 40X
### subset
str(data_cgMLST)
data_cgMLST$targeted_depth
data_cgMLST_40X = subset(data_cgMLST,data_cgMLST$targeted_depth %in% c("Dr100-Dk75","Dr90-Dk68","Dr80-Dk60","Dr70-Dk53","Dr60-Dk45","Dr50-Dk38","Dr40-Dk31"))
dim(data_cgMLST_40X)
### pass from long to short dataframe
data_cgMLST_40X_short <- dcast(data_cgMLST_40X, formula = sample+reference_strain+targeted_depth~workflow, value.var = "precision")
### check variables
str(data_cgMLST_40X_short)
### check dimensions
dim(data_cgMLST_40X_short)
### all ATCC included >= 40X
table_ATCC_40=ddply(data_cgMLST_40X_short, .(targeted_depth), summarize, 
   BIGSdb_mean=mean(BIGSdb), 
   BIGSdb_sd=sd(BIGSdb), 
   INNUENDO_mean=mean(INNUENDO), 
   INNUENDO_sd=sd(INNUENDO),
   GENPAT_mean=mean(GENPAT), 
   GENPAT_sd=sd(GENPAT),
   SeqSphere_mean=mean(SeqSphere), 
   SeqSphere_sd=sd(SeqSphere),
   BioNumerics_mean=mean(BioNumerics), 
   BioNumerics_sd=sd(BioNumerics),
   MentaLiST_mean=mean(MentaLiST), 
   MentaLiST_sd=sd(MentaLiST))
write.table(table_ATCC, file = "Precision_ATCC_40X.csv", sep = ",", col.names = NA)
### synthesis
summary(data_cgMLST_40X_short$BIGSdb)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.66   99.71  100.00   99.90  100.00  100.00
summary(data_cgMLST_40X_short$INNUENDO)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.66   99.71   99.83   99.79   99.83   99.83
summary(data_cgMLST_40X_short$GENPAT)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.66   99.71   99.83   99.79   99.83   99.83
summary(data_cgMLST_40X_short$SeqSphere)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.66   99.71  100.00   99.90  100.00  100.00
summary(data_cgMLST_40X_short$BioNumerics)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 97.65   97.83  100.00   99.25  100.00  100.00
summary(data_cgMLST_40X_short$MentaLiST)
# => Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# => 99.26   99.71   99.94   99.86  100.00  100.00

#### Additional non parametric tests between strains (i.e. datarownameNArow: all cgMLST workflows excepted MentaLIST which does not present assembly output) ####

# switch from long to short dataframe for Wilcoxon tests
str(datarownameNArow)
dim(datarownameNArow)
# => [1] 2016   62
datarownameNArow_subset = subset(datarownameNArow,datarownameNArow$targeted_depth %in% c("Dr100-Dk75", "Dr90-Dk68", "Dr80-Dk60", "Dr70-Dk53", "Dr60-Dk45", "Dr50-Dk38", "Dr40-Dk31", "Dr30-Dk23"))
dim(datarownameNArow_subset)
# => [1] 1680   62

## differences of GC% (GC)
datarownameNArow_subset_short_GC <- dcast(datarownameNArow_subset, formula = workflow+sample_origin+DNA_extraction_replicate+sequencing_replicate+targeted_depth~reference_strain, value.var = "GC")
str(datarownameNArow_subset_short_GC)
dim(datarownameNArow_subset_short_GC)
# => [1] 560   8
mean(datarownameNArow_subset_short_GC$ATCC19114, na.rm = TRUE)
# => [1] 38.08109
sd(datarownameNArow_subset_short_GC$ATCC19114, na.rm = TRUE)
# => [1] 0.007570248
mean(datarownameNArow_subset_short_GC$ATCC19115, na.rm = TRUE)
# => [1] 37.8798
sd(datarownameNArow_subset_short_GC$ATCC19115, na.rm = TRUE)
# => [1] 0.006617372
mean(datarownameNArow_subset_short_GC$ATCCBAA679, na.rm = TRUE)
# => 37.86525
sd(datarownameNArow_subset_short_GC$ATCCBAA679, na.rm = TRUE)
# => [1] 0.006919548
wilcox.test(datarownameNArow_subset_short_GC$ATCC19114, datarownameNArow_subset_short_GC$ATCC19115, paired = FALSE, na.rm = TRUE)
# => W = 313600, p-value < 2.2e-16
wilcox.test(datarownameNArow_subset_short_GC$ATCC19114, datarownameNArow_subset_short_GC$ATCCBAA679, paired = FALSE, na.rm = TRUE)
# => W = 313600, p-value < 2.2e-16
wilcox.test(datarownameNArow_subset_short_GC$ATCC19115, datarownameNArow_subset_short_GC$ATCCBAA679, paired = FALSE, na.rm = TRUE)
# => W = 285422, p-value < 2.2e-16

## differences of duplication ratio (DR)
datarownameNArow_subset_short_DR <- dcast(datarownameNArow_subset, formula = workflow+sample_origin+DNA_extraction_replicate+sequencing_replicate+targeted_depth~reference_strain, value.var = "DR")
str(datarownameNArow_subset_short_DR)
dim(datarownameNArow_subset_short_DR)
# => [1] 560   8
mean(datarownameNArow_subset_short_DR$ATCC19114, na.rm = TRUE)
# => [1] 1.000163
sd(datarownameNArow_subset_short_DR$ATCC19114, na.rm = TRUE)
# => [1] 0.0003692387
mean(datarownameNArow_subset_short_DR$ATCC19115, na.rm = TRUE)
# => [1] 1.00018
sd(datarownameNArow_subset_short_DR$ATCC19115, na.rm = TRUE)
# => [1] 0.0003940159
mean(datarownameNArow_subset_short_DR$ATCCBAA679, na.rm = TRUE)
# => 1.000061
sd(datarownameNArow_subset_short_DR$ATCCBAA679, na.rm = TRUE)
# => [1] 0.0008991506
wilcox.test(datarownameNArow_subset_short_DR$ATCC19114, datarownameNArow_subset_short_DR$ATCC19115, paired = FALSE, na.rm = TRUE)
# => W = 154469, p-value = 0.5078
wilcox.test(datarownameNArow_subset_short_DR$ATCC19114, datarownameNArow_subset_short_DR$ATCCBAA679, paired = FALSE, na.rm = TRUE)
# => W = 313600, W = 162778, p-value = 0.07424
wilcox.test(datarownameNArow_subset_short_DR$ATCC19115, datarownameNArow_subset_short_DR$ATCCBAA679, paired = FALSE, na.rm = TRUE)
# => W = 165075, p-value = 0.01535

