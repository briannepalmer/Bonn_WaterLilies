library(data.table)
library(tidyverse)
library(plyr)

count.18 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/feature-table_18S.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a") 

colnames(count.18)<-as.character(unlist(count.18[2,]))
count.18<-count.18[-c(1,2),]
colnames(count.18)[1]<-"Feature.ID"
head(count.18[1:2,]) # check dataframe structure
x<-dim(count.18)[2];x # number of columns
# Convert class in count1 to numerics
count.18<-count.18
x<-dim(count.18)[2];x # number of columns
count.18[2:x] <- lapply(count.18[2:x], function(x) as.numeric(as.character(x)))

seq_total<-apply(count.18[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count.18[2:x]>0) # OTUs per sample
OTU_single<-colSums(count.18[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count.18[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count.18[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count.18)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count.18[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count.18[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count.18)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname.18 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/taxonomy_18S.tsv", header=TRUE, row.names=NULL)

taxname.18 <- separate(data = taxname.18, col = Taxon, into = c("kingdom", "sub_kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

counts_wtax <-join(count.no1, taxname.18, by="Feature.ID", type="left", match="first")
euk <- counts_wtax %>% filter(kingdom %in% c("Eukaryota", "Eukaryota:mito", "Eukaryota:plas"))
# save as excel and then remove __ in excel because it is easier, in excel, I laso had to find all the bacteria_x, archae_x, and terrabacteria and move the columns so the phylums matched with the others and had to change phyla so they mtches (i.e. actinobacteria and actinobacteriota)

fwrite(counts_wtax,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization /Dataframes_to_merge/experiment7_18S.csv")
##################################################################
####################################################################

counts_wtax <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization /Dataframes_to_merge/experiment7_18S.csv")

# the metadata 
SampleID <- c("NuBotR1", "NuBotR2", "NuBotR3", "NuTopR1", "NuTopR2", "NuTopR3", "NyBotR1", "NyBotR2", "NyBotR3", "NyTopR1", "NyTopR2", "NyTopR3")
Location <-c("buried", "buried", "buried", "floating", "floating", "floating", "buried", "buried", "buried", "floating", "floating", "floating")
plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")
treatment <- c("NuBot", "NuBot", "NuBot", "NuTop", "NuTop", "NuTop", "NyBot", "NyBot", "NyBot", "NyTop", "NyTop", "NyTop")

metadata <- data.frame(SampleID, Location, plant, treatment)

x <- seq(1:4161)
otu.names <- as.vector(paste0("ASV", x))
counts_wtax$otu <- otu.names

tax.only <- counts_wtax[,c(1, 14:20)]
count.only <- counts_wtax[,c(2:13)]

otu.table <- counts_wtax[,c(22, 2:13)]
otu.table <- otu.table[,-1]
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- counts_wtax$otu

tax.table <- counts_wtax[,c(22,14:20)]
tax.table <- as.matrix(tax.table)
rownames(tax.table) <- counts_wtax$otu
tax.table <- tax.table[,-1]

# make OTUs columns 
transformed <- t(count.only)
colnames(transformed) <- tax.only$Feature.ID
transformed[is.na(transformed)] <-0
transformed[c(1:4161)] <- sapply(transformed[c(1:4161)],as.numeric)
transformed <- as.matrix(transformed)
transformed <- matrix(
  as.numeric(transformed), ncol = 4161) 

# make otu plot 
library(vegan)
otu.bray <- vegdist(transformed, method = "bray")

library(ecodist)
pcoaVS <- pco(otu.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero

# Make dataframe for plot with vectors and metadata

pcoa1 = pcoaVS$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2 = pcoaVS$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1 = pcoa1, pcoa2 = pcoa2)
MDS = cbind(MDS, metadata)

set.seed(12345)
ggplot(MDS, aes(x = pcoa1, y = pcoa2, color = Location, shape = plant)) + geom_point(size =3) + 
  theme_bw()

adonis2(transformed ~ Location, metadata) # 0.024
adonis2(transformed ~ plant, metadata) # 0.61
adonis2(transformed[c(1:3),] ~ SampleID, metadata %>% filter(treatment == "NuBot"))
adonis2(transformed[c(4:6),] ~ SampleID, metadata %>% filter(treatment == "NuTop"))
adonis2(transformed[c(7:9),] ~ SampleID, metadata %>% filter(treatment == "NyBot"))
adonis2(transformed[c(10:12),] ~ SampleID, metadata %>% filter(treatment == "NyTop"))

Nubot.t <- transformed[c(1:3),]

library(ampvis2)

rownames(metadata) <- metadata$SampleID
ps <- phyloseq(otu_table(otu.table, taxa_are_rows = TRUE), tax_table(tax.table), sample_data(metadata))
ps.eu <- subset_taxa(ps, sub_kingdom %in% c("Archaeplastida:plas", "Stramenopiles:plas", "Opisthokonta", "Hacrobia:plas"))

otu.eu <- data.frame(ps.eu@otu_table)
tax.eu <- data.frame(ps.eu@tax_table)


d <- amp_load(otutable = otu.eu,
              metadata = metadata,
              taxonomy = tax.eu)

amp_heatmap(data = d, 
            group = "SampleID",
            tax_show = 25,
            tax_aggregate = "Genus",
            plot_colorscale = "log10",
            plot_na = TRUE)


amp_heatmap (data = d, 
             tax_aggregate = "Genus",
             tax_add = "Phylum",
             group_by = c("Location", "plant"),
             tax_show = 25,
             plot_colorscale = "log10", facet_by = "SampleID")

amp_heatmap(d,
            group_by = "Location",
            facet_by = "plant",
            tax_aggregate = "Genus",
            tax_add = "Phylum",
            tax_show = 20,
            color_vector = c("white", "darkred"),
            plot_colorscale = "sqrt",
            plot_values = FALSE) +
  theme(axis.text.x = element_text(angle = 45, size=12, vjust = 1),
        axis.text.y = element_text(size=12),
        legend.position="right")



amp_heatmap(d,
            group_by = c("Location"),
            tax_aggregate = "Genus",
            plot_functions = TRUE,
            tax_show = 25,
            functions = c("Filamentous", "Chemoautotroph_mixotroph", "AOB", "NOB", "Aerobic heterotroph", "Nitrite reduction", "Sulfate reduction", "Fermentation", "Acetogen", "Short-chain fatty acids", "Sugars", "Proteins_amino acids"))

# remove ASVs that have no information 
pivot <- counts_wtax %>% filter(kingdom != "Unassigned", kingdom != "Alveolata") %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")

# if there are blanks, these should be NA
pivot[pivot == ""] <- NA   

all.data <- merge(metadata, pivot, by = "SampleID")
