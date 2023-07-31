library(data.table)
library(tidyverse)
library(plyr)

count.ITS <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/feature-table_ITS.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

count.ITS.sw <- read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/PondWater_PondSoil_GinkgoAquarium_ITS/demux/Brianne_SoilWater_ITS/ASVs/feature-table.tsv", , header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

colnames(count.ITS)<-as.character(unlist(count.ITS[2,]))
count.ITS<-count.ITS[-c(1,2),]
colnames(count.ITS)[1]<-"Feature.ID"
head(count.ITS[1:2,]) # check dataframe structure
x<-dim(count.ITS)[2];x # number of columns
# Convert class in count1 to numerics
count.ITS<-count.ITS
x<-dim(count.ITS)[2];x # number of columns
count.ITS[2:x] <- lapply(count.ITS[2:x], function(x) as.numeric(as.character(x)))

colnames(count.ITS.sw)<-as.character(unlist(count.ITS.sw[2,]))
count.ITS.sw<-count.ITS.sw[-c(1,2),]
colnames(count.ITS.sw)[1]<-"Feature.ID"
head(count.ITS.sw.sw[1:2,]) # check dataframe structure
x<-dim(count.ITS.sw)[2];x # number of columns
# Convert class in count1 to numerics
count.ITS.sw<-count.ITS.sw
x<-dim(count.ITS.sw)[2];x # number of columns
count.ITS.sw[2:x] <- lapply(count.ITS.sw[2:x], function(x) as.numeric(as.character(x)))

seq_total<-apply(count.ITS[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count.ITS[2:x]>0) # OTUs per sample
OTU_single<-colSums(count.ITS[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count.ITS[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count.ITS[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

seq_total<-apply(count.ITS.sw[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count.ITS.sw[2:x]>0) # OTUs per sample
OTU_single<-colSums(count.ITS.sw[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count.ITS.sw[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count.ITS.sw[2:x]>2) # Number of OTUs with >2 seqs
dataset_info.sw<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count.ITS.sw)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count.ITS[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count.ITS[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count.ITS)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

rowsum<-apply(count.ITS.sw[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1.sw = count.ITS.sw[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count.ITS.sw)[1] - dim(count.no1.sw)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname.ITS <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/taxonomy_ITS.tsv", header=TRUE, row.names=NULL)

taxname.ITS.sw <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/PondWater_PondSoil_GinkgoAquarium_ITS/demux/Brianne_SoilWater_ITS/tax_dir/taxonomy.tsv", header=TRUE, row.names=NULL)

taxname.ITS <- separate(data = taxname.ITS, col = Taxon, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

taxname.ITS.sw <- separate(data = taxname.ITS.sw, col = Taxon, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")


counts_wtax <-join(count.no1, taxname.ITS, by="Feature.ID", type="left", match="first")
counts_wtax.sw <-join(count.no1.sw, taxname.ITS.sw, by="Feature.ID", type="left", match="first")

fwrite(counts_wtax,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/experiment7_ITS.csv")
fwrite(counts_wtax.sw,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/experiment7_ITS_sw.csv")

##################################################################
####################################################################

counts_wtax <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization /Dataframes_to_merge/experiment7_ITS.csv")

metadata<-fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization /NymphaeaNuphar_18S/sample-metadata.tsv")
metadata$Location <-c("bottom", "bottom", "bottom", "top", "top", "top", "bottom", "bottom", "bottom", "top", "top", "top")
metadata$plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")

x <- seq(1:630)
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
tax.table[tax.table == ""] <- NA   

# make OTUs columns 
transformed <- t(count.only)
colnames(transformed) <- tax.only$Feature.ID
transformed[is.na(transformed)] <-0
transformed[c(1:630)] <- sapply(transformed[c(1:630)],as.numeric)
transformed <- as.matrix(transformed)
transformed <- matrix(
  as.numeric(transformed), ncol = 630) 

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

library(ampvis2)

d <- amp_load(otutable = otu.table,
              metadata = metadata,
              taxonomy = tax.table)

amp_heatmap(data = d, 
            group = "Location",
            tax_show = 10,
            tax_aggregate = "Class",
            plot_colorscale = "log10",
            plot_na = TRUE)


amp_heatmap (data = d, 
             tax_aggregate = "Family",
             group_by = c("Location", "plant"),
             tax_show = 25,
             plot_colorscale = "log10")

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
