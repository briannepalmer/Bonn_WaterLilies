# Import data from biom convert output


## NEED TO DO THIS STEP IN QIIME2 

#biom convert \
#-i /Users/briannepalmer/Downloads/NymphaeaNuphar_ITSRaw/ASVs/feature-table.biom \
#-o /Users/briannepalmer/Downloads/NymphaeaNuphar_ITSRaw/ASVs/feature-table.tsv  --to-tsv 
# count<-read.delim("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/ASVs/feature-table.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

# Re-classify and re-format

colnames(count)<-as.character(unlist(count[2,]))
count1<-count[-c(1,2),]
colnames(count1)[1]<-"Feature.ID"
head(count1[1:2,]) # check dataframe structure
x<-dim(count1)[2];x # number of columns
# Convert class in count1 to numerics
count2<-count1
x<-dim(count2)[2];x # number of columns
count2[2:x] <- lapply(count2[2:x], function(x) as.numeric(as.character(x)))

# Get base stats and generate a table called "dataset_info"

seq_total<-apply(count2[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count2[2:x]>0) # OTUs per sample
OTU_single<-colSums(count2[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count2[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count2[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count2)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count2[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count2[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count2)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname<-read.delim("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/ASVs/tax_dir/taxonomy.tsv", header=TRUE, row.names=NULL)
library(tidyverse)
taxname <- separate(data = taxname, col = Taxon, into = c("kingdom", "sub_kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

library(plyr) # Load plyr library into R environment
counts_wtax<-join(count.no1, taxname, by="Feature.ID", type="left")
counts_wtax$ASV <- c(paste("ASV", seq(c(1:4161))))

# Write output table with OTU or ASV cluster information and taxonomy IDs:
write.table(counts_wtax, file="/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/ASVs/OutputTable_wtax.txt", quote=FALSE, sep="\t", row.names=FALSE)

library(data.table)
metadata<-fread("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/sample-metadata.tsv")
metadata$Location <-c("bottom", "bottom", "bottom", "top", "top", "top", "bottom", "bottom", "bottom", "top", "top", "top")
metadata$plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")

# make new columns with count and relative abundance 

counts_wtax$n <- rowSums(counts_wtax[ , c(2:13)], na.rm=TRUE)
counts_wtax <- counts_wtax %>% mutate(rel = n / sum(n)) 

tax.only <- counts_wtax[c(1, 14:23),]
count.only <- counts_wtax[c(2:13, 23)]

# make relative abundance df
count.only[-13] <-sweep(count.only[-13], 2, colSums(count.only[,-13]), `/`)

# make OTUs columns 
transformed <- t(count.only)
colnames(transformed) <- transformed[13,]
transformed <- transformed[-13,]
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

adonis2(transformed ~ Location, metadata) # 0.017
adonis2(transformed ~ plant, metadata) # 0.67

# pivot 

pivot <- counts_wtax %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.phylum <- pivot %>% filter(kingdom == c("Eukaryota","Eukaryota:plas", "Eukaryota:mito")) %>% group_by(SampleID, phylum) %>% summarize(value = sum(value)) %>% mutate(rel = value / sum(value)) 

phylum <- pivot.phylum %>% dplyr::select(SampleID, phylum, rel)
phylum <- phylum %>% pivot_wider(names_from = phylum, values_from = rel)
phylum[c(2:9)] <- sapply(phylum[c(2:9)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = Fungi:Streptophyta.mito,names_to = "phylum", values_to = "count")

ggplot(phylum.ggplot, aes(fill = phylum, y = count, x = Location)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  ))


pivot <- counts_wtax %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.phylum <- pivot %>% filter(kingdom == c("Bacteria")) %>% group_by(SampleID, phylum) %>% summarize(value = sum(value)) %>% mutate(rel = value / sum(value)) 

phylum <- pivot.phylum %>% dplyr::select(SampleID, phylum, rel)
phylum <- phylum %>% pivot_wider(names_from = phylum, values_from = rel)
phylum[c(2:22)] <- sapply(phylum[c(2:22)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = Acidobacteria:Verrucomicrobia,names_to = "phylum", values_to = "count")

ggplot(phylum.ggplot, aes(fill = phylum, y = count, x = Location)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  ))



