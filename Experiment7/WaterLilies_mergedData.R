### LOAD LIBRARIES ####
library(phyloseq)
library(plyr)
library(data.table)
library(tidyverse)
library(vegan)
library(ecodist)
library(patchwork)
library(RColorBrewer)
library(ape)
library(microbiome)


### IMPORT DATA ####
# (this has already been done so skip to ### LOAD DATA)

# combine datasets from 16S, 18S and ITS from water lily experiment 

# import 16S data `

count.16 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/feature-table_16S.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")
count.18 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/feature-table_18S.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")
count.ITS <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/feature-table_ITS.tsv", header=FALSE, row.names=NULL, stringsAsFactors=FALSE, na.strings = "n/a")

# Re-classify and re-format

colnames(count.16)<-as.character(unlist(count.16[2,]))
count.16<-count.16[-c(1,2),]
colnames(count.16)[1]<-"Feature.ID"
head(count.16[1:2,]) # check dataframe structure
x<-dim(count.16)[2];x # number of columns
# Convert class in count1 to numerics
count.16<-count.16
x<-dim(count.16)[2];x # number of columns
count.16[2:x] <- lapply(count.16[2:x], function(x) as.numeric(as.character(x)))

colnames(count.18)<-as.character(unlist(count.18[2,]))
count.18<-count.18[-c(1,2),]
colnames(count.18)[1]<-"Feature.ID"
head(count.18[1:2,]) # check dataframe structure
x<-dim(count.18)[2];x # number of columns
# Convert class in count1 to numerics
count.18<-count.18
x<-dim(count.18)[2];x # number of columns
count.18[2:x] <- lapply(count.18[2:x], function(x) as.numeric(as.character(x)))

colnames(count.ITS)<-as.character(unlist(count.ITS[2,]))
count.ITS<-count.ITS[-c(1,2),]
colnames(count.ITS)[1]<-"Feature.ID"
head(count.ITS[1:2,]) # check dataframe structure
x<-dim(count.ITS)[2];x # number of columns
# Convert class in count1 to numerics
count.ITS<-count.ITS
x<-dim(count.ITS)[2];x # number of columns
count.ITS[2:x] <- lapply(count.ITS[2:x], function(x) as.numeric(as.character(x)))

count.all <- rbind(count.16, count.18, count.ITS)

# Get base stats and generate a table called "dataset_info"

seq_total<-apply(count.all[2:x],2,sum) #number of sequences per sample
OTU_count<-colSums(count.all[2:x]>0) # OTUs per sample
OTU_single<-colSums(count.all[2:x]==1) # OTUs with only 1 seq
OTU_double<-colSums(count.all[2:x]==2) # OTUs that have only 2 seqs
OTU_true<-colSums(count.all[2:x]>2) # Number of OTUs with >2 seqs
dataset_info<-data.frame(seq_total,OTU_count,OTU_single,OTU_double,OTU_true)

# Totals for whole dataset:
sum(seq_total) # total number of sequences over whole dataset
dim(count.all)[1] # total unique OTUs over whole dataset

# Option to remove global singletons
rowsum<-apply(count.all[2:x],1,sum) # Number of sequences per OTU (for whole dataset)
count.no1 = count.all[ rowsum>1, ] # Remove rowsums =< 1 (only one sequence for 1 OTU)
dim(count.all)[1] - dim(count.no1)[1] # Total number of global singletons removed from data
## Switch values above to remove global doubletons as well!

# Import taxonomy file and join with count data:
taxname.16 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/taxonomy_16S.tsv", header=TRUE, row.names=NULL)
taxname.16 <- taxname.16 %>% dplyr::rename(Consensus = Confidence)
taxname.18 <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/taxonomy_18S.tsv", header=TRUE, row.names=NULL)
taxname.ITS <-read.delim("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/taxonomy_ITS.tsv", header=TRUE, row.names=NULL)

taxname.all <- rbind(taxname.16, taxname.18, taxname.ITS)
taxname.all <- separate(data = taxname.all, col = Taxon, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";")


counts_wtax <-join(count.no1, taxname.all, by="Feature.ID", type="left", match="first")
# save as excel and then remove __ in excel because it is easier, in excel, I laso had to find all the bacteria_x, archae_x, and terrabacteria and move the columns so the phylums matched with the others and had to change phyla so they mtches (i.e. actinobacteria and actinobacteriota)
fwrite(counts_wtax,file = "/Users/briannepalmer/Desktop/Dataframes_to_merge/merged.csv")

### LOAD DATA ####

counts_wtax <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/merged_edited.csv")

metadata<-fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/NymphaeaNuphar_18S/sample-metadata.tsv")
metadata$Location <-c("bottom", "bottom", "bottom", "top", "top", "top", "bottom", "bottom", "bottom", "top", "top", "top")
metadata$plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")
metadata$treatment <- c("NuBot", "NuBot", "NuBot", "NuTop", "NuTop", "NuTop", "NyBot", "NyBot", "NyBot", "NyTop", "NyTop", "NyTop")

# make new columns with count and relative abundance 

counts_wtax$n <- rowSums(counts_wtax[ , c(2:13)], na.rm=TRUE)
counts_wtax <- counts_wtax %>% mutate(rel = n / sum(n)) 
x <- seq(1:8786)
otu.names <- as.vector(paste0("ASV", x))
counts_wtax$otu <- otu.names

tax.only <- counts_wtax[,c(1, 14:20)]
count.only <- counts_wtax[,c(2:13)]

otu.table <- counts_wtax[,c(24, 2:13)]
otu.table <- otu.table[,-1]
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- counts_wtax$otu

fwrite(otu.table, file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesOTU.csv")

tax.table <- counts_wtax[,c(24,14:20)]
tax.table <- as.matrix(tax.table)
rownames(tax.table) <- counts_wtax$otu
tax.table <- tax.table[,-1]

fwrite(tax.table, file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesTAX.csv")

# make relative abundance df
count.only <-sweep(count.only, 2, colSums(count.only), `/`)
count.only[c(1:12)] <- sapply(count.only[c(1:12)], as.numeric)

# make OTUs columns 
transformed <- t(count.only)
colnames(transformed) <- tax.only$Feature.ID
transformed[is.na(transformed)] <-0
transformed[c(1:8786)] <- sapply(transformed[c(1:8786)],as.numeric)
transformed <- as.matrix(transformed)
transformed <- matrix(
  as.numeric(transformed), ncol = 8786) 

# prokaryote only 
prokaryote <- counts_wtax %>% filter(kingdom == "Bacteria")
transformed.p <- t(prokaryote[,c(2:13)])
colnames(transformed.p) <- prokaryote$Feature.ID
transformed.p[is.na(transformed.p)] <-0
transformed.p <- as.matrix(transformed.p)

# eukaryote only 
eukaryote <- counts_wtax %>% filter(kingdom == "Eukaryota")
transformed.e <- t(eukaryote[,c(2:13)])
colnames(transformed.e) <- eukaryote$Feature.ID
transformed.e[is.na(transformed.e)] <-0
transformed.e <- as.matrix(transformed.e)

# archaea only 
archaea <- counts_wtax %>% filter(kingdom == "Archaea")
transformed.a <- t(archaea[,c(2:13)])
colnames(transformed.a) <- archaea$Feature.ID
transformed.a[is.na(transformed.a)] <-0
transformed.a <- as.matrix(transformed.a)

### MAKE PCOA PLOTS ####
otu.bray.p <- vegdist(transformed.p, method = "bray")
otu.bray.e <- vegdist(transformed.e, method = "bray")

pcoaVS.p <- ecodist::pco(otu.bray.p, negvals = "zero", dround = 0)
pcoaVS.e <- ecodist::pco(otu.bray.e, negvals = "zero", dround = 0)

# Make dataframe for plot with vectors and metadata

pcoa1.p = pcoaVS.p$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2.p = pcoaVS.p$vectors[,2] # adds the points of the NMDS 2 dimmension
pcoa1.e = pcoaVS.e$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2.e = pcoaVS.e$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1.p = pcoa1.p, pcoa2.p = pcoa2.p, pcoa1.e = pcoa1.e, pcoa2.e = pcoa2.e)
MDS = cbind(MDS, metadata)

set.seed(12345)

p <- ggplot(MDS, aes(x = pcoa1.p, y = pcoa2.p, color = Location, shape = plant)) + geom_point(size =3) + theme_bw() + ggtitle("A) Prokaryotes") + theme(legend.position = "none") + labs(x = "PCOA1", y = "PCOA2")
e <- ggplot(MDS, aes(x = pcoa1.e, y = pcoa2.e, color = Location, shape = plant)) + geom_point(size =3) + theme_bw()+ ggtitle("B) Eukaryotes")+ labs(x = "PCOA1", y = "PCOA2")


p + e

adonis2(transformed.p ~ Location, metadata) # 0.024
adonis2(transformed.e ~ plant, metadata) # 0.61

### MAKE PHYLUM REL ABUND PLOTS ####

pivot <- counts_wtax %>% filter(kingdom != "Unassigned", kingdom != "Alveolata") %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot[pivot == ""] <- NA   
pivot.phylum <- pivot %>% group_by(SampleID, kingdom, phylum) %>% dplyr::summarize(value = sum(value))%>% mutate(rel = value / sum(value)) 


phylum <- pivot.phylum %>% dplyr::select(SampleID, kingdom, phylum, rel)
phylum <- phylum %>% pivot_wider(names_from = phylum, values_from = rel)
phylum[c(3:65)] <- sapply(phylum[c(3:65)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = Crenarchaeota:Trebouxiophyceae,names_to = "phylum", values_to = "rel") 

nb.cols <- 50
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

p <- ggplot(phylum.ggplot %>% filter(kingdom == "Bacteria", rel >0), aes(fill = phylum, y = rel, x = treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + theme_bw() + scale_fill_manual(values = mycolors) + ggtitle("A) Prokaryotic Phyla")

nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
e <- ggplot(phylum.ggplot %>% filter(kingdom == "Eukaryota", rel > 0), aes(fill = phylum, y = rel, x = treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + ggtitle("B) Eukaryotic Phyla") + theme_bw()+ scale_fill_manual(values = mycolors) 

nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
a <- ggplot(phylum.ggplot %>% filter(kingdom == "Archaea", rel >0 ), aes(fill = phylum, y = rel, x = treatment)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + ggtitle("C) Archaea Phyla") + theme_bw()+ scale_fill_manual(values = mycolors) 

p+e+a

count.matrix <- as.matrix(count.only)
heatmap(count.matrix, scale = "row")

### heat map and plot for class ####
pivot.class <- pivot %>% group_by(SampleID, class) %>% dplyr::summarize(value = sum(value))%>% mutate(rel = value / sum(value)) 

class <- pivot.class %>% dplyr::select(SampleID, class, value)
class <- class %>% pivot_wider(names_from = class, values_from = value)
class[c(2:147)] <- sapply(class[c(2:147)],as.numeric)
class[is.na(class)] <- 0

class <- as.data.frame(c(metadata, class))

# class without NAs
class.ggplot <- pivot_longer(class, cols = Abditibacteria:NA., names_to = "class", values_to = "count") 

class.order <- class.ggplot %>% group_by(class) %>% dplyr::summarize(count = sum(count)) %>% mutate(rel = count/sum(count)) %>% arrange(desc(rel))
top20class <- class.order[c(1:20),]

class.20 <- class.ggplot %>% filter(class %in% top20class$class)

class.20 <- class.20 %>% pivot_wider(names_from = class, values_from = count)

class.20.matrix <- as.matrix(class.20[c(12:31)])
rownames(class.20.matrix) <- class.20$X.SampleID

pdf("experiment7_heatmap_class", width = 8, height = 7)
heatmap(class.20.matrix, Rowv = NA,
        Colv = NA, margins =c(18,14))
dev.off()

### phylum anova ####

phylum.anova <- phylum.ggplot %>% group_by(phylum) %>%
  do(broom::tidy(TukeyHSD(aov(rel ~ Location, data = .)))) %>%
  ungroup

phylum.anova.sig <- phylum.anova %>% filter(adj.p.value < 0.05) # nothing 

### class anova ####
class.anova <- class.ggplot %>% group_by(class) %>%
  do(broom::tidy(TukeyHSD(aov(count ~ Location*plant, data = .)))) %>%
  ungroup

class.anova.sig <- class.anova %>% filter(adj.p.value < 0.05) 

### LOAD OTU AND TAX TABLES ####
otu.table <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesOTU.csv")
tax.table <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesTAX.csv")
### AMPVIS2 ####

d <- amp_load(otutable = otu.table,
              metadata = metadata,
              taxonomy = tax.table)

amp_heatmap(data = d, 
            group = "Location",
            tax_show = 25,
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
            tax_show = 50,
            functions = c("Filamentous", "Chemoautotroph_mixotroph", "AOB", "NOB", "Aerobic heterotroph", "Nitrite reduction", "Sulfate reduction", "Fermentation", "Acetogen", "Short-chain fatty acids", "Sugars", "Proteins_amino acids"))

### PHYLOSEQ ####

otu.table <- otu_table(otu.table, taxa_are_rows = TRUE)
tax.table <- tax_table(tax.table)
sam.data <- as.matrix(metadata)
sam.data <- sample_data(metadata)
rownames(sam.data) <- metadata$`#SampleID`

ps <- phyloseq(otu.table, tax.table, sam.data)
random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps = merge_phyloseq(ps,random_tree)
ps <- subset_taxa(ps, kingdom !="Unassigned")
ps <- subset_taxa(ps, phylum != "NA")

plot_richness(ps, x="Location", measures=c("Shannon", "Simpson")) + geom_boxplot(aes(fill = Location)) + theme_bw()



### richness and diversity ####

pivot <- counts_wtax %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot[pivot == ""] <- NA   
pivot.ASV <- pivot %>% group_by(SampleID, otu) %>% dplyr::summarize(value = sum(value)) 


ASV <- pivot.ASV %>% dplyr::select(SampleID, otu, value)
ASV <- ASV %>% pivot_wider(names_from = otu, values_from = value)
ASV[c(2:8787)] <- sapply(ASV[c(2:8787)],as.numeric)
ASV[is.na(ASV)] <- 0

ASV <- as.data.frame(c(metadata, ASV))
ASV.all <- pivot_longer(ASV, cols = ASV1:ASV999, names_to = "ASV", values_to = "value") 
ASV.all <- ASV.all %>% pivot_wider(names_from = ASV, values_from = value)

nutop.rich <- specnumber((ASV.all %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)])
nubot.rich <- specnumber((ASV.all %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)])
nytop.rich <- specnumber((ASV.all %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)])
nybot.rich <- specnumber((ASV.all %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)])

nutop.div <- vegan::diversity((ASV.all %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nubot.div <- vegan::diversity((ASV.all %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nytop.div <- vegan::diversity((ASV.all %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)], index = "shannon")
nybot.div <- vegan::diversity((ASV.all %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)],index = "shannon")

richness <- c(nubot.rich, nutop.rich, nybot.rich, nytop.rich)
diversity <- c(nubot.div, nutop.div, nybot.div, nytop.div)

rich.div <- data.frame(ASV.all$Location, ASV.all$plant, ASV.all$treatment, richness, diversity)

rich.lm <- lm(richness ~ ASV.all.Location, data = rich.div)
rich.aov <- aov(rich.lm)
summary(rich.aov)# there is a difference in OTU richness 

div.lm <- lm(diversity ~ ASV.all.Location, data = rich.div)
div.aov <- aov(div.lm)
summary(div.aov)# there is no difference in diversity 

richness.plot <- ggplot(data = rich.div, aes(x = ASV.all.Location, y = richness, fill = ASV.all.Location)) + geom_boxplot() + theme_bw() + theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "# of ASVs", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + facet_wrap(ASV.all.plant ~ .)+ theme_bw() + theme(legend.position = "none") + ggtitle("A) Richness")

diversity.plot <- ggplot(data = rich.div, aes(x = ASV.all.Location, y = diversity, fill = ASV.all.Location)) + geom_boxplot()+ theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "Shannon-Weiner Diversity", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ facet_wrap(ASV.all.plant ~ .) + theme_bw() + theme(legend.position = "none") + ggtitle("B) Diversity")
richness.plot + diversity.plot

### richness and diversity for bacteria only ####

pivot.bact <- counts_wtax %>% filter(kingdom == "Bacteria") %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.bact[pivot.bact == ""] <- NA   
pivot.ASV.bact <- pivot.bact %>% group_by(SampleID, otu) %>% dplyr::summarize(value = sum(value)) 


ASV.bact <- pivot.ASV.bact %>% dplyr::select(SampleID, otu, value)
ASV.bact <- ASV.bact %>% pivot_wider(names_from = otu, values_from = value)
ASV.bact[c(2:6184)] <- sapply(ASV.bact[c(2:6184)],as.numeric)
ASV.bact[is.na(ASV.bact)] <- 0

ASV.bact <- as.data.frame(c(metadata, ASV.bact))
ASV.all.bact <- pivot_longer(ASV.bact, cols = ASV10:ASV999, names_to = "ASV", values_to = "value") 
ASV.all.bact <- ASV.all.bact %>% pivot_wider(names_from = ASV, values_from = value)

nutop.rich <- specnumber((ASV.all.bact %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)])
nubot.rich <- specnumber((ASV.all.bact  %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)])
nytop.rich <- specnumber((ASV.all.bact  %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)])
nybot.rich <- specnumber((ASV.all.bact  %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)])

nutop.div <- vegan::diversity((ASV.all.bact %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nubot.div <- vegan::diversity((ASV.all.bact  %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nytop.div <- vegan::diversity((ASV.all.bact  %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)], index = "shannon")
nybot.div <- vegan::diversity((ASV.all.bact  %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)],index = "shannon")

richness <- c(nubot.rich, nutop.rich, nybot.rich, nytop.rich)
diversity <- c(nubot.div, nutop.div, nybot.div, nytop.div)

rich.div <- data.frame(ASV.all.bact$Location, ASV.all.bact$plant, ASV.all.bact$treatment, richness, diversity)

rich.lm <- lm(richness ~ ASV.all.bact$Location, data = rich.div)
rich.aov <- aov(rich.lm)
summary(rich.aov)# there is a difference in OTU richness 

div.lm <- lm(diversity ~ ASV.all.bact$Location, data = rich.div)
div.aov <- aov(div.lm)
summary(div.aov)# there is no difference in diversity 

richness.plot <- ggplot(data = rich.div, aes(x = ASV.all.bact.Location, y = richness, fill = ASV.all.bact.Location)) + geom_boxplot() + theme_bw() + theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "# of ASVs", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + facet_wrap(ASV.all.bact.plant ~ .)+ theme_bw() + theme(legend.position = "none") + ggtitle("A) Richness")

diversity.plot <- ggplot(data = rich.div, aes(x = ASV.all.bact.Location, y = diversity, fill = ASV.all.bact.Location)) + geom_boxplot()+ theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "Shannon-Weiner Diversity", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ facet_wrap(ASV.all.bact.plant ~ .) + theme_bw() + theme(legend.position = "none") + ggtitle("B) Diversity")

richness.plot + diversity.plot

### richness and diversity for eukaryota only ####

pivot.euk <- counts_wtax %>% filter(kingdom == "Eukaryota") %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.euk[pivot.euk == ""] <- NA   
pivot.ASV.euk <- pivot.euk %>% group_by(SampleID, otu) %>% dplyr::summarize(value = sum(value)) 


ASV.euk <- pivot.ASV.euk %>% dplyr::select(SampleID, otu, value)
ASV.euk <- ASV.euk %>% pivot_wider(names_from = otu, values_from = value)
ASV.euk[c(2:596)] <- sapply(ASV.euk[c(2:596)],as.numeric)
ASV.euk[is.na(ASV.euk)] <- 0

ASV.euk <- as.data.frame(c(metadata, ASV.euk))
ASV.all.euk <- pivot_longer(ASV.euk, cols = ASV1:ASV8786, names_to = "ASV", values_to = "value") 
ASV.all.euk <- ASV.all.euk %>% pivot_wider(names_from = ASV, values_from = value)

nutop.rich <- specnumber((ASV.all.euk %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)])
nubot.rich <- specnumber((ASV.all.euk %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)])
nytop.rich <- specnumber((ASV.all.euk %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)])
nybot.rich <- specnumber((ASV.all.euk %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)])

nutop.div <- vegan::diversity((ASV.all.euk %>% filter(Location == "top", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nubot.div <- vegan::diversity((ASV.all.euk %>% filter(Location == "bottom", plant == "Nuphar"))[-c(1:11)], index = "shannon")
nytop.div <- vegan::diversity((ASV.all.euk %>% filter(Location == "top", plant == "Nymphaea"))[-c(1:11)], index = "shannon")
nybot.div <- vegan::diversity((ASV.all.euk %>% filter(Location == "bottom", plant == "Nymphaea"))[-c(1:11)],index = "shannon")

richness <- c(nubot.rich, nutop.rich, nybot.rich, nytop.rich)
diversity <- c(nubot.div, nutop.div, nybot.div, nytop.div)

rich.div <- data.frame(ASV.all.euk$Location, ASV.all.euk$plant, richness, diversity)

rich.lm <- lm(richness ~ ASV.all.euk.Location, data = rich.div)
rich.aov <- aov(rich.lm)
summary(rich.aov)# there is a difference in OTU richness 

div.lm <- lm(diversity ~ ASV.all.euk.Location, data = rich.div)
div.aov <- aov(div.lm)
summary(div.aov)# there is no difference in diversity 

richness.plot <- ggplot(data = rich.div, aes(x = ASV.all.euk.Location, y = richness, fill = ASV.all.euk.Location)) + geom_boxplot() + theme_bw() + theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "# of ASVs", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + facet_wrap(ASV.all.euk.plant ~ .)+ theme_bw() + theme(legend.position = "none") + ggtitle("A) Richness")

diversity.plot <- ggplot(data = rich.div, aes(x = ASV.all.euk.Location, y = diversity, fill = ASV.all.euk.Location)) + geom_boxplot()+ theme(legend.position = "none") + theme(text=element_text(size=20)) + labs(y = "Shannon-Weiner Diversity", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+ facet_wrap(ASV.all.euk.plant ~ .) + theme_bw() + theme(legend.position = "none") + ggtitle("B) Diversity")

richness.plot + diversity.plot

### indicator species ####

library(indicspecies)

# indicators with all OTUS ####
pc = read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesOTU.csv", header = TRUE)
pc$otu <- counts_wtax$otu
# find OTUS that are actually assigned to the kingdom level from the previously created ps 

otu.assigned <- rownames(ps@otu_table)
pc <- pc %>% filter(otu %in% otu.assigned)
otu.names.ass <- pc$otu
pc<- pc[,-13]
pc = t(pc)
pc <- as.matrix(pc)
colnames(pc) <- otu.names.ass
pc <- as.data.frame(pc)
pc$Location <- as.vector(metadata$Location)
pc$plant <- metadata$plant
pc <- pc %>% dplyr::select(Location, plant, ASV1:ASV6370)

abund = pc[,3:ncol(pc)]
location = pc$Location

# indicator species
set.seed(12345)
inv = multipatt(abund, location,func = "r.g", control = how(nperm=9999), duleg = TRUE)
summary(inv, indvalcomp = TRUE)

bottom.inv <- c("ASV4158",
                "ASV5701",
                "ASV3596",
                "ASV3512",
                "ASV3546",
                "ASV5513",
                "ASV1712",
                "ASV2629",
                "ASV2650",
                "ASV3159",
                "ASV3369",
                "ASV2135",
                "ASV4540",
                "ASV4594",
                "ASV4162",
                "ASV4432",
                "ASV5333",
                "ASV4249",
                "ASV4170",
                "ASV1696",
                "ASV4177",
                "ASV1764",
                "ASV6088",
                "ASV5491",
                "ASV4125",
                "ASV6013",
                "ASV947",
                "ASV99",
                "ASV3377",
                "ASV110",
                "ASV176",
                "ASV4464",
                "ASV4095",
                "ASV746",
                "ASV3591",
                "ASV2802",
                "ASV756",
                "ASV3580",
                "ASV958",
                "ASV3027",
                "ASV5358",
                "ASV750",
                "ASV624",
                "ASV100",
                "ASV5199",
                "ASV6091",
                "ASV6016",
                "ASV951",
                "ASV3167",
                "ASV1293",
                "ASV952",
                "ASV5895",
                "ASV5844",
                "ASV5838",
                "ASV760",
                "ASV803",
                "ASV6274",
                "ASV1010",
                "ASV3579",
                "ASV2783",
                "ASV749",
                "ASV702",
                "ASV342",
                "ASV5139",
                "ASV634",
                "ASV2964",
                "ASV260",
                "ASV3368",
                "ASV2712",
                "ASV690",
                "ASV801",
                "ASV4216",
                "ASV3075",
                "ASV1007",
                "ASV4446",
                "ASV23",
                "ASV3590",
                "ASV779",
                "ASV641",
                "ASV4541",
                "ASV3157",
                "ASV5354",
                "ASV5338",
                "ASV3367",
                "ASV61",
                "ASV3177",
                "ASV700",
                "ASV5494",
                "ASV3156",
                "ASV713",
                "ASV5261",
                "ASV4191",
                "ASV5865",
                "ASV984",
                "ASV5372",
                "ASV781",
                "ASV3155",
                "ASV6011",
                "ASV6086",
                "ASV2939",
                "ASV778",
                "ASV5497",
                "ASV708",
                "ASV4707",
                "ASV2062",
                "ASV2115",
                "ASV5528",
                "ASV6028",
                "ASV4827",
                "ASV980",
                "ASV1672",
                "ASV707",
                "ASV992",
                "ASV1657",
                "ASV3182",
                "ASV6015",
                "ASV684",
                "ASV6090",
                "ASV761",
                "ASV3392",
                "ASV4440",
                "ASV4487",
                "ASV4231",
                "ASV4185",
                "ASV941",
                "ASV4599",
                "ASV618",
                "ASV688",
                "ASV3179",
                "ASV5866",
                "ASV1210",
                "ASV5204",
                "ASV963",
                "ASV5688",
                "ASV3356",
                "ASV5642",
                "ASV5882",
                "ASV3158",
                "ASV639",
                "ASV752",
                "ASV966",
                "ASV622",
                "ASV6017",
                "ASV955",
                "ASV3178",
                "ASV6014",
                "ASV5860",
                "ASV764",
                "ASV3358",
                "ASV5845",
                "ASV938",
                "ASV6089",
                "ASV5207",
                "ASV736",
                "ASV4482",
                "ASV266",
                "ASV3146",
                "ASV4228",
                "ASV3152",
                "ASV347",
                "ASV2134",
                "ASV3168",
                "ASV3160",
                "ASV1213",
                "ASV4463",
                "ASV3573",
                "ASV3653",
                "ASV1261",
                "ASV982",
                "ASV734",
                "ASV3685",
                "ASV4200",
                "ASV3683",
                "ASV935",
                "ASV17",
                "ASV5206",
                "ASV2202",
                "ASV2308",
                "ASV2216",
                "ASV2321",
                "ASV4593",
                "ASV3572",
                "ASV5036",
                "ASV3737",
                "ASV6018",
                "ASV274",
                "ASV355",
                "ASV5489",
                "ASV2314",
                "ASV2208",
                "ASV5332",
                "ASV4465",
                "ASV4207")

top.inv <- c("ASV4189",
             "ASV4445",
             "ASV6341",
             "ASV4064",
             "ASV257",
             "ASV338",
             "ASV4078",
             "ASV2196",
             "ASV2136",
             "ASV2510",
             "ASV402",
             "ASV4174",
             "ASV4429",
             "ASV1556",
             "ASV1207",
             "ASV751",
             "ASV1209",
             "ASV1256",
             "ASV1258",
             "ASV427",
             "ASV395",
             "ASV1582",
             "ASV418",
             "ASV1550",
             "ASV392",
             "ASV954",
             "ASV2656",
             "ASV754",
             "ASV934",
             "ASV732",
             "ASV549",
             "ASV2108",
             "ASV2055",
             "ASV2119",
             "ASV1053",
             "ASV1593",
             "ASV953",
             "ASV2693",
             "ASV2066",
             "ASV1551",
             "ASV927",
             "ASV930",
             "ASV725",
             "ASV1581",
             "ASV729",
             "ASV1968",
             "ASV1954",
             "ASV358",
             "ASV929",
             "ASV565",
             "ASV727",
             "ASV1587",
             "ASV1260",
             "ASV4058",
             "ASV4428",
             "ASV1211",
             "ASV1100",
             "ASV1099",
             "ASV278",
             "ASV4075",
             "ASV1970",
             "ASV1991",
             "ASV4090",
             "ASV337",
             "ASV3352",
             "ASV3652",
             "ASV731",
             "ASV4173",
             "ASV256",
             "ASV4179",
             "ASV4434",
             "ASV4431",
             "ASV4176",
             "ASV4054",
             "ASV2047",
             "ASV2007",
             "ASV4071",
             "ASV1051",
             "ASV4616",
             "ASV4615",
             "ASV1054",
             "ASV2509",
             "ASV1934",
             "ASV1992",
             "ASV4057",
             "ASV4074",
             "ASV4052",
             "ASV4069",
             "ASV3603",
             "ASV2110",
             "ASV2057")

### TOP 500 Only ####
Top500 <- top_taxa(ps, n = 500)

pc = read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/waterliliesOTU.csv", header = TRUE)
pc$otu <- otu.names
pc <- pc %>% filter(otu %in% Top500)
otu.names <- pc$otu
pc<- pc[,-13]
pc = t(pc)
pc <- as.matrix(pc)
colnames(pc) <- otu.names
pc <- as.data.frame(pc)
pc$Location <- as.vector(metadata$Location)
pc$plant <- metadata$plant
pc <- pc %>% dplyr::select(Location, plant, ASV15:ASV6339)

abund = pc[,3:ncol(pc)]
location = pc$Location

set.seed(12345)
inv = multipatt(abund, location,func = "r.g", control = how(nperm=9999), duleg = TRUE)
summary(inv, indvalcomp = TRUE)


bottom.inv <- c("ASV3512",
                "ASV3546",
                "ASV2629",
                "ASV2650",
                "ASV4432",
                "ASV5333",
                "ASV4177",
                "ASV6088",
                "ASV5491",
                "ASV6013",
                "ASV947",
                "ASV99",
                "ASV176",
                "ASV4464",
                "ASV746",
                "ASV756",
                "ASV958",
                "ASV750",
                "ASV624",
                "ASV951",
                "ASV952",
                "ASV760",
                "ASV749",
                "ASV702",
                "ASV342",
                "ASV260",
                "ASV690",
                "ASV4446",
                "ASV3590",
                "ASV3177",
                "ASV700",
                "ASV5494",
                "ASV4191",
                "ASV5865",
                "ASV6011",
                "ASV6086",
                "ASV4707",
                "ASV5528",
                "ASV684",
                "ASV6090",
                "ASV761",
                "ASV4440",
                "ASV4185",
                "ASV941",
                "ASV4599",
                "ASV618",
                "ASV688",
                "ASV1210",
                "ASV963",
                "ASV752",
                "ASV966",
                "ASV622",
                "ASV955",
                "ASV6014",
                "ASV764",
                "ASV938",
                "ASV6089",
                "ASV736",
                "ASV266",
                "ASV347",
                "ASV4463",
                "ASV3653",
                "ASV1261",
                "ASV734",
                "ASV4200",
                "ASV3683",
                "ASV935",
                "ASV2202",
                "ASV2308",
                "ASV4593",
                "ASV3572",
                "ASV5036",
                "ASV274",
                "ASV355",
                "ASV5489",
                "ASV5332",
                "ASV4465",
                "ASV4207")
top.inv <- c("ASV4189",
             "ASV4445",
             "ASV257",
             "ASV338",
             "ASV4174",
             "ASV4429",
             "ASV1207",
             "ASV751",
             "ASV1209",
             "ASV1256",
             "ASV1258",
             "ASV395",
             "ASV1582",
             "ASV392",
             "ASV954",
             "ASV754",
             "ASV934",
             "ASV732",
             "ASV549",
             "ASV2108",
             "ASV2055",
             "ASV1053",
             "ASV953",
             "ASV927",
             "ASV930",
             "ASV725",
             "ASV1581",
             "ASV729",
             "ASV1968",
             "ASV1954",
             "ASV358",
             "ASV929",
             "ASV727",
             "ASV1260",
             "ASV4058",
             "ASV4428",
             "ASV1211",
             "ASV1100",
             "ASV1099",
             "ASV278",
             "ASV4075",
             "ASV1970",
             "ASV1991",
             "ASV4090",
             "ASV337",
             "ASV3652",
             "ASV731",
             "ASV4173",
             "ASV256",
             "ASV4179",
             "ASV4434",
             "ASV4431",
             "ASV4176",
             "ASV4054",
             "ASV2047",
             "ASV2007",
             "ASV4071",
             "ASV1051",
             "ASV1934",
             "ASV1992",
             "ASV4057",
             "ASV4074",
             "ASV4052",
             "ASV4069",
             "ASV3603",
             "ASV2110",
             "ASV2057")
#####
bottom.ind <- counts_wtax %>% filter(otu %in% bottom.inv)
top.ind <- counts_wtax %>% filter(otu %in% top.inv)

pivot.bottom <- bottom.ind %>% dplyr::select(NuBotR1:species,otu) %>%  pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.bottom <- pivot.bottom %>% mutate(rel = value/sum(value)) 
pivot.bottom[pivot.bottom == ""] <- NA
metadata <- metadata %>% dplyr::rename(SampleID = `#SampleID`)
bottom.md <- merge(pivot.bottom, metadata, by = "SampleID")


pivot.top <- top.ind %>% dplyr::select(NuBotR1:species,otu) %>%  pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.top <- pivot.top %>% mutate(rel = value/sum(value)) 
pivot.top[pivot.top == ""] <- NA
top.md <- merge(pivot.top, metadata, by = "SampleID")

all.ind <- rbind(top.md, bottom.md)

ggplot(all.ind, aes(x = Location, y = rel, fill = phylum)) + geom_bar(stat = "identity") + facet_wrap(.~phylum, scales = "free") + theme_bw() + theme(legend.position = "none")

ggplot(all.ind %>% filter(phylum == "Acidobacteriota"), aes(x = Location, y = rel, fill = class)) + geom_boxplot() + facet_wrap(family~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Actinobacteriota"), aes(x = Location, y = rel, fill = class)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Ascomycota"), aes(x = Location, y = rel, fill = class)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Bacteroidota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Basidiomycota"), aes(x = Location, y = rel, fill = kingdom)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Bdellovibrionota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Campilobacterota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Chloroflexi"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(family~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Chlorophyta"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(family~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Cyanobacteria"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Deinococcota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Desulfobacterota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Firmicutes"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Gemmatimonadota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Halobacterota"), aes(x = Location, y = rel, fill = kingdom)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Myxococcota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Nitrospirae"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Patescibacteria"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Planctomycetota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(class == "Alphaproteobacteria"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(class == "Deltaproteobacteria"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(class == "Gammaproteobacteria"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()

ggplot(all.ind %>% filter(phylum == "Verrucomicrobiota"), aes(x = Location, y = rel, fill = order)) + geom_boxplot() + facet_wrap(genus~otu, scales = "free")+ theme_bw()



d <- amp_load(otutable = ps@otu_table,
              metadata = metadata,
              taxonomy = ps@tax_table)

tax_vector <- c("Acidobacteriota")
Acidobacteriota <- amp_filter_taxa(d, tax_vector = tax_vector)

p1 <- amp_heatmap(data = Acidobacteriota, 
            group = "Location",
            tax_show = 30,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)


tax_vector <- c("Actinobacteriota")
Actinobacteriota <- amp_filter_taxa(d, tax_vector = tax_vector)

p2 <- amp_heatmap(data = Actinobacteriota, 
            group = "Location",
            tax_show = 30,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)


tax_vector <- c("Ascomycota")
Ascomycota <- amp_filter_taxa(d, tax_vector = tax_vector)

p3 <- amp_heatmap(data = Ascomycota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)


tax_vector <- c("Bacteroidota")
Bacteroidota <- amp_filter_taxa(d, tax_vector = tax_vector)

p4 <- amp_heatmap(data = Bacteroidota, 
            group = "Location",
            tax_show = 40,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Basidiomycota")
Basidiomycota <- amp_filter_taxa(d, tax_vector = tax_vector)

p5 <- amp_heatmap(data = Basidiomycota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)


tax_vector <- c("Bdellovibrionota")
Bdellovibrionota <- amp_filter_taxa(d, tax_vector = tax_vector)

p6 <- amp_heatmap(data = Bdellovibrionota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Campilobacterota")
Campilobacterota <- amp_filter_taxa(d, tax_vector = tax_vector)

p7 <- amp_heatmap(data = Campilobacterota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Chloroflexi")
Chloroflexi <- amp_filter_taxa(d, tax_vector = tax_vector)

p8 <- amp_heatmap(data = Chloroflexi, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Chlorophyta")
Chlorophyta <- amp_filter_taxa(d, tax_vector = tax_vector)

p9 <- amp_heatmap(data = Chlorophyta, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Cyanobacteria")
Cyanobacteria <- amp_filter_taxa(d, tax_vector = tax_vector)

p10 <- amp_heatmap(data = Cyanobacteria, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Deinococcota")
Deinococcota <- amp_filter_taxa(d, tax_vector = tax_vector)

p11 <- amp_heatmap(data = Deinococcota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Desulfobacterota")
Desulfobacterota <- amp_filter_taxa(d, tax_vector = tax_vector)

p12 <- amp_heatmap(data = Desulfobacterota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Firmicutes")
Firmicutes <- amp_filter_taxa(d, tax_vector = tax_vector)

p13 <- amp_heatmap(data = Firmicutes, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Gemmatimonadota")
Gemmatimonadota <- amp_filter_taxa(d, tax_vector = tax_vector)

p14 <- amp_heatmap(data = Gemmatimonadota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Halobacterota")
Halobacterota <- amp_filter_taxa(d, tax_vector = tax_vector)

p15 <- amp_heatmap(data = Halobacterota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Myxococcota")
Myxococcota <- amp_filter_taxa(d, tax_vector = tax_vector)

p16 <- amp_heatmap(data = Myxococcota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Nitrospirae")
Nitrospirae <- amp_filter_taxa(d, tax_vector = tax_vector)

p17<- amp_heatmap(data = Nitrospirae, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Patescibacteria")
Patescibacteria <- amp_filter_taxa(d, tax_vector = tax_vector)

p18 <- amp_heatmap(data = Patescibacteria, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Planctomycetota")
Planctomycetota <- amp_filter_taxa(d, tax_vector = tax_vector)

p19<- amp_heatmap(data = Planctomycetota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("Proteobacteria")
Proteobacteria <- amp_filter_taxa(d, tax_vector = tax_vector)

p20 <- amp_heatmap(data = Proteobacteria, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

tax_vector <- c("SAR324_clade(Marine_group_B)")
SAR324 <- amp_filter_taxa(d, tax_vector = tax_vector)

p21 <- amp_heatmap(data = SAR324, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_na = TRUE)

tax_vector <- c("Verrucomicrobiota")
Verrucomicrobiota <- amp_filter_taxa(d, tax_vector = tax_vector)

p22<- amp_heatmap(data = Verrucomicrobiota, 
            group = "Location",
            tax_show = 25,
            tax_aggregate = "OTU",
            plot_colorscale = "log10",
            plot_na = TRUE)

library(patchwork)
pdf("Desktop/asvheatmap.pdf", width = 20, height =60)
(p1 + p2)/(p3 +p4)/(p5+p6)/(p7+p8)/(p9+p10)/(p11+p12)/(p13+p14)/(p15+p16)/(p17+p18)/(p19+p20)/(p21+p22)
dev.off()

### genera relative abundance plots #### 
pivot.genus <- counts_wtax %>% filter(kingdom != "Unassigned", genus != "NA", genus != "uncultured") %>% dplyr::select(NuBotR1:genus) 
pivot.genus <- pivot.genus %>% mutate(value = rowSums(pivot.genus[,1:12])) %>% mutate(rel = value/sum(value)) %>% group_by(class, genus) %>% summarize(rel = sum(rel)) %>% filter(rel > 0.01)

genera <- pivot.genus$genus

genus <- counts_wtax %>% dplyr::select(NuBotR1:genus) %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID", values_to = "value") %>% group_by(SampleID, class, genus) %>% filter(genus %in% genera) 

metadata <- metadata %>% rename(SampleID = `#SampleID`)
genus.md <- merge(genus, metadata, by = "SampleID")

ggplot(genus.md, aes(fill = genus, y = value, x = Location)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + theme_bw() + ggtitle("Genera")

genus.anova <- genus.md %>% group_by(genus) %>%
  do(broom::tidy(TukeyHSD(aov(value ~ Location, data = .)))) %>%
  ungroup

genus.anova.sig <- genus.anova %>% filter(adj.p.value < 0.05) # nothing 

#### Bug Base ####

bugbase <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/BugBaseWaterLilies.csv")

library(ggpubr)

biofilm <-ggboxplot(bugbase %>% filter(Function == c("FormsBiofilm")), 
                    x = "Treatment", y = "Rel", fill = "Treatment", palette = "jco", 
                    facet.by = "Function",short.panel.labs = FALSE) + 
  stat_compare_means(label.x = 1.5, label.y = 0.05) + theme(legend.position = "none") + labs(x= "Location", y = "Relative Abundance", title ="C")

compare_means(Rel ~ Treatment, data = bugbase, 
              group.by = "Function")

aerobic <- ggboxplot(bugbase %>% filter(Function == c("Aerobic")), x = "Treatment", y = "Rel",
                 fill = "Treatment", palette = "jco",
                 facet.by = "Function", short.panel.labs = FALSE) + theme(legend.position = "none") + labs(x= "Location", y = "Relative Abundance", title = "A") + stat_compare_means(label.x = 1.5, label.y = 0.0082)

anaerobic <-  ggboxplot(bugbase %>% filter(Function == c("Anaerobic")), x = "Treatment", y = "Rel",
                    fill = "Treatment", palette = "jco",
                    facet.by = "Function", short.panel.labs = FALSE) + theme(legend.position = "none") + labs(x= "Location", y = "Relative Abundance", title = "") + stat_compare_means(label.x = 1.3, label.y = 0.014)

facultative <-  ggboxplot(bugbase %>% filter(Function == c("FacultiveAnaerobic")), x = "Treatment", y = "Rel",
                        fill = "Treatment", palette = "jco",
                        facet.by = "Function", short.panel.labs = FALSE) + theme(legend.position = "none") + labs(x= "Location", y = "Relative Abundance", title = "") + stat_compare_means(label.x = 1.5, label.y = 0.021)


pal <- palette(c("darkgoldenrod1", "orchid"))
gram <- ggplot(bugbase %>% filter(Function == c("GramNegative", "GramPositive")), aes(x = Treatment, y = Rel/3, fill= Function)) + geom_bar(stat = "identity", position = "stack") + theme_bw()+ labs(x= "Location", y = "Relative Abundance", title = "D") + scale_fill_manual(values = pal)

stress <- ggboxplot(bugbase %>% filter(Function == c("StressTolerant")), x = "Treatment", y = "Rel",fill = "Treatment", palette = "jco", facet.by = "Function",short.panel.labs = FALSE) + stat_compare_means(label.x = 1.5, label.y = 0.04)+ theme(legend.position = "none") + labs(x= "Location", y = "Relative Abundance", title = "B")

pdf("Desktop/bugbase.pdf", height = 8, width = 12)
(aerobic + anaerobic +facultative)/(stress + biofilm + gram)
dev.off()

### Venn Diagram ####
library(ggVennDiagram)
pivot.otu <- counts_wtax %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
pivot.otu[pivot.otu == ""] <- NA   
pivot.otu <- pivot.otu %>% group_by(SampleID, otu) %>% dplyr::summarize(value = sum(value)) 
OTU <- pivot.otu %>% group_by(SampleID, otu) %>% summarize(n = sum(value))%>% mutate(rel = n / sum(n))

OTU[is.na(OTU)] <- NA
OTU <- OTU %>% pivot_wider(names_from = otu, values_from = rel)
OTU[c(3:8788)] <- sapply(OTU[c(3:8788)],as.numeric)
OTU[is.na(OTU)] <- 0

OTU <- merge(metadata, OTU, by = "SampleID")
OTU.ggplot <- pivot_longer(OTU, cols = ASV1:ASV999, names_to = "OTU", values_to = "count")

y <- OTU.ggplot %>% dplyr::select(Location, OTU, count)

Top = y %>% filter(Location == "top")  %>% filter(count >0) %>% dplyr::select(OTU) #7956
Bottom = y %>% filter(Location == "bottom")  %>% filter(count >0) %>% dplyr::select(OTU) #14133

y <- rbind(Top, Bottom)

Location <- c(rep("Top", 7956), rep("Bottom", 14133))

y$Location <- Location

x <- list(Top = y$OTU[1:7956], Bottom = y$OTU[7957:22089])

# load Venn diagram package
library(ggvenn)
smallpal2 <- palette(c("#0073C2FF", "#EFC000FF" , "#f08080", "#CD534CFF"))

ggvenn(
  x, 
  fill_color = smallpal2,
  stroke_size = 0.5, set_name_size = 4, show_percentage = TRUE
)
