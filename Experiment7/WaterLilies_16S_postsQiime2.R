library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")
devtools::install_github("jbisanz/qiime2R")# change the ref argument to get other versions

### Load Libraries ####
library(data.table)
library(dada2)
library(beepr)
library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)
library(qiime2R)
theme_set(theme_bw())

SVs<-read_qza("/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/demux/table-dada2.qza")
SVs$data[1:5,1:5]

test <- read_qza("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/tax_vsearch.qza")

metadata<-fread("/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/sample-metadata.txt")
rownames(metadata) <- metadata$`#SampleID`
taxonomy<-read_qza("/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/demux/taxonomy.qza")
head(taxonomy$data)

physeq<-qza_to_phyloseq(
  features="/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/demux/table-dada2.qza",
  tree="/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/demux/rooted-tree-dada2.qza",
  taxonomy="/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/demux/taxonomy.qza",
  metadata = "/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/sample-metadata.txt")

physeq

plot_richness(physeq, x="Plant", measures=c("Shannon", "Simpson"))
plot_heatmap(physeq, method = "NMDS", distance = "bray", sample.label = "Location", taxa.label = "Family")

# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(physeq, function(x) x/sum(x))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Location", shape = "Plant",title="Bray NMDS")



library(ecodist)
library(vegan)
ord.pcoa.bray <- ordinate(ps.prop, method="PCoA", distance="bray")

ord.pcoa.bray$values

pcoa1 = ord.pcoa.bray$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2 = ord.pcoa.bray$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1 = pcoa1, pcoa2 = pcoa2)
MDS = cbind(MDS, metadata)

ggplot(MDS, aes(x = pcoa1, y = pcoa2, color = Location, shape = Plant)) + geom_point(size =3) + 
  theme_bw()

taxa_names(physeq)
n_seqs <- seq(ntaxa(physeq))
len_n_seqs <- nchar(max(n_seqs))
taxa_names(physeq) <- paste("Seq", formatC(n_seqs, 
                                       width = len_n_seqs, 
                                       flag = "0"), sep = "_")
taxa_names(physeq)

otu.only <- as.data.frame(physeq@otu_table)
taxa.only <- as.data.frame(physeq@tax_table)
otu.t <- t(otu.only)
adonis2(otu.t ~ Location, metadata) # p = 0.039
adonis2(otu.t ~ Plant, metadata) # n.s.

merged <- merge(taxa.only, otu.only,
                by = 'row.names', all = TRUE)
merged <- merged %>% dplyr::rename(OTU = Row.names)

fwrite(merged, "/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/waterlilie_16s.csv")

bacteria <- fread("/Users/briannepalmer/Downloads/NymphaeaNuphar_16SRaw/waterlilie_16s.csv")
bacteria.t <- bacteria %>% pivot_longer(cols = starts_with("N"), names_to = "SampleID")

bacteria.t <- bacteria.t %>% filter(Kingdom == "d__Bacteria")
bacteria.gp <- bacteria.t %>% group_by(SampleID, OTU) %>% summarize(n = sum(value))  
bacteria.gp <- bacteria.gp %>% mutate(rel = n / sum(n))

bacteria.otu <- bacteria.gp %>% pivot_wider(names_from =OTU, values_from = rel)
bacteria.otu[c(3:3971)] <- sapply(bacteria.otu[c(3:3971)],as.numeric)
bacteria.otu[c(3:3971)][is.na(bacteria.otu[c(3:3971)])] <-0
bacteria.otu <- bacteria.otu %>% group_by(SampleID) %>% summarize_all(sum)
otu.only <- bacteria.otu[3:3971]


phylum <-  bacteria.t %>% group_by(SampleID, Phylum) %>% summarize(n = sum(value))%>% mutate(rel = n / sum(n))
phylum <- phylum %>% dplyr::select(SampleID, Phylum, rel)
phylum[phylum == ""] <- NA  
phylum <- phylum %>% pivot_wider(names_from = Phylum, values_from = rel)
phylum[c(2:44)] <- sapply(phylum[c(2:44)],as.numeric)
phylum[is.na(phylum)] <- 0

phylum <- as.data.frame(c(metadata, phylum))
phylum.ggplot <- pivot_longer(phylum, cols = NA.:Zixibacteria,names_to = "phylum", values_to = "count")


ggplot(phylum.ggplot, aes(fill = phylum, y = count, x = Location)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance") +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + facet_grid(. ~ Plant)

phylum.anova <- phylum.ggplot %>% group_by(phylum) %>%
  do(broom::tidy(TukeyHSD(aov(count ~ Location, data = .)))) %>%
  ungroup

phylum.anova.treatment <- phylum.anova %>% filter(adj.p.value < 0.05)

bacteria.all <- as.data.frame(c(metadata, bacteria.otu))

library(vegan)
NupharBottom.rich <- specnumber((bacteria.all %>% filter(Location == "bottom", Plant == "Nuphar"))[-c(1:11)])
NymphaeaBottom.rich <- specnumber((bacteria.all %>% filter(Location == "bottom", Plant == "Nymphaea"))[-c(1:11)])
NupharTop.rich <- specnumber((bacteria.all %>% filter(Location == "top", Plant == "Nuphar"))[-c(1:11)])
NymphaeaTop.rich <- specnumber((bacteria.all %>% filter(Location == "top", Plant == "Nymphaea"))[-c(1:11)])


NupharBottom.div <- vegan::diversity((bacteria.all %>% filter(Location == "bottom", Plant == "Nuphar"))[-c(1:11)])
NymphaeaBottom.div <- vegan::diversity((bacteria.all %>% filter(Location == "bottom", Plant == "Nymphaea"))[-c(1:11)])
NupharTop.div <- vegan::diversity((bacteria.all %>% filter(Location == "top", Plant == "Nuphar"))[-c(1:11)])
NymphaeaTop.div <- vegan::diversity((bacteria.all %>% filter(Location == "top", Plant == "Nymphaea"))[-c(1:11)])


richness <- c(NupharBottom.rich, NymphaeaBottom.rich, NupharTop.rich, NymphaeaTop.rich)
diversity <- c(NupharBottom.div, NymphaeaBottom.div, NupharTop.div, NymphaeaTop.div)

# make new data frame 

rich.div <- data.frame(bacteria.all$Plant, bacteria.all$SampleID, bacteria.all$Location, richness, diversity)

rich.lm <- lm(richness ~ bacteria.all.Plant, data = rich.div)
rich.aov <- aov(rich.lm)
summary(rich.aov) #there is a difference in richness between plants but not location 

div.lm <- lm(diversity ~ bacteria.all.Location, data = rich.div)
div.aov <- aov(div.lm)
summary(div.aov) # divesity did not vary 


richness.plot <- ggplot(data = rich.div, aes(x = bacteria.all.Location, y = richness, fill = bacteria.all.Plant)) + geom_boxplot() + theme_bw()  + theme(text=element_text(size=20)) + labs(y = "# of OTUs", x = "Location") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
richness.plot

diversiy.plot <- ggplot(data = rich.div, aes(x = bacteria.all.Location, y = diversity, fill = bacteria.all.Plant)) + geom_boxplot() + theme(text=element_text(size=20)) + labs(y = "Shannon-Weiner Diversity", x = "Treatment") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

diversiy.plot


