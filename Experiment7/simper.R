library(tidyverse)
library(data.table)
library(vegan)
library(ggvegan)

S16 <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwater_03april.csv")

metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")

sample <- metadata$SampleID
metadata <- as.matrix(metadata)
rownames(metadata) <- sample
metadata <- as.data.frame(metadata)

x <- seq(1:10318)
otu.names <- as.vector(paste0("ASV", x))
S16$ASV <- otu.names

otu.table <- S16[,c(9:26)]
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- S16$ASV
otu.table <- as.data.frame(otu.table) # make matrix again
otu.table[is.na(otu.table)] <-0
otu.table$sums <- rowSums(otu.table) # this is an optional step, I usually run it as a test to ensure that there is a value for every column, there should be no rowSums = 0 
otu.table <- otu.table %>% filter(sums > 0)
otu.table <- otu.table[-19]

as <- anosim(
  t(otu.table),
  metadata$Treatment,
  permutations = 999,
  distance = "bray",
  strata = NULL,
  parallel = getOption("mc.cores")
)

summary(as)
as16S <- autoplot(as, notch = FALSE, title = "16S (Bacteria)") + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


ITS <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwaterITS.csv")

metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")

sample <- metadata$SampleID
metadata <- as.matrix(metadata)
rownames(metadata) <- sample
metadata <- as.data.frame(metadata)

x <- seq(1:1288)
otu.names <- as.vector(paste0("ASV", x))
ITS$ASV <- otu.names

otu.table <- ITS[,c(9:26)]
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- ITS$ASV
otu.table <- as.data.frame(otu.table) # make matrix again
otu.table[is.na(otu.table)] <-0
otu.table$sums <- rowSums(otu.table) # this is an optional step, I usually run it as a test to ensure that there is a value for every column, there should be no rowSums = 0 
otu.table <- otu.table %>% filter(sums > 0)
otu.table <- otu.table[-19]

as <- anosim(
  t(otu.table),
  metadata$Treatment,
  permutations = 999,
  distance = "bray",
  strata = NULL,
  parallel = getOption("mc.cores")
)

summary(as)

asITS <- autoplot(as, notch = FALSE, title = "ITS (Fungi)") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

S18 <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_18S.csv")
S18 <- S18 %>% filter(kingdom %in% c("Eukaryota", "Eukaryota:mito", "Eukaryota:plas"))

metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")
metadata <- metadata %>% filter(Type == "Plant")

sample <- metadata$SampleID
metadata <- as.matrix(metadata)
rownames(metadata) <- sample
metadata <- as.data.frame(metadata)

x <- seq(1:54)
otu.names <- as.vector(paste0("ASV", x))
S18$ASV <- otu.names

otu.table <- S18[,c(2:13)]
otu.table <- as.matrix(otu.table)
rownames(otu.table) <- S18$ASV
otu.table <- as.data.frame(otu.table) # make matrix again
otu.table[is.na(otu.table)] <-0
otu.table$sums <- rowSums(otu.table) # this is an optional step, I usually run it as a test to ensure that there is a value for every column, there should be no rowSums = 0 
otu.table <- otu.table %>% filter(sums > 0)
otu.table <- otu.table[-13]

as <- anosim(
  t(otu.table),
  metadata$Treatment,
  permutations = 999,
  distance = "bray",
  strata = NULL,
  parallel = getOption("mc.cores")
)

summary(as)
as18S <- autoplot(as, notch = FALSE, title = "18s (Eukaryotes)") + theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
plot(as)

library(patchwork)
png("figures/Experiment7/June_Revision_Figures /anosim.png", res = 300, units = "cm", height = 12, width = 48)
as16S + asITS + as18S
dev.off()
