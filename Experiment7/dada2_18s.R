# https://benjjneb.github.io/dada2/tutorial.html

library("devtools")

### Load Libraries ####
library(data.table)
library(dada2)
library(beepr)

### Set Path ####
path <- "/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/demux"
list.files(path)

### Forward and reverse gz files have the format samplename_1 and samplename_2 ####
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


### Extract sample names, assuming filenames have format: SAMPLENAME_XXX.gz####
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list.files(path)

###plot the read quality, these look good! Q score is above 30####
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

### Filter and Trim ####
# Place filtered files in filtered/ subdirectory

filtFs <- file.path("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/demux", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("/Users/briannepalmer/Downloads/NymphaeaNuphar_18S/demux", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names
list.files(path)


# We’ll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and maxEE=2. The maxEE parameter sets the maximum number of “expected errors” allowed in a read, which is a better filter than simply averaging quality scores.

# this step took about 2 hours to run for 30 samples 
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=200,
                     maxN=0, maxEE=2, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE, trimLeft =c(18,19));beep(sound = 3) # On Windows set multithread=FALSE
head(out)

### dereplication -- this reduces the downstrean computationsl time.  ####

derepF1 <- derepFastq(filtFs, verbose=TRUE); beep(sound = 3) 
derepR1 <- derepFastq(filtRs, verbose=TRUE); beep(sound = 3) 

### Learn the error rates ####

errF <- learnErrors(filtFs, multithread = TRUE); beep(sound = 3) 
errR <- learnErrors(filtRs, multithread = TRUE); beep(sound = 3)
plotErrors(errF, nominalQ=TRUE)

### Sample inference ####
dadaFs <- dada(filtFs, err=errF, multithread = TRUE); beep(sound = 3)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE); beep(sound = 3)
dadaFs[[1]]

### merged paired reads ####

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE); beep(sound = 3)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

### Construct sequence table #### 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

### Remove chimeras ####

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

### Assign taxonomy ####

## HERE is where you can add the CyanoSeq database to better classify the cyano sequences ##
# put the fasta file with the silva classification in the path folder

to_split <- seq(1, ncol(seqtab), by = 50000)
to_split2 <- c(to_split[2:length(to_split)]-1, ncol(seqtab))

taxtab = NULL

for(i in 1:length(to_split)){
  seqtab2 <- seqtab[, to_split[i]:to_split2[i]]
  taxtab2 <- assignTaxonomy(seqtab2, refFasta = ref_fasta, multithread = TRUE)
  if(!is.null(ref_fasta_spp)){taxtab2 <- addSpecies(taxtab2, refFasta = ref_fasta_spp, verbose = TRUE)}
  if(!is.null(taxtab)){taxtab <- rbind(taxtab, taxtab2)}
  if(is.null(taxtab)){taxtab <- taxtab2}
}
# save files in case phylogeny does not run
saveRDS(taxtab, paste(output_path, 'temp', 'taxtab.rds', sep = '/'))

taxa <- assignTaxonomy(seqtab.nochim, "/Users/briannepalmer/Downloads/SILVA_132_SSURef_tax_silva.fasta.gz"); beep(sound = 3) #13:47 

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

### Handoff to phyloseq ####

library(phyloseq)
library(Biostrings)
library(ggplot2)
library(tidyverse)
theme_set(theme_bw())

# We can construct a simple sample data.frame from the information encoded in the filenames. Usually this step would instead involve reading the sample data in from a file.

# change this to match your samples 


meta <- read.csv("/Users/briannepalmer/Downloads/metadata.csv", sep = ";")
meta <- meta %>% filter(SampleID != c("B1", "B2", "B3"))


samples.out <- rownames(seqtab.nochim)
metadata <- meta

samdf <- metadata
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# save.image("~/experiment 3 13Okober2022.RData")
# load("~/dada2 experiment 3 13Okober2022.RData"); beep(sound = 3) 


## make OTU table for other analyses ####
taxa_names(ps) <- paste0("SV", seq(ntaxa(ps)))
# Extract abundance matrix from the phyloseq object
OTU1 = as(otu_table(ps), "matrix")
# transpose if necessary
#if(taxa_are_rows(ps)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)


# make tax table based on relative abundance ####
TGroup <- tax_glom(ps, taxrank = "Genus")


# write OTU and TAX tables ####
fwrite(OTUdf, "Desktop/otu_cy.txt", sep = ",")
write.csv(ps@tax_table, "Desktop/tax_cy.csv")


plot_richness(ps, x="Treatment", measures=c("Shannon", "Simpson"))
# Transform data to proportions as appropriate for Bray-Curtis distances
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Treatment", title="Bray NMDS")

library(ecodist)
library(vegan)
#make bray-curtis matrix
otu.bray <- vegdist(ps.prop@otu_table, method = "bray")# dissimilarity matrix using bray-curtis distance indices
#Calculate PCoA values 
pcoaVS <- pco(otu.bray, negvals = "zero", dround = 0) # if negvals = 0 sets all negative eigenvalues to zero

pcoa1 = pcoaVS$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2 = pcoaVS$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1 = pcoa1, pcoa2 = pcoa2)
MDS = cbind(MDS, metadata)

smallpal2 <- palette(c("#0073C2FF", "#EFC000FF" , "#f08080", "#CD534CFF"))
smallpal2 <- colorRampPalette(smallpal2)(5)
ggplot(MDS, aes(x = pcoa1, y = pcoa2, color = Treatment.Name)) + geom_point(size =3) + 
  theme_bw() + stat_ellipse(aes(color = Treatment.Name)) + scale_color_manual(values = smallpal2)
adonis(ps.prop@otu_table~ Treatment.Name, metadata) # p = 0.001

library(pairwiseAdonis)
set.seed(12345)
pairwise.adonis(ps.prop@otu_table, metadata$Treatment.Name, p.adjust.m = "BH") 

top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)




