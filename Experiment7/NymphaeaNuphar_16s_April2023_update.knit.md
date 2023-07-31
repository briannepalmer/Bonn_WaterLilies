
<!-- rnb-text-begin -->

---
title: "Nymphaea Nuphar 16S "
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook.
When you execute code within the notebook, the results appear beneath the code.

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*.

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file).

The preview shows you a rendered HTML copy of the contents of the editor.
Consequently, unlike *Knit*, *Preview* does not run any R code chunks.
Instead, the output of the chunk when it was last run in the editor is displayed.


<!-- rnb-text-end -->



<!-- rnb-text-begin -->


This document contains all the code for Experiment 7 -- Nymphaea Nuphar Buried and Floating.
The primary file used here contains the 16S sequences for the water lily leaves and for the soil and water.


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubGlicmFyeShwaHlsb3NlcSlcbmxpYnJhcnkocGx5cilcbmxpYnJhcnkoZGF0YS50YWJsZSlcbmxpYnJhcnkodGlkeXZlcnNlKVxubGlicmFyeSh2ZWdhbilcbmxpYnJhcnkoZWNvZGlzdClcbmxpYnJhcnkoUkNvbG9yQnJld2VyKVxubGlicmFyeShhcGUpXG5saWJyYXJ5KGFtcHZpczIpXG5saWJyYXJ5KGdncHVicilcbmxpYnJhcnkoZ2d2ZW5uKVxubGlicmFyeShnZ3BhdHRlcm4pXG5saWJyYXJ5KHJzdGF0aXgpXG5saWJyYXJ5KGdncmVwZWwpXG5saWJyYXJ5KGZ1bnJhcilcbmxpYnJhcnkoREVTZXEyKVxuYGBgIn0= -->

```r
library(phyloseq)
library(plyr)
library(data.table)
library(tidyverse)
library(vegan)
library(ecodist)
library(RColorBrewer)
library(ape)
library(ampvis2)
library(ggpubr)
library(ggvenn)
library(ggpattern)
library(rstatix)
library(ggrepel)
library(funrar)
library(DESeq2)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuUzE2IDwtIGZyZWFkKFwiL1VzZXJzL2JyaWFubmVwYWxtZXIvTGlicmFyeS9DbG91ZFN0b3JhZ2UvT25lRHJpdmUtUGVyc29uYWwvTGVhZkZvc3NpbGl6YXRpb24vTGVhZiBGb3NzaWxpemF0aW9uIFIvRXhwZXJpbWVudCA3L2RhdGEvRXhwZXJpbWVudDcvd2F0ZXJsaWxpZXNfd19zb2lsd2F0ZXJfMDNhcHJpbC5jc3ZcIilcbmBgYCJ9 -->

```r
S16 <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwater_03april.csv")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWV0YWRhdGEgPC0gZnJlYWQoXCIvVXNlcnMvYnJpYW5uZXBhbG1lci9MaWJyYXJ5L0Nsb3VkU3RvcmFnZS9PbmVEcml2ZS1QZXJzb25hbC9MZWFmRm9zc2lsaXphdGlvbi9MZWFmIEZvc3NpbGl6YXRpb24gUi9FeHBlcmltZW50IDcvZGF0YS9FeHBlcmltZW50Ny9leHBlcmltZW50N19tZXRhZGF0YS5jc3ZcIilcbmBgYCJ9 -->

```r
metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxudHJhbnNmb3JtZWQuMTYgPC0gdChTMTZbLGMoOToyNildKVxudHJhbnNmb3JtZWQuMTZbaXMubmEodHJhbnNmb3JtZWQuMTYpXSA8LTBcbnRyYW5zZm9ybWVkLjE2W2MoMToxODU3MjQpXSA8LSBzYXBwbHkodHJhbnNmb3JtZWQuMTZbYygxOjE4NTcyNCldLGFzLm51bWVyaWMpXG5cblMxNi5yZWwgPC0gbWFrZV9yZWxhdGl2ZSh0cmFuc2Zvcm1lZC4xNilcblxub3R1LmJyYXkuMTYgPC0gdmVnZGlzdChTMTYucmVsLCBtZXRob2QgPSBcImJyYXlcIilcbnBjb2FWUy4xNiA8LSBlY29kaXN0OjpwY28ob3R1LmJyYXkuMTYsIG5lZ3ZhbHMgPSBcInplcm9cIiwgZHJvdW5kID0gMClcblxucGNvYTEuMTYgPSBwY29hVlMuMTYkdmVjdG9yc1ssMV0gIyBhZGRzIHRoZSBwb2ludHMgb2YgdGhlIE5NRFMgMSBkaW1tZW5zaW9uXG5wY29hMi4xNiA9IHBjb2FWUy4xNiR2ZWN0b3JzWywyXSAjIGFkZHMgdGhlIHBvaW50cyBvZiB0aGUgTk1EUyAyIGRpbW1lbnNpb25cblxuTURTID0gZGF0YS5mcmFtZShwY29hMS4xNiA9IHBjb2ExLjE2LCBwY29hMi4xNiA9IHBjb2EyLjE2KVxuTURTID0gY2JpbmQoTURTLCBtZXRhZGF0YSlcbmBgYCJ9 -->

```r
transformed.16 <- t(S16[,c(9:26)])
transformed.16[is.na(transformed.16)] <-0
transformed.16[c(1:185724)] <- sapply(transformed.16[c(1:185724)],as.numeric)

S16.rel <- make_relative(transformed.16)

otu.bray.16 <- vegdist(S16.rel, method = "bray")
pcoaVS.16 <- ecodist::pco(otu.bray.16, negvals = "zero", dround = 0)

pcoa1.16 = pcoaVS.16$vectors[,1] # adds the points of the NMDS 1 dimmension
pcoa2.16 = pcoaVS.16$vectors[,2] # adds the points of the NMDS 2 dimmension

MDS = data.frame(pcoa1.16 = pcoa1.16, pcoa2.16 = pcoa2.16)
MDS = cbind(MDS, metadata)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


RELATIVE ABUNDANCE OF PHYLA


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Are there phyla that differ between treatments?


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Make table based on ASVs to then use as OTU and Tax Tables in phyloseq


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Plot richness and diversity


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Make a new dataframe with the ASVs from the pivot dataframe (leaves only for the venn diagram)


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBsZWF2ZXMub25seSA8LSBTMTYgJT4lIGZpbHRlcihraW5nZG9tICE9IFwiVW5hc3NpZ25lZFwiKSAlPiUgcGl2b3RfbG9uZ2VyKGNvbHMgPSBjKE51Qm90UjE6TnlUb3BSMyksIG5hbWVzX3RvID0gXCJTYW1wbGVJRFwiKVxuIyBcbiMgbGVhdmVzLm9ubHkgPC0gbGVhdmVzLm9ubHkgJT4lIGdyb3VwX2J5KFNhbXBsZUlELCBraW5nZG9tLCBwaHlsdW0sIGNsYXNzLCBvcmRlciwgZmFtaWx5LCBnZW51cywgc3BlY2llcywgQVNWKSAlPiUgZHBseXI6OnN1bW1hcml6ZSh2YWx1ZSA9IHN1bSh2YWx1ZSkpXG4jIFxuIyBBU1YgPC0gbGVhdmVzLm9ubHkgJT4lIGRwbHlyOjpzZWxlY3QoU2FtcGxlSUQsIEFTViwgdmFsdWUpICNcbiMgQVNWIDwtIEFTViAlPiUgcGl2b3Rfd2lkZXIobmFtZXNfZnJvbSA9IEFTViwgdmFsdWVzX2Zyb20gPSB2YWx1ZSlcbiMgQVNWW2MoOToxMDMyMyldIDwtIHNhcHBseShBU1ZbYyg5OjEwMzIzKV0sYXMubnVtZXJpYylcbiMgQVNWW2MoOToxMDMyMyldW2lzLm5hKEFTVltjKDk6MTAzMjMpXSldIDwtIDBcbiMgQVNWIDwtIG1lcmdlKG1ldGFkYXRhLCBBU1YpXG4jIEFTVi5ncCA8LSBBU1YgJT4lIGdyb3VwX2J5KFNhbXBsZUlELCBMb2NhdGlvbiwgUGxhbnQsIFRyZWF0bWVudCkgJT4lIHN1bW1hcml6ZV9hdCh2YXJzKEFTVjIzNjY6QVNWOTE1KSwgc3VtKVxuYGBgIn0= -->

```r
# leaves.only <- S16 %>% filter(kingdom != "Unassigned") %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
# 
# leaves.only <- leaves.only %>% group_by(SampleID, kingdom, phylum, class, order, family, genus, species, ASV) %>% dplyr::summarize(value = sum(value))
# 
# ASV <- leaves.only %>% dplyr::select(SampleID, ASV, value) #
# ASV <- ASV %>% pivot_wider(names_from = ASV, values_from = value)
# ASV[c(9:10323)] <- sapply(ASV[c(9:10323)],as.numeric)
# ASV[c(9:10323)][is.na(ASV[c(9:10323)])] <- 0
# ASV <- merge(metadata, ASV)
# ASV.gp <- ASV %>% group_by(SampleID, Location, Plant, Treatment) %>% summarize_at(vars(ASV2366:ASV915), sum)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Make heatmaps with ampvis


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2VudXMgPC0gbWVyZ2UobWV0YWRhdGEsIHBpdm90LjE2UywgYnkgPSBcIlNhbXBsZUlEXCIpXG5cbmdlbnVzLmFub3ZhIDwtIGdlbnVzICU+JSBncm91cF9ieShnZW51cykgJT4lXG4gIGRvKGJyb29tOjp0aWR5KFR1a2V5SFNEKGFvdih2YWx1ZSB+IExvY2F0aW9uLCBkYXRhID0gLikpKSkgJT4lXG4gIHVuZ3JvdXBcblxuZ2VudXMuYW5vdmEuc2lnIDwtZ2VudXMuYW5vdmEgJT4lIGZpbHRlcihhZGoucC52YWx1ZSA8IDAuMDUpIFxuXG5nZW51cy5hbm92YS5zaWdcbmBgYCJ9 -->

```r
genus <- merge(metadata, pivot.16S, by = "SampleID")

genus.anova <- genus %>% group_by(genus) %>%
  do(broom::tidy(TukeyHSD(aov(value ~ Location, data = .)))) %>%
  ungroup

genus.anova.sig <-genus.anova %>% filter(adj.p.value < 0.05) 

genus.anova.sig
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZ2dwbG90KGdlbnVzICU+JSBmaWx0ZXIoZ2VudXMgJWluJSBnZW51cy5hbm92YS5zaWckZ2VudXMpLCBhZXMoeCA9IGdlbnVzLCB5ID0gdmFsdWUsIGZpbGwgPSBMb2NhdGlvbikpICsgZ2VvbV9iYXIoc3RhdCA9IFwiaWRlbnRpdHlcIiwgcG9zaXRpb24gPSBcImRvZGdlXCIpIFxuXG5gYGAifQ== -->

```r
ggplot(genus %>% filter(genus %in% genus.anova.sig$genus), aes(x = genus, y = value, fill = Location)) + geom_bar(stat = "identity", position = "dodge") 

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGRmKFwiZmlndXJlcy9FeHBlcmltZW50Ny9OeW1waGFlYU51cGhhcl9GaWd1cmVzX0FwcmlsMjAyMy9oZWF0bWFwXzE2Uy5wZGZcIiwgd2lkdGggPSAyMCwgaGVpZ2h0ID0gMjApXG5nZW51cy5oZWF0bWFwLnBcbmRldi5vZmYoKVxuYGBgIn0= -->

```r
pdf("figures/Experiment7/NymphaeaNuphar_Figures_April2023/heatmap_16S.pdf", width = 20, height = 20)
genus.heatmap.p
dev.off()
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Determine differences in richness and diversity


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Look for ASVs that are overrepresented in floating or buried leaves using Deseq


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBGaXJzdCwgd2Ugd2lsbCBjb252ZXJ0IHRoZSBzZXF1ZW5jZSBjb3VudCBkYXRhIGZyb20gb3VyIHBoeWxvc2VxIG9iamVjdCBpbnRvIHRoZSBwcm9wZXIgZm9ybWF0IGZvciBERVNlcTJcbm1iLmRkcyA8LSBwaHlsb3NlcV90b19kZXNlcTIobWIsIH4gTG9jYXRpb24pXG5cbiMgTmV4dCB3ZSBuZWVkIHRvIGVzdGltYXRlIHRoZSBzaXplIGZhY3RvcnMgZm9yIG91ciBzZXF1ZW5jZXNcbmdtX21lYW4gPSBmdW5jdGlvbih4LCBuYS5ybT1UUlVFKXtcbiAgICBleHAoc3VtKGxvZyh4W3ggPiAwXSksIG5hLnJtPW5hLnJtKSAvIGxlbmd0aCh4KSlcbn1cbmdlb01lYW5zID0gYXBwbHkoY291bnRzKG1iLmRkcyksIDEsIGdtX21lYW4pXG5tYi5kZHMgPSBlc3RpbWF0ZVNpemVGYWN0b3JzKG1iLmRkcywgZ2VvTWVhbnMgPSBnZW9NZWFucylcblxuIyBGaW5hbGx5LCB3ZSBjYW4gcnVuIHRoZSBjb3JlIERFU2VxMiBhbGdvcml0aG1cbm1iLmRkcyA9IERFU2VxKG1iLmRkcywgZml0VHlwZT1cImxvY2FsXCIpXG5gYGAifQ== -->

```r
# First, we will convert the sequence count data from our phyloseq object into the proper format for DESeq2
mb.dds <- phyloseq_to_deseq2(mb, ~ Location)

# Next we need to estimate the size factors for our sequences
gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(mb.dds), 1, gm_mean)
mb.dds = estimateSizeFactors(mb.dds, geoMeans = geoMeans)

# Finally, we can run the core DESeq2 algorithm
mb.dds = DESeq(mb.dds, fitType="local")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuIyBDcmVhdGUgYSBuZXcgZGF0YWZyYW1lIGZyb20gdGhlIHJlc3VsdHMgb2Ygb3VyIERFU2VxMiBydW5cbm1iLnJlcyA9IHJlc3VsdHMobWIuZGRzKVxuXG4jIFJlb3JkZXIgdGhlIHNlcXVlbmNlcyBieSB0aGVpciBhZGp1c3RlZCBwLXZhbHVlc1xubWIucmVzID0gbWIucmVzW29yZGVyKG1iLnJlcyRwYWRqLCBuYS5sYXN0PU5BKSwgXVxuXG4jIFNldCB5b3VyIGFscGhhIGZvciB0ZXN0aW5nIHNpZ25pZmljYW5jZSwgYW5kIGZpbHRlciBvdXQgbm9uLXNpZ25pZmljYW50IHJlc3VsdHNcbmFscGhhID0gMC4wNVxubWIuc2lndGFiID0gbWIucmVzWyhtYi5yZXMkcGFkaiA8IGFscGhhKSwgXVxuXG5tYi5yZXMgPSBjYmluZChhcyhtYi5yZXMsIFwiZGF0YS5mcmFtZVwiKSwgYXModGF4X3RhYmxlKG1iKVtyb3duYW1lcyhtYi5yZXMpLCBdLCBcIm1hdHJpeFwiKSlcblxuIyBBZGQgdGF4b25vbXkgaW5mb3JtYXRpb24gdG8gZWFjaCBzZXF1ZW5jZVxubWIuc2lndGFiID0gY2JpbmQoYXMobWIuc2lndGFiLCBcImRhdGEuZnJhbWVcIiksIGFzKHRheF90YWJsZShtYilbcm93bmFtZXMobWIuc2lndGFiKSwgXSwgXCJtYXRyaXhcIikpXG5oZWFkKG1iLnNpZ3RhYilcblxubWIucmVzLmZ1bmd1cyA8LSBmcmVhZChcImV4cG9ydCAvZGVzZXFmdW5ndXMuY3N2XCIpXG5tYi5zaWd0YWIuZnVuZ3VzIDwtIG1iLnJlcy5mdW5ndXMgJT4lIGZpbHRlcihwYWRqID4gMC4wNSlcblxubWIuc2lndGFiLmFsbCA8LSByYmluZChtYi5zaWd0YWIsIG1iLnNpZ3RhYi5mdW5ndXMpXG5cbmBgYCJ9 -->

```r
# Create a new dataframe from the results of our DESeq2 run
mb.res = results(mb.dds)

# Reorder the sequences by their adjusted p-values
mb.res = mb.res[order(mb.res$padj, na.last=NA), ]

# Set your alpha for testing significance, and filter out non-significant results
alpha = 0.05
mb.sigtab = mb.res[(mb.res$padj < alpha), ]

mb.res = cbind(as(mb.res, "data.frame"), as(tax_table(mb)[rownames(mb.res), ], "matrix"))

# Add taxonomy information to each sequence
mb.sigtab = cbind(as(mb.sigtab, "data.frame"), as(tax_table(mb)[rownames(mb.sigtab), ], "matrix"))
head(mb.sigtab)

mb.res.fungus <- fread("export /deseqfungus.csv")
mb.sigtab.fungus <- mb.res.fungus %>% filter(padj > 0.05)

mb.sigtab.all <- rbind(mb.sigtab, mb.sigtab.fungus)

```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->





<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxucGRmKFwiZmlndXJlcy9FeHBlcmltZW50Ny9OeW1waGFlYU51cGhhcl9GaWd1cmVzX0FwcmlsMjAyMy92b2xjYW5vcGxvdF9hbGwucGRmXCIpXG52b2xjYW5vLnBsb3RcbmRldi5vZmYoKVxuYGBgIn0= -->

```r
pdf("figures/Experiment7/NymphaeaNuphar_Figures_April2023/volcanoplot_all.pdf")
volcano.plot
dev.off()
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Bugbase for phenotye prediction


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Functional Predictions with PiCrust


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxua2VnZ19icml0ZV9tYXBwaW5nIDwtIHJlYWQuY3N2KFwiL1VzZXJzL2JyaWFubmVwYWxtZXIvTGlicmFyeS9DbG91ZFN0b3JhZ2UvT25lRHJpdmUtUGVyc29uYWwvTGVhZkZvc3NpbGl6YXRpb24vRXhwZXJpbWVudDdfTnltcGhhZWFOdXBoYXIvcGljcnVzdF9rZWdnX2VkaXQuY3N2XCIpXG5gYGAifQ== -->

```r
kegg_brite_mapping <- read.csv("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Experiment7_NymphaeaNuphar/picrust_kegg_edit.csv")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuZGF0YSA8LSByZWFkLnRhYmxlKFwiL1VzZXJzL2JyaWFubmVwYWxtZXIvTGlicmFyeS9DbG91ZFN0b3JhZ2UvT25lRHJpdmUtUGVyc29uYWwvRE5BX0RhdGFiYXNlcy9waWNydXN0Mi0yLjUuMS9waWNydXN0Ml9vdXRfcGlwZWxpbmUvS09fbWV0YWdlbm9tZV9vdXQvcHJlZF9tZXRhZ2Vub21lX3Vuc3RyYXQudHN2Lmd6XCIsIGhlYWRlcj1UUlVFLCBzZXA9XCJcXHRcIiwgcm93Lm5hbWVzPTEpXG5kYXRhJHBhdGh3YXkgPC0gcm93bmFtZXMoZGF0YSlcbnRlc3QgPC0gbWVyZ2UoZGF0YSwga2VnZ19icml0ZV9tYXBwaW5nLCBieSA9IFwicGF0aHdheVwiKVxuI3Rlc3QgPC0gdGVzdCAlPiUgc2VwYXJhdGUobWV0YWRhdGFfS0VHR19QYXRod2F5cywgYyhcImxldmVsMVwiLCBcImxldmVsMlwiLCBcImxldmVsM1wiLCBcImxldmVsNFwiLCBcImxldmVsNVwiLCBcImxldmVsNlwiLCBcImxldmVsN1wiLCBcImxldmVsOFwiLCBcImxldmVsOVwiLCBcImxldmVsMTBcIiwgXCJsZXZlbDExXCIsIFwibGV2ZWwxMlwiLCBcImxldmVsMTNcIiwgXCJsZXZlbDE0XCIpLCBzZXAgPSBcIjtcIilcbnBpY3J1c3QgPC0gdGVzdCAlPiUgcGl2b3RfbG9uZ2VyKGNvbHMgPSBjKE51Qm90UjE6TnlUb3BSMyksIG5hbWVzX3RvID0gXCJTYW1wbGVJRFwiKVxuYGBgIn0= -->

```r
data <- read.table("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/DNA_Databases/picrust2-2.5.1/picrust2_out_pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", header=TRUE, sep="\t", row.names=1)
data$pathway <- rownames(data)
test <- merge(data, kegg_brite_mapping, by = "pathway")
#test <- test %>% separate(metadata_KEGG_Pathways, c("level1", "level2", "level3", "level4", "level5", "level6", "level7", "level8", "level9", "level10", "level11", "level12", "level13", "level14"), sep = ";")
picrust <- test %>% pivot_longer(cols = c(NuBotR1:NyTopR3), names_to = "SampleID")
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxuU2FtcGxlSUQgPC0gYyhcIk51Qm90UjFcIiwgXCJOdUJvdFIyXCIsIFwiTnVCb3RSM1wiLCBcIk51VG9wUjFcIiwgXCJOdVRvcFIyXCIsIFwiTnVUb3BSM1wiLCBcIk55Qm90UjFcIiwgXCJOeUJvdFIyXCIsIFwiTnlCb3RSM1wiLCBcIk55VG9wUjFcIiwgXCJOeVRvcFIyXCIsIFwiTnlUb3BSM1wiKVxuTG9jYXRpb24gPC1jKFwiYnVyaWVkXCIsIFwiYnVyaWVkXCIsIFwiYnVyaWVkXCIsIFwiZmxvYXRpbmdcIiwgXCJmbG9hdGluZ1wiLCBcImZsb2F0aW5nXCIsIFwiYnVyaWVkXCIsIFwiYnVyaWVkXCIsIFwiYnVyaWVkXCIsIFwiZmxvYXRpbmdcIiwgXCJmbG9hdGluZ1wiLCBcImZsb2F0aW5nXCIpXG5wbGFudCA8LSBjKFwiTnVwaGFyXCIsIFwiTnVwaGFyXCIsIFwiTnVwaGFyXCIsIFwiTnVwaGFyXCIsIFwiTnVwaGFyXCIsIFwiTnVwaGFyXCIsIFwiTnltcGhhZWFcIiwgXCJOeW1waGFlYVwiLCBcIk55bXBoYWVhXCIsIFwiTnltcGhhZWFcIiwgXCJOeW1waGFlYVwiLCBcIk55bXBoYWVhXCIpXG50cmVhdG1lbnQgPC0gYyhcIk51Qm90XCIsIFwiTnVCb3RcIiwgXCJOdUJvdFwiLCBcIk51VG9wXCIsIFwiTnVUb3BcIiwgXCJOdVRvcFwiLCBcIk55Qm90XCIsIFwiTnlCb3RcIiwgXCJOeUJvdFwiLCBcIk55VG9wXCIsIFwiTnlUb3BcIiwgXCJOeVRvcFwiKVxuXG5tZXRhZGF0YSA8LSBkYXRhLmZyYW1lKFNhbXBsZUlELCBMb2NhdGlvbiwgcGxhbnQsIHRyZWF0bWVudClcbmBgYCJ9 -->

```r
SampleID <- c("NuBotR1", "NuBotR2", "NuBotR3", "NuTopR1", "NuTopR2", "NuTopR3", "NyBotR1", "NyBotR2", "NyBotR3", "NyTopR1", "NyTopR2", "NyTopR3")
Location <-c("buried", "buried", "buried", "floating", "floating", "floating", "buried", "buried", "buried", "floating", "floating", "floating")
plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")
treatment <- c("NuBot", "NuBot", "NuBot", "NuTop", "NuTop", "NuTop", "NyBot", "NyBot", "NyBot", "NyTop", "NyTop", "NyTop")

metadata <- data.frame(SampleID, Location, plant, treatment)
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->


<!-- rnb-source-begin eyJkYXRhIjoiYGBgclxubWluZXJhbC5hYnNvcnB0aW9uIDwtIGdncGxvdChkYXRhID0gcGljcnVzdCAlPiUgZmlsdGVyKGxldmVsMyA9PSBcIk1pbmVyYWwgYWJzb3JwdGlvblwiKSwgYWVzKHggPSBMb2NhdGlvbiwgeSA9IHZhbHVlLCBmaWxsID0gTG9jYXRpb24pKSArIGdlb21fYmFyKHN0YXQgPSBcImlkZW50aXR5XCIsIHBvc2l0aW9uID0gXCJkb2RnZVwiKSAgKyB0aGVtZV9idygpKyBzY2FsZV9maWxsX21hbnVhbCh2YWx1ZXMgPSBjKCcjODgyMjU1JywgJyM4OENDRUUnKSkgKyB0aGVtZSh0ZXh0PWVsZW1lbnRfdGV4dChzaXplPTIwKSwgbGVnZW5kLnBvc2l0aW9uID0gXCJub25lXCIpICsgbGFicyh5ID0gXCIjIEdlbmVzXCIsIHggPSBcIkxvY2F0aW9uXCIpICsgZmFjZXRfd3JhcCgufmxldmVsMywgc2NhbGVzID0gXCJmcmVlXCIpXG5cbm1pbmVyYWwuYWJzb3JwdGlvblxuYGBgIn0= -->

```r
mineral.absorption <- ggplot(data = picrust %>% filter(level3 == "Mineral absorption"), aes(x = Location, y = value, fill = Location)) + geom_bar(stat = "identity", position = "dodge")  + theme_bw()+ scale_fill_manual(values = c('#882255', '#88CCEE')) + theme(text=element_text(size=20), legend.position = "none") + labs(y = "# Genes", x = "Location") + facet_wrap(.~level3, scales = "free")

mineral.absorption
```

<!-- rnb-source-end -->

<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->


Combine graphs into prokaryotic figure 


<!-- rnb-text-end -->


<!-- rnb-chunk-begin -->



<!-- rnb-chunk-end -->


<!-- rnb-text-begin -->



<!-- rnb-text-end -->

