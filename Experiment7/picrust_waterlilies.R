BiocManager::install("ALDEx2")
browseVignettes("ALDEx2")

library("data.table")   # Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("tidyverse")
library("R.utils")

picrust2 = "/Users/briannepalmer/Downloads/picrust2-2.5.1/picrust2_out_stratified"
list.files(picrust2, recursive = TRUE)

p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

mapfile = "/Users/briannepalmer/Downloads/picrust2-2.5.1/picrust2/default_files/description_mapfiles/"
list.files(mapfile, recursive = TRUE)

mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")

mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")


p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)


metadata<-fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/NymphaeaNuphar_18S/sample-metadata.tsv")
metadata$Location <-c("bottom", "bottom", "bottom", "top", "top", "top", "bottom", "bottom", "bottom", "top", "top", "top")
metadata$plant <- c("Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nuphar", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea", "Nymphaea")
sample_data <- sample_data(metadata)

# EC
set.seed(12345)
system.time({
  aldex2_EC1 = aldex(p2EC, sample_data$Location, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO, sample_data$Location, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})

# Pathway
set.seed(12345)
system.time({
  aldex2_PW1 = aldex(p2PW, sample_data$Location, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


par(mfrow = c(2,2))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")

# since the effect size is bigger than 1, there is a difference in the functions for some of the samples 

png("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/picrust_waterlilies/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())


png("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/picrust_waterlilies/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())


df_EC1 = aldex2_EC1 %>% tibble::rownames_to_column(var = "EC") %>% 
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO1 = aldex2_KO1 %>% tibble::rownames_to_column(var = "KO") %>% 
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)

write.table(df_EC1, file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/picrust_waterlilies/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_KO1, file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/picrust_waterlilies/ALDEx2_picrust2_KO_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.table(df_PW1, file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/picrust_waterlilies/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

BiocManager::install("microbiomeMarker")
library(microbiomeMarker)
BiocManager::install("microbiome")
library(microbiome)
sam_tab <- "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/NymphaeaNuphar_18S/sample-metadata.csv"

library(ape)

ps <- import_picrust2(
  p2_PW,
  sam_tab,
  trait = "PATHWAY")

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))
plot(random_tree)

ps = merge_phyloseq(ps,random_tree)
ps@sam_data <- sample_data(metadata)
rownames(ps@sam_data) <- ps@sam_data$SampleID
colnames(ps@otu_table) <- ps@sam_data$SampleID

# https://rstudio-pubs-static.s3.amazonaws.com/566863_687400bd7e8742568e73bf167fc42d3d.html

wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps))
GP1 = prune_taxa(wh0, ps)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(GP1, "NMDS", "bray")
plot_ordination(GP1, GP.ord , type="samples", color = "Location") + theme_bw()

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)

pw <- plot_ordination(GP1, ordu, color="Location", shape="plant") + theme_bw() + labs(x = "PCOA1", y = "PCOA2")  + geom_point(size =3) + ggtitle("PiCRUST2 Pathways")

library(vegan)
metadata <- as(sample_data(GP1), "data.frame")
adonis2(phyloseq::distance(GP1, method="bray") ~ Location,
       data = metadata) # 0.352

adonis2(phyloseq::distance(GP1, method="bray") ~ plant,
        data = metadata) # 0.37

pw <- t(p2PW)
pw <- cbind(metadata, pw)
pw.t <- pw %>% pivot_longer(cols = `1CMET2-PWY`:`VALSYN-PWY`, names_to = "pathway", values_to = "sum")

pw.anova <- pw.t %>% group_by(pathway) %>%
  do(broom::tidy(TukeyHSD(aov(sum ~ Location, data = .)))) %>%
  ungroup

pw.anova.sig <- pw.anova %>% filter(adj.p.value < 0.05) # nothing 

ggplot(pw.t %>% filter(pathway %in% pw.anova.sig$pathway), aes(fill = pathway, y = sum, x = Location)) + geom_bar(position = "fill", stat = "identity")  + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + labs(y = "Relative Abundance")  +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = .5,
    hjust = 1
  )) + theme_bw()

ggplot(pw.t %>% filter(pathway %in% pw.anova.sig$pathway), aes(y = sum, x = Location, fill = Location)) + geom_boxplot() + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + theme_bw() + facet_wrap(.~pathway, scales = "free") + theme(legend.position = "none")

ggplot(pw.t %>% filter(pathway %in% c("DENITRIFICATION-PWY", "METH-ACETATE-PWY", "METHANOGENESIS-PWY", "PWY-6830", "PWY-7084")), aes(y = sum, x = Location, fill = Location)) + geom_boxplot() + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + theme_bw() + facet_wrap(.~pathway, scales = "free") + theme(legend.position = "none") 


ec <- t(p2EC)
ec <- cbind(metadata, ec)
ec.t <- ec %>% pivot_longer(cols = `EC:1.1.1.1`:`EC:6.6.1.2`, names_to = "function", values_to = "sum")
ec.t <- merge(ec.t, mapEC, by = "function")

ec.anova <- ec.t %>% group_by(description) %>%
  do(broom::tidy(TukeyHSD(aov(sum ~ Location, data = .)))) %>%
  ungroup

ec.anova.sig <- ec.anova %>% filter(adj.p.value < 0.05)




ggplot(pw.t %>% filter(pathway %in% pw.anova.sig$pathway), aes(y = sum, x = Location, fill = Location)) + geom_boxplot() + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + theme_bw() + facet_wrap(.~pathway, scales = "free") + theme(legend.position = "none")

ggplot(pw.t %>% filter(pathway %in% c("DENITRIFICATION-PWY", "METH-ACETATE-PWY", "METHANOGENESIS-PWY", "PWY-6830", "PWY-7084")), aes(y = sum, x = Location, fill = Location)) + geom_boxplot() + theme(axis.text.x = element_text(size = 10)) + theme(legend.position = "right") + theme_bw() + facet_wrap(.~pathway, scales = "free") + theme(legend.position = "none") 



###############################################

ps <- import_picrust2(
  p2_EC,
  sam_tab,
  trait = "EC")

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps = merge_phyloseq(ps,random_tree)
ps@sam_data <- sample_data(metadata)
rownames(ps@sam_data) <- ps@sam_data$SampleID
colnames(ps@otu_table) <- ps@sam_data$SampleID

wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps))
GP1 = prune_taxa(wh0, ps)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(GP1, "NMDS", "bray")
plot_ordination(GP1, GP.ord , type="samples", color = "Location") + theme_bw()

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)
ec <- plot_ordination(GP1, ordu, color="Location", shape="plant") + theme_bw() + labs(x = "PCOA1", y = "PCOA2")  + geom_point(size =3) + ggtitle("PiCRUST2 EC")

metadata <- as(sample_data(GP1), "data.frame")
adonis2(phyloseq::distance(GP1, method="bray") ~ Location,
        data = metadata) # 0.4

adonis2(phyloseq::distance(GP1, method="bray") ~ plant,
        data = metadata) # 0.357

ps <- import_picrust2(
  p2_KO,
  sam_tab,
  trait = "KO")

random_tree = rtree(ntaxa(ps), rooted=TRUE, tip.label=taxa_names(ps))

ps = merge_phyloseq(ps,random_tree)
ps@sam_data <- sample_data(metadata)
rownames(ps@sam_data) <- ps@sam_data$SampleID
colnames(ps@otu_table) <- ps@sam_data$SampleID

wh0 = genefilter_sample(ps, filterfun_sample(function(x) x > 5), A=0.5*nsamples(ps))
GP1 = prune_taxa(wh0, ps)
GP1 = transform_sample_counts(GP1, function(x) 1E6 * x/sum(x))

GP.ord <- ordinate(GP1, "NMDS", "bray")
plot_ordination(GP1, GP.ord , type="samples", color = "Location") + theme_bw()

ordu = ordinate(GP1, "PCoA", "unifrac", weighted=TRUE)
KO <- plot_ordination(GP1, ordu, color="Location", shape="plant") + theme_bw() + labs(x = "PCOA1", y = "PCOA2")  + geom_point(size =3) + ggtitle("PiCRUST2 KO")

metadata <- as(sample_data(GP1), "data.frame")
adonis2(phyloseq::distance(GP1, method="bray") ~ Location,
        data = metadata) # 0.386

adonis2(phyloseq::distance(GP1, method="bray") ~ plant,
        data = metadata) # 0.372

library(patchwork)
KO + ec + pw

