### water lilies venn diagram with soil and water 

library(tidyverse)
library(data.table)

### Bacteria ####

S16 <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwater_03april.csv")

metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")

S16[,c(9:26)][is.na(S16[,c(9:26)])] <- 0


S16 <- S16 %>% filter(kingdom != "Unassigned") %>% pivot_longer(cols = c(BS1:NyTopR3), names_to = "SampleID")

S16 <- S16 %>% group_by(SampleID, kingdom, phylum, class, order, family, genus) %>% dplyr::summarize(value = sum(value))

genus <- S16 %>% dplyr::select(SampleID, genus, value) %>% filter(genus != "")
genus <- genus %>% pivot_wider(names_from = genus, values_from = value)
genus[c(7:905)] <- sapply(genus[c(7:905)],as.numeric)
genus[c(7:905)][is.na(genus[c(7:905)])] <- 0
genus <- merge(metadata, genus)
genus.gp <- genus %>% group_by(SampleID, Location, Plant, Treatment) %>% summarize_at(vars(g__Bathyarchaeia:g__Zixibacteria), sum)


genus.long <- pivot_longer(genus.gp, cols = g__Bathyarchaeia:g__Zixibacteria, names_to = "genus", values_to = "count")

y <- genus.long %>% dplyr::select(Treatment, genus, count)

Sunken = y %>% filter(Location == "Sunken")  %>% filter(count >0) %>% dplyr::select(genus) #1959
Floating = y %>% filter(Location == "Floating")  %>% filter(count >0) %>% dplyr::select(genus) #1216
Soil = y %>% filter(Location == "Soil")  %>% filter(count >0) %>% dplyr::select(genus) #1508
Water = y %>% filter(Location == "Water")  %>% filter(count >0) %>% dplyr::select(genus) #505


y <- rbind(Water, Soil)

Treatment <- c(rep("Water", 505), rep("Soil", 1508))

y$Treatment <- Treatment

x <- list(Water = y$genus[1:505], Soil = y$genus[506:2014])

library(ggvenn)

venn.soil.sunken.bacteria <- ggvenn(
  x,
  fill_color = c('#882255', '#117733'),
  stroke_size = 0.5, set_name_size = 3, show_percentage = TRUE, text_size = 4, auto_scale = TRUE) + theme(text = element_text(size = 12)) + ggtitle ("A")
venn.soil.sunken.bacteria

venn.soil.water.bacteriua <- ggvenn(
  x,
  fill_color = c('#882255', '#117733'),
  stroke_size = 0.5, set_name_size = 3, show_percentage = TRUE, text_size = 4, auto_scale = FALSE) + theme(text = element_text(size = 12)) + ggtitle ("Bacteria in Pond Mud and Water")
venn.soil.water.bacteriua


y <- rbind(Floating, Water)

Treatment <- c(rep("Floating", 1216), rep("Water", 505))

y$Treatment <- Treatment

x <- list(Floating = y$genus[1:1216], Water = y$genus[1217:1722])

venn.water.floating.bacteria <- ggvenn(
  x,
  fill_color = c('#88CCEE','#332288'),
  stroke_size = 0.5, set_name_size = 3, show_percentage = TRUE, text_size = 4, auto_scale = TRUE) + theme(text = element_text(size = 12)) + ggtitle ("B")
venn.water.floating.bacteria


#### Fungi ####

ITS <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwaterITS.csv")

metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")

ITS[,c(9:26)][is.na(ITS[,c(9:26)])] <- 0

ITS <- ITS %>% filter(kingdom != "Unassigned") %>% pivot_longer(cols = c(BS1:NyTopR3), names_to = "SampleID")

ITS <- ITS %>% group_by(SampleID, kingdom, phylum, class, order, family, genus) %>% dplyr::summarize(value = sum(value))

genus <- ITS %>% dplyr::select(SampleID, genus, value) %>% filter(genus != "")
genus <- genus %>% pivot_wider(names_from = genus, values_from = value)
genus[c(7:139)] <- sapply(genus[c(7:139)],as.numeric)
genus[c(7:139)][is.na(genus[c(7:139)])] <- 0
genus <- merge(metadata, genus)
genus.gp <- genus %>% group_by(SampleID, Location, Plant, Treatment) %>% summarize_at(vars(g__Scrippsiella:g__Rozellomycota_gen_Incertae_sedis), sum)

genus.long <- pivot_longer(genus.gp, cols = g__Scrippsiella:g__Rozellomycota_gen_Incertae_sedis, names_to = "genus", values_to = "count")

y <- genus.long %>% dplyr::select(Location, genus, count)

Sunken = y %>% filter(Location == "Sunken")  %>% filter(count >0) %>% dplyr::select(genus) #96
Floating = y %>% filter(Location == "Floating")  %>% filter(count >0) %>% dplyr::select(genus) #96
Soil = y %>% filter(Location == "Soil")  %>% filter(count >0) %>% dplyr::select(genus) #75
Water = y %>% filter(Location == "Water")  %>% filter(count >0) %>% dplyr::select(genus) #175


y <- rbind(Water, Soil)

Treatment <- c(rep("Water", 175), rep("Soil", 75))

y$Treatment <- Treatment

x <- list(Water = y$genus[1:175], Soil = y$genus[176:251])


venn.soilsunken.fungi <- ggvenn(
  x,
  fill_color = c('#882255', '#117733'),
  stroke_size = 0.5, set_name_size = 3, show_percentage = TRUE, text_size = 4, auto_scale = TRUE) + theme(text = element_text(size = 12)) + ggtitle ("C")
venn.soilsunken.fungi


venn.soil.water.fungi <- ggvenn(
  x,
  fill_color = c('#88CCEE','#332288'),
  stroke_size = 0.5, set_name_size = 3, show_percentage = TRUE, text_size = 4, auto_scale = FALSE) + theme(text = element_text(size = 12)) + ggtitle ("Fungi in Pond Mud and Water")
venn.soil.water.fungi


y <- rbind(Floating, Water)

Treatment <- c(rep("Floating", 96), rep("Water", 175))

y$Treatment <- Treatment

x <- list(Floating = y$genus[1:96], Water = y$genus[97:271])

venn.water.floating.fungi <- ggvenn(
  x,
  fill_color = c('#88CCEE','#332288'),
  stroke_size = 0.5, set_name_size = 2.5, show_percentage = TRUE, text_size = 4, auto_scale = TRUE) + theme(text = element_text(size = 12)) + ggtitle ("D")
venn.water.floating.fungi

library(patchwork)
png("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/figures/Experiment7/NymphaeaNuphar_Figures_April2023/soilwater_venn.png", width = 39, height = 24, units = 'cm', res = 300)
(venn.soil.sunken.bacteria + venn.water.floating.bacteria)/(venn.soilsunken.fungi + venn.water.floating.fungi)
dev.off()

png("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/figures/Experiment7/NymphaeaNuphar_Figures_April2023/soilwater_venn_nobiofilm.png", width = 39, height = 24, units = 'cm', res = 300)
(venn.soil.water.bacteriua + venn.soil.water.fungi)
dev.off()

#find ASVs that are the shared in all treatments 

y <- genus.long %>% dplyr::select(Location, genus, count)

asv.shared <- y %>% dplyr::group_by(genus, Location) %>%
  mutate(row = row_number())  %>% dplyr::summarize(row = sum(row)) %>% pivot_wider(names_from = genus, values_from = row)
# remove any column with NA 
asv.shared <- asv.shared %>% select_if(~ !any(is.na(.))) 


#What are the ASVs that are present in all substrates?

shared.asvs <- y %>% filter(genus %in% colnames(asv.shared))
shared.asvs$present <- +(shared.asvs$count >=1 &!is.na(shared.asvs$count))

shared.asvs.phylum <- shared.asvs %>% filter(present > 0) %>% group_by(genus) %>% tally()
shared.asvs.phylum <- shared.asvs.phylum %>% group_by(genus) %>% tally() 



