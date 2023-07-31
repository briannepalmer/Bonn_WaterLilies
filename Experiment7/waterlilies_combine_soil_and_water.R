library(tidyverse)
library(phyloseq)

leaves <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/nymphaeanuphar_16S_leavesonly.csv")
soilwater <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_16S_sw.csv")

# remove the leaves from the soil and water file 

soilwater <- soilwater[,-c(8:19)]
metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/experiment7_metadata.csv")

leaves.t <- leaves %>% pivot_longer(cols = NuBotR1:NyTopR3, names_to = "SampleID")
soilwater.t <- soilwater %>% pivot_longer(cols = BS1:BW3, names_to = "SampleID")

all.data <- rbind(leaves.t, soilwater.t)
all.data <- merge(all.data, metadata, by = "SampleID")

all.data.wide <- all.data %>% select(SampleID:species, value) %>% pivot_wider(names_from = SampleID, values_from = value)

fwrite(all.data,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwater_and_metadata_03april.csv")

fwrite(all.data.wide,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization R/Experiment 7/data/Experiment7/waterlilies_w_soilwater_03april.csv")


########################################################

leaves <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/experiment7_ITS.csv")
soilwater <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Dataframes_to_merge/experiment7_ITS_sw.csv")
metadata <- fread("/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/experiment7_metadata.csv")

leaves.t <- leaves %>% pivot_longer(cols = NuBotR1:NyTopR3, names_to = "SampleID")
soilwater.t <- soilwater %>% pivot_longer(cols = BS1:BW3, names_to = "SampleID")

all.data <- rbind(leaves.t, soilwater.t)
all.data <- merge(all.data, metadata, by = "SampleID")

all.data.wide <- all.data %>% select(SampleID:species, value) %>% pivot_wider(names_from = SampleID, values_from = value)

fwrite(all.data,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/data/waterlilies_w_soilwater_and_metadataITS.csv")

fwrite(all.data.wide,file = "/Users/briannepalmer/Library/CloudStorage/OneDrive-Personal/LeafFossilization/Leaf Fossilization Experiments_3_4_7_8/Experiment 7/data/waterlilies_w_soilwaterITS.csv")


