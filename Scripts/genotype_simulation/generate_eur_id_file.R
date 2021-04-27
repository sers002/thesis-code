library(data.table)

pop.1000g <- fread(here::here("Data/1000g_metadata/20130606_g1k.ped"))
# table(pop.1000g$Population)

## Population values:
## https://www.internationalgenome.org/category/population/

eur.pop <- c("CEU", "TSI", "GBR", "IBS")
eas.pop <- c("CHB", "JPT", "CHS", "CDX", "KHV")

pop.1000g <- pop.1000g[Population %in% eur.pop]
cat(pop.1000g$`Individual ID`, file = "./Output/1000g_eur_ids.txt", sep = "\n")
