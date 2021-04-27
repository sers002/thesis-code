library(data.table)

pop.1000g <- fread(here::here("Data/1000g_metadata/20130606_g1k.ped"))
# table(pop.1000g$Population)

## Population values:
## https://www.internationalgenome.org/category/population/

afr.pop <- c("YRI")

pop.1000g <- pop.1000g[Population %in% afr.pop]
cat(pop.1000g$`Individual ID`, file = "./Output/1000g_yri_ids.txt", sep = "\n")
