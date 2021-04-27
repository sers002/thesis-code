library(data.table)

source(here::here("Scripts/eqtls/dataprep_and_exploratory/1._process_eqtl_data.R"))

## check orginal vcf
ref_data = fread(
  here::here("Output/genotype_simulation/1000g_chr19_eur_subset.vcf.gz"), 
  skip = 100, 
  select = 1:5
)

length(unique(ref_data$ID))

ref_data = ref_data[nchar(REF) == 1]
ref_data = ref_data[nchar(ALT) == 1]

ref_data = ref_data[!duplicated(ID)]


file_match = merge(ss, ref_data, by.x = "snp", by.y = "ID")

length(unique(ref_data$ID))
length(unique(ss$snp))
length(unique(file_match$snp))

file_match[, orig1 := paste0(a1, a2)]
file_match[, new1 := paste0(REF, ALT)]
file_match[, new2 := paste0(ALT, REF)]

nrow(file_match[(orig1 != new1) & (orig1 != new2)])
file_match = file_match[!((orig1 != new1) & (orig1 != new2))]

length(unique(file_match$sid))
length(unique(file_match$pid))

file_match = file_match[!duplicated(file_match[, list(pid = pid, b = abs(b))])]

nrow(file_match[(orig1 == new1)])
nrow(file_match[(orig1 == new2)])

######################################################################
# Check simulated genotype

ref_data <- fread(
  here::here("Output/genotype_simulation/1000g_chr19_eur_simulated.controls.gen"),
  select = 1:5
)


length(unique(ref_data$V2))

ref_data = ref_data[nchar(V4) == 1]
ref_data = ref_data[nchar(V5) == 1]

ref_data = ref_data[!duplicated(V2)]


file_match = merge(ss, ref_data, by.x = "snp", by.y = "V2")

length(unique(ref_data$V2))
length(unique(ss$snp))
length(unique(file_match$snp))

file_match[, orig1 := paste0(a1, a2)]
file_match[, new1 := paste0(V4, V5)]
file_match[, new2 := paste0(V5, V4)]

nrow(file_match[(orig1 != new1) & (orig1 != new2)])
file_match = file_match[!((orig1 != new1) & (orig1 != new2))]

length(unique(file_match$sid))
length(unique(file_match$pid))

file_match = file_match[!duplicated(file_match[, list(pid = pid, b = abs(b))])]

nrow(file_match[(orig1 == new1)])
nrow(file_match[(orig1 == new2)])
