library(data.table)

load(file = here::here("Output/sojo_output/all_betas.RData"))

source(here::here("Scripts/helper_functions.R"))

## Get a vector of directory names that contain the results, e.g.:
##   matrixEQTL/rep_sims/n100, matrixEQTL/rep_sims/n500, etcs
n_dirs <- list.dirs("Output/matrixEQTL/rep_sims", full.names = TRUE, recursive = FALSE)

## For each directory, list filenames containing "cis"
cis_files <- lapply(
  n_dirs,
  list.files,
  pattern = "cis",
  full.names = TRUE
)

comp_time.dt <- data.table()

## For each file, calculate the power based on the results
pwr.cis <- lapply(
  cis_files,
  function(file_list) {
    lapply(
      file_list,
      function(filename) {
        message(filename)

        begin_time <- proc.time()

        power_results <- calculate_power_egenes(
          fread(
            filename,
            select = c("SNP", "gene", "beta", "t-stat", "p-value"),
            col.names = c("snps", "gene", "beta", "statistic", "pvalue")
          ),
          n_filtered = 382146,
          n_all      = 11102563
        )


        end_time <- proc.time()

        comp_time.dt <<- rbind(
          comp_time.dt,
          data.table(
            file = filename,
            time = end_time["elapsed"] - begin_time["elapsed"]
          )
        )

        power_results
      }
    )
  }
)

# p-value thresholds
# setNames(lapply(pwr.cis, function(l) mean(unlist(lapply(l, function(l2) unique(l2$p_threshold) )))), gsub(".*n_([0-9]+)$", "\\1", n_dirs))

## Each entry of the pwr.cis list is the directory name, so rename it to just
## the sample size (e.g. 100, 250, etc)
names(pwr.cis) <- gsub(".*n_([0-9]+)$", "\\1", n_dirs)

## Row bind datasets together
power.dt  <- rbindlist(lapply(pwr.cis, rbindlist), idcol = "n")

## Grouping by sample size, effect size and MAF, calculate average power
power.agg <- power.dt[, list(
    power_filtered_fdr = mean(power_filtered_fdr),
    power_all_fdr      = mean(power_all_fdr),
    power_filtered_bon = mean(power_filtered_bon),
    power_all_bon      = mean(power_all_bon)
  ),
  by = c("n", "egene")
]

power.agg$n <- as.numeric(as.character(power.agg$n))

power.agg <- power.agg[order(n, egene)]

fwrite(power.agg, file = "Output/power_egenes/matrixEQTL.csv")
