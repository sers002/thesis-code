library(callr)
library(magrittr)
library(data.table)

## Vector of directory names (rep_sims/n_100, rep_sims/n_250, etc)
dirs <- list.dirs(here::here("Output/genotype_simulation/rep_sims"), recursive = FALSE)

## For each directory, get filenames within that
simulation_files <- lapply(
  dirs,
  list.files,
  pattern = "n[0-9]+_[0-9]+.controls.haps.gz",
  full.names = TRUE
) %>%
  lapply(function(x) x[1:min(5, length(x))]) %>%
  unlist()

simulation_files <- na.omit(simulation_files)

comp_time.dt <- data.table()

## Loop over simulated datasets, using a new R session for each one (each new R
## session will source the file "one_loop.R" with the given filename as an argument)
for (filename in simulation_files) {
  begin_time <- proc.time()

  rscript(
    here::here("Scripts/full_simulations/empirical/analyse_one_hapgen_dataset.R"),
    cmdargs = filename,
    stderr = "test_output.txt",
    spinner = TRUE,
    echo = TRUE
  )

  end_time <- proc.time()

  comp_time.dt <- rbind(
    comp_time.dt,
    data.table(
      file = filename,
      time = end_time["elapsed"] - begin_time["elapsed"]
    )
  )
}

# fwrite(comp_time.dt, file = "Output/matrixEQTL/rep_sims/time_taken_matrixeqtl.csv.gz")
