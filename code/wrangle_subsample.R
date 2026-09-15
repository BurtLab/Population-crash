library(tidyverse)

##### Wrangle subsampled Tajima's D outputs into a single CSV #####
#
# Combines the per-seed "{model}_{seed}_subsample_tajimasD.csv" files
# produced by subsample_tajimasD.py (columns: Time, N, Tajimas_D) into one
# combined file with columns: Time, N, Seed, Model, Value.

process_subsample_tajimasD <- function(input_dir, output_dir) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  files <- list.files(input_dir, pattern = "_subsample_tajimasD\\.csv$", full.names = TRUE)

  combined <- files %>%
    map_dfr(function(file) {
      seed_num <- str_extract(basename(file), "\\d+(?=_subsample_tajimasD\\.csv)") %>% as.integer()
      model_type <- str_extract(basename(file), "^[^_]+")

      read_csv(file, show_col_types = FALSE) %>%
        mutate(Seed = seed_num, Model = model_type) %>%
        rename(Value = Tajimas_D) %>%
        select(Time, N, Seed, Model, Value)
    })

  write_csv(combined, file.path(output_dir, "TajimasD_subsample.csv"))
}
