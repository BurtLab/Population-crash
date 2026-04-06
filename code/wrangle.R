library(tidyverse)
library(reticulate)  # for Python integration
np <- import("numpy")


# Define statistic name mappings
stat_specs <- tribble(
  ~orig_name,                     ~output_name,
  "Number of mutations",          "Number_of_mutations",
  "Density of segregating sites", "Density_of_segregating_sites",
  "Nucleotide diversity (pi)",    "Nucleotide_diversity_pi",
  "Tajima's D",                   "Tajimas_D"
)


##### Main function to process files in a directory #####

process_directory_stats <- function(input_dir, output_dir) {

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # List all CSV files
  all_csvs <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
  
  ld_files <- all_csvs[str_detect(basename(all_csvs), "_ld_summary")]
  main_csvs <- all_csvs[!str_detect(basename(all_csvs), "_ld_summary")]
  
  # 1. Process Main CSV files
  if(length(main_csvs) > 0) {
    main_csvs %>%
      map_dfr(function(file) {
        seed_num <- str_extract(basename(file), "\\d+(?=\\.csv)") %>% as.integer()
        # Model detection: capture string before the first underscore
        model_type <- str_extract(basename(file), "^[^_]+")
        
        read_csv(file, show_col_types = FALSE) %>%
          rename(Time = 1) %>%
          mutate(Seed = seed_num, Model = model_type)
      }) %>%
      pivot_longer(-c(Time, Seed, Model), names_to = "Statistic", values_to = "Value") %>%
      filter(Statistic %in% stat_specs$orig_name) %>%
      left_join(stat_specs, by = c("Statistic" = "orig_name")) %>%
      select(-Statistic) %>%
      # Output separate CSVs per statistic (Time, Seed, Model, Value)
      group_by(output_name) %>%
      group_split() %>%
      walk(~{
        stat_name <- first(.x$output_name)
        write_csv(select(.x, -output_name), file.path(output_dir, paste0(stat_name, ".csv")))
      })
  }
  
  # 2. Process LD Summary files
  if(length(ld_files) > 0) {
    ld_data <- ld_files %>%
      map_dfr(function(file) {
        # Filename expected format: model_seed_ld_summary.csv
        seed_num <- str_extract(basename(file), "\\d+(?=_ld_summary)") %>% as.integer()
        # Model detection
        model_type <- str_extract(basename(file), "^[^_]+")
        
        read_csv(file, show_col_types = FALSE) %>%
          mutate(Seed = seed_num, Model = model_type)
      })
    
    # Save combined LD summary
    write_csv(ld_data, file.path(output_dir, "LD.csv"))
  }


  # 3. Process Ne files
  # Save the Ne values into a single CSV
  ne_files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)

  if(length(ne_files) > 0) {
    ne_data <- ne_files %>%
      map_dfr(function(file) {
        # Extract seed number and model type from filename
        seed_num <- str_extract(basename(file), "\\d+(?=\\.txt)") %>% as.integer()
        model_type <- str_extract(basename(file), "^[^_]+")
        
        # Read the single number from the file
        ne_value <- scan(file, quiet = TRUE, what = numeric())
        
        # Create data frame
        tibble(
          Seed = seed_num,
          Model = model_type,
          Ne = ne_value
        )
      })
    
    # Save Ne data as CSV
    write_csv(ne_data, file.path(output_dir, "Ne.csv"))
  }
}
  