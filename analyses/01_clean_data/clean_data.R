# Libraries ---------------------------------------------------------------
library(here)

library(tidyr)
library(dplyr)

library(RSnetwork)

# Get data folders ---------------------------------------------------------------
andean_folder <- here("data/ANDEAN_frugivory")

out_folder <- here("outputs/01_clean_data")

all_files_andean <- list.files(andean_folder)

# Read traits tables ----------
plant_traits_file <- file.path(andean_folder,
                               grep("Plant_traits", all_files_andean, value = TRUE))
plant_traits_all <- read.csv(plant_traits_file,
                             sep = "\t")

animal_traits_file <- file.path(andean_folder,
                                grep("Bird_traits", all_files_andean, value = TRUE))
animal_traits_all <- read.csv(animal_traits_file,
                              sep = "\t")

# Format traits table ----------
plant_traits_all <- plant_traits_all |>
  rename("plant_species" = "Species")

animal_traits_all <- animal_traits_all |>
  rename("animal_species" = "Species")

# List network files ----------
network_files <- grep("NW_", all_files_andean,
                      value = TRUE)

for (i in 1:length(network_files)) {
  filename_i <- network_files[i]

  interactions <- read.csv(file.path(andean_folder,
                                     filename_i),
                           sep = "\t")

  # Transform matrix to dataframe ----------
  colnames(interactions)[1] <- "plant_species"
  interactions <- interactions |> pivot_longer(cols = 2:ncol(interactions),
                                               names_to = "animal_species",
                                               values_to = "frequency")

  # Select traits subset ----------
  country <- gsub(filename_i,
                  pattern = "(^NW_)([[:alpha:]]+).*(\\.txt$)",
                  replacement = "\\2")

  animal_traits <- animal_traits_all |>
    filter(Country == country) |>
    filter(animal_species %in% interactions$animal_species)

  plant_traits <- plant_traits_all |>
    filter(Country == country) |>
    filter(plant_species %in% interactions$plant_species)

  # Create codes ----------
  plant_codes <- create_code(interactions$plant_species,
                             name = "plant_species")
  animal_codes <- create_code(interactions$animal_species,
                              name = "animal_species",
                              code_pattern = "B")

  ## Add codes to interactions (remove whole name) -----
  interactions <- interactions |>
    left_join(animal_codes)

  interactions <- interactions |>
    left_join(plant_codes)

  interactions <- interactions |>
    select(plant_species, plant_species_code,
           animal_species, animal_species_code,
           frequency, everything())

  ## Add codes to traits -----
  animal_traits <- animal_traits |>
    left_join(animal_codes) |>
    relocate(animal_species_code, .before = animal_species) |>
    select(-Country)

  plant_traits <- plant_traits |>
    left_join(plant_codes) |>
    relocate(plant_species_code, .before = plant_species) |>
    select(-Country)

  file_id <- gsub(filename_i,
                  pattern = "(^NW_)(.*)(\\.txt$)",
                  replacement = "\\2")

  subfolder <- file.path(out_folder, paste("ANDEAN", file_id, sep = "_"))

  if(!dir.exists(subfolder)) {
    dir.create(subfolder, recursive = TRUE)
  }

  write.csv(interactions, file.path(subfolder, "interactions.csv"),
            row.names = FALSE)
  write.csv(plant_traits, file.path(subfolder, "plant_traits.csv"),
            row.names = FALSE)
  write.csv(animal_traits, file.path(subfolder, "animal_traits.csv"),
            row.names = FALSE)
}
