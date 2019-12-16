##########################################
####  Select species for modeling
##########################################
#### | Project name: Bioclimatic envelopes
#### | Script type: Data processing
#### | What it does: Description
#### | Date created: August 21, 2019.
#### | Creator: ------
#### | Contact: ---------
##########################################


# Script setup ------------------------------------------------------------


pacman::p_load(Rahat, tidyverse)

# Load data ---------------------------------------------------------------


# Load files in the folder

species_selection_filename <- "Projects/Bioclim_envelopes/Data/Species_selection_data.csv" %>% 
  milkunize2("archive")

if (!file.exists(species_selection_filename))
{
  all_files_raw <- "Projects/Bioclim_envelopes/Data/Species_selection_biomes" %>% 
    milkunize2("archive") %>% 
    list.files(full.names = TRUE)
  
  all_files_loaded <- all_files_raw %>% 
    map(data.table::fread) %>% 
    reduce(rbind)
  
  write_csv(all_files_loaded, species_selection_filename)
  
} else {
  all_files_loaded <- data.table::fread(species_selection_filename)
}

#### Load species list file (contains all species of fish, amphibians, and mammals)
species_list <- read_csv(milkunize2("Data_RAW/Species_ranges/Fish_amphibians_mammals/Species_list.csv", "archive"))  
# Load species ranges -----------------------------------------------------


ranges_amphibians <- "Data_RAW/Species_ranges/Amphibians_IUCN/AMPHIBIANS.shp" %>% 
  milkunize2("archive") %>% 
  st_read() %>% 
  transmute(
    species = str_replace(binomial, " ", "_"),
    group = "Amphibians"
  )

ranges_mammals <- "Data_RAW/Species_ranges/Mammals_IUCN/TERRESTRIAL_MAMMALS.shp" %>% 
  milkunize2("archive") %>% 
  st_read() %>% 
  transmute(
    species = str_replace(binomial, " ", "_"),
    group = "Mammals"
  )


ranges_fish_pt1 <- "Data_RAW/Species_ranges/Fish_IUCN/FW_FISH_PART1.shp" %>% 
  milkunize2("archive") %>% 
  st_read()

ranges_fish_pt2 <- "Data_RAW/Species_ranges/Fish_IUCN/FW_FISH_PART2.shp" %>% 
  milkunize2("archive") %>% 
  st_read()

ranges_fish <- rbind(ranges_fish_pt1, ranges_fish_pt2) %>% 
  transmute(
    species = str_replace(binomial, " ", "_"),
    group = "Fish"
  )


species_all <- rbind(ranges_fish,
                     ranges_amphibians,
                     ranges_mammals)
#################
# Filter species that appear in more than 2 realms, and that can have more than x presences/absences

####
# Subset species with criteria
species_selected <- all_files_loaded %>% 
  filter(n_realms >= 2, n_realms_abs >= 2, train_1 >= 100, test_1 >= 100, train_0 >= 10000, test_0 >= 10000) %>% 
  inner_join(species_list, by = "species") %>% 
  distinct(species, .keep_all = TRUE)

# This is how many species satisfy species selection criteria, per taxonomic group
species_selected %>% 
  add_count(group) %>% 
  select(group, n) %>% 
  distinct(group, .keep_all = TRUE)

#### Exclude species that have all models done

# Load species that have all 72 models done
# note: this could also contain species with models that haven't converged for example. Check.

species_finished_list <- "Projects/Bioclim_envelopes/Data/Species_list_finished.csv" %>% 
  milkunize2("archive") %>% 
  read_csv()

species_finished_list %>% 
  pull(Species)

# Exclude finished species
species_selected_excluded <- species_selected %>% 
  filter(species %notin% pull(species_finished_list, Species)) %>% 
  pull(species)



species_selected_excluded %>% 
  select(species, group)


#### Save geopackage file with ~1900 species to be used for further modeling
# Run those models, and out of those, select 100 species per group with 72 finished models.

species_to_model_filename <- "Projects/Bioclim_envelopes/Data/Species_to_model_all.gpkg" %>% 
  milkunize2("archive")

species_to_model <- species_all %>% 
  filter(species %in% species_selected_excluded) 

st_write(species_to_model, species_to_model_filename)
