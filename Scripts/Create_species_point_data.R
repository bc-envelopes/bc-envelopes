##########################################
####  Run bioclim models - updated
##########################################
#### | Project name: Bioclim envelopes
#### | Script type: Data processing/modeling
#### | What it does: Integrated script that creates species data, calculates hypervolume, and runs models
#### | Date created: March 04, 2019.
#### | Creator: -----
#### | Contact: -----
##########################################


# Script setup ------------------------------------------------------------

# Load packages
pacman::p_load(Rahat, biomod3, raster, tidyverse, hypervolume, sf, sp, fasterize, janitor, glue, reshape2)

source(milkunize2("Projects/Bioclim_envelopes/R/hypervolume_overlap_percent.R"))

source(milkunize2("Projects/Bioclim_envelopes/R/niche_overlap.R"))

source(milkunize2("Projects/Bioclim_envelopes/R/modeling_functions.R"))




# Load data ---------------------------------------------------------------

# Load csv that contains list of 300 species to model (100 species of fish, amphibians, and mammals).
# Species list is taken from manuscript Supplementary Info
# processing_species <- "Projects/Bioclim_envelopes/Data/Species_list.csv" %>%
#   milkunize2() %>%
#   read_csv()

species_to_model_filename <- "Projects/Bioclim_envelopes/Data/Species_to_model_all.gpkg" %>% 
  milkunize2("archive")

species_to_model_sf <- st_read(species_to_model_filename)

# "/home/mirza/milkunarc/Projects/Bioclim_envelopes/Data/"

# Load global realms raster. This will be used to limit presence/absence selection (see methods).
# Spatial resolution is 5 arc-minutes.
biomes_raster <- "Projects/amphibians.pressures/Data_raw/Rasters/Biomes_5m.tif" %>%
  milkunize2() %>%
  raster()

# Load global basins raster. This is used to delinate species data subsetting regions
basins_5m <- raster("/vol/milkunarc/vbarbarossa/xGlobio-aqua/pcrglobwb_hydrography/basins_5min.tif") %>%
  setMinMax()

# Load raster predictors at 5 arc degree spatial resolution -- set with 4 layers
bioclim_stack <- "Data_RAW/Current_climate/5m" %>%
  milkunize2("archive") %>%
  list.files(pattern = "tif", full.names = TRUE) %>%
  stack()

# Change layer names
names(bioclim_stack) <- names(bioclim_stack) %>%
  str_replace("bio10_", "bio_")

# Load normalized raster predictors at 5 arc degree spatial resolution
# To be used for niche overlap stuff
bioclim_stack_norm <- "Data_RAW/Current_climate/5m/Normalized" %>%
  milkunize2("archive") %>%
  list.files(pattern = "tif", full.names = TRUE) %>%
  stack()


# Check how many species
species_to_model_sf %>%
  st_set_geometry(NULL) %>%
  distinct(species, .keep_all = TRUE) %>%
  tabyl(group)


# Create species data -----------------------------------------------------

# Get array index from slurm scheduler
i <- as.numeric(commandArgs(trailingOnly = TRUE))


species_to_model_string <- as.character(unique(species_to_model_sf$species))



  species_name <- species_to_model_string[i]
  
  cat(glue("Creating data for {species_name}, index {i}. {length(species_to_model_string) - i} more species to go."), "\n")
  
  
  filename_sp_pts <- milkunize2(glue("Projects/Bioclim_envelopes/Data/Species_point_data_all/{species_name}_pts.gpkg"), "archive")
  
  if (!file.exists(filename_sp_pts))
  {
    # species_to_model_sf
    
    input_species_range <- species_to_model_sf %>%
      filter(species == species_name)
    
    
    taxonomic_group <- unique(as.character(input_species_range$group))
    
    range_sp <- fasterize(input_species_range, biomes_raster)
    
    ####
    #### Create pseudoabsences ####
    # print("Calculating species mask")
    
    if (taxonomic_group != "Fish")
    {
      #select code of realms overlaping with species range
      mask_realms <- biomes_raster %>%
        mask(range_sp) %>%
        unique()
      
      # Create realms raster which keep the values of original realms,
      realms_species <- biomes_raster
      realms_species[!realms_species %in% mask_realms] <- NA
      realms_species[mask(biomes_raster, range_sp)] <- NA
      
      # Sample 100000 pseudoabsences
      pas <- sampleRandom(realms_species, 100000, na.rm = TRUE, sp = TRUE) %>%
        st_as_sf()
      
      #### Create species presence/absence data
      
      # Get presence data from rasterized species range
      species_presences <- biomes_raster %>%
        mask(range_sp) %>%
        rasterToPoints(spatial = TRUE) %>%
        st_as_sf() %>%
        rename(
          area_id = 1
        )
      
    }
    
    if (taxonomic_group == "Fish")
    {
      mask_realms <- basins_5m %>%
        mask(range_sp) %>%
        unique()
      
      # Create realms raster which keep the values of original realms,
      realms_species <- basins_5m
      realms_species[!realms_species %in% mask_realms] <- NA
      realms_species[mask(basins_5m, range_sp)] <- NA
      
      # Sample 100000 pseudoabsences
      pas <- sampleRandom(realms_species, 100000, na.rm = TRUE, sp = TRUE) %>%
        st_as_sf()
      
      #### Create species presence/absence data
      
      # Get presence data from rasterized species range
      species_presences <- basins_5m %>%
        mask(range_sp) %>%
        rasterToPoints(spatial = TRUE) %>%
        st_as_sf() %>%
        rename(
          area_id = 1
        )
      
    }
    
    # Get realm with the highest number of points (might not be the correct approach).
    
    # If there are 2 realms, realm with the largest number of points will be for training data, other one for model testing
    # If there are more than 2 realms, merge 3 realm with the 1st
    
    # To deal with this problem, realm that is second when it comes to point quanity will be "test", and all other (even if n > 2) will be train
    
    points_freq <- species_presences %>%
      tabyl(1) %>%
      arrange(desc(n))
    
    my_row <- nrow(points_freq)
    my_vector <- rep(c("train", "test"), my_row)
    
    my_vector <- my_vector[1:my_row]
    
    points_freq <- points_freq %>%
      transmute(
        area_id = .[,1],
        n,
        Set = my_vector)
    
    points_presence <- species_presences %>%
      inner_join(points_freq, by = "area_id") %>%
      mutate(
        PA = 1
      )
    
    
    processing_species_sf <- pas %>%
      rename(area_id = 1) %>%
      inner_join(points_freq, by = "area_id") %>%
      mutate(
        PA = 0
      ) %>%
      rbind(points_presence) %>%
      mutate(binomial = species_name)
    
    st_write(processing_species_sf, filename_sp_pts)
    
  }
  
  

###



cat(glue("All done for {species_name}."), "\n")
