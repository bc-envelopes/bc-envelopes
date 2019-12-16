##########################################
####  Run bioclim models - updated
##########################################
#### | Project name: Bioclimatic envelopes
#### | Script type: Data processing/modeling
#### | What it does: Integrated script that creates species data, calculates hypervolume, and runs models
#### | Date created: August 21, 2019.
#### | Creator: ------
#### | Contact: ------
##########################################


# Script setup ------------------------------------------------------------

# Load packages
pacman::p_load(Rahat, biomod3, raster, tidyverse, hypervolume,
               sf, sp, fasterize, janitor, glue, reshape2)

source(milkunize2("Projects/Bioclim_envelopes/R/hypervolume_overlap_percent.R"))

source(milkunize2("Projects/Bioclim_envelopes/R/niche_overlap.R"))

source(milkunize2("Projects/Bioclim_envelopes/R/modeling_functions.R"))




# Load data ---------------------------------------------------------------

# Load csv that contains list of 300 species to model (100 species of fish, amphibians, and mammals). 

processing_species_files <- "Projects/Bioclim_envelopes/Data/Species_point_data_all" %>%
  milkunize2("archive") %>%
  list.files(full.names = TRUE)

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


# Do stuff ----------------------------------------------------------------

# Get index from slurm array job
i <- as.numeric(commandArgs(trailingOnly = TRUE))


species_list <- "Projects/Bioclim_envelopes/Data/All_species_list.csv" %>%
  milkunize2("archive") %>%
  read_csv()

input_species_range <- st_read(processing_species_files[i])


# Get species name
species_name <- as.character(unique(input_species_range$binomial))
# Get taxonomic group
taxonomic_group <- species_list %>% 
  filter(species == species_name) %>% 
  pull(group)

cat(glue("Creating data for {species_name}. {length(processing_species_files) - i} more species to go."), "\n")


# Niche overlap -----------------------------------------------------------


#### Get predictor set
# clim_var_sets <- c("2var", "4var", "10var", "allvar")
clim_var_sets <- c("2var", "4var", "10var")

bioclim_2var <- paste0("bio_", str_pad(c(1, 12), width = 2, pad = "0", side = "left"))
bioclim_4var <- paste0("bio_", str_pad(c(1, 4, 12, 15), width = 2, pad = "0", side = "left"))
bioclim_10var <- paste0("bio_", str_pad(c(2:4, 8:9, 13:15, 18:19), width = 2, pad = "0", side = "left"))
bioclim_allvar <- paste0("bio_", str_pad(1:19, width = 2, pad = "0", side = "left"))

# Set folder for species


species_data_path <- glue("Projects/Bioclim_envelopes/Output_2019/{taxonomic_group}/{species_name}") %>% 
  milkunize2("archive") 

dir.create(species_data_path, recursive = TRUE)
dir.create(glue("{species_data_path}/Outputs"), recursive = TRUE)


# Count the number of files in the folder
files_number <- species_data_path %>% 
  list.files(recursive = TRUE, pattern = "assessment*.*.csv$") %>% 
  length()


if (files_number == 48)
{
  print("All done")
} else {
  # Niche overlap calc
  overlap_filename <- glue("{species_data_path}/Outputs/{species_name}_overlap.csv")
  
  
  if (!file.exists(overlap_filename))
  {
    
    cat(glue("Calculating niche overlap for {species_name}."), "\n")
    
    # Extract values of rescaled environmental variables for locations of species data 
    processing_species_sp <- as(input_species_range, "Spatial")
    
    bioclim_stack_norm <- subset(bioclim_stack_norm, bioclim_4var)
    
    species_bioclim_raw <- raster::extract(bioclim_stack_norm, processing_species_sp, sp = TRUE)
    
    
    
    species_bioclim <- species_bioclim_raw %>% 
      as.data.frame() %>%  
      # rename(Realm = Realm_5m) %>% 
      dplyr::select(-c(coords.x1, coords.x2))
    
    try(
      {
        niche_overlap_ex <- niche_overlap(species_bioclim, var_set = "4var", mode = "ex")
        niche_overlap_cv <- niche_overlap(species_bioclim, var_set = "4var", mode = "cv")
      }
    )
    
    try(
      {
        overlap_stats_ex <- niche_overlap_ex %>% 
          transmute(
            index,
            value,
            species = species_name,
            var_set = "4var",
            mode = "ex"
          )
        
        overlap_stats_cv <- niche_overlap_cv %>% 
          transmute(
            index,
            value,
            species = species_name,
            var_set = "4var",
            mode = "cv"
          )
        overlap_stats <- rbind(overlap_stats_ex, overlap_stats_cv)
        write_csv(overlap_stats, overlap_filename)
      }
    )
  } else {
    print("Niche file exists")
  }
  # Main loop ---------------------------------------------------------------
  for (predictor_set in c("2var", "4var", "10var", "allvar"))
  # for (predictor_set in c("2var", "4var"))
  {
    for (mod_eval in c("ex", "cv"))
    {
      # for (pseudo_method in "PA1")
      for (pseudo_method in c("PA1", "PA2", "PA3"))
      {
        # for (predictor_set in c("allvar"))
        cat(glue("{predictor_set}, {mod_eval}, {pseudo_method}"), "\n")
        
        scenario_name <- glue("{predictor_set}_{pseudo_method}_{mod_eval}")
        
        model_assessment_name <- glue("{species_data_path}/Outputs/{species_name}_{scenario_name}_assessment.csv")
        model_assessment_ensemble_name <- glue("{species_data_path}/Outputs/{species_name}_{scenario_name}_assessment_ensemble.csv")
        
        if (file.exists(model_assessment_ensemble_name))
        {
          cat("File exists", "\n")
          next()
        }
        
        setwd(species_data_path)
        
        fitted_model <- run_model(data_sf = input_species_range, 
                                  evaluation_mode = mod_eval,
                                  variable_set = predictor_set,
                                  pseudoabs_method = pseudo_method,
                                  climate_stack = bioclim_stack)
        
        ####
        model_evaluation <- evaluate_model(model = fitted_model,
                                           evaluation_mode = mod_eval,
                                           variable_set = predictor_set,
                                           pseudoabs_method = pseudo_method)
        
        
        write_csv(model_evaluation, model_assessment_name)
        
        # Ensemble model
        ensemble_model_out <- BIOMOD_EnsembleModeling(modeling.output = fitted_model,
                                                      chosen.models = "all",
                                                      em.by = "PA_dataset+repet",
                                                      eval.metric = c("CSI", "TSS"), #metric used to scale the ensamble
                                                      eval.metric.quality.threshold = c(-1, -1),
                                                      models.eval.meth = c("CSI", "TSS", "ROC"),
                                                      prob.mean = TRUE,
                                                      prob.mean.weight = TRUE)
        
        ####
        ensemble_model_evaluation <- evaluate_ensemble(model = fitted_model, 
                                                       ensemble_model = ensemble_model_out, 
                                                       evaluation_mode = mod_eval,
                                                       variable_set = predictor_set,
                                                       pseudoabs_method = pseudo_method)
        
        
        write_csv(ensemble_model_evaluation, model_assessment_ensemble_name)  
        
      }
    }
  }
  
  glue("{getwd()}/{str_replace(species_name, \"_\", \".\")}") %>% 
    list.files(recursive = TRUE, full.names = TRUE) %>% 
    file.remove()
  
}


