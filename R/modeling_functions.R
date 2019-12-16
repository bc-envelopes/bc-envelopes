#' Title
#'
#' @param data simple features species
#' @param evaluation_mode 
#' @param variable_set 
#' @param pseudoabs_method 
#'
#' @return
#' @export
#'
#' @examples
run_model <- function(data_sf, evaluation_mode, variable_set, pseudoabs_method,
                      climate_normalized, climate_stack, which_model = "all")
{
  
  assertthat::assert_that(
    !missing(data_sf),
    evaluation_mode %in% c("ex", "cv"),
    variable_set %in% c("2var", "4var", "10var", "allvar"),
    pseudoabs_method %in% c("PA1", "PA2", "PA3"),
    which_model %in% c("all", "maxent_bckg")
  )
  # 
  # cat(glue("Calculating overlap for {species_name}, {variable_set}iables, for {evaluation_mode}"), "\n")
  
  species_name <- unique(data_sf$binomial)
  
  
  
  
  #### Create dataset for model fitting
  # Subset presences and absences, sample x number of absences and merge back together
  # Subset data for training
  species_presences_train <- data_sf %>% 
    filter(PA == 1) %>% 
    filter(Set == "train")
  
  species_absences_train <- data_sf %>% 
    filter(PA == 0) %>% 
    filter(Set == "train")
  # Subset x number of pseudo-absences
  # Define number of pseudoabsences based on set
  if (pseudoabs_method == "PA1")
  {
    pa_number <- 1000
  }
  if (pseudoabs_method == "PA2")
  {
    pa_number <- 10000
  }
  if (pseudoabs_method == "PA3")
  {
    pa_number <- nrow(species_presences_train)
  }
  
  set.seed(666)
  species_absences_sampled <- sample_n(species_absences_train, pa_number, replace = FALSE)
  # species_absences_sampled <- sample_n(species_absences_train, pa_number, replace = FALSE)
  
  
  if (evaluation_mode == "cv")
  {
    
    species_data_training <- rbind(species_presences_train, 
                                   species_absences_sampled)
    
    set.seed(666)
    
    my_sample <- sample(1:nrow(species_data_training), round(0.8 * nrow(species_data_training)), replace = FALSE)
    
    training_data <- species_data_training[my_sample, ]
    testing_data  <- species_data_training[-my_sample, ]
    
  }
  # If we assess the model performance using independent data, then training data is the full dataset with sampled absences, 
  # and testing data is the predetermined test dataset
  if (evaluation_mode == "ex")
  {
    species_presences_test <- data_sf %>% 
      filter(PA == 1) %>% 
      filter(Set == "test")
    
    species_absences_test <- data_sf %>% 
      filter(PA == 0) %>% 
      filter(Set == "test")
    
    set.seed(666)
    species_absences_test_sampled <- sample_n(species_absences_test, pa_number, replace = FALSE)
    
    
    training_data <- rbind(species_presences_train, 
                           species_absences_sampled)
    
    testing_data  <- rbind(species_presences_test, 
                           species_absences_test_sampled)
  }
  
  # Extract values of rescaled environmental variables for locations of species data 
  # processing_species_sp <- as(training_data, "Spatial")
  # species_bioclim_raw <- raster::extract(climate_normalized, processing_species_sp, sp = TRUE)
  # 
  # species_bioclim <- species_bioclim_raw %>% 
  #   as.data.frame() %>%  
  #   # rename(Realm = Realm_5m) %>% 
  #   dplyr::select(-c(coords.x1, coords.x2))
  # 
  
  # Calculate niche overlap
  
  # overlap_filename <- glue("{species_data_path}{species_name}_{scenario_name}_overlap.csv")
  
  # if (!file.exists(overlap_filename))
  # {
  #   try(
  #     overlap_2var <- niche_overlap(species_bioclim, var_set = variable_set, mode = evaluation_mode)  
  #   )
  #   try(overlap_stats <- overlap_2var %>% 
  #         transmute(
  #           index,
  #           value,
  #           species = species_name,
  #           PA_set,
  #           var_set = variable_set,
  #           mode = evaluation_mode
  #         )
  #   )
  #   
  #   try(write_csv(overlap_stats, overlap_filename)
  #   )
  # }
  
  
  # Select climate layers to subset according to the predictors set     
  layers_to_subset <- paste0("bioclim_", variable_set)
  
  bioclim_stack_subset <- subset(climate_stack, get(layers_to_subset))
  print(nlayers(bioclim_stack_subset))
  
  
  
  # BIOMOD modeling section -------------------------------------------------
  # MaxEnt
  maxent_params <- list(path_to_maxent.jar = "/vol/milkunB/mcengic/MAXENT/maxent.jar",
                        memory_allocated = 2048)
  #### Still to add more modeling tecniques; check from meeting notes what have we decided upon
  biomod_options <- BIOMOD_ModelingOptions(MAXENT.Phillips = maxent_params)
  
  # if (!file.exists(model_assessment_name))
  # {
  # Prepare data for modeling
  model_data <- BIOMOD_FormatingData(resp.var = as.numeric(training_data$PA),
                                     resp.xy = st_coordinates(training_data),
                                     eval.resp.var = as.numeric(testing_data$PA),
                                     eval.resp.xy = st_coordinates(testing_data),
                                     expl.var = bioclim_stack_subset,
                                     resp.name = as.character(species_name),
                                     na.rm = TRUE,
                                     PA.strategy = "user.defined")
  
  print(model_data)
  # Fit model
  if (which_model == "all")
  {
    # CSI method was modified so it calculates Soerensen index
    fitted_models <- BIOMOD_Modeling(model_data,
                                     models = c("GLM", "GAM", "GBM", "MARS", "CTA", "RF", "MAXENT.Phillips"),
                                     # models = c("GLM", "MAXENT.Phillips"),
                                     models.options = biomod_options,
                                     NbRunEval = 1,
                                     DataSplit = 100,
                                     Prevalence = 0.5,
                                     models.eval.meth = c("CSI", "TSS", "ROC"),
                                     SaveObj = TRUE)  
  }
  
  # 
  if (which_model == "maxent_bckg")
  {
    fitted_models <- BIOMOD_Modeling(model_data,
                                     models = "MAXENT.Phillips",
                                     models.options = biomod_options,
                                     NbRunEval = 1,
                                     DataSplit = 100,
                                     Prevalence = 0.5,
                                     models.eval.meth = c("CSI", "TSS", "ROC"),
                                     SaveObj = TRUE)
  }
  
  return(fitted_models)
  
  # }  
  
}


evaluate_model <- function(model, pseudoabs_method,
                           variable_set,
                           evaluation_mode)
{
  assertthat::assert_that(
    !missing(model),
    evaluation_mode %in% c("ex", "cv"),
    variable_set %in% c("2var", "4var", "10var", "allvar"),
    pseudoabs_method %in% c("PA1", "PA2", "PA3")
  )
  
  model_assessment <- model %>%
    biomod3::get_evaluations() %>%
    as.data.frame() %>%
    rownames_to_column("Metric") %>%
    dplyr::select(-contains("Cutoff"), -contains("Sensitivity"), -contains("Specificity")) %>% 
    gather("Variable", "value", 2:ncol(.)) %>%
    transmute(
      Metric,
      Algorithm = str_remove(Variable, ".Full.AllData") %>%
        str_remove("Testing.data.") %>% 
        str_remove("Evaluating.data."),
      value,
      evaluation_set = ifelse(str_detect(Variable, "Testing"), "Fit", "Eval"),
      Species = str_replace(model@sp.name, "[.]", "_"),
      PA_set = pseudoabs_method,
      var_set = variable_set,
      mode = evaluation_mode)
  
  return(model_assessment)
  
}

####
evaluate_ensemble <- function(model, ensemble_model, evaluation_mode, variable_set, pseudoabs_method)
{
  
  ensemble_assessment_testing <- ensemble_model %>%
    biomod3::get_evaluations() %>%
    as.data.frame() %>%
    rownames_to_column("Metric") %>%
    dplyr::select(-contains("Cutoff"), -contains("Sensitivity"), -contains("Specificity")) %>%
    # dplyr::select(Metric, contains("Testing")) %>%
    gather("Variable", "value", 2:ncol(.)) %>%
    mutate(
      evaluation_set = ifelse(str_detect(Variable, "Testing"), "Fit", "Eval"),
      Variable = Variable %>%
        str_remove_all("_mergedAlgo_Full_AllData.Testing.data") %>%
        str_remove("_mergedAlgo_Full_AllData.Evaluating.data")) %>%
    separate(Variable, into = c("Species", "Method"), sep = "_") %>%
    filter(evaluation_set == "Fit") %>% 
    transmute(
      Metric,
      Algorithm = Method,
      value,
      evaluation_set,
      Species = str_replace(Species, "[.]", "_"),
      PA_set = pseudoabs_method,
      var_set = variable_set,
      mode = evaluation_mode
    )
  
  
  
  ## As the evaluating ensemble doesn't work when using an independent dataset, we evaluate the extrapolationTSS manually
  data_eval <- cbind(get_formal_data(model,'eval.resp.var'),
                     get_formal_data(model,'eval.expl.var'))
  names(data_eval)[1] <- as.character(model@sp.name)
  
  
  ## Retrieve evaluation metrics, take the average over pseudo-absence sets  
  ensemble_assessment <- as.data.frame(biomod2::evaluate(model = ensemble_model,
                                                         data = data_eval,
                                                         stat = c("CSI", 'ROC','TSS')))
  
  ensemble_assessment_eval <- ensemble_assessment %>%
    rownames_to_column("Metric") %>%
    dplyr::select(-contains("Cutoff"), -contains("Sensitivity"), -contains("Specificity")) %>%
    gather("Variable", "value", 2:ncol(.)) %>%
    mutate(
      Variable = Variable %>% str_remove("_mergedAlgo_Full_AllData.Evaluating.data")
    ) %>%
    separate(Variable, into = c("Species", "Method"), sep = "_") %>%
    transmute(
      Metric,
      Algorithm = Method,
      value,
      evaluation_set = "Eval",
      Species = str_replace(Species, "[.]", "_"),
      PA_set = pseudoabs_method,
      var_set = variable_set,
      mode = evaluation_mode)
  
  ensemble_assessment_combined <- rbind(ensemble_assessment_eval, ensemble_assessment_testing)
  
  
  return(ensemble_assessment_combined)
}
