# ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

library(readxl)
library(stringr)

pipeline <- raveio::pipeline("fragility_explorer", paths = "./modules/")

parameters <- list(
  t_window = c(250,200,150,100),
  t_step = c(125,100,75,50),
  nSearch = c(100,200)
)

# # set settings
# t_window <- 100
# t_step <- 50
# nSearch <- 200

# # choose where to export results
# export_path <- paste0("/Volumes/bigbrain/Fragility_Results/",t_window,"-",t_step,"_lambda_n_",nSearch)
# export_path
# read patient_data file

# define which patients to process from patient_data file
#pts <- dipsaus::parse_svec("1-35,37-42,50-60,65-67,75-76,121,125,127,135,157-159")
pts <- 1:length(patient_key$subject)
#pts <- 80:99
# check which patients are being run
patient_key$subject[pts]

# check pt metadata if needed
# for(i in pts){
#   subject_code <- patient_key$subject_code[i]
#   project <- patient_key$project_name[i]
#   subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
#   print(subject$reference_names)
# }
for(t in 1:4) {
  for(n in 1:2) {
    t_window <- parameters$t_window[t]
    t_step <- parameters$t_step[t]
    nSearch <- parameters$nSearch[n]

    # choose where to export results
    print(export_path)

    for(i in pts){

      subject_code <- patient_key$subject_code[i]
      project <- patient_key$project_name[i]
      electrodes <- dipsaus::parse_svec(patient_key$load_electrodes[i])
      display <- electrodes # display all electrodes

      epoch_name <- patient_key$epoch_file_name[i]
      reference_name <- patient_key$reference_name[i]
      condition <- patient_key$condition[i]

      soz <- dipsaus::parse_svec(patient_key$SOZ[i])
      sozc <- electrodes[!(electrodes%in%soz)]

      if(!all(c(soz,sozc) %in% electrodes)){
        warning("Not all electrodes specified in soz are loaded! Will omit unloaded electrodes from soz.")
        soz <- soz[soz %in% electrodes]
        sozc <- sozc[sozc %in% electrodes]
      }

      subject_check <- raveio::validate_subject(paste0(project,"/",subject_code),
                                                method = "basic", verbose = FALSE)

      if(!subject_check$paths$data_path$valid){
        stop("Subject data path is not valid!")
      }

      print(paste0("starting pipeline for pt: ", subject_code, ", ", condition))

      # Fragility ----------------------------------------------
      # Add `path` to force using devel pipeline
      fragility_pipeline <- raveio::pipeline("fragility_explorer", paths = "./modules/")
      #fragility_pipeline$target_table
      #raveio::pipeline_visualize(fragility_pipeline$pipeline_path)
      #fragility_pipeline$get_settings()

      # set subject object from rave
      subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
      elec_list <- subject$get_electrode_table()

      # create export directory for this subject
      export <- file.path(export_path, subject_code)
      raveio::dir_create2(export)

      trial_num <- match(condition, subject$get_epoch(epoch_name)$table$Condition)

      fragility_pipeline$set_settings(
        project_name = project,
        subject_code = subject_code,
        epoch_name = epoch_name,
        epoch_time_window = c(-10,20),
        reference_name = reference_name,
        load_electrodes = electrodes,
        display_electrodes = display,
        condition = condition,
        trial_num = trial_num,
        t_window = t_window,
        t_step = t_step,
        sz_onset = 0,
        lambda = NULL,
        nSearch = nSearch,
        fs_new = 1000,
        soz = soz,
        sozc = sozc
      )

      # for saving output files ---------------------------------------
      # env <- fragility_pipeline$load_shared()
      source("./modules/fragility_explorer/R/shared-functions.R")
      source("./modules/fragility_explorer/R/shared-plots.R")

      tryCatch(
        error = function(e){
          # if there are any errors it will output this error file instead
          if (file.exists(export)) {
            file.create(file.path(export, paste0(subject_code,"_",fragility_pipeline$get_settings("condition"),"_ERROR")))
          }
        },{
          results <- c(fragility_pipeline$run(c("repository", "adj_frag_info")))

          # force evaluation
          #env <- c(fragility_pipeline$eval(c("repository", "adj_frag_info")), shortcut = TRUE)
          #results <- list(repository = env[[1]]$repository, adj_frag_info = env[[1]]$adj_frag_info)

          output_files(results$repository, results$adj_frag_info$epoch,
                       results$adj_frag_info$fragres,
                       fragility_pipeline$get_settings(), export)
        })
    }

  }
}

# DIPSAUS DEBUG START
# load_electrodes <- dipsaus::parse_svec(1-4,7-24,26-36,42-43,46-54,56-69,72-95)
# repository <- results$repository
# adj_frag_info <- results$adj_frag_info
# trial_num = 1
# t_window = 250
# t_step = 125
# display_electrodes <- c(33,34,62:69,88:91)
