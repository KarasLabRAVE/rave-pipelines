# ravedash::debug_modules(module_root = rstudioapi::getActiveProject())

library(readxl)
library(stringr)

pipeline <- raveio::pipeline("fragility_explorer", paths = "./modules/")

# choose where to export results
export_path <- "/karaslab/Results_FragilityEEGDataset"

# define which patients to process from patient_data file
pts <- dipsaus::parse_svec("1-35,37-42,50-60,65-67,75-76,121,125,127,135,157-159")

# read patient_data file
patient_key <- read.csv("/Volumes/bigbrain/Multipatient/patient_data_all.csv")

# check which patients are being run
patient_key$subject[pts]

# check pt metadata if needed
# for(i in pts){
#   subject_code <- patient_key$subject_code[i]
#   project <- patient_key$project_name[i]
#   subject <- raveio::as_rave_subject(paste0(project,"/",subject_code))
#   print(subject$reference_names)
# }

# set settings
t_window <- 250
t_step <- 125

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
  subject_check$paths$data_path$valid

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
    lambda = FALSE,
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
      results <- c(fragility_pipeline$run(c("repository", "adj_frag_info","quantiles")))

      # force evaluation
      #env <- c(fragility_pipeline$eval(c("repository", "adj_frag_info","quantiles")), shortcut = TRUE)
      #results <- list(repository = env[[1]]$repository, adj_frag_info = env[[1]]$adj_frag_info, quantiles = env[[1]]$quantiles)

      moving_average_width <- 10

      # save unranked results
      output_files(results$repository, results$adj_frag_info$frag, results$quantiles,
                   fragility_pipeline$get_settings(),export,"norank", moving_average_width)

      # save ranked results
      output_files(results$repository, results$adj_frag_info$frag_ranked, results$quantiles,
                   fragility_pipeline$get_settings(),export,"ranked", moving_average_width)

      # save R2 results
      output_R2(results$repository, results$adj_frag_info$R2, results$adj_frag_info$lambdas,
                fragility_pipeline$get_settings(),export)

      # print results to pdf
      pdf_path <- file.path(export, paste0(subject_code,'_',fragility_pipeline$get_settings("condition"),"_reconstruction.pdf"))
      grDevices::pdf(pdf_path, width = 12, height = 7)
      par(mfrow=c(2,1),mar=rep(2,4))

      # voltage reconstruction example using electrode 1
      g <- do.call(voltage_recon_plot, c(results[1:2],
                                         list(fragility_pipeline$get_settings("t_window"),
                                              fragility_pipeline$get_settings("t_step"),
                                              fragility_pipeline$get_settings("trial_num"),
                                              timepoints = 1:1000,
                                              elec_num = 1,
                                              lambda = fragility_pipeline$get_settings("lambda"))
      ))
      print(g)

      grDevices::dev.off()
    })
}

# DIPSAUS DEBUG START
# repository <- results$repository
# adj_frag_info <- results$adj_frag_info
# trial_num = 3
# t_window = 250
# t_step = 125
# timepoints = 1:1000
# elec_num = 1
# percentile = 0.1
# display_electrodes <- c(33,34,62:69,88:91)
