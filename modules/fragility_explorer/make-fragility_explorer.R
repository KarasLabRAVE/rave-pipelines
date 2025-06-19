library(targets)
library(raveio)
source("common.R", local = TRUE, chdir = TRUE)
._._env_._. <- environment()
._._env_._.$pipeline <- pipeline_from_path(".")
lapply(sort(list.files(
  "R/", ignore.case = TRUE,
  pattern = "^shared-.*\\.R", 
  full.names = TRUE
)), function(f) {
  source(f, local = ._._env_._., chdir = TRUE)
})
targets::tar_option_set(envir = ._._env_._.)
rm(._._env_._.)
...targets <- list(`__Check_settings_file` = targets::tar_target_raw("settings_path", 
    "settings.yaml", format = "file"), `__Load_settings` = targets::tar_target_raw("settings", 
    quote({
        yaml::read_yaml(settings_path)
    }), deps = "settings_path", cue = targets::tar_cue("always")), 
    input_fs_new = targets::tar_target_raw("fs_new", quote({
        settings[["fs_new"]]
    }), deps = "settings"), input_nSearch = targets::tar_target_raw("nSearch", 
        quote({
            settings[["nSearch"]]
        }), deps = "settings"), input_project_name = targets::tar_target_raw("project_name", 
        quote({
            settings[["project_name"]]
        }), deps = "settings"), input_subject_code = targets::tar_target_raw("subject_code", 
        quote({
            settings[["subject_code"]]
        }), deps = "settings"), input_epoch_name = targets::tar_target_raw("epoch_name", 
        quote({
            settings[["epoch_name"]]
        }), deps = "settings"), input_epoch_time_window = targets::tar_target_raw("epoch_time_window", 
        quote({
            settings[["epoch_time_window"]]
        }), deps = "settings"), input_reference_name = targets::tar_target_raw("reference_name", 
        quote({
            settings[["reference_name"]]
        }), deps = "settings"), input_load_electrodes = targets::tar_target_raw("load_electrodes", 
        quote({
            settings[["load_electrodes"]]
        }), deps = "settings"), input_display_electrodes = targets::tar_target_raw("display_electrodes", 
        quote({
            settings[["display_electrodes"]]
        }), deps = "settings"), input_condition = targets::tar_target_raw("condition", 
        quote({
            settings[["condition"]]
        }), deps = "settings"), input_trial_num = targets::tar_target_raw("trial_num", 
        quote({
            settings[["trial_num"]]
        }), deps = "settings"), input_t_window = targets::tar_target_raw("t_window", 
        quote({
            settings[["t_window"]]
        }), deps = "settings"), input_t_step = targets::tar_target_raw("t_step", 
        quote({
            settings[["t_step"]]
        }), deps = "settings"), input_sz_onset = targets::tar_target_raw("sz_onset", 
        quote({
            settings[["sz_onset"]]
        }), deps = "settings"), input_lambda = targets::tar_target_raw("lambda", 
        quote({
            settings[["lambda"]]
        }), deps = "settings"), input_soz = targets::tar_target_raw("soz", 
        quote({
            settings[["soz"]]
        }), deps = "settings"), input_sozc = targets::tar_target_raw("sozc", 
        quote({
            settings[["sozc"]]
        }), deps = "settings"), load_subject = targets::tar_target_raw(name = "subject", 
        command = quote({
            .__target_expr__. <- quote({
                subject <- raveio::RAVESubject$new(project_name = project_name, 
                  subject_code = subject_code, strict = TRUE)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(subject)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "subject", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = "rave-subject", 
            target_export = "subject", target_expr = quote({
                {
                  subject <- raveio::RAVESubject$new(project_name = project_name, 
                    subject_code = subject_code, strict = TRUE)
                }
                subject
            }), target_depends = c("project_name", "subject_code"
            )), deps = c("project_name", "subject_code"), cue = targets::tar_cue("always"), 
        pattern = NULL, iteration = "list"), load_electrodes = targets::tar_target_raw(name = "loading_elec", 
        command = quote({
            .__target_expr__. <- quote({
                loading_elec <- dipsaus::parse_svec(load_electrodes)
                loading_elec <- subject$electrodes[subject$electrodes %in% 
                  loading_elec]
                if (!length(loading_elec)) {
                  stop("No valid electrode to load!")
                }
            })
            tryCatch({
                eval(.__target_expr__.)
                return(loading_elec)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "loading_elec", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "loading_elec", target_expr = quote({
                {
                  loading_elec <- dipsaus::parse_svec(load_electrodes)
                  loading_elec <- subject$electrodes[subject$electrodes %in% 
                    loading_elec]
                  if (!length(loading_elec)) {
                    stop("No valid electrode to load!")
                  }
                }
                loading_elec
            }), target_depends = c("load_electrodes", "subject"
            )), deps = c("load_electrodes", "subject"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), display_electrodes = targets::tar_target_raw(name = "displayed_elec", 
        command = quote({
            .__target_expr__. <- quote({
                displayed_elec <- dipsaus::parse_svec(display_electrodes)
                displayed_elec <- subject$electrodes[subject$electrodes %in% 
                  displayed_elec]
                if (!length(displayed_elec)) {
                  stop("No valid electrode to load!")
                }
            })
            tryCatch({
                eval(.__target_expr__.)
                return(displayed_elec)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "displayed_elec", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "displayed_elec", target_expr = quote({
                {
                  displayed_elec <- dipsaus::parse_svec(display_electrodes)
                  displayed_elec <- subject$electrodes[subject$electrodes %in% 
                    displayed_elec]
                  if (!length(displayed_elec)) {
                    stop("No valid electrode to load!")
                  }
                }
                displayed_elec
            }), target_depends = c("display_electrodes", "subject"
            )), deps = c("display_electrodes", "subject"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), load_voltage = targets::tar_target_raw(name = "repository", 
        command = quote({
            .__target_expr__. <- quote({
                repository <- raveio::prepare_subject_voltage_with_epoch(subject = subject, 
                  electrodes = loading_elec, epoch_name = epoch_name, 
                  reference_name = reference_name, time_windows = epoch_time_window)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(repository)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "repository", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "repository", target_expr = quote({
                {
                  repository <- raveio::prepare_subject_voltage_with_epoch(subject = subject, 
                    electrodes = loading_elec, epoch_name = epoch_name, 
                    reference_name = reference_name, time_windows = epoch_time_window)
                }
                repository
            }), target_depends = c("subject", "loading_elec", 
            "epoch_name", "reference_name", "epoch_time_window"
            )), deps = c("subject", "loading_elec", "epoch_name", 
        "reference_name", "epoch_time_window"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"), find_adj_and_frag = targets::tar_target_raw(name = "adj_frag_info", 
        command = quote({
            .__target_expr__. <- quote({
                adj_frag_info <- calc_adj_frag(repository = repository, 
                  trial_num = trial_num, t_window = t_window, 
                  t_step = t_step, soz = soz, sozc = sozc, lambda = lambda, 
                  nSearch = nSearch, fs_new = fs_new)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(adj_frag_info)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "adj_frag_info", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "adj_frag_info", target_expr = quote({
                {
                  adj_frag_info <- calc_adj_frag(repository = repository, 
                    trial_num = trial_num, t_window = t_window, 
                    t_step = t_step, soz = soz, sozc = sozc, 
                    lambda = lambda, nSearch = nSearch, fs_new = fs_new)
                }
                adj_frag_info
            }), target_depends = c("repository", "trial_num", 
            "t_window", "t_step", "soz", "sozc", "lambda", "nSearch", 
            "fs_new")), deps = c("repository", "trial_num", "t_window", 
        "t_step", "soz", "sozc", "lambda", "nSearch", "fs_new"
        ), cue = targets::tar_cue("thorough"), pattern = NULL, 
        iteration = "list"), find_quantiles = targets::tar_target_raw(name = "quantiles", 
        command = quote({
            .__target_expr__. <- quote({
                quantiles <- frag_quantile(repository = repository, 
                  f = adj_frag_info$frag, t_window = t_window, 
                  t_step = t_step, soz = soz, sozc = sozc, fs_new = fs_new)
            })
            tryCatch({
                eval(.__target_expr__.)
                return(quantiles)
            }, error = function(e) {
                asNamespace("raveio")$resolve_pipeline_error(name = "quantiles", 
                  condition = e, expr = .__target_expr__.)
            })
        }), format = asNamespace("raveio")$target_format_dynamic(name = NULL, 
            target_export = "quantiles", target_expr = quote({
                {
                  quantiles <- frag_quantile(repository = repository, 
                    f = adj_frag_info$frag, t_window = t_window, 
                    t_step = t_step, soz = soz, sozc = sozc, 
                    fs_new = fs_new)
                }
                quantiles
            }), target_depends = c("repository", "adj_frag_info", 
            "t_window", "t_step", "soz", "sozc", "fs_new")), 
        deps = c("repository", "adj_frag_info", "t_window", "t_step", 
        "soz", "sozc", "fs_new"), cue = targets::tar_cue("thorough"), 
        pattern = NULL, iteration = "list"))
