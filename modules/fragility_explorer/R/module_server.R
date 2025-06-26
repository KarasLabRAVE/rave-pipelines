
module_server <- function(input, output, session, ...){


  # Local reactive values, used to store reactive event triggers
  local_reactives <- shiny::reactiveValues(
    update_outputs = NULL
  )

  # Local non-reactive values, used to store static variables
  local_data <- dipsaus::fastmap2()

  # get server tools to tweek
  server_tools <- get_default_handlers(session = session)

  # # Run analysis once the following input IDs are changed
  # # This is used by auto-recalculation feature
  # server_tools$run_analysis_onchange(
  #   component_container$get_input_ids(c(
  #     "electrode_text", "baseline_choices",
  #     "analysis_ranges", "condition_groups"
  #   ))
  # )

  # Register event: main pipeline need to run
  shiny::bindEvent(
    ravedash::safe_observe({

      # Run analysis button is clicked

      # # Invalidate previous results (stop them because they are no longer needed)
      # if(!is.null(local_data$results)) {
      #   local_data$results$invalidate()
      #   ravedash::logger("Invalidating previous run", level = "trace")
      # }
      #

      # Collect input data
      # settings <- component_container$collect_settings(ids = c(
      #   "display_electrodes", "trial_num", "t_window", "t_step", "sz_onset"
      # ))

      display <- dipsaus::parse_svec(input$display_electrodes)

      trial_num <- which(input$condition == component_container$data$trial_choices)
      t_step <- input$t_window * (as.numeric(gsub("%","",input$t_step_percentage)) / 100)
      if (input$lambda == "") {
        lambda <- NULL
      } else {
        lambda <- as.numeric(input$lambda)
      }

      electrodes <- component_container$data$repository$electrode_table$Electrode

      sozc <- electrodes[!(electrodes%in%dipsaus::parse_svec(input$soz))]

      if(input$fs_new_TF) {
        fs_new <- input$fs_new
      } else {
        fs_new <- NULL
      }

      pipeline$set_settings(
        display_electrodes = display,
        condition = input$condition,
        trial_num = trial_num,
        t_window = input$t_window,
        t_step = t_step,
        sz_onset = input$sz_onset,
        lambda = lambda,
        nSearch = as.numeric(input$nSearch),
        fs_new = fs_new,
        soz = dipsaus::parse_svec(input$soz),
        sozc = sozc
      )

      #print(pipeline$get_settings())

      #' Run pipeline without blocking the main session
      #' The trick to speed up is to set
      #' `async=TRUE` will run the pipeline in the background
      #' `shortcut=TRUE` will ignore the dependencies and directly run `names`
      #' `names` are the target nodes to run
      #' `scheduler="none"` will try to avoid starting any schedulers and
      #' run targets sequentially. Combined with `callr_function=NULL`,
      #' scheduler's overhead can be removed.
      #' `type="smart"` will start `future` plan in the background, allowing
      #' multicore calculation
      results <- pipeline$run(
        as_promise = TRUE,
        scheduler = "none",
        # type = "smart",
        # callr_function = NULL,
        # progress_title = "Calculating in progress",
        # async = TRUE,
        # check_interval = 0.1,
        # shortcut = TRUE,
        names = c(
          "adj_frag_info",
          "quantiles"
        )
      )


      local_data$results <- results
      ravedash::logger("Scheduled: ", pipeline$pipeline_name,
                       level = 'debug', reset_timer = TRUE)

      results$promise$then(
        onFulfilled = function(...){
          ravedash::logger("Fulfilled: ", pipeline$pipeline_name,
                           level = 'debug')
          shidashi::clear_notifications(class = "pipeline-error")
          local_reactives$update_outputs <- Sys.time()
          return(TRUE)
        },
        onRejected = function(e, ...){
          msg <- paste(e$message, collapse = "\n")
          if(inherits(e, "error")){
            ravedash::logger(msg, level = 'error')
            ravedash::logger(traceback(e), level = 'error', .sep = "\n")
            shidashi::show_notification(
              message = msg,
              title = "Error while running pipeline", type = "danger",
              autohide = FALSE, close = TRUE, class = "pipeline-error"
            )
          }
          return(msg)
        }
      )

      return()

    }),
    server_tools$run_analysis_flag(),
    ignoreNULL = TRUE, ignoreInit = TRUE
  )


  # (Optional) check whether the loaded data is valid
  shiny::bindEvent(
    ravedash::safe_observe({
      loaded_flag <- ravedash::watch_data_loaded()
      if(!loaded_flag){ return() }
      new_repository <- pipeline$read("repository")
      if(!inherits(new_repository, "rave_prepare_subject_voltage_with_epoch")){
        ravedash::logger("Repository read from the pipeline, but it is not an instance of `rave_prepare_subject_voltage_with_epoch`. Abort initialization", level = "warning")
        return()
      }
      ravedash::logger("Repository read from the pipeline; initializing the module UI", level = "debug")

      # check if the repository has the same subject as current one
      old_repository <- component_container$data$repository
      if(inherits(old_repository, "rave_prepare_subject_voltage_with_epoch")){

        if( !attr(loaded_flag, "force") &&
            identical(old_repository$signature, new_repository$signature) ){
          ravedash::logger("The repository data remain unchanged ({new_repository$subject$subject_id}), skip initialization", level = "debug", use_glue = TRUE)
          return()
        }
      }


      # Reset preset UI & data
      component_container$reset_data()
      component_container$data$repository <- new_repository
      component_container$initialize_with_new_data()

      # construct initialization data that can be reused elsewhere
      epoch_table <- new_repository$epoch_table
      component_container$data$trial_choices <- sprintf("%s (%d)", epoch_table$Condition, epoch_table$Trial)

      # TODO: reset UIs to default
      # new_repository <- pipeline$read("repository")

      # Update the trial selection input
      shiny::updateSelectInput(
        session = session,
        inputId = "condition",
        choices = component_container$data$trial_choices,
        selected = pipeline$get_settings("condition")
      )

      shiny::updateTextInput(
        session = session,
        inputId = "soz",
        value = dipsaus::deparse_svec(pipeline$get_settings("soz"))
      )

      shiny::updateTextInput(
        session = session,
        inputId = "display_electrodes",
        value = dipsaus::deparse_svec(component_container$data$repository$electrode_list)
      )

      shiny::updateSliderInput(
        session = session,
        inputId = "sz_onset",
        min = component_container$data$repository$time_windows[[1]][1],
        max = component_container$data$repository$time_windows[[1]][2],
        value = 0
      )

      shiny::updateTextInput(
        session = session,
        inputId = "fs_new",
        value = component_container$data$repository$sample_rate
      )

      # shiny::updateSliderInput(
      #   session = session,
      #   inputId = "threshold_limits",
      #   min = component_container$data$repository$time_windows[[1]][1],
      #   max = component_container$data$repository$time_windows[[1]][2],
      #   value = c(0,component_container$data$repository$time_windows[[1]][2])
      # )

      # Reset outputs
      # shidashi::reset_output("collapse_over_trial")

    }, priority = 1001),
    ravedash::watch_data_loaded(),
    ignoreNULL = FALSE,
    ignoreInit = FALSE
  )





  # Register outputs
  ravedash::register_output(
    outputId = "v_plot",
    output_type = "image",
    title = "Voltage Reconstruction",
    shiny::renderPlot({
      shiny::validate(
        shiny::need(
          length(local_reactives$update_outputs) &&
            !isFALSE(local_reactives$update_outputs),
          message = "Please run the module first"
        )
      )
      shiny::validate(
        shiny::need(
          isTRUE(local_data$results$valid),
          message = "One or more errors while executing pipeline. Please check the notification."
        )
      )

      results <- pipeline$read(var_names = c("repository","adj_frag_info"))
      display_electrodes <- dipsaus::parse_svec(input$display_electrodes)

      voltage_plot(results$repository, results$adj_frag_info, display_electrodes,
                   pipeline$get_settings("soz"), pipeline$get_settings("sozc"), input$sepsoz)

      # do.call(voltage_recon_plot, c(results[1:2],
      #                               list(pipeline$get_settings("t_window"),
      #                                    pipeline$get_settings("t_step"),
      #                                    pipeline$get_settings("trial_num"),
      #                                    timepoints = 1:1000,
      #                                    elec_num = 1,
      #                                    lambda = pipeline$get_settings("lambda"))

      # ))
    })
  )

  ravedash::register_output(
    outputId = "f_plot",
    output_type = "image",
    title = "Fragility Map",
    shiny::renderPlot({
      shiny::validate(
        shiny::need(
          length(local_reactives$update_outputs) &&
            !isFALSE(local_reactives$update_outputs),
          message = "Please run the module first"
        )
      )
      shiny::validate(
        shiny::need(
          isTRUE(local_data$results$valid),
          message = "One or more errors while executing pipeline. Please check the notification."
        )
      )

      results <- pipeline$read(var_names = c("repository","adj_frag_info"))

      fragility_plot(results$repository, results$adj_frag_info, pipeline$get_settings(),
                     dipsaus::parse_svec(input$display_electrodes),
                     input$ranked, input$sepsoz, input$thresholding, input$buckets)
    })
  )

  ravedash::register_output(
    outputId = "mean_f_plot",
    output_type = "image",
    title = "Average Fragility Over Time",
    shiny::renderPlot({
      shiny::validate(
        shiny::need(
          length(local_reactives$update_outputs) &&
            !isFALSE(local_reactives$update_outputs),
          message = "Please run the module first"
        )
      )
      shiny::validate(
        shiny::need(
          isTRUE(local_data$results$valid),
          message = "One or more errors while executing pipeline. Please check the notification."
        )
      )

      results <- pipeline$read(var_names = c("repository","adj_frag_info"))

      avg_f_over_time_plot(results$repository, results$adj_frag_info, pipeline$get_settings(),input$moving_avg_width)
    })
  )

  ravedash::register_output(
    outputId = "q_plot",
    output_type = "image",
    title = "Quantile Plot",
    shiny::renderPlot({
      shiny::validate(
        shiny::need(
          length(local_reactives$update_outputs) &&
            !isFALSE(local_reactives$update_outputs),
          message = "Please run the module first"
        )
      )
      shiny::validate(
        shiny::need(
          isTRUE(local_data$results$valid),
          message = "One or more errors while executing pipeline. Please check the notification."
        )
      )

      results <- pipeline$read(var_names = c("repository","quantiles"))

      quantiles_plot(results$repository, results$quantiles, pipeline$get_settings(), input$thresholding, input$buckets)
    })
  )
}
