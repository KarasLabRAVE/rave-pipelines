

module_html <- function(){

  shiny::fluidPage(
    shiny::fluidRow(

      shiny::column(
        width = 3L,
        shiny::div(
          # class = "row fancy-scroll-y stretch-inner-height",
          class = "row screen-height overflow-y-scroll",
          shiny::column(
            width = 12L,

            ravedash::input_card(
              title = "Fragility Parameters",

              ravedash::flex_group_box(
                title = "Basic",
                # parameters needed: sz_onset, display_electrodes, trial_num
                shidashi::flex_container(
                  shidashi::flex_item(

                    # Bin the time using `t_window` with min=100ms, max=1000ms, step=100ms
                    shiny::sliderInput(
                      inputId = ns("t_window"),
                      label = "Time Window Size (ms)",
                      min = 50, max = 1000, step = 50, value = 250,
                      width = "100%", post = " ms"
                    )
                  ),

                  # new line
                  shidashi::flex_break(),

                  shidashi::flex_item(
                    shiny::selectInput(
                      inputId = ns("t_step_percentage"),
                      label = "Time Step relative to window size",
                      choices = c("100%", "50%", "25%", "10%"),
                      selected = "50%",
                      width = "100%"
                    )
                  ),
                  shidashi::flex_item(
                    shiny::selectInput(
                      inputId = ns("nSearch"),
                      label = "# of lambdas to search through",
                      choices = 100:200,
                      width = "100%"
                    )
                  )
                ),


                # new line
                shidashi::flex_break(),

                shidashi::flex_container(
                  shidashi::flex_item(
                    shiny::textInput(
                      inputId = ns("soz"),
                      label = "SOZ Electrodes",
                      value = "",
                      width = "100%"
                    )
                  ),
                  shidashi::flex_item(
                    shiny::selectInput(
                      inputId = ns("condition"),
                      label = "Choose trial",
                      choices = 1:4,
                      width = "100%"
                    )
                  )
                ),
                shidashi::flex_container(
                  shidashi::flex_item(
                    shiny::checkboxInput(
                      inputId = ns("fs_new_TF"),
                      label = "Resample?",
                      width = "100%"
                    )
                  ),
                  shidashi::flex_item(
                    shiny::textInput(
                      inputId = ns("fs_new"),
                      label = "Resample rate:",
                      width = "100%"
                    )
                  )
                )
              ),
              ravedash::flex_group_box(
                title = "Display",
                  shidashi::flex_item(
                    shiny::textInput(
                      inputId = ns("display_electrodes"),
                      label = "Electrodes to display",
                      value = "",
                      width = "100%"
                    )
                  ),
                  shidashi::flex_item(
                    shiny::sliderInput(
                      inputId = ns("sz_onset"),
                      label = "Seizure Onset Marker",
                      min = -10, max = 10, value = 0, width = "100%"
                    )
                  ),
                shidashi::flex_item(
                  shiny::checkboxInput(
                    inputId = ns("sepsoz"),
                    label = "Separate SOZ electrodes? (displayed in red)",
                    value = T
                  )
                )
              ),
              ravedash::flex_group_box(
                title = "Thresholding",
                shidashi::flex_item(
                  shiny::checkboxInput(
                    inputId = ns("thresholding"),
                    label = "Threshold?",
                    value = F
                  )
                ),
                shidashi::flex_item(
                  shiny::textInput(
                    inputId = ns("buckets"),
                    label = "Threshold limits (must start with 0 and end with 1)",
                    value = "0, 0.25, 0.5, 0.9, 1"
                  )
                ),
                shidashi::flex_item(
                  shiny::numericInput(
                    inputId = ns("moving_avg_width"),
                    label = "Moving average width (set at 1 for no moving average)",
                    value = 10
                  )
                )
              ),
              ravedash::flex_group_box(
                title = "Advanced",
                shidashi::flex_item(
                  shiny::textInput(
                    inputId = ns("lambda"),
                    label = "Choose lambda (leave blank for auto-calculation)",
                    value = "",
                    width = "100%"
                  )
                ),
                shidashi::flex_item(
                  shiny::checkboxInput(
                    inputId = ns("ranked"),
                    label = "Enhance contrast by ranking Fragility?",
                    value = F,
                    width = "100%"
                  )
                )
                # shidashi::flex_item(
                #   shiny::numericInput(
                #     inputId = ns("threshold"),
                #     label = "Threshold value",
                #     min = -1, max = 1, value = 0.5, width = "100%"
                #   )
                # )
              )
            )
            # electrode_selector$ui_func(),
            #
            # comp_condition_groups$ui_func(),
            #
            # baseline_choices$ui_func(),
            #
            # comp_analysis_ranges$ui_func()

          )
        )
      ),

      shiny::column(
        width = 9L,
        shiny::div(
          class = "row screen-height overflow-y-scroll output-wrapper",
          shiny::column(
            width = 12L,
            ravedash::output_card(
              'Voltage Viewer',
              class_body = "no-padding fill-width height-450 min-height-450 resize-vertical",
              shiny::div(
                class = 'position-relative fill',
                shiny::plotOutput(ns("v_plot"), width = '100%', height = "100%")
              )
            ),
            ravedash::output_card(
              'Fragility Map',
              class_body = "no-padding fill-width height-450 min-height-450 resize-vertical",
              shiny::div(
                class = 'position-relative fill',
                shiny::plotOutput(ns("f_plot"), width = '100%', height = "100%")
              )
            ),
            ravedash::output_card(
              'Average Fragility Over Time',
              class_body = "no-padding fill-width height-450 min-height-450 resize-vertical",
              shiny::div(
                class = 'position-relative fill',
                shiny::plotOutput(ns("mean_f_plot"), width = '100%', height = "100%")
              )
            ),
            ravedash::output_card(
              'Quantiles',
              class_body = "no-padding fill-width height-450 min-height-450 resize-vertical",
              shiny::div(
                class = 'position-relative fill',
                shiny::plotOutput(ns("q_plot"), width = '100%', height = "100%")
              )
            )
          )
        )
      )

    )
  )
}
