voltage_plot <- function(repository, adj_frag_info, display_electrodes, soz, sozc, sepsoz){

  display_electrodes_i <- match(display_electrodes, repository$electrode_list)
  vmat <- adj_frag_info$voltage[,display_electrodes_i]
  times <- repository$voltage$dimnames$Time

  colorelec <- ifelse(display_electrodes%in%soz,"red","black")

  vmat <- vmat[seq(1,nrow(vmat),4),] # compress by factor of 4 to save plotting time
  times <- times[seq(1,length(times),4)]

  t <- dim(vmat)[1] # of timepoints
  N <- dim(vmat)[2] # of electrodes
  z <- (vmat - mean(vmat)) / sd(vmat) # z score

  display_electrode_names <- repository$electrode_table$Label[display_electrodes_i]

  if (sepsoz){
    # create voltage plot with soz electrodes separated from sozc electrodes
    soz_i <- match(soz,display_electrodes)
    soz_i <- soz_i[!is.na(soz_i)]
    sozc_i <- match(sozc,display_electrodes)
    sozc_i <- sozc_i[!is.na(sozc_i)]
    z <- z[,c(soz_i,sozc_i)]
    colorelec <- colorelec[c(soz_i,sozc_i)]
    display_electrode_names <- display_electrode_names[c(soz_i,sozc_i)]
  }

  z <- z[,ncol(z):1] # reverse matrix for display

  par(mar=c(3,6,0,0))
  rutabaga::plot_clean(repository$voltage$dimnames$Time, -3:N*3)
  axis(1)
  abline(h=0:(N-1)*3, col = "skyblue") # add  horizontal lines

  for(ii in 1:ncol(vmat)) {
    lines(x = times, y = z[,ii]+(ii-1)*3) # plot lines
  }

  #axis(2, at=0:(N-1)*3, label = repository$electrode_table$Label[display_electrodes_i], las = 1, tcl=0) # electrode labels
  axis(2, at=0:(N-1)*3, label = FALSE, tcl=0)

  max_labels <- floor(par("pin")[2] / 0.1) # 0.1 inches per label
  spacing <- ceiling(N / max_labels)
  label_idx <- seq(1, N, by = spacing)
  label_positions <- (label_idx - 1) * 3
  label_text <- display_electrode_names[label_idx]
  label_colors <- colorelec[label_idx]

  text(x = par("usr")[1] - 0.5,  # left of plot
       y = label_positions,
       labels = rev(label_text), # reversed for display
       col = rev(label_colors), # reversed for display
       xpd = TRUE,
       adj = 1,
       cex = 0.8) # electrode labels
}

fragility_plot <- function(repository, adj_frag_info, pipeline_settings, display_electrodes, ranked, sepsoz, thresholding, buckets) {

  subject_code <- pipeline_settings$subject_code
  t_window <- pipeline_settings$t_window
  t_step <- pipeline_settings$t_step
  sz_num <- pipeline_settings$condition
  soz <- pipeline_settings$soz
  sozc <- pipeline_settings$sozc
  sz_onset <- pipeline_settings$sz_onset

  sozNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% soz]
  sozcNames <- repository$electrode_table$Label[repository$electrode_table$Electrode %in% sozc]
  #soz_i <- match(soz,repository$electrode_list)
  #sozc_i <- match(sozc,repository$electrode_list)

  if(ranked){
    note <- "ranked"
    f <- adj_frag_info$frag_ranked
  } else {
    note <- "unranked"
    f <- adj_frag_info$frag
  }

  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- floor((n_tps - t_window) / t_step) + 1
  epoch_time_window <- repository$time_windows[[1]]
  fs <- round(repository$sample_rate,-1)

  # make soz electrodes separate color for heatmap labels
  colorelec <- ifelse(display_electrodes%in%soz,"red","black")
  display_i <- match(display_electrodes,repository$electrode_list)

  if(any(repository$electrode_table$Label == "NoLabel")) {
    #elec_names <- repository$electrode_table$Electrode[match(c(soz,sozc), repository$electrode_table$Electrode)]
    #elec_names <- as.character(elec_names)

    display_electrode_names <- repository$electrode_table$Electrode[match(display_electrodes,repository$electrode_table$Electrode)]
    display_electrode_names <- as.character(display_electrode_names)
    names(colorelec) <- as.character(display_electrodes)
  } else {
    #elec_names <- repository$electrode_table$Label[match(c(soz,sozc), repository$electrode_table$Electrode)]

    display_electrode_names <- repository$electrode_table$Label[match(display_electrodes,repository$electrode_table$Electrode)]
    names(colorelec) <- dimnames(f)$Electrode[display_i]
  }

  f_display <- f[display_electrode_names,]

  if (sepsoz){
    # create fragility map with soz electrodes separated from sozc electrodes
    soz_i <- match(sozNames,attr(f_display, "dimnames")$Electrode)
    soz_i <- soz_i[!is.na(soz_i)]
    sozc_i <- match(sozcNames,attr(f_display, "dimnames")$Electrode)
    sozc_i <- sozc_i[!is.na(sozc_i)]
    f_display <- f_display[c(soz_i,sozc_i),]
    colorelec <- colorelec[c(soz_i,sozc_i)]
    display_electrode_names <- display_electrode_names[c(soz_i,sozc_i)]
  }

  # convert time windows to seconds
  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]

  # raw fragility map
  fplot <- expand.grid(Time = stimes, Electrode = display_electrode_names)

  fplot$Value <- c(t(f_display))

  titlepng=paste(subject_code,as.character(sz_num),note,"fragility",sep=" ")

  # only display selected electrodes

  if (thresholding) {
    fplot$Value <- threshold_buckets(as.matrix(fplot$Value),as.numeric(unlist(strsplit(buckets, ","))))
  }

  fplot$Electrode <- factor(fplot$Electrode, levels = rev(unique(fplot$Electrode))) # reversed for display

  electrodes <- rev(unique(fplot$Electrode))
  N <- length(electrodes)
  max_labels <- 50
  spacing <- max(1, ceiling(N / max_labels))
  keep_idx <- seq(1, N, by = spacing)
  visible_labels <- electrodes[keep_idx]

  y_labels <- ifelse(electrodes %in% visible_labels, levels(electrodes), "")
  y_colors <- ifelse(electrodes %in% visible_labels, rev(colorelec), "transparent")

  ggplot(fplot, aes(x = Time, y = Electrode, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Electrode") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
    geom_vline(xintercept = sz_onset, color = "blue") +
    scale_y_discrete(breaks = electrodes, labels = y_labels) +
    theme_minimal() +
    theme(
      axis.text.y = ggtext::element_markdown(size = 10,colour=rev(colorelec)), # reversed for display
    )
}

avg_f_over_time_plot <- function(repository, adj_frag_info, pipeline_settings, moving_avg_width = 1) {

  subject_code <- pipeline_settings$subject_code
  t_window <- pipeline_settings$t_window
  t_step <- pipeline_settings$t_step
  sz_num <- pipeline_settings$condition
  soz <- pipeline_settings$soz
  sozc <- pipeline_settings$sozc
  sz_onset <- pipeline_settings$sz_onset
  epoch_time_window <- repository$time_windows[[1]]

  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- floor((n_tps - t_window) / t_step) + 1
  epoch_time_window <- repository$time_windows[[1]]
  fs <- round(repository$sample_rate,-1)

  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]

  # average f over time windows
  mean_f <- mean_f_calc(repository,adj_frag_info$frag,soz,sozc)

  if (moving_avg_width > 1) {
    mean_f <- lapply(mean_f, function(x) {
      x_ma <- moving_average(x,moving_avg_width)
      return(x_ma)
    })
  }

  df <- data.frame(stimes,mean_f)

  titlepng=paste(subject_code,"Seizure",as.character(sz_num),"Average Fragility Over Time",sep=" ")

  ggplot(df, aes(stimes)) +
    geom_line(aes(y=mean_f_soz, color = "SOZ +/- sem")) +
    geom_line(aes(y=mean_f_sozc, color = "SOZc +/- sem")) +
    geom_ribbon(aes(ymin = mean_f_soz - se_f_soz, ymax = mean_f_soz + se_f_soz), fill = "indianred3", alpha = 0.7) +
    geom_ribbon(aes(ymin = mean_f_sozc - se_f_sozc, ymax = mean_f_sozc + se_f_sozc), fill = "grey30", alpha = 0.6) +
    geom_vline(xintercept = sz_onset, color = "blue") +
    scale_color_manual(values = c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")) +
    labs(x = "Time", y = "Average Fragility", color = "Legend") +
    ggtitle(titlepng)
}

quantiles_plot <- function(repository, quantile_results, pipeline_settings, thresholding, buckets) {

  subject_code <- pipeline_settings$subject_code
  sz_num <- pipeline_settings$condition
  sz_onset <- pipeline_settings$sz_onset

  # quantile map
  titlepng=paste(subject_code,as.character(sz_num),"Quantiles",sep=" ")

  if (thresholding) {
    quantile_results$q_plot$Value <- threshold_buckets(as.matrix(quantile_results$q_plot$Value),as.numeric(unlist(strsplit(buckets, ","))))
  }

  ggplot(quantile_results$q_plot, aes(x = Time, y = Stats, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Statistic") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
    theme_minimal() +
    geom_vline(xintercept = sz_onset, color = "blue") +
    theme(
      axis.text.y = element_text(size = 10),     # Adjust depending on electrodes
    )
}

output_files <- function(repository, epoch, fragres, pipeline_settings, export) {
  raveio::dir_create2(paste0(export))

  # save raw fragility csvs
  sz_num <- pipeline_settings$condition

  raveio::safe_write_csv(
    fragres@frag,
    file.path(export, paste0(subject_code, "_", sz_num,"_fragility.csv"))
  )

  # save fragility heatmaps
  sozIndex <- Epoch::metaData(epoch)$sozNames
  EZFragility::plotFragHeatmap(fragres, groupIndex = sozIndex)
  ggsave(paste0(export,"/",subject_code,"_",sz_num,"_map.png"))

  # save quantiles csv
  fragstats <- fragStat(fragres, groupIndex = sozIndex)
  raveio::safe_write_csv(
    fragstats@qmatrix,
    file.path(export, paste0("/",subject_code,"_",sz_num,"_quantile.csv"))
  )

  # save quantiles map
  EZFragility::plotFragQuantile(fragres, groupIndex = sozIndex)
  ggsave(paste0(export,"/",subject_code,"_",sz_num,"_qmap.png"))

  # save mean fragility csv
  meandata <- cbind(fragstats@meanGroup, fragstats@meanRef, fragstats@sdGroup, fragstats@sdRef, fragres@startTimes)
  colnames(meandata) <- c("mean_f_group","mean_f_ref","sd_f_group","sd_f_ref","time")
  raveio::safe_write_csv(
    meandata,
    file.path(export, paste0("/",subject_code,"_",sz_num,"_meandata.csv"))
  )

  # save mean fragility plot
  EZFragility::plotFragDistribution(fragres, groupIndex = sozIndex)
  ggsave(paste0(export,"/",subject_code,"_",sz_num,"_meanplot.png"))

  # save R2
  raveio::safe_write_csv(
    rbind(fragres$R2, fragres$lambdas),
    file.path(export, paste0("/",subject_code,"_",sz_num,"_R2.csv"))
  )
}

output_files_old <- function(repository, f, quantile_results, pipeline_settings, export, note, moving_avg_width = NULL) {

  raveio::dir_create2(paste0(export,"/",note))

  subject_code <- pipeline_settings$subject_code
  sz_num <- pipeline_settings$condition
  t_window <- pipeline_settings$t_window
  t_step <- pipeline_settings$t_step
  soz <- pipeline_settings$soz
  sozc <- pipeline_settings$sozc
  sz_onset <- pipeline_settings$sz_onset
  epoch_time_window <- pipeline_settings$epoch_time_window

  # save fragility matrix results to csv
  raveio::safe_write_csv(
    f,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_fragility_",note,".csv"))
  )

  # fragility heatmap
  colorelec <- rep("black",length(c(soz,sozc)))
  colorelec[1:length(soz)]="red"

  titlepng=paste(subject_code,as.character(sz_num),note,sep=" ")

  ggplot(quantile_results$fplot, aes(x = Time, y = Electrode, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Electrode") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5)+  #
    geom_vline(xintercept = sz_onset, color = "black", linewidth = 1.5, linetype = 5) +
    theme_minimal() +
    theme(
      axis.text.y = ggtext::element_markdown(size = 5,colour=colorelec),     # Adjust depending on electrodes
    )
  img <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_map_",note,".png")
  ggsave(img)

  # quantile
  raveio::safe_write_csv(
    quantile_results$q_matrix,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_quantile_",note,".csv"))
  )

  # quantile map
  titlepng=paste(subject_code,as.character(sz_num),"Quantiles",note,sep=" ")

  ggplot(quantile_results$q_plot, aes(x = Time, y = Stats, fill = Value)) +
    geom_tile() +
    ggtitle(titlepng)+
    labs(x = "Time (s)", y = "Statistic") +
    scale_fill_gradient2(low="navy", mid="white", high="red",midpoint=0.5) +  #
    geom_vline(xintercept = sz_onset, color = "black", linewidth = 1.5, linetype = 5) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),     # Adjust depending on electrodes
    )
  q_image <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_qmap_",note,".png")
  ggsave(q_image)

  # average f over time windows
  mean_f <- mean_f_calc(repository,f,soz,sozc)

  stimes <- as.numeric(attr(quantile_results$q_matrix,"dimnames")$Time)

  if (!is.null(moving_avg_width)) {
    mean_f <- raveio::lapply_async(mean_f, function(x) {
      x_ma <- moving_average(x,moving_avg_width)
      return(x_ma)
    })
  }

  df <- data.frame(stimes,mean_f)

  titlepng=paste(subject_code,"Seizure",as.character(sz_num),"Avg Fragility Over Time",note,sep=" ")

  ggplot(df, aes(stimes)) +
    geom_line(aes(y=mean_f_soz, color = "SOZ +/- sem")) +
    geom_line(aes(y=mean_f_sozc, color = "SOZc +/- sem")) +
    geom_ribbon(aes(ymin = mean_f_soz - se_f_soz, ymax = mean_f_soz + se_f_soz), fill = "indianred3", alpha = 0.7) +
    geom_ribbon(aes(ymin = mean_f_sozc - se_f_sozc, ymax = mean_f_sozc + se_f_sozc), fill = "grey30", alpha = 0.6) +
    geom_vline(xintercept = sz_onset, color = "black", linewidth = 1.5, linetype = 5) +
    scale_color_manual(values = c("SOZ +/- sem" = "red", "SOZc +/- sem" = "black")) +
    labs(x = "Time", y = "Average Fragility", color = "Legend") +
    ggtitle(titlepng)

  mean_plot <- paste0(export,"/",note,"/",subject_code,"_",sz_num,"_meanplot_",note,".png")
  ggsave(mean_plot)

  line_plot_df <- data.frame(mean_f,attr(quantile_results$q_matrix, "dimnames")$Time)
  names(line_plot_df)[5] <- "time"

  raveio::safe_write_csv(
    line_plot_df,
    file.path(export, paste0(note,"/",subject_code, "_", sz_num,"_meandata_",note,".csv"))
  )
}

output_R2_old <- function(repository, R2, lambdas, pipeline_settings, export, fs_new = NULL) {

  subject_code <- pipeline_settings$subject_code
  sz_num <- pipeline_settings$condition
  t_window <- pipeline_settings$t_window
  t_step <- pipeline_settings$t_step
  soz <- pipeline_settings$soz
  sozc <- pipeline_settings$sozc
  sz_onset <- pipeline_settings$sz_onset
  epoch_time_window <- pipeline_settings$epoch_time_window

  # save R2 to csv with optimal lambdas
  raveio::safe_write_csv(
    rbind(R2, lambdas),
    file.path(export, paste0("/",subject_code,"_",sz_num,"_R2.csv"))
  )

  n_tps <- length(repository$voltage$dimnames$Time)
  n_elec <- length(repository$voltage$dimnames$Electrode)
  n_steps <- dim(R2)[2]
  epoch_time_window <- repository$time_windows[[1]]
  fs <- round(repository$sample_rate,-1)
  if(any(repository$electrode_table$Label == "NoLabel")) {
    elec_names <- repository$electrode_table$Electrode[match(c(soz,sozc), repository$electrode_table$Electrode)]
    elec_names <- as.character(elec_names)
  } else {
    elec_names <- repository$electrode_table$Label[match(c(soz,sozc), repository$electrode_table$Electrode)]
  }

  # create fragility map with soz electrodes separated from sozc electrodes

  if(!is.null(fs_new)){
    fs <- fs_new
  }

  stimes <- (seq_len(n_steps)-1)*t_step/fs+epoch_time_window[1]
  R2plot <- expand.grid(Time = stimes, Electrode = elec_names)
  R2plot$Value <- c(t(R2))
  R2plot$lambda <- lambdas

  colorelec <- rep("black",length(c(soz,sozc)))
  colorelec[1:length(soz)]="red"

  df <- data.frame(stimes,lambdas)

  titlepng=paste(subject_code, "Seizure",as.character(sz_num),"R2",sep=" ")

  ggplot() +
    geom_tile(data = R2plot, aes(x = Time, y = Electrode, fill = Value)) +
    geom_line(data = df, aes(x = stimes, y=lambdas, color = "Optimal lambda")) +
    ggtitle(titlepng)+
    geom_vline(xintercept = sz_onset, color = "blue") +
    labs(x = "Time (s)", y = "Electrode", color = "Legend") +
    scale_fill_gradient2(low="red4", mid="red", high="white",midpoint=0) +  #
    scale_color_manual(values = c("Optimal lambda" = "black")) +
    theme_minimal()
  img <- paste0(export, "/",subject_code,"_",sz_num,"_R2.png")
  ggsave(img)

  titlepng=paste(subject_code,"Seizure",as.character(sz_num),"Optimal Lambda Over Time",sep=" ")

  ggplot(df, aes(stimes)) +
    geom_point(aes(y=lambdas)) +
    geom_line(aes(y=lambdas)) +
    geom_vline(xintercept = sz_onset, color = "blue") +
    labs(x = "Time", y = "Optimal lambda value") +
    ggtitle(titlepng)

  lambda_plot <- paste0(export,"/",subject_code,"_",sz_num,"_lambda_plot.png")
  ggsave(lambda_plot)
}
