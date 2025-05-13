library(tidyverse)
library(glue)
library(ggcorrplot)
library(corrplot)
library(pracma)
library(stats)
library(h2o)
library(signal)

# define which patients to read from patient_data file
#pts <- dipsaus::parse_svec("1-35,37-42,50-60,65-67,75-76,121,125,127,135,157-159") # old patient data file
pts <- dipsaus::parse_svec("1-40,47-71") # Fragility and NIH pts only
pts <- dipsaus::parse_svec("1-40,47-92,95-99") # Fragility, NIH, and Karas pts
pts <- dipsaus::parse_svec("72-92") # Karas pts only

# read patient_data file
patient_key <-
    read.csv("/Volumes/bigbrain/Multipatient/patient_data_FINAL.csv", header = TRUE, stringsAsFactors = FALSE)|>
    mutate(
      condition = gsub("\\s.*", "", condition)
    )|>
    tibble()

# check which patients are being run
patient_key <- patient_key[pts,]
patient_key$subject_code

# choose folder with results for analysis
folder <- "/Volumes/bigbrain/Fragility_Results/250-125_lambda_n_100"

# specify norank
note <- "norank"
csv_files <- NULL

for (subject_code in unique(patient_key$subject_code)) {
  directory <- paste0(folder,"/",subject_code,"/",note)
  if (file.exists(directory)) {
    csv_files <- append(csv_files,list.files(directory, pattern = paste0(note,"\\.csv$"), full.names = TRUE))
  } else {
    print(paste0("patient ", subject_code, " error"))
  }
}
csv_files

quant_csv_files <- csv_files[grepl("quantile",csv_files)]
quant_csv_files

# Get list of CSV files in each directory
# csv_files_list <- lapply(directories, function(directory) {
#   list.files(directory, pattern = "\\.csv$", full.names = TRUE)
# })

# Initialize an empty list to store dataframes
all_dataframes <- list()
all_dataframes_annotation <- list()

idx <- 1
for (file in quant_csv_files) {
  # Extract patient code and condition from file name
  file_name <- basename(file)
  patient_code <- gsub("_.*", "", file_name)
  cnd <- gsub(".*_", "", gsub("_quantile_norank.csv", "", file_name))

  # convert "seizure" in file name to "sz" condition name
  cnd <- gsub("seizure","sz",cnd)

  ## only read the data if the patient code and condition exist in the patient key
  exist <- patient_key|>dplyr::filter(subject_code==patient_code, condition==cnd)
  if(nrow(exist)==0){
    next
  }


  # Read CSV file into a dataframe
  df <- tibble(read.csv(file))|>dplyr::select(-X)

  # Add patient code, condition, and frequency as columns
  all_dataframes_annotation[[idx]] <-
  tibble(
    patient_code = patient_code,
    condition = cnd
  )
  # Store dataframe in the list
  all_dataframes[[idx]] <- df

  idx <- idx + 1
}

all_dataframes_annotation <- do.call(rbind, all_dataframes_annotation)

non_exist_patient <-
  patient_key|>
  dplyr::select(subject_code, condition)|>
  unique()|>
  anti_join(all_dataframes_annotation,
            by = c("subject_code" = "patient_code",
                   "condition" = "condition"))

if(nrow(non_exist_patient)!=0){
  warning("Patient key do not match")
  patient_key_orig <- patient_key
  patient_key <- dplyr::filter(patient_key, !(subject_code %in% non_exist_patient$subject_code & condition %in% non_exist_patient$condition))
}



all_dataframes_annotation2 <-
  all_dataframes_annotation|>
  left_join(
    patient_key,
    by = c("patient_code" = "subject_code",
          "condition" = "condition")
          )

## times
times <- lapply(all_dataframes, function(x) as.numeric(gsub("^\\.","\\-",gsub("^X", "", colnames(x)))))

# identify patient with lowest sample rate
# sapply(times,length)
# minNumReadings_i <- which.min(sapply(times, length))
#minNumReadings <- all_dataframes[[minNumReadings_i]]

minNumReadings <- min(sapply(times,length))

# Initialize patient_data dataframe
patient_data <- list()
for(i in seq_len(nrow(all_dataframes_annotation2))){
      patient <- all_dataframes_annotation2$patient_code[i]
      condition <- all_dataframes_annotation2$condition[i]
      # truncate data to match lowest sample rate patient
      #print(i)
      #data_over_time_per_elec <- all_dataframes[[i]][, colnames(minNumReadings)]

      mult <- ncol(all_dataframes[[i]])/minNumReadings

      idx <- round(seq(1,ncol(all_dataframes[[i]]),mult))

      data_over_time_per_elec <- all_dataframes[[i]][, idx]

      colnames(data_over_time_per_elec) <- paste0("X", seq_len(ncol(data_over_time_per_elec)))

      # Create the new column with labels
      data_over_time_per_elec <- data_over_time_per_elec %>%
        mutate(label = c(paste0("resect_", seq(10, 100, by=10)), paste0("other_", seq(10, 100, by=10))))

      # Reshape data to a long format
      data_long <- data_over_time_per_elec %>%
        mutate(id = row_number()) %>%
        pivot_longer(cols = -c(id, label), names_to = "variable", values_to = "value")

      # Separate resect and other data
      resect_data <- data_long %>% dplyr::filter(grepl("resect", label))
      other_data <- data_long %>% dplyr::filter(grepl("other", label))

      # Ensure ordering for subtraction
      resect_data <- resect_data %>% arrange(id, label)
      other_data <- other_data %>% arrange(id, label)

      # changed from 0 to 1
      # 1 to 10 is good
      resect_data <- resect_data %>%
        mutate(value = scales::rescale(value, to = c(1, 10)))

      other_data <- other_data %>%
        mutate(value = scales::rescale(value, to = c(1, 10)))

      divide_data <- resect_data %>%
        mutate(value = value / other_data$value) %>%
        mutate(label = sub("divide_", "", label))

      divide_data <- resect_data %>%
        mutate(value = value / other_data$value) %>%
        mutate(label = sub("divide_", "", label)) %>%
        mutate(value = if_else(is.finite(value), value, 0))

      divide_data <- divide_data %>%
        mutate(across(where(is.numeric), ~ if_else(is.finite(.), ., 0)))

      # Remove old rows and combine the new results
      final_data <- data_long %>%
        dplyr::filter(!label %in% c(resect_data$label, other_data$label)) %>%
        bind_rows(divide_data) %>%
        arrange(id) %>%
        pivot_wider(names_from = variable, values_from = value, id_cols = c(id, label))

      # Optionally drop the 'id' column if it was just for tracking
      data_over_time_per_elec <- final_data %>%
        dplyr::select(-id)

      dt <- tibble(
        patient_code = patient,
        condition = condition,
      )|>
      bind_cols(data_over_time_per_elec)

      patient_data <- c(patient_data, list(dt))
}

patient_data <- do.call(rbind, patient_data)

Positive <- c("subpt01", "subpt2", "subpt3", "subpt8", "subpt11", "subpt13",
              "subpt17", "subpt15", "subpt16","subummc002","subummc005","subummc009",
              "subjh105", "subNIH032", "subNIH070", "subHUP070", "subHUP134", "subHUP163")
Negative <- c("subpt6", "subpt7", "subpt10", "subpt12", "subpt14", "subjh101", "subjh103",
              "subNIH016", "subNIH041", "subNIH046")

# Create outcome column
patient_data$outcome <- ifelse(patient_data$patient_code %in% Positive, 1,
                               ifelse(patient_data$patient_code %in% Negative, 0, NA))

## make the table wide format
#value_columns <- paste0("X", seq_len(length(minNumReadings)))
value_columns <- paste0("X", seq_len(minNumReadings))

long_data <- patient_data %>%
  pivot_longer(
    cols = starts_with("X"),
    names_to = "Measurement_ID",
    values_to = "Measurement_Value"
  )

final_data <- long_data %>%
  pivot_wider(
    id_cols = c("patient_code", "condition", "outcome"),
    names_from = "label",
    values_from = "Measurement_Value",
    values_fn = list(Measurement_Value = mean)
  )

# final_data_scaled <- final_data %>%
#   mutate(across(4:ncol(final_data), scale))
#
# remove_attributes <- function(x) {
#   attributes(x) <- NULL
#   x
# }
#
# # Assuming your data is in a tibble called final_data_scaled
# final_data_scaled <- final_data_scaled %>%
#   mutate(across(where(is.numeric), remove_attributes))
#
final_data<- final_data %>%
  dplyr::select_if(~ !any(is.na(.)))

final_data <- final_data %>%
  dplyr::select_if(~ !any(is.infinite(.)))

dt_final <- final_data

# Helper Functions
convert_range_to_vector <- function(range_string) {
  # Remove any whitespace
  range_string <- gsub("\\s", "", range_string)

  ## If empty, return empty vector
  if(range_string==""){
    return(c())
  }

  # Split the string by comma
  ranges <- strsplit(range_string, ",")[[1]]

  # Initialize an empty vector to store individual numbers
  numbers <- c()

  # Loop through each range
  for (range in ranges) {
    # Split the range by hyphen
    range_parts <- strsplit(range, "-")[[1]]

    # If only one part exists, it's a single number
    if (length(range_parts) == 1) {
      number <- as.integer(range_parts)
      if (is.na(number)) {
        warning("Invalid number detected: ", range)
        next  # Skip to the next range
      }
      numbers <- c(numbers, number)
    } else if (length(range_parts) == 2) {
      # Convert range parts to numbers
      start <- as.integer(range_parts[1])
      end <- as.integer(range_parts[2])

      # Check for NA or NaN values
      if (is.na(start) || is.na(end)) {
        warning("Invalid range detected: ", range)
        next  # Skip to the next range
      }

      # Add numbers within the range to the vector
      numbers <- c(numbers, start:end)
    } else {
      warning("Invalid range format: ", range)
    }
  }

  # if(length(numbers)==0){
  #   stop("stop")
  # }
  # if(length(numbers)==1&&is.na(numbers)){
  #   stop("stop")
  # }
  return(numbers)
}

nanpow2db <- function(y){
  # Power to dB conversion, setting negatives and zeros to NaN
  #
  # params:
  #         y: power --required
  #
  # returns:
  #         ydB: dB (with 0s and negatives set to 0)

  # Convert negative values to 0
  y <- pmax(y, 0)

  if(length(y)==1){
    if(y==0){
      return(NaN)
    } else {
      ydB <- 10*log10(y)
    }
  } else {
    y[y==0] <- NaN
    ydB <- 10*log10(y)
  }
  return(ydB)
}

