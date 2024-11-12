# define which patients to read from patient_data file
pts <- dipsaus::parse_svec("1-35,37-42,50-60,65-67,75-76,121,125,127,135,157-159")

# read patient_data file
patient_key <- read.csv("/Volumes/bigbrain/Multipatient/patient_data_all.csv")

# check which patients are being run
patient_key$subject[pts]

# choose folder with results for analysis
folder <- "/Volumes/bigbrain/Fragility2024/Results/250-125"

csv_files <- NULL

for (subject_code in unique(patient_key$subject_code[pts])) {
  directory <- paste0(folder,"/",subject_code)
  if (file.exists(directory)) {
    csv_files <- append(csv_files,list.files(directory, pattern = paste0("R2","\\.csv$"), full.names = TRUE))
  } else {
    print(paste0("patient ", subject_code, " error"))
  }
}
csv_files

# Initialize an empty list to store dataframes
dataframes <- list()

# Loop through each CSV file
for (file in csv_files) {
  # Extract patient code and condition from file name
  file_name <- basename(file)
  patient_code <- gsub("_.*", "", file_name)
  condition <- strsplit(file_name, "_")[[1]][2]
  #condition <- gsub("seizure","sz",condition)

  # Read CSV file into a dataframe
  df <- read.csv(file)
  df$patient_code <- patient_code
  df$condition <- condition

  # remove lambda row and index row
  df <- df[1:(dim(df)[1]-1),2:dim(df)[2]]
  #df <- df[1:(dim(df)[1]-1),]

  # Store dataframe in the list
  dataframes[[file_name]] <- df
}

patient_data <- NULL

# Loop through each unique subject in patient key
for (patient in unique(patient_key$subject_code[pts])) {
  temp_patient <- patient_key[patient_key$subject_code == patient, ]
  for (condition in temp_patient$condition) {
    temp_condition <- temp_patient[temp_patient$condition == condition, ]
    subject <- temp_condition$subject_code
    condition <- gsub("\\s.*", "", condition)

    # Find dataframe corresponding to the subject and condition
    R2 <- NULL
    for (df_name in names(dataframes)) {
      df <- dataframes[[df_name]]
      if (subject == df$patient_code[1] && condition == df$condition[1]) {
        df <- as.matrix(dplyr::select(df,!(patient_code | condition)))
        R2 <- mean(apply(df,2,quantile,0.1))
        break
      }
    }

    #R2 <- data_over_time_per_elec[c("patient_code", "condition", "time", "R2")]
    patient_data <- append(patient_data, R2)
  }
}

names(patient_data) <- paste0(patient_key$subject[pts]," ",patient_key$condition[pts])
patient_data

mean(patient_data)
sd(patient_data)

