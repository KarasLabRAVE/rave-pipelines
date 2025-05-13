library('pROC')
library('ggplot2')
library('tidyr')
library('dplyr')
library("pracma")

# define which patients to read from patient_data file
pts <- dipsaus::parse_svec("1-40,47-71") # Fragility and NIH pts only
#pts <- dipsaus::parse_svec("1-40,47-92,95-99") # Fragility, NIH, and Karas pts
#pts <- dipsaus::parse_svec("1-12,16-21,26-34,39,47-53,55,61-92,95-99") # Fragility, NIH, and Karas pts with avg R2 > 0.8
#pts <- dipsaus::parse_svec("1-12,16-21,26-34,39,47-53,55,61-71") # Fragility and NIH pts with avg R2 >0.8 only
#pts <- dipsaus::parse_svec("72-92,95-99") # Karas pts only

# read patient_data file
patient_key <- read.csv("/Volumes/bigbrain/Multipatient/patient_data_FINAL.csv")

# check which patients are being run
patient_key$subject[pts]

# choose folder with results for analysis
folder <- "/Volumes/bigbrain/EZFragility_Results/250-125_lambda_n_100"

note <- "norank"
file_identifiers <- paste0(patient_key$subject[pts],"_", patient_key$condition[pts])

csv_files <- NULL

for (subject_code in unique(patient_key$subject_code[pts])) {
  directory <- paste0(folder,"/",subject_code,"/",note)
  if (file.exists(directory)) {
    for (filename in file_identifiers) {
      csv_files <- append(csv_files,list.files(directory, pattern = paste0(filename,"_meandata_",note,"\\.csv$"), full.names = TRUE))
    }
  } else {
    print(paste0("patient ", subject_code, " error"))
  }
}
csv_files

mean_csv_files <- csv_files[grepl("meandata",csv_files)]
mean_csv_files

# Initialize an empty list to store dataframes
dataframes <- list()

time_window <- c(0,20)

# Loop through each CSV file
for (file in mean_csv_files) {
  # Extract patient code and condition from file name
  file_name <- basename(file)
  patient_code <- gsub("_.*", "", file_name)
  condition <- strsplit(file_name, "_")[[1]][2]
  condition <- gsub("seizure","sz",condition)

  # Read CSV file into a dataframe
  df <- read.csv(file)
  df$patient_code <- patient_code
  df$condition <- condition

  df <- subset(df,between(df$time,time_window[1],time_window[2]))

  # Store dataframe in the list
  dataframes[[file_name]] <- df
}

# Set time step from multitaper
time_step <- 0.125

patient_data <- data.frame(patient_code = character(),
                           condition = character(),
                           time = numeric(),
                           other_power = numeric(),
                           soz_power = numeric())

# Loop through each unique subject in patient key
for (patient in unique(patient_key$subject_code[pts])) {
  temp_patient <- patient_key[patient_key$subject_code == patient, ]
  for (condition in temp_patient$condition) {
    temp_condition <- temp_patient[temp_patient$condition == condition, ]
    subject <- temp_condition$subject_code
    condition <- gsub("\\s.*", "", condition)

    # Find dataframe corresponding to the subject and condition
    data_over_time_per_elec <- NULL
    for (df_name in names(dataframes)) {
      df <- dataframes[[df_name]]
      if (subject == df$patient_code[1] && condition == df$condition[1]) {
        data_over_time_per_elec <- df
        break
      }
    }

    data_over_time_per_elec <- data_over_time_per_elec[c("patient_code", "condition", "time", "mean_f_sozc", "mean_f_soz")]
    patient_data <- rbind(patient_data, data_over_time_per_elec)
  }
}

colnames(patient_data) <- c("patient_code", "condition", "time", "other_fragility", "soz_fragility")


#Positive <- c("subpt01", "subpt2", "subpt3", "subpt8", "subpt11", "subpt13", "subpt15", "subpt16", "subpt17", "subummc002", "subummc005", "subummc009", "subjh105")
#Negative <- c("subpt6", "subpt7", "subpt10", "subpt12", "subpt14", "subjh101", "subjh103")

Positive <- unique(patient_key$subject_code[which(patient_key$outcome == "S")])
Negative <- unique(patient_key$subject_code[which(patient_key$outcome == "F")])

# Create outcome column
patient_data$outcome <- ifelse(patient_data$patient_code %in% Positive, 1,
                               ifelse(patient_data$patient_code %in% Negative, 0, NA))


data <- patient_data %>%
  group_by(patient_code, condition) %>%
  arrange(time) %>%
  summarise(
    auc_other_fragility = trapz(time, other_fragility),
    auc_soz_fragility = trapz(time, soz_fragility),
    outcome = first(outcome),
    .groups = 'drop'
  )


# separate seizure free and not seizure free
seizure_free <- data[data$outcome==1, ]
seizure_free$auc_other_fragility <- as.numeric(seizure_free$auc_other_fragility)
seizure_free$auc_soz_fragility <- as.numeric(seizure_free$auc_soz_fragility)

not_seizure_free <- data[data$outcome==0, ]
not_seizure_free$auc_other_fragility <- as.numeric(not_seizure_free$auc_other_fragility)
not_seizure_free$auc_soz_fragility <- as.numeric(not_seizure_free$auc_soz_fragility)


# Calculate the minimum and maximum values for y-axis
# min_value <- min(seizure_free$power, not_seizure_free$power)
min_value <- 0
max_value <- max(max( (mean(seizure_free$auc_other_fragility) + quantile(seizure_free$auc_other_fragility, 0.75)), ((mean(seizure_free$auc_soz_fragility) + quantile(seizure_free$auc_soz_fragility, 0.75)))),
                 max( (mean(not_seizure_free$auc_other_fragility) + quantile(not_seizure_free$auc_other_fragility, 0.75)), (mean(not_seizure_free$auc_soz_fragility) + quantile(not_seizure_free$auc_soz_fragility, 0.75))))
# max_value <- max(seizure_free$power, not_seizure_free$power)

# Filter data for resect 0 and 1 groups
SzF_soz <- seizure_free$auc_soz_fragility
SzF_other <- seizure_free$auc_other_fragility

# Function to calculate mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Summary statistics

summary_stats_SzF_soz <- c(
  N = length(SzF_soz),
  Min = format(min(SzF_soz), digits = 3),
  Max = format(max(SzF_soz), digits = 3),
  Mean = format(mean(SzF_soz), digits = 3),
  Range = format(diff(range(SzF_soz)), digits = 3),
  Median = format(median(SzF_soz), digits = 3),
  Mode = format(Mode(SzF_soz), digits = 3),
  SD = format(sd(SzF_soz), digits = 3)
)

summary_stats_SzF_other <- c(
  N = length(SzF_other),
  Min = format(min(SzF_other), digits = 3),
  Max = format(max(SzF_other), digits = 3),
  Mean = format(mean(SzF_other), digits = 3),
  Range = format(diff(range(SzF_other)), digits = 3),
  Median = format(median(SzF_other), digits = 3),
  Mode = format(Mode(SzF_other), digits = 3),
  SD = format(sd(SzF_other), digits = 3)
)

# Statistical test (two-sample t-test)
t_test_result1 <- t.test(SzF_soz, SzF_other, conf.level = 0.99)

# Print summary statistics for resect EZ group
print("Seizure Free: Summary statistics for resect EZ:")
print(summary_stats_SzF_soz)

# Print summary statistics for resect not EZ group
print("Seizure Free: Summary statistics for resect Not EZ:")
print(summary_stats_SzF_other)

# Print t-test result
print("Seizure Free: T-test result:")
print(t_test_result1)

NSzF_soz <- not_seizure_free$auc_soz_fragility
NSzF_other <- not_seizure_free$auc_other_fragility

# Summary statistics

summary_stats_NSzF_soz <- c(
  N = length(NSzF_soz),
  Min = format(min(NSzF_soz), digits = 3),
  Max = format(max(NSzF_soz), digits = 3),
  Mean = format(mean(NSzF_soz), digits = 3),
  Range = format(diff(range(NSzF_soz)), digits = 3),
  Median = format(median(NSzF_soz), digits = 3),
  Mode = format(Mode(NSzF_soz), digits = 3),
  SD = format(sd(NSzF_soz), digits = 3)
)

summary_stats_NSzF_other <- c(
  N = length(NSzF_other),
  Min = format(min(NSzF_other), digits = 3),
  Max = format(max(NSzF_other), digits = 3),
  Mean = format(mean(NSzF_other), digits = 3),
  Range = format(diff(range(NSzF_other)), digits = 3),
  Median = format(median(NSzF_other), digits = 3),
  Mode = format(Mode(NSzF_other), digits = 3),
  SD = format(sd(NSzF_other), digits = 3)
)

# Statistical test (two-sample t-test)
t_test_result2 <- t.test(NSzF_soz, NSzF_other, conf.level = 0.99)

# Print summary statistics for resect "0" group
print("Not Seizure Free: Summary statistics for resect Not EZ:")
print(summary_stats_NSzF_other)

# Print summary statistics for resect "1" group
print("Not Seizure Free: Summary statistics for resect EZ:")
print(summary_stats_NSzF_soz)

# Print t-test result
print("Not Seizure Free: T-test result:")
print(t_test_result2)

# Create box plot
# Define the directory where you want to save the plots
output_dir <- "/Volumes/bigbrain/EZFragility_Results/Data_Analysis"

# Add a new 'group' column
seizure_free$group <- "Seizure Free"
not_seizure_free$group <- "Not Seizure Free"

# Combine the two datasets
combined_data <- bind_rows(seizure_free, not_seizure_free)

# Reshape the combined data to long format
long_data <- pivot_longer(combined_data,
                          cols = c("auc_other_fragility", "auc_soz_fragility"),
                          names_to = "Label",
                          values_to = "auc_values")

long_data$Label <- recode(long_data$Label,
                          "auc_other_fragility" = "Not EZ",
                          "auc_soz_fragility" = "EZ")

long_data <- long_data %>%
  mutate(group = factor(group, levels = c("Seizure Free", "Not Seizure Free")))

plot <- ggplot(long_data, aes(x = group, y = auc_values, fill = Label)) +
  geom_boxplot() +
  labs(title = paste0("Comparison of AUC: "),
       x = "Patient Outcome",
       y = paste0("AUC ")) +
  scale_fill_manual(values = c("Not EZ" = "lightblue", "EZ" = "salmon")) +  # Updated color mapping with new labels
  theme_minimal() +
  theme(
    text = element_text(size = 8),  # This changes global text size
    axis.title = element_text(size = 8),  # This changes axis titles size
    plot.title = element_text(size = 8, face = "bold")) +  # This changes plot title size and makes it bold
  #scale_y_continuous(labels = scales::label_number(scale = 1, accuracy = 0.1), breaks = waiver()) +
  coord_cartesian(ylim = c(min_value, max_value))

# Display the plot
print(plot)

# Save the plot as an image
# ggsave(file.path(output_dir, paste0(note, "_seizure_free_auc_plot.png")), plot = plot, width = 85, height = 34, units = "cm")

# Perform more pairwise testing
t_test_result3 <- t.test(SzF_soz, NSzF_soz)
# Print t-test result
print("Seizure Free EZ vs Not Seizure Free EZ: T-test result:")
print(t_test_result3)


# Perform more pairwise testing
t_test_result4 <- t.test(SzF_other, NSzF_other)
# Print t-test result
print("Seizure Free Not EZ vs Not Seizure Free Not EZ: T-test result:")
print(t_test_result4)

p.adjust(c(t_test_result1$p.value, t_test_result2$p.value, t_test_result3$p.value, t_test_result4$p.value), method = "bonferroni")

summary_stats_SzF_soz
summary_stats_SzF_other
summary_stats_NSzF_soz
summary_stats_NSzF_other

# Helper Functions
convert_range_to_vector <- function(range_string) {
  # Remove any whitespace
  range_string <- gsub("\\s", "", range_string)

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
