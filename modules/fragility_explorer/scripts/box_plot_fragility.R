library('pROC')
library('ggplot2')
library('tidyr')
library('dplyr')
library("pracma")

# define which patients to read from patient_data file
pts <- dipsaus::parse_svec("1-35,37-42,50-60,65-67,75-76,121,125,127,135,157-159")

# read patient_data file
patient_key <- read.csv("/Volumes/bigbrain/Multipatient/patient_data_all.csv")

# check which patients are being run
patient_key$subject[pts]

# choose folder with results for analysis
folder <- "/Volumes/bigbrain/Fragility2024/Results/250-125"

note <- "norank"
csv_files <- NULL

for (subject_code in unique(patient_key$subject_code[pts])) {
  directory <- paste0(folder,"/",subject_code,"/",note)
  if (file.exists(directory)) {
    csv_files <- append(csv_files,list.files(directory, pattern = paste0(note,"\\.csv$"), full.names = TRUE))
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
group_0 <- seizure_free$auc_other_fragility
group_1 <- seizure_free$auc_soz_fragility

# Function to calculate mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Summary statistics
summary_stats_0 <- c(
  N = length(group_0),
  Min = format(min(group_0), digits = 3),
  Max = format(max(group_0), digits = 3),
  Mean = format(mean(group_0), digits = 3),
  Range = format(diff(range(group_0)), digits = 3),
  Median = format(median(group_0), digits = 3),
  Mode = format(Mode(group_0), digits = 3),
  SD = format(sd(group_0), digits = 3)
)

summary_stats_1 <- c(
  N = length(group_1),
  Min = format(min(group_1), digits = 3),
  Max = format(max(group_1), digits = 3),
  Mean = format(mean(group_1), digits = 3),
  Range = format(diff(range(group_1)), digits = 3),
  Median = format(median(group_1), digits = 3),
  Mode = format(Mode(group_1), digits = 3),
  SD = format(sd(group_1), digits = 3)
)

# Statistical test (two-sample t-test)
t_test_result <- t.test(group_0, group_1)

# Print summary statistics for resect "0" group
print("Seizure Free: Summary statistics for resect Not EZ:")
print(summary_stats_0)

# Print summary statistics for resect "1" group
print("Seizure Free: Summary statistics for resect EZ:")
print(summary_stats_1)

# Print t-test result
print("Seizure Free: T-test result:")
print(t_test_result)

group_0 <- not_seizure_free$auc_other_fragility
group_1 <- not_seizure_free$auc_soz_fragility

# Summary statistics
summary_stats_0 <- c(
  N = length(group_0),
  Min = format(min(group_0), digits = 3),
  Max = format(max(group_0), digits = 3),
  Mean = format(mean(group_0), digits = 3),
  Range = format(diff(range(group_0)), digits = 3),
  Median = format(median(group_0), digits = 3),
  Mode = format(Mode(group_0), digits = 3),
  SD = format(sd(group_0), digits = 3)
)

summary_stats_1 <- c(
  N = length(group_1),
  Min = format(min(group_1), digits = 3),
  Max = format(max(group_1), digits = 3),
  Mean = format(mean(group_1), digits = 3),
  Range = format(diff(range(group_1)), digits = 3),
  Median = format(median(group_1), digits = 3),
  Mode = format(Mode(group_1), digits = 3),
  SD = format(sd(group_1), digits = 3)
)

# Statistical test (two-sample t-test)
t_test_result <- t.test(group_0, group_1)

# Print summary statistics for resect "0" group
print("Not Seizure Free: Summary statistics for resect Not EZ:")
print(summary_stats_0)

# Print summary statistics for resect "1" group
print("Not Seizure Free: Summary statistics for resect EZ:")
print(summary_stats_1)

# Print t-test result
print("Not Seizure Free: T-test result:")
print(t_test_result)

# Create box plot
# Define the directory where you want to save the plots
output_dir <- "/Users/ozhou/Library/CloudStorage/OneDrive-TexasA&MUniversity/Karas Lab"

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
ggsave(file.path(output_dir, paste0(note, "_seizure_free_auc_plot.png")), plot = plot, width = 85, height = 34, units = "cm")

# Perform more pairwise testing
group_0 <- seizure_free$auc_soz_fragility
group_1 <- not_seizure_free$auc_soz_fragility
t_test_result <- t.test(group_0, group_1)
# Print t-test result
print("Seizure Free EZ vs Not Seizure Free EZ: T-test result:")
print(t_test_result)


# Perform more pairwise testing
group_0 <- seizure_free$auc_other_fragility
group_1 <- not_seizure_free$auc_other_fragility
t_test_result <- t.test(group_0, group_1)
# Print t-test result
print("Seizure Free Not EZ vs Not Seizure Free Not EZ: T-test result:")
print(t_test_result)

p.adjust(c(5.279e-08, 0.03786, 0.8105, 0.00257), method = "bonferroni")

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
