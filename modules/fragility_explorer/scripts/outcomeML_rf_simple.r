source("./modules/fragility_explorer/Outcome_Classification_ML/read_data.R")
library(glue)
library(h2o)
h2o.shutdown()
h2o.init(max_mem_size = "10g")

## leave one out
predictions <- list()
dt_final$outcome <- as.factor(dt_final$outcome)

# patient <- unique(dt_weight$patient_code)[1]
for (patient in unique(dt_final$patient_code)) {
    message(patient)
    dt_train <- dt_final|>
        dplyr::filter(patient_code!=patient)
    dt_test <- dt_final|>
        dplyr::filter(patient_code==patient)

    dt_train_h2o <- as.h2o(dt_train)
    dt_test_h2o <- as.h2o(dt_test)

    columns <- colnames(dt_train_h2o[, -c(1:3)])


    rf_model <- h2o.randomForest(x = columns, y = "outcome",
                                 training_frame = dt_train_h2o, nfolds = 0,
                                 seed = 1234)

    rf_probabilities <- predict(rf_model, newdata = dt_test_h2o, type = "prob")
    preds <- as.data.frame(rf_probabilities)$p1
    dt_test[["preds"]] <- preds
    predictions[[patient]] <- dt_test|>
      select(patient_code, condition, outcome, preds)
}

predictions <- do.call(rbind, predictions)

# Compute statistics for ensemble predictions
predictions_voted <- predictions|>
  group_by(patient_code)|>
  summarise(preds = mean(preds),
            outcome = unique(outcome))


## Performance metrics
library(pROC)
roc_score=roc(predictions_voted$outcome, predictions_voted$preds)
coordinates <- coords(roc_score, x = "all", input = "threshold", ret = "all")
thresholds <- coordinates[,"threshold"]
ppv <- coordinates[,"ppv"]
npv <- coordinates[,"npv"]
tp <- coordinates[,"tp"]
fp <- coordinates[,"fp"]
fpr <- coordinates[,"fpr"]
tpr <- coordinates[,"tpr"]


## max number of true positives when false positives are 0
max(tp[which(fp<4)])

roc_obj <- roc(predictions_voted$outcome, predictions_voted$preds)

roc_data <- data.frame(
  TPR = roc_obj$sensitivities,
  FPR = roc_obj$specificities,
  Thresholds = roc_obj$thresholds
)

roc_plot <- ggplot(roc_data, aes(x = 1 - FPR, y = TPR)) +  # FPR is 1 - Specificity
  geom_line(color = "blue", size = 1.2) +
  geom_area(alpha = 0.2) +  # Add area under the curve
  labs(title = "ROC Curve",
       x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)") +
  theme_minimal()

auc_value <- auc(roc_obj)
roc_plot +
  annotate("text", x = 0.2, y = 0.1, label = paste("AUC =", round(auc_value, 3)), size = 5, color = "red")



# # Function to permeate importance of frequency bands
# perm_importance <- function(ensemble, data, base_models, target) {
#   # Calculate the original performance
#   original_perf <- h2o.performance(ensemble, newdata = data)
#   original_score <- original_perf@metrics$MSE
#
#   importances <- sapply(seq_along(base_models), function(i) {
#     # Exclude one model and retrain the ensemble
#     subset_models <- base_models[-i]
#     temp_ensemble <- h2o.stackedEnsemble(
#       x = names(data)[!names(data) %in% target],
#       y = target,
#       training_frame = data,
#       base_models = subset_models
#     )
#
#     # Calculate performance with one model omitted
#     temp_perf <- h2o.performance(temp_ensemble, newdata = data)
#     temp_score <- temp_perf@metrics$MSE
#
#     # Difference (higher indicates more importance)
#     original_score - temp_score
#   })
#
#   names(importances) <- sapply(base_models, function(model) model@model_id)
#   return(importances)
# }
