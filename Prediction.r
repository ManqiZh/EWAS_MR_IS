# Set working directory
setwd("../UKB/IS-EWAS&MR/XGBoost")
workdir <- getwd()


# Load packages
pkgs <- c("dplyr", "lattice", "Matrix", "data.table", "survival")
inst <- lapply(pkgs, require, character.only = TRUE)


# 1 Prepare Data----
data <- readRDS("../UKB/IS-EWAS&MR/ISdata_final.RDS")

data_Initial <- data %>% dplyr::select(Status_IS, Time_IS, dplyr::all_of(variable_names))
data_XGB <- data %>% dplyr::select(Status_IS, Time_IS, dplyr::all_of(variable_names_XGB))

surv_cols <- c("Status_IS", "Time_IS")

feature_cols_initial <- colnames(data_Initial)[3:(ncol(data_Initial))]
design_matrix_initial <- model.matrix(
    ~.,
    data = data_Initial[, feature_cols_initial, with = FALSE]
)[, -1]

feature_cols_XGB <- colnames(data_XGB)[3:(ncol(data_XGB))]
design_matrix_XGB <- model.matrix(
    ~.,
    data = data_XGB[, feature_cols_XGB, with = FALSE]
)[, -1]


# 2 Predictions----

## 2.1 Define formulas for the Cox models----
FML_initial <- as.formula(
    paste(
        "survival::Surv(Time_IS, Status_IS == 1) ~",
        paste(feature_cols_initial, collapse = " + ")
    )
)
FML_XGB <- as.formula(
    paste(
        "survival::Surv(Time_IS, Status_IS == 1) ~",
        paste(feature_cols_XGB, collapse = " + ")
    )
)


## 2.2 Fit Cox models----
fit_initial <- survival::coxph(
    formula = FML_initial,
    data = data_Initial,
    ties = "efron",
    x = TRUE
)
fit_XGB <- survival::coxph(
    formula = FML_XGB,
    data = data_XGB,
    ties = "efron",
    x = TRUE
)


## 2.3 Obtain the beta estimates----
betaest_initial <- as.matrix(coefficients(fit_initial))
betaest_XGB <- as.matrix(coefficients(fit_XGB))


## 2.4 Generate the marker matrix----
score_initial <- as.matrix(design_matrix_initial) %*% betaest_initial
score_XGB <- as.matrix(design_matrix_XGB) %*% betaest_XGB


## 2.5 Generate Predicted Survival Probabilities----
predict_initial <- predict(fit_initial, type = "risk")
predict_XGB <- predict(fit_XGB, type = "risk")
