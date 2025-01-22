# Set working directory
setwd("../UKB/IS-EWAS&MR")
workdir <- getwd()


# Load packages
pkgs <- c("dplyr", "tidyr", "readxl", "survival", "grprep")
inst <- lapply(pkgs, require, character.only = TRUE)


# Load data
data <- readRDS("../UKB/IS-EWAS&MR/ISdata_final.RDS")


# 1 EWAS analysis----

## 1.1 Define functions for EWAS analysis----
check_data_and_variable <- function(data, x) {
    if (!is.data.frame(data)) { # nolint
        stop("Provided 'data' is not a data frame.") # nolint
    } # nolint
    if (!x %in% names(data)) {
        stop(paste0("Variable '", x, "' not found in 'data'.")) # nolint
    } # nolint
}

uni_cox_model <- function(data, x) {
    check_data_and_variable(data, x)

    # Fit a Cox proportional hazards model
    FML <- as.formula(paste0(
        "survival::Surv(Time_IS, Status_IS == 1) ~ ",
        x,
        " + n_21022_0_0 + n_31_0_0 + n_54_0_0 + Ethnic"
    ))
    coxphModel <- survival::coxph(FML, data = data)
    SumModel <- summary(coxphModel)
    PH_test <- survival::cox.zph(coxphModel)

    # Extract the results
    results <- list()
    variable_PH_p_value <- PH_test$table[1, "p"]

    if (is.factor(data[[x]]) && length(levels(data[[x]])) > 1) {
        levels_x <- levels(data[[x]])

        # Extract results for each level of the factor variable
        for (i in 2:length(levels_x)) {
            level <- levels_x[i]
            coef_name <- paste0(x, level)
            results[[coef_name]] <- list(
                Coef = SumModel$coefficients[coef_name, "coef"],
                HR = exp(SumModel$coefficients[coef_name, "coef"]),
                CI_Lower = exp(confint(coxphModel, level = 0.95)[coef_name, 1]),
                CI_Upper = exp(confint(coxphModel, level = 0.95)[coef_name, 2]),
                P_Value = SumModel$coefficients[coef_name, "Pr(>|z|)"],
                PH_P_Value = variable_PH_p_value
            )
        }
    } else {
        results[[x]] <- list(
            Coef = SumModel$coefficients[1],
            HR = exp(SumModel$coefficients[1]),
            CI_Lower = exp(confint(coxphModel, level = 0.95)[1, 1]),
            CI_Upper = exp(confint(coxphModel, level = 0.95)[1, 2]),
            P_Value = SumModel$coefficients[1, "Pr(>|z|)"],
            PH_P_Value = variable_PH_p_value
        )
    }

    results_df <- dplyr::bind_rows(
        lapply(results, dplyr::bind_rows),
        .id = "Variable_Level"
    )
    results_df$FDR_P_Value <- p.adjust(results_df$P_Value, method = "fdr")
    results_df$Holm_P_Value <- p.adjust(results_df$P_Value, method = "holm")
    results_df$Hochberg_P_Value <- p.adjust(
        results_df$P_Value,
        method = "hochberg"
    )
    results_df$Types_of_Variables <- ifelse(
        is.factor(data[[x]]), "Factor", "Continuous"
    )
    return(results_df)
}

process_variables <- function(data, variable_names) {
    results_df <- lapply(variable_names, function(x) {
        tryCatch(
            {
                uni_cox_model(data, x)
            },
            error = function(e) {
                message("Error with variable ", x, ": ", e$message)
                NULL
            }
        )
    })
    do.call(rbind, results_df)
}


## 1.2 Formal analysis----
all_results <- process_variables(data, variable_names)


## 1.3 Multivariable analysis----
# Fit a Cox proportional hazards model
formula <- as.formula(paste(
    "survival::Surv(Time_IS, Status_IS) ~",
    paste(variable_names, collapse = " + ")
))
cox_model <- survival::coxph(formula, data = data)

# Extract the results
summary_cox <- summary(cox_model)
coef <- summary_cox$coefficients
conf_int <- summary_cox$conf.int


# 2 Collinearity detection----

## Data reduction
data_selected <- data %>% dplyr::select(dplyr::all_of(var_names))

## Transform the data into numeric format,
## as R::cor() requires 'x' to be numeric
data_selected_numeric <- data_selected
for (col in names(data_selected_numeric)) {
    if (!is.numeric(data_selected_numeric[[col]])) {
        data_selected_numeric[[col]] <- as.numeric(
            as.character(data_selected_numeric[[col]])
        )
    }
}

## Detect collinearity
collinearity <- caret::findCorrelation(
    x = cor(data_selected_numeric),
    cutoff = 0.9,
    verbose = TRUE,
    names = TRUE
)


# 3 Group LASSO-related Algorithms----

## 3.1 Data Preparation----
# Load the data
y <- cbind(time = data$Time_IS, status = as.character(data$Status_IS))
y <- as.data.frame(y)
y$status <- as.factor(y$status)
y <- as.matrix(y)
X <- data %>% dplyr::select(dplyr::all_of(var_names))

# Identify the type of each variable
continuous_vars <- sapply(X, is.numeric)
continuous_var_names <- names(X)[continuous_vars]

binary_vars <- sapply(X, function(x) length(unique(x)) == 2)
binary_var_names <- names(X)[binary_vars]

categorical_vars <- !continuous_vars & !binary_vars
categorical_var_names <- names(X)[categorical_vars]

# Select continuous and binary variables
x_continuous_binary <- X[, continuous_vars | binary_vars]

# One-hot encode the categorical variables
x_categorical <- X[, categorical_vars]
x_one_hot <- dummy::dummy(x_categorical, int = FALSE)

# Combine the continuous and binary variables with the one-hot encoded variables
x_final <- dplyr::bind_cols(x_continuous_binary, as.data.frame(X_one_hot))

# Generate group names
continuous_binary_names <- colnames(x_continuous_binary)
one_hot_names <- colnames(X_one_hot)
modified_names <- sub("_[0-9]+$", "", one_hot_names) # Remove the suffixes
group <- c(continuous_binary_names, modified_names)

# Preparation for the grpreg function
group <- factor(group)
# R::grpreg requires 'x' to be a matrix and numeric
x_final_matrix <- as.matrix(x_final)
x_final_matrix <- apply(x_final_matrix, 2, as.numeric)


## 3.2 grSCAD----
set.seed(123)
cvfit_grSCAD_gamma <- grpreg::cv.grpsurv(
    X = x_final_matrix,
    y = y,
    group = group,
    penalty = "grSCAD",
    nlambda = 100,
    log.lambda = TRUE,
    alpha = 1,
    eps = 1e-4,
    max.iter = 10000,
    tau = 1 / 3,
    nfolds = 10,
    seed = 123,
    se = "quick",
    warn = TRUE,
    returnY = TRUE,
    trace = TRUE
)

summary(cvfit_grSCAD_gamma)
cvfit_grSCAD_gamma$lambda.min
coefficients_grSCAD_gamma <- coef(cvfit_grSCAD_gamma)
selected_variables_grSCAD_gamma <- which(coefficients_grSCAD_gamma != 0)


## 3.3 grLasso----
set.seed(123)
cvfit_grLasso <- grpreg::cv.grpsurv(
    X = x_final_matrix,
    y = y,
    group = group,
    penalty = "grLasso",
    nlambda = 100,
    log.lambda = TRUE,
    alpha = 1,
    eps = 1e-4,
    max.iter = 10000,
    tau = 1 / 3,
    nfolds = 10,
    seed = 123,
    se = "quick",
    warn = TRUE,
    returnY = TRUE,
    trace = TRUE
)

summary(cvfit_grLasso)
cvfit_grLasso$lambda.min
coefficients_grLasso <- coef(cvfit_grLasso)
selected_variables_grLasso <- which(coefficients_grLasso != 0)
