# Set working directory
setwd("../UKB/IS-EWAS&MR/XGBoost")
workdir <- getwd()


# Load packages
pkgs <- c(
    "dplyr", "lattice", "Matrix", "data.table", "survival", "mlsurvlrnrs",
    "xgboost", "Ckmeans.1d.dp", "mlexperiments", "ParBayesianOptimization",
)
inst <- lapply(pkgs, require, character.only = TRUE)


# 1 Preprocessing----

## 1.1 Import and Prepare Data----
data <- readRDS("../UKB/IS-EWAS&MR/ISdata_final.RDS")
data_211 <- data %>% dplyr::select(
    Status_IS, Time_IS, dplyr::all_of(variable_names)
)

surv_cols <- c("Status_IS", "Time_IS")
feature_cols <- colnames(data_211)[3:(ncol(data_211))]


## 1.2 General Configurations----
seed <- 123

# Set the number of cores to use intelligently
if (isTRUE(as.logical(Sys.getenv("_R_CHECK_LIMIT_CORES_")))) {
    # on cran
    ncores <- 2L
} else {
    ncores <- ifelse(
        test = parallel::detectCores() > 4,
        yes = 4L,
        no = ifelse(
            test = parallel::detectCores() < 2L,
            yes = 1L,
            no = parallel::detectCores()
        )
    )
}

options("mlexperiments.bayesian.max_init" = 10L)
options("mlexperiments.optim.xgb.nrounds" = 1000L)
options("mlexperiments.optim.xgb.early_stopping_rounds" = 100L)


## 1.3 Generate Training- and Test Data----

# Stratify the data into 4 clusters
split_vector <- splitTools::multi_strata(
    df = data_211[
        ,
        .SD,
        .SDcols = surv_cols
    ],
    strategy = "kmeans",
    k = 4L
)

# Split the data into training and test sets
data_split <- splitTools::partition(
    y = split_vector,
    p = c(train = 0.7, test = 0.3),
    type = "stratified",
    n_bins = 10L,
    split_into_list = TRUE,
    use_names = TRUE,
    shuffle = FALSE,
    seed = seed
)

# Prepare training data
train_x <- model.matrix(
    ~ -1 + ., # -1 to remove intercept
    data_211[
        data_split$train,
        .SD,
        .SDcols = feature_cols
    ]
)
train_y <- survival::Surv(
    event = (data_211[data_split$train, get("Status_IS")] |>
        as.character() |>
        as.integer()),
    time = data_211[data_split$train, get("Time_IS")],
    type = "right"
)

split_vector_train <- splitTools::multi_strata(
    df = data_211[
        data_split$train,
        .SD,
        .SDcols = surv_cols
    ],
    strategy = "kmeans",
    k = 4L
)

# Prepare test data
test_x <- model.matrix(
    ~ -1 + ., # -1 to remove intercept
    data_211[
        data_split$test,
        .SD,
        .SDcols = feature_cols
    ]
)
test_y <- survival::Surv(
    event = (data_211[data_split$test, get("Status_IS")] |>
        as.character() |>
        as.integer()),
    time = data_211[data_split$test, get("Time_IS")],
    type = "right"
)


## 1.4 Generate Training Data Folds----
fold_list <- splitTools::create_folds(
    y = split_vector_train,
    k = 5L,
    type = "stratified",
    n_bins = 10L,
    m_rep = 1L,
    use_names = TRUE,
    invert = FALSE,
    shuffle = FALSE,
    seed = seed
)


# 2 Experiments----

## 2.1 Prepare Experiments----
# Required learner arguments, not optimized
learner_args <- list(
    objective = "survival:cox",
    eval_metric = "cox-nloglik"
)

# Set arguments for predict function and performance metric,
# required for mlexperiments::MLCrossValidation and mlexperiments::MLNestedCV
predict_args <- NULL
performance_metric <- c_index
performance_metric_args <- NULL
return_models <- TRUE

# Required for grid search and initialization of bayesian optimization
parameter_grid <- expand.grid(
    subsample = seq(0.6, 1, .1),
    colsample_bytree = seq(0.6, 1, .1),
    min_child_weight = seq(1, 5, 1),
    learning_rate = seq(0.01, 0.3, 0.01),
    max_depth = seq(4, 8, 1)
)

# Reduce to a maximum of 10 rows to speed up the process
if (nrow(parameter_grid) > 10) {
    set.seed(123)
    sample_rows <- sample(seq_len(nrow(parameter_grid)), 10, FALSE)
    parameter_grid <- kdry::mlh_subset(parameter_grid, sample_rows)
}


## 2.2 Nested Cross Validation: Inner Grid Search----
validator <- mlexperiments::MLNestedCV$new(
    learner = LearnerSurvXgboostCox$new(
        metric_optimization_higher_better = FALSE
    ),
    strategy = "grid",
    fold_list = fold_list,
    k_tuning = 5L,
    ncores = ncores,
    seed = seed
)

validator$parameter_grid <- parameter_grid
validator$learner_args <- learner_args
validator$split_type <- "stratified"
validator$split_vector <- split_vector_train

validator$predict_args <- predict_args
validator$performance_metric <- performance_metric
validator$performance_metric_args <- performance_metric_args
validator$return_models <- return_models

validator$set_data(
    x = train_x,
    y = train_y
)

validator_results <- validator$execute()


## 2.3 Holdout Test Dataset Performance----
### Predict Outcome in Holdout Test Dataset
preds_xgboost <- mlexperiments::predictions(
    object = validator,
    newdata = test_x
)

### Evaluate Performance on Holdout Test Dataset
perf_xgboost <- mlexperiments::performance(
    object = validator,
    prediction_results = preds_xgboost,
    y_ground_truth = test_y
)


# 3 XGBoost Survival Model----

## 3.1 XGBoost requires the time variable to be negative for censored observations----
train_y_matrix <- ifelse(
    data_211[data_split$train, get("Status_IS")] == 1,
    data_211[data_split$train, get("Time_IS")],
    -data_211[data_split$train, get("Time_IS")]
)
train_y_matrix <- matrix(train_y_matrix, ncol = 1)

test_y_matrix <- ifelse(
    data_211[data_split$test, get("Status_IS")] == 1,
    data_211[data_split$test, get("Time_IS")],
    -data_211[data_split$test, get("Time_IS")]
)
test_y_matrix <- matrix(test_y_matrix, ncol = 1)


## 3.2 Group the data and the label in a list----
dtrain <- xgboost::xgb.DMatrix(
    data = train_x,
    label = train_y_matrix,
    missing = NA,
    silent = FALSE,
    nthread = 4
)
dtest <- xgboost::xgb.DMatrix(
    data = test_x,
    label = test_y_matrix,
    missing = NA,
    silent = FALSE,
    nthread = 4
)


## 3.3 Measure learning progress----
set.seed(123)
bst <- xgboost::xgboost(
    data = train_x,
    label = train_y_matrix,
    missing = NA,
    weight = NULL,
    params = list(
        # General Parameters
        booster = "gbtree",
        nthread = 4, # Default to maximum number of threads available
        ## Parameters for Tree Booster
        eta = 0.01, # Range: [0,1]. Default: 0.3. Control the learning rate. Scale the contribution of each tree by a factor of 0 < eta < 1 when it is added to the current approximation. Used to prevent overfitting by making the boosting process more conservative. Lower value for eta implies larger value for nrounds : low eta value means model more robust to overfitting but slower to compute.
        gamma = 0, # Range: [0,∞]. Default: 0. Minimum loss reduction required to make a further partition on a leaf node of the tree. The larger gamma is, the more conservative the algorithm will be.
        max.depth = 6, # Range: [0,∞]. Default: 6. Maximum depth of a tree. Increasing this value will make the model more complex and more likely to overfit. Beware that XGBoost aggressively consumes memory when training a deep tree. The larger, the more conservative the algorithm will be.
        min_child_weight = 1, # Range: [0,∞]. Default: 1. Minimum sum of instance weight (hessian) needed in a child. If the tree partition step results in a leaf node with the sum of instance weight less than min_child_weight, then the building process will give up further partitioning. The larger min_child_weight is, the more conservative the algorithm will be.
        subsample = 0.8, # Range: (0,1]. Default: 1. subsample ratio of the training instance. Setting it to 0.5 means that xgboost randomly collected half of the data instances to grow trees and this will prevent overfitting. It makes computation shorter (because less data to analyse). It is advised to use this parameter with eta and increase nrounds.
        sampling_method = "uniform", # Default: "uniform". "uniform" or "gradient_based". The method to use to sample the training instances. "uniform" means that each training instance has an equal probability of being selected. Typically set subsample >= 0.5 for good results. "gradient_based" means that the selection probability for each training instance is proportional to the regularized absolute value of gradients. subsample may be set to as low as 0.1 without loss of model accuracy. Note that this sampling method is only supported when tree_method is set to hist and the device is cuda; other tree methods only support uniform sampling.
        colsample_bytree = 0.7, # Range: (0,1]. Default: 1. This is a family of parameters for subsampling of columns (colsample_bytree, colsample_bylevel, colsample_bynode).
        lambda = 1, # Range: [0,∞]. Default: 1. L2 regularization term on weights. Increasing this value will make model more conservative.
        alpha = 0, # Range: [0,∞]. Default: 0. L1 regularization term on weights. Increasing this value will make model more conservative.
        tree_method = "auto", # Default: "auto". The tree construction algorithm used in XGBoost. Choices: auto, exact, approx, hist.
        refresh_leaf = 1, # Default: 1. This is a parameter of the refresh updater. When this flag is 1, tree leafs as well as tree nodes’ stats are updated. When it is 0, only node stats are updated.
        process_type = "default", # Choices: default, update. A type of boosting process to run.

        # Learning Task Parameters
        objective = "survival:cox",
        eval_metric = "cox-nloglik"
    ),
    nrounds = 10000,
    verbose = 2,
    print_every_n = 1L,
    early_stopping_rounds = 100,
    maximize = TRUE,
    save_period = 100,
    save_name = "xgboost.model.%04d",
    xgb_model = NULL
)


# 4 View feature importance/influence from the learnt model----
importance_matrix <- xgboost::xgb.importance(
    feature_names = NULL,
    model = bst,
    trees = NULL,
    data = NULL,
    label = NULL,
    target = NULL
)

xgboost::xgb.ggplot.importance(
    importance_matrix = importance_matrix,
    top_n = 100,
    measure = "Gain",
    rel_to_first = FALSE, # Whether importance values should be represented as relative to the highest ranked feature
    n_clusters = c(1:10)
)
