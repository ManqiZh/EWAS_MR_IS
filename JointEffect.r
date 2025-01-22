# Set working directory
setwd("../UKB/IS-EWAS&MR")
workdir <- getwd()


# Load packages
pkgs <- c("dplyr", "tidyr", "caret", "data.table", "survival")
inst <- lapply(pkgs, require, character.only = TRUE)


# 1 Prepare data----

## 1.1 Load data----
data <- readRDS("../UKB/IS-EWAS&MR/ISdata_final.RDS")
data <- data %>% dplyr::select(
    n_eid, Status_IS, Time_IS, dplyr::all_of(variable_names)
)

Covariate <- data %>% dplyr::select(
    n_eid, n_21022_0_0, n_31_0_0, n_54_0_0, Ethnic
)
covariate.names <- paste(
    colnames(Covariate)[-1],
    collapse = "+"
) # Exclude n_eid
Label_IS <- data %>% dplyr::select(n_eid, Status_IS, Time_IS)


## 1.2  Domain----

### 1.2.1 Lifestyle----
Domain <- data %>% dplyr::select(dplyr::all_of(FieldID))
EWAS_Results <- data.table::as.data.table(EWAS_Results)


## 1.3 Combine data----
df <- dplyr::inner_join(
    dplyr::inner_join(
        Label_IS, Covariate,
        by = "n_eid"
    ), Domain,
    by = "n_eid"
)


## 1.4 Convert factor to dummy variables----
dmy <- caret::dummyVars(
    ~.,
    data = df[, 8:ncol(df)],
    sep = ".",
    levelsOnly = FALSE,
    fullRank = TRUE
) # Create a converter dmy for the factor variables

b <- as.data.frame(
    predict(dmy, df[, 8:ncol(df)])
) # Apply the converter to the data

b$n_eid <- df$n_eid
b <- dplyr::select(b, n_eid, dplyr::everything())

for (i in 2:ncol(b)) {
    b[[i]] <- as.factor(b[[i]])
}


## 1.5 Adjust HR<1 values to make them risk factors----
protectEXPOS <- EWAS_Results$X[
    EWAS_Results$HR < 1 & EWAS_Results$Significance == 1
]
protectEXPOS_NoSig <- EWAS_Results$X[
    EWAS_Results$HR < 1 & EWAS_Results$Significance == 0
]

# For factor variables
cols_to_recode <- which(colnames(b) %in% protectEXPOS)
cols_to_recode_NoSig <- which(colnames(b) %in% protectEXPOS_NoSig)

for (i in cols_to_recode) {
    b[[i]] <- dplyr::recode_factor(b[[i]], `1` = "0", `0` = "1")
}
for (i in cols_to_recode_NoSig) {
    b[[i]] <- dplyr::recode_factor(b[[i]], `1` = "0", `0` = "1")
}

# For numeric variables
b <- df %>% dplyr::select(-c(2:7))
cols_to_recode <- which(colnames(b) %in% protectEXPOS)
cols_to_recode_NoSig <- which(colnames(b) %in% protectEXPOS_NoSig)

for (i in cols_to_recode) {
    b[[i]] <- 1 - b[[i]]
}
for (i in cols_to_recode_NoSig) {
    b[[i]] <- 1 - b[[i]]
}


# 2 Construct the Cox model and generate the beta value, on which the weighted scores depend----
df <- dplyr::inner_join(
    dplyr::inner_join(b, Covariate, by = "n_eid"), Label_IS,
    by = "n_eid"
)

variable_names <- paste(colnames(df)[2:(ncol(b) + 1)], collapse = "+")
FML_Weighted <- as.formula(
    paste0(
        "survival::Surv(Time_IS, Status_IS == 1) ~ ",
        variable_names,
        " + ", covariate.names
    )
)
fit_Weighted <- survival::coxph(FML_Weighted, data = df, ties = "efron")
model_summary_Weighted <- summary(fit_Weighted)
coef_Weighted <- as.data.frame(model_summary_Weighted$coefficients[, 1])


# 3 Unweighted scores----

## 3.1 Calculate unweighted scores----

# For factor variables
for (i in 2:ncol(b)) {
    b[[i]] <- as.integer(as.character(b[[i]]))
} # Convert to integer to calculate the sum

# For numeric variables
for (i in 2:ncol(b)) {
    b[[i]] <- as.numeric(b[[i]])
} # Convert to num to calculate the sum

# Continue
b_modified <- b
for (i in cols_to_recode_NoSig) {
    b_modified[[i]] <- b[[i]] * 0.01
}
b$Score_Unweighted <- rowSums(b_modified[, 2:ncol(b_modified)], na.rm = TRUE)


## 3.2 Stratify the unweighted scores based on quantiles of the scores----
b$Score_Unweighted_Strat <- cut(
    x = b$Score_Unweighted,
    breaks = quantile(
        x = b$Score_Unweighted,
        probs = seq(0, 1, length.out = 4)
    ),
    labels = NULL,
    include.lowest = TRUE # include the lowest value of the interval
)
table(b$Score_Unweighted_Strat)
prop.table(table(b$Score_Unweighted_Strat))


## 3.3 Independent domain Cox model for unweighted scores----
df <- dplyr::inner_join(
    dplyr::inner_join(b, Covariate, by = "n_eid"), Label_IS,
    by = "n_eid"
)

FML_Unweighted <- as.formula(
    paste0(
        "survival::Surv(Time_IS, Status_IS == 1) ~ ",
        "Score_Unweighted_Strat",
        " + ", covariate.names
    )
)
fit_Unweighted <- survival::coxph(
    FML_Unweighted,
    data = df,
    ties = "efron"
)
summary(fit_Unweighted)


# 4 Weighted scores----

## 4.1 Calculate weighted scores----
b2 <- b[, 2:(ncol(b) - 2)] # Exclude the unweighted score
b2 <- b[, 2:ncol(b)]

b2 <- as.matrix(b2)
b2 <- apply(b2, 2, as.numeric)
b2 <- b2 %*% diag(beta) # Multiply the beta values to the variables
b2 <- as.data.frame(b2)
b2$n_eid <- b$n_eid
b2 <- dplyr::select(b2, n_eid, dplyr::everything())

b2$Score_Weighted <- rowSums(b2[2:ncol(b2)], na.rm = TRUE)


## 4.2 Stratify the weighted scores based on quantiles of the scores----
b2$Score_Weighted_Strat <- cut(
    b2$Score_Weighted,
    breaks = quantile(
        x = b2$Score_Weighted,
        probs = seq(0, 1, length.out = 4)
    ),
    labels = NULL,
    include.lowest = TRUE
)
table(b2$Score_Weighted_Strat)
prop.table(table(b2$Score_Weighted_Strat))


## 4.3 Independent domain Cox model for weighted scores----
df <- dplyr::inner_join(
    dplyr::inner_join(b2, Covariate, by = "n_eid"),
    Label_IS,
    by = "n_eid"
)

FML_Weighted_Domain <- as.formula(
    paste0(
        "survival::Surv(Time_IS, Status_IS == 1) ~ ",
        "Score_Weighted_Strat",
        " + ", covariate.names
    )
) # Regular

fit_Weighted_Domain <- survival::coxph(
    FML_Weighted_Domain,
    data = df,
    ties = "efron"
)
summary(fit_Weighted_Domain)


## 4.4 Adjusted All-Domains Model for weighted scores----

# Combine data
data_origin <- dplyr::inner_join(Label_IS, Covariate, by = "n_eid")
data_origin$Score_Weighted_Strat_Lifestyle <- b2$Score_Weighted_Strat

# Develop the formula
df <- data_origin
fit_Weighted <- survival::coxph(
    survival::Surv(Time_IS, Status_IS == 1) ~ .,
    data = df,
    ties = "efron"
)
sum_fit_Weighted <- summary(fit_Weighted)


# 5 Unweighted Score-Adjusted All-Domains Model----

## 5.1 Calculate unweighted scores----
b$Score_Unweighted_Strat_Lifestyle <- cut(
    x = b$Score_Unweighted,
    breaks = quantile(
        x = b$Score_Unweighted,
        probs = seq(0, 1, length.out = 4)
    ),
    labels = NULL,
    include.lowest = TRUE # include the lowest value of the interval
)
b_Lifestyle <- b


## 5.2 Cox model----
# Combine data
df <- b_Lifestyle %>%
    dplyr::left_join(b_HMH, by = "n_eid") %>%
    dplyr::left_join(b_Psychosocial, by = "n_eid") %>%
    dplyr::left_join(b_SES, by = "n_eid") %>%
    dplyr::left_join(b_Physic, by = "n_eid") %>%
    dplyr::left_join(Covariate, by = "n_eid") %>%
    dplyr::left_join(Label_IS, by = "n_eid")

# Develop the formula
FML_Unweighted <- as.formula(
    paste0(
        "survival::Surv(Time_IS, Status_IS == 1) ~ ",
        "Score_Unweighted_Strat_Lifestyle + Score_Unweighted_Strat_HMH + Score_Unweighted_Strat_Psychosocial + Score_Unweighted_Strat_SES + Score_Unweighted_Strat_Physic",
        " + ", covariate.names
    )
)
fit_Unweighted <- survival::coxph(
    FML_Unweighted,
    data = df,
    ties = "efron"
)
summary(fit_Unweighted)
