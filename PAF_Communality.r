# Set working directory
setwd("../UKB/IS-EWAS&MR")
workdir <- getwd()


# Load packages
pkgs <- c(
    "dplyr", "tidyr", "caret", "data.table",
    "survival", "stdReg", "forcats", "psych"
)
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
covariate.names <- paste(colnames(Covariate)[-1], collapse = "+")
Label_IS <- data %>% dplyr::select(n_eid, Status_IS, Time_IS)


## 1.2 Weighted scores----

### 1.2.1 Calculate weighted scores----
b2$Score_Weighted_Strat_Lifestyle <- cut(
    b2$Score_Weighted,
    breaks = quantile(
        x = b2$Score_Weighted,
        probs = seq(0, 1, length.out = 4)
    ),
    labels = NULL,
    include.lowest = TRUE
)
b2_Lifestyle <- b2


### 1.2.2 Combine weighted data----
df <- b2_Lifestyle %>%
    dplyr::left_join(b2_HMH, by = "n_eid") %>%
    dplyr::left_join(b2_Psychosocial, by = "n_eid") %>%
    dplyr::left_join(b2_SES, by = "n_eid") %>%
    dplyr::left_join(b2_Physic, by = "n_eid") %>%
    dplyr::left_join(Covariate, by = "n_eid") %>%
    dplyr::left_join(Label_IS, by = "n_eid")


### 1.2.3 Weighted data ready for PAF calculation----
b2 <- b2 %>%
    dplyr::select(
        n_eid,
        Score_Weighted_Strat_Lifestyle,
        Score_Weighted_Strat_HMH,
        Score_Weighted_Strat_Psychosocial,
        Score_Weighted_Strat_SES,
        Score_Weighted_Strat_Physic
    )
df <- dplyr::inner_join(
    dplyr::inner_join(Label_IS, Covariate, by = "n_eid"), b2,
    by = "n_eid"
)


## 1.3 Unweighted scores----

### 1.3.1 Unweighted data ready for PAF calculation----
b2 <- b2 %>%
    dplyr::select(
        n_eid,
        Score_Unweighted_Strat_Lifestyle,
        Score_Unweighted_Strat_HMH,
        Score_Unweighted_Strat_Psychosocial,
        Score_Unweighted_Strat_SES,
        Score_Unweighted_Strat_Physic
    )

df <- dplyr::inner_join(
    dplyr::inner_join(Label_IS, Covariate, by = "n_eid"), b2,
    by = "n_eid"
)


# 2 Convert the detrimental to beneficial and moderate ones----
df2 <- filter(df, df$Time_IS >= 2)

## 2.1 Convert to numeric levels----

### 2.1.1 Weighted scores----
unique(df2$Score_Weighted_Strat_Lifestyle)
df2$Score_Weighted_Strat_Lifestyle <- factor(
    df2$Score_Weighted_Strat_Lifestyle,
    levels = c("[2,3]", "(3,3.27]", "(3.27,5.37]"),
    labels = c(1, 1, 2)
)


### 2.1.2 Unweighted scores----
unique(df2$Score_Unweighted_Strat_Lifestyle)
df2$Score_Unweighted_Strat_Lifestyle <- factor(
    df2$Score_Unweighted_Strat_Lifestyle,
    levels = c("[13.2,21.2]", "(21.2,22.2]", "(22.2,30.2]"),
    labels = c(1, 1, 2)
)


## 2.2 Convert to integer----
for (i in c(8:12)) {
    df2[[i]] <- as.integer(df2[[i]])
}


# 3 Reduce the detrimental and moderate to beneficial ones----
df2 <- filter(df, df$Time_IS >= 2)

## 3.1 Convert to numeric levels----

### 3.1.1 Weighted scores----
unique(df2$Score_Weighted_Strat_Lifestyle)
df2$Score_Weighted_Strat_Lifestyle <- factor(
    df2$Score_Weighted_Strat_Lifestyle,
    levels = c("[2,3]", "(3,3.27]", "(3.27,5.37]"),
    labels = c(1, 2, 2)
)


### 3.1.2 Unweighted scores----
unique(df2$Score_Unweighted_Strat_Lifestyle)
df2$Score_Unweighted_Strat_Lifestyle <- factor(
    df2$Score_Unweighted_Strat_Lifestyle,
    levels = c("[13.2,21.2]", "(21.2,22.2]", "(22.2,30.2]"),
    labels = c(1, 2, 2)
)


## 3.2 Convert to integer----
for (i in c(8:12)) {
    df2[[i]] <- as.integer(df2[[i]])
}


# 4 PAF calculation----
# Function to calculate the attributable fraction (PAF)
calculate_PAF <- function(est) {
    p <- est[1]
    p0 <- est[2]
    PAF <- 1 - p0 / p
    return(PAF)
}

# Function to fit a univariate GLM model and calculate various statistics
fit_univariate_glm <-
    function(predictor) {
        formula <- as.formula(
            paste("Status_IS ~ ", predictor, "+", covariate.names)
        )
        glm_model <- glm(formula, data = df2, family = "binomial")

        # Standardized regression
        fit_std <- stdReg::stdGlm(
            fit = glm_model,
            data = df2,
            X = predictor,
            x = c(2, 1), # an vector containing the specific values of X at which to estimate the standardized mean. If X is binary (0/1) or a factor, then x defaults to all values of X . If X is numeric, then x defaults to the mean of X . If x is set to NA , then X is not altered. This produces an estimate of the marginal mean E(Y)=E\{E(Y|X,Z)\} .
            case.control = FALSE
        )
        # c(NA,0) is the value of the OR>1 variable, NA is the factual distribution, and 0 is the counterfactual distribution; The specific value of X to be written in x=c() is not necessarily 0

        # Calculate PAF and confidence intervals
        AF_estimate <- calculate_PAF(fit_std$est)
        conf_intervals <- confint(object = fit_std, fun = calculate_PAF, level = 0.95)

        # Calculate OR (95% CI) and its P-values
        glm_summary <- summary(glm_model)
        OR <- round(exp(coef(glm_model)[2]), 4)
        SE <- coef(glm_summary)[2, 2]
        CI_lower <- round(exp(coef(glm_model)[2] - 1.96 * SE), 4)
        CI_upper <- round(exp(coef(glm_model)[2] + 1.96 * SE), 4)
        P_values <- coef(glm_summary)[2, 4]

        # Create a data frame with the results
        results <- data.frame(
            "predictors" = predictor,
            "OR" = OR,
            "OR 95% CI_Lower" = CI_lower,
            "OR 95 %CI_Upper" = CI_upper,
            "PAF" = AF_estimate,
            "PAF 95% CI_Lower" = conf_intervals[1],
            "PAF 95 %CI_Upper" = conf_intervals[2],
            "P Value" = P_values
        )

        return(results)
    }

variable_names <- colnames(df2)[c(8:12)]

# Initialize and populate the Uni_glm list
# with the results of fit_univariate_glm
uni_glm <- lapply(variable_names, fit_univariate_glm)
results <- do.call(rbind, uni_glm)


# 5 Communality & Weighted PAF Calculation----

## 5.1 Generate a correlation matrix----
# The tetrachoric correlation is the inferred Pearson Correlation from a two x two table with the assumption of bivariate normality.
# The polychoric correlation generalizes this to the n x m table.
# Particularly important when doing Item Response Theory or converting comorbidity statistics using normal theory to correlations.
# Input may be a 2 x 2 table of cell frequencies, a vector of cell frequencies, or a data.frame or matrix of dichotomous data (for tetrachoric).

correlation <- psych::tetrachoric(df2[, 8:12], na.rm = TRUE)

cor <- correlation$rho # Correlation coefficient matrix
psych::corPlot(
    cor,
    numbers = TRUE, # Display the numeric value of the correlations
    colors = TRUE, # colors=FALSE will use a grey scale
    n = 51, # The number of levels of shading to use
    main = NULL, # The title of the plot
    zlim = c(-1, 1), # The range of values to color.
    show.legend = TRUE,
    labels = NULL, # if NULL, use column and row names
    n.legend = 10, # How many categories should be labelled in the legend?
    keep.par = TRUE, # restore the graphic parameters when exiting
    select = NULL, # Select the subset of variables to plot
    pval = NULL, # scale the numbers by their pvals, categorizing them based upon the values of cuts
    digits = 2,
    trailing = TRUE, # Show trailing zeros
    cuts = c(.001, .01), # Scale the numbers by the categories defined by pval < cuts
    scale = TRUE, # Should the size of the numbers be scaled by the significance level?
    upper = TRUE,
    diag = TRUE,
    symmetric = TRUE,
    stars = TRUE,
    adjust = "holm",
    xaxis = 1,
    xlas = 0, # Orientation of the x axis labels
    ylas = 2,
    ysrt = 0, #  Rotation of y labels in degrees
    xsrt = 0,
    alpha = .75, # The degree of transparency (0 = completely, 1= not)
    min.length = NULL, # If not NULL, then the maximum number of characters to use in row/column labels
    sort = FALSE, # If true, then sort the variables using the iclust algorithm
    n.obs = NULL #  If you want to show "stars" for symmetric input matrices (i.e. correlations), specify the number of observations
)
dev.off()


## 5.2 Eigenvalue decomposition----
# According to the Kaiser criterion, retain eigenvalues>1
ev <- eigen(cor)
u <- as.matrix(ev$vectors)[, which(ev$values > 1)]


## 5.3 Communality Calculation----
u$communality <- rowSums(u[, c("V1", "V2")]^2)
u$communality <- rowSums(u^2) # When there is only one principal component

u$predictors <- colnames(df2)[8:12]


## 5.4 Weighted PAF Calculation----
communality_data <- u %>% select(communality, predictors)
results <- results %>% dplyr::left_join(communality_data, by = "predictors")

overall_PAF <- 1 - cumprod(1 - (1 - results$Communality) * results$PAF) %>% tail(1)
results <- results %>% dplyr::mutate(weighted_PAF = (PAF / sum(PAF)) * overall_PAF)
