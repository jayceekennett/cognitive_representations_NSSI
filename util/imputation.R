# KENNETT (2026) SENSITIVITY ANALYSIS ----
# MULTIPLE IMPUTATION OF TWO INCOMPLETE B-COPE ITEMS USING MICE 

### load packages ----
library(mice)
library(lattice)

# PREP DATA ----

### declare variables for imputation----
imputed_variables <- c("Q12_27", "Q12_28")

# eight B-IPQ predictors used in the primary regression analyses
independent_vars <- c("consequences",
                      "timeline",
                      "pers_control",
                      "treat_control",
                      "identity",
                      "concern",
                      "coherence",
                      "emot_rep")

### set outcome vars ----
dependent_vars <- c("SD",
                    "AC",
                    "DEN",
                    "SU",
                    "UES",
                    "UIS",
                    "BD",
                    "VEN",
                    "PF",
                    "PLA",
                    "HUM",
                    "ACC",
                    "REL",
                    "SB",
                    "PSYCH_OPEN",
                    "HELP_SEEK",
                    "STIG_INDIFF")

cope_subscales <- c("SD",
                    "AC",
                    "DEN",
                    "SU",
                    "UES",
                    "UIS",
                    "BD",
                    "VEN",
                    "PF",
                    "PLA",
                    "HUM",
                    "ACC",
                    "REL",
                    "SB")

### set sig threshold ----
p_value_threshold <- 0.05
 
### load original data ----
df_imp <- read.csv("data/s1_data.csv", stringsAsFactors = FALSE)

# REPLICATE PRIMARY ANALYSIS ----

### prep data for replicaton of primary analysis ----
prepare_analysis_data <- function(df) {
  
  # match the original analysis:
  # convert columns 1-8 to numeric.
  df[1:8] <- lapply(df[1:8], as.numeric)
  return(df)
}

primary_df <- prepare_analysis_data(df_imp)

### replicate primary analysis correlation ----
spearman_test <- function(data, x, y) {
  
  complete_rows <- complete.cases(data[, c(x, y)])
  
  # extract only complete pairs 
  x_values <- data[[x]][complete_rows]
  y_values <- data[[y]][complete_rows]
 
  if (length(x_values) < 3) {
    return(NULL)
  }
  # run Spearman's rank correlation
  suppressWarnings(cor.test(x_values, y_values, method = "spearman", exact = FALSE))
}

### get primary correlation matrix ----
primary_correlation_matrix <- cor(
  primary_df[, c(independent_vars, dependent_vars)],
  method = "spearman",
  use = "pairwise.complete.obs")

### get MLR model formulas ----
get_primary_formulas <- function(data, dependent_vars,independent_vars,
    threshold = 0.05) {
  formulas <- list()
  
  for (dv in dependent_vars) {
    significant_ivs <- character(0)
    
    for (iv in independent_vars) {
      result <- spearman_test(data = data, x = dv, y = iv)
      
      # if the test could be calculated and p < .05 retain the predictor 
      if (!is.null(result) && !is.na(result$p.value) && result$p.value < threshold) {
        
        significant_ivs <- c(significant_ivs, iv)}}
    
    # only create a model if at least one predictor was selected 
    if (length(significant_ivs) > 0) {
      
      formulas[[dv]] <- reformulate(
        termlabels = significant_ivs,
        response = dv)}}
  return(formulas)
}

primary_formulas <- get_primary_formulas(
  data = primary_df,
  dependent_vars = dependent_vars,
  independent_vars = independent_vars,
  threshold = p_value_threshold)

### declare MLR formulas being used ----
cat("\nPRIMARY MODEL SPECIFICATIONS\n")

for (model_name in names(primary_formulas)) {
  cat(model_name,
      ": ", paste(deparse(primary_formulas[[model_name]]), collapse = ""),
      "\n", sep = "")
}

### replicate primary regression analysis ----
primary_models <- lapply(primary_formulas, function(formula) {
    lm(formula, data = primary_df)})


# PREPARE DATA FOR MICE -----

### create imputation df ----
dat <- df_imp[, 1:60]

### conv to numeric ----
dat[] <- lapply(dat, as.numeric)

### set missing variable to ordered ----
dat$Q12_27 <- ordered(dat$Q12_27)
dat$Q12_28 <- ordered(dat$Q12_28)

### inspect missing data ----
cat("\nMISSING DATA PATTERN\n")
print(md.pattern(dat, plot = FALSE))

### check missingness specifically for the two target variables ----
missingness_table <- table(
  Q12_27_missing = is.na(dat$Q12_27),
  Q12_28_missing = is.na(dat$Q12_28))

### specify imputation method ----
meth <- make.method(dat)
meth[] <- ""

### set polr for missing items ----
meth["Q12_27"] <- "polr"
meth["Q12_28"] <- "polr"

### create predictor matrix ----
pred <- make.predictorMatrix(dat)
pred[,] <- 0

### specify predictors for each item ----
# Q12_27 predictors 
pred["Q12_27", c("Q12_2","Q12_10","Q12_14","Q12_21","Q12_23")] <- 1

# Q12_28 predictors 
pred["Q12_28", c("Q12_4","Q12_11","Q12_26")] <- 1

# check no variable is predicting itself 
diag(pred) <- 0

### check no predictors have missing data ----
other_missing <- setdiff(names(which(colSums(is.na(dat)) > 0)), imputed_variables)

cat("\nOTHER VARIABLES WITH MISSING DATA\n")
if (length(other_missing) == 0) {
  cat("None.\n")
} else {
  print(other_missing)
}

invalid_predictors_Q12_27 <- intersect(other_missing,
                                       names(which(pred["Q12_27", ] == 1)))

invalid_predictors_Q12_28 <- intersect(other_missing,
                                       names(which(pred["Q12_28", ] == 1)))


# RUN MICE ----

### set seed ----
set.seed(7224)

imp <- mice(
  dat,
  method = meth,
  predictorMatrix = pred,
  m = 40,
  maxit = 20,
  seed = 7224,
  printFlag = TRUE)

### check events log ----
cat("\nMICE LOGGED EVENTS\n")

if (is.null(imp$loggedEvents)) {
  cat("No logged events.\n")
} else {
  print(imp$loggedEvents)
}

# RUN DIAGNOSTICS ----

#dir.create("imputation_diagnostics", recursive = TRUE,showWarnings = FALSE)

### convergence plot ----
#png("imputation_diagnostics/convergence_plot.png", width = 1800,height = 1000,res = 300)
#plot(imp)
#dev.off()

### stripplot ----
#png("imputation_diagnostics/strip_plot_Q12_27.png", width = 1800, height = 1000, res = 300)
#stripplot(imp,Q12_27 ~ .imp,pch = 20)
#dev.off()

#png("imputation_diagnostics/strip_plot_Q12_28.png", width = 1800, height = 1000, res = 300)
#stripplot(imp,Q12_28 ~ .imp,pch = 20)
#dev.off()

### extract imputed values ----
Q12_27_imputed_values <- imp$imp$Q12_27
Q12_28_imputed_values <- imp$imp$Q12_28

#write.csv(Q12_27_imputed_values,"imputation_diagnostics/Q12_27_imputed_values.csv",
  #row.names = TRUE)

#write.csv(Q12_28_imputed_values,"imputation_diagnostics/Q12_28_imputed_values.csv",
  #row.names = TRUE)

### compare observed/imputed distribution ----
compare_imputed_distribution <- function(imp, original_data, variable) {
  
  # define response categories 
  categories <- levels(original_data[[variable]])
  
 # observed distribution 
  observed <- original_data[[variable]][!is.na(original_data[[variable]])]
  observed_counts <- table(factor(observed, levels = categories))
  observed_prop <- prop.table(observed_counts)

  # imputed only distributions 
  missing_rows <- which(
    is.na(original_data[[variable]]))

  imputed_props <- sapply(seq_len(imp$m),function(i) {
      completed <- complete(imp, i)
      imputed_values <- completed[missing_rows, variable]
      counts <- table(factor(imputed_values, levels = categories))
      as.numeric(prop.table(counts))})
 
  result <- data.frame(category = categories,
                        observed_proportion = as.numeric(observed_prop),
                        mean_imputed_proportion = rowMeans(imputed_props),
                        min_imputed_proportion = apply(imputed_props, 1, min),
                        max_imputed_proportion = apply(imputed_props, 1, max))
  return(result)
  }

distribution_Q12_27 <- compare_imputed_distribution(
  imp = imp,
  original_data = dat,
  variable = "Q12_27")

distribution_Q12_28 <- compare_imputed_distribution(
  imp = imp,
  original_data = dat,
  variable = "Q12_28")

# print comparisons
cat("\nQ12_27 DISTRIBUTION COMPARISON\n")
print(distribution_Q12_27)

cat("\nQ12_28 DISTRIBUTION COMPARISON\n")
print(distribution_Q12_28)

### save to csv ----
#write.csv(distribution_Q12_27, "imputation_diagnostics/distribution_Q12_27.csv",
  #row.names = FALSE)

#write.csv(distribution_Q12_28, "imputation_diagnostics/distribution_Q12_28.csv",
  #row.names = FALSE)

### create imputed datasets -----
completed_datasets <- lapply(seq_len(imp$m), function(i) {
    completed <- complete(imp, i)
    completed$Q12_27 <- as.numeric(as.character(completed$Q12_27))
    completed$Q12_28 <- as.numeric(as.character(completed$Q12_28))
    
    # replace incomplete raw items with imputation 
    df_complete <- df_imp
    df_complete$Q12_27 <- completed$Q12_27
    df_complete$Q12_28 <- completed$Q12_28

    # reconstruct subscales 
    df_complete$REL <- rowSums(df_complete[, c("Q12_27", "Q12_25")], na.rm = FALSE)
    df_complete$SD <- rowSums(df_complete[, c("Q12_28","Q12_1")], na.rm = FALSE)

    # match primary-analysis numeric conversion 
    df_complete[1:8] <- lapply(df_complete[1:8], as.numeric)
    
    return(df_complete)
  })

### write to csv ----
#dir.create("data/imputed_data", recursive = TRUE, showWarnings = FALSE)
#for (i in seq_along(completed_datasets)) {
#  write.csv(completed_datasets[[i]],
#    paste0("data/imputed_data/imputed_dataset_", i, ".csv"), row.names = FALSE)
#}

# RUN CORRELATION ANALYSIS----
correlation_vars <- append(independent_vars, dependent_vars)
correlation_matrices <- lapply(completed_datasets,function(df){
    cor(df[, correlation_vars], method = "spearman", use = "pairwise.complete.obs")
  })

### stack matrices ----
cor_array <- simplify2array(correlation_matrices)

### calculate summary statistics ----
mean_cor <- apply(cor_array, c(1,2), mean)
sd_cor <- apply(cor_array, c(1,2), sd)
min_cor <- apply(cor_array, c(1,2), min)
max_cor <- apply(cor_array, c(1,2), max)

# 2.5th percentile
lower_cor <- apply(cor_array, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)

# 97.5th percentile
upper_cor <- apply(cor_array, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)

# get absolute mean difference
difference <- mean_cor  - primary_correlation_matrix
abs_difference <- abs(difference)
summary(abs_difference)
which(abs_difference == max(abs_difference), arr.ind = TRUE)
original_cor <- primary_correlation_matrix

## present differences ----
cor_table <- expand.grid(
  var1 = rownames(original_cor),
  var2 = colnames(original_cor),
  stringsAsFactors = FALSE)

cor_table$complete_case <- as.vector(original_cor)
cor_table$mean_imputed <- as.vector(mean_cor)
cor_table$sd_imputed <- as.vector(sd_cor)
cor_table$lower_2.5 <- as.vector(lower_cor)
cor_table$upper_97.5 <- as.vector(upper_cor)
cor_table$absolute_difference <- as.vector(abs_difference)

keep <- upper.tri(original_cor)
cor_table <- data.frame(
  var1 = rownames(original_cor)[row(original_cor)[keep]],
  var2 = colnames(original_cor)[col(original_cor)[keep]],
  complete_case = original_cor[keep],
  mean_imputed = mean_cor[keep],
  sd_imputed = sd_cor[keep],
  lower_2.5 = lower_cor[keep],
  upper_97.5 = upper_cor[keep])

cor_table$difference <- (cor_table$mean_imputed -cor_table$complete_case)
cor_table$absolute_difference <- abs(cor_table$difference)

# sort changes 
cor_table <- cor_table[order(cor_table$absolute_difference, decreasing = TRUE),]

# round numeric columns
numeric_cols <- sapply(cor_table, is.numeric)
cor_table[numeric_cols] <- round(cor_table[numeric_cols], 3)
rownames(cor_table) <- NULL

### write to csv ----
#write.csv(cor_table, "imputation_results/cor_imputation_results.csv", row.names = FALSE)

models_by_analysis <- lapply(names(primary_formulas),
  function(model_name) {
    
    formula <- primary_formulas[[model_name]]
    models <- lapply(completed_datasets, function(data) {
        lm(formula,data = data)})
    
    return(models)
  })

# RUN REGRESSIONS ----

# restore model names
names(models_by_analysis) <- names(primary_formulas)

### check model formulas are identical ----
for (model_name in names(models_by_analysis)) {
  
  formulas <- vapply(models_by_analysis[[model_name]],
    function(model) {
      
      paste(deparse(formula(model)), collapse = "")},
    FUN.VALUE = character(1))
  
  if (length(unique(formulas)) != 1) {
    stop(paste0("Model formulas differ across imputations for ", model_name,
                ". Pooling has been stopped."))}}

### convert objects to mira ----
pooled_models <- lapply(models_by_analysis, function(models) {
    mira_object <- mice::as.mira(models)
    mice::pool(mira_object)})

pooled_results <- lapply(pooled_models, summary, conf.int = TRUE)

### combine pooled results ----
pooled_results_df <- do.call(
  rbind, lapply(names(pooled_models), function(model_name) {
      result <- summary(pooled_models[[model_name]], conf.int = TRUE)
      result$model <- model_name
      
      return(result)
    }))

rownames(pooled_results_df) <- NULL

# put model name first 
pooled_results_df <- pooled_results_df[
  c("model", setdiff(names(pooled_results_df), "model"))]

### save pooled results ----
#dir.create("imputation_results", recursive = TRUE,showWarnings = FALSE)

#write.csv(pooled_results_df, "imputation_results/pooled_imputation_results.csv",
          #row.names = FALSE)

