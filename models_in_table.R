library(lavaan)
library(tidyverse)
library(MASS)

# ==============================================================================
# Reproducible lavaan code for all models in Tables 2-5 of
# "Model Identification Made Clear: From Conceptual Foundations to CFA
# Identification Rules"
#
# IMPORTANT NOTE FOR THE REVIEWER:
# Across all models, the factor scale is set by fixing the loading of the
# first indicator of each factor to 1 (e.g., F =~ 1*RSES_1 + ...).
# This is one of the two standard ways to identify the factor scale
# (the other is fixing the factor variance to 1 with F ~~ 1*F).
# It is NOT correct to apply BOTH constraints at the same time; doing so
# over-identifies the scale and can make an unidentified model appear
# identified. All models below use ONLY the marker-loading constraint.
#
# If the reviewer failed to fix the factor scale (or, equivalently, fixed
# both a loading and the factor variance), even identified models can appear
# unidentified and the warning "Could not compute standard errors!" will appear
# for the wrong reason. The models in Table 2 are unidentified because of
# correlated-error confounding (Theorem 1 / Theorem 2), not because of scale
# indeterminacy. Fixing the loading of x1 only sets the scale; it does not
# resolve the correlated-error problem in Table 2 model (a) or (c).
# ==============================================================================

# ------------------------------------------------------------------------------
# Load the RSES dataset (used for Tables 2, 3, and 5)
# ------------------------------------------------------------------------------
rses_data <- read.csv("dataset RSES Zhang et al 2023.csv") %>%
  filter(RSES_condition == "Original")

# ------------------------------------------------------------------------------
# Helper function: fit a model, print key info, and collect results
# ------------------------------------------------------------------------------
fit_and_summarize <- function(model_string, data, label, data_source) {
  cat("\n\n============================================================\n")
  cat(label, "\n")
  cat("Data source:", data_source, "\n")
  cat("============================================================\n")

  captured_warnings <- character()
  add_warning <- function(w) {
    captured_warnings <<- c(captured_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }

  fit <- tryCatch(
    withCallingHandlers(cfa(model_string, data = data, std.lv = FALSE),
                        warning = add_warning),
    error = function(e) {
      cat("\nERROR: model could not be estimated.\n")
      cat(conditionMessage(e), "\n")
      return(NULL)
    }
  )

  if (is.null(fit)) {
    return(data.frame(
      Model = label,
      Data = data_source,
      npar = NA,
      df = NA,
      Warning = "Estimation error",
      stringsAsFactors = FALSE
    ))
  }

  df_val <- fit@test[[1]]$df
  npar_val <- fit@loglik$npar

  cat("\nNumber of free parameters:", npar_val, "\n")
  cat("Degrees of freedom:", df_val, "\n")

  cat("\n--- Model summary ---\n")
  withCallingHandlers(
    print(summary(fit, fit.measures = TRUE, standardized = TRUE)),
    warning = add_warning
  )

  # Check for identification problems. lavaan may signal this either as a
  # warning or by returning NA standard errors / a non-positive-definite
  # information matrix.
  id_warning <- "None"
  if (any(grepl("information matrix could not be inverted", captured_warnings))) {
    id_warning <- "Information matrix not inverted (identification warning)"
  } else {
    vcov_mat <- tryCatch(lavInspect(fit, "vcov"), error = function(e) NULL)
    if (is.null(vcov_mat) || any(is.na(vcov_mat)) ||
        any(eigen(vcov_mat, only.values = TRUE)$values <= 1e-10)) {
      id_warning <- "Singular or non-positive-definite information matrix (identification warning)"
    }
  }

  if (id_warning != "None") {
    cat("\n*** IDENTIFICATION WARNING ***\n")
    cat(id_warning, "\n")
  } else if (!fit@optim$converged) {
    id_warning <- "Model did not converge"
    cat("\n*** CONVERGENCE WARNING ***\n")
    cat(id_warning, "\n")
  }

  return(data.frame(
    Model = label,
    Data = data_source,
    npar = npar_val,
    df = df_val,
    Warning = id_warning,
    stringsAsFactors = FALSE
  ))
}

# Store results for the final summary table
results_list <- list()


# ==============================================================================
# Table 2: Examples of One-Factor Model With Correlated Errors
# (RSES data; items 1-6 used for models c-e)
# ==============================================================================

# -----------------------------------------------------------------------------
# Table 2 Model (a): Unidentified (df = 1)
# Single factor F on items 1-5. Item 4's residual is correlated with the
# residual of every other item, violating Theorem 1 (At-Least-One-
# Uncorrelated-Error Rule).
# -----------------------------------------------------------------------------
model_2a <- '
F =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5
RSES_4 ~~ RSES_1 + RSES_2 + RSES_3 + RSES_5
'
results_list[["2a"]] <- fit_and_summarize(
  model_2a, rses_data,
  "Table 2 Model (a): Unidentified",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 2 Model (b): Just-identified (df = 0)
# Satisfies Theorem 3 (One-Factor-Plus-Correlated-Errors Rule):
#   (1) RSES_1, RSES_2, RSES_3 have mutually uncorrelated residuals;
#   (2) every other item has at least one residual correlation fixed to 0
#       with an item whose loading is identified.
# -----------------------------------------------------------------------------
model_2b <- '
F =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5
RSES_1 ~~ RSES_4 + RSES_5
RSES_3 ~~ RSES_4 + RSES_5
RSES_4 ~~ RSES_5
'
results_list[["2b"]] <- fit_and_summarize(
  model_2b, rses_data,
  "Table 2 Model (b): Just-identified",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 2 Model (c): Unidentified CTCU model with two method effects (df = 3)
# Positively worded items: RSES_1, RSES_3, RSES_4
# Negatively worded items: RSES_2, RSES_5, RSES_6
# All residuals within each wording group are correlated; no cross-group
# correlated residuals. This is globally unidentified by Theorem 2.
# -----------------------------------------------------------------------------
model_2c <- '
F =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6
# positively worded items
RSES_1 ~~ RSES_3 + RSES_4
RSES_3 ~~ RSES_4
# negatively worded items
RSES_2 ~~ RSES_5 + RSES_6
RSES_5 ~~ RSES_6
'
results_list[["2c"]] <- fit_and_summarize(
  model_2c, rses_data,
  "Table 2 Model (c): Unidentified CTCU (two method effects)",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 2 Model (d): Over-identified CTCU model (df = 4)
# Same as Model (c) except the correlated residual between RSES_1 and RSES_3
# is removed. The model now satisfies Theorem 3.
# -----------------------------------------------------------------------------
model_2d <- '
F =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6
# positively worded items (one pair removed)
RSES_1 ~~ RSES_4
RSES_3 ~~ RSES_4
# negatively worded items
RSES_2 ~~ RSES_5 + RSES_6
RSES_5 ~~ RSES_6
'
results_list[["2d"]] <- fit_and_summarize(
  model_2d, rses_data,
  "Table 2 Model (d): Over-identified CTCU",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 2 Model (e): Over-identified CTCU model (df = 6)
# Correlated residuals are added only among the negatively worded items.
# The positively worded items serve as the reference method effect.
# -----------------------------------------------------------------------------
model_2e <- '
F =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6
# no correlated residuals among positively worded items
# correlated residuals only among negatively worded items
RSES_2 ~~ RSES_5 + RSES_6
RSES_5 ~~ RSES_6
'
results_list[["2e"]] <- fit_and_summarize(
  model_2e, rses_data,
  "Table 2 Model (e): Over-identified CTCU (one method effect)",
  "RSES Original"
)


# ==============================================================================
# Table 3: Examples of Correlated-Factor Model with Correlated Errors
#            but No Cross-Loadings (RSES data; items 1-8)
# ==============================================================================

# -----------------------------------------------------------------------------
# Table 3 Model (a): Over-identified (df = 15)
# Two identified one-factor models combined by allowing F1 and F2 to correlate
# (Theorem 4). Within F1, item 4 is correlated with items 1 and 2. Within F2,
# item 8 is correlated with items 6 and 7.
# -----------------------------------------------------------------------------
model_3a <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4
F2 =~ 1*RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ F2
RSES_4 ~~ RSES_1 + RSES_2
RSES_8 ~~ RSES_6 + RSES_7
'
results_list[["3a"]] <- fit_and_summarize(
  model_3a, rses_data,
  "Table 3 Model (a): Over-identified (Model-Combining Rule)",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 3 Model (b): Just-identified (df = 0)
# On top of Model (a), all cross-factor correlated residuals are added
# EXCEPT between RSES_1 and RSES_5, so that the factor covariance can be
# identified (Theorem 5).
# -----------------------------------------------------------------------------
model_3b <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4
F2 =~ 1*RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ F2
# within-factor correlated residuals (same as Model a)
RSES_4 ~~ RSES_1 + RSES_2
RSES_8 ~~ RSES_5 + RSES_6
# cross-factor correlated residuals (all pairs except RSES_1--RSES_5)
RSES_1 ~~ RSES_6 + RSES_7 + RSES_8
RSES_2 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
RSES_3 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
RSES_4 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
'
results_list[["3b"]] <- fit_and_summarize(
  model_3b, rses_data,
  "Table 3 Model (b): Just-identified (Cross-Factor-Correlated-Error Rule)",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 3 Model (c): Just-identified (df = 0)
# Modified Model (b): factor correlation fixed to 0, but the previously
# omitted cross-factor correlated residual (RSES_1--RSES_5) is added.
# -----------------------------------------------------------------------------
model_3c <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4
F2 =~ 1*RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ 0*F2
# within-factor correlated residuals (same as Model a)
RSES_4 ~~ RSES_1 + RSES_2
RSES_8 ~~ RSES_5 + RSES_6
# cross-factor correlated residuals (all 16 pairs)
RSES_1 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
RSES_2 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
RSES_3 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
RSES_4 ~~ RSES_5 + RSES_6 + RSES_7 + RSES_8
'
results_list[["3c"]] <- fit_and_summarize(
  model_3c, rses_data,
  "Table 3 Model (c): Just-identified (orthogonal factors)",
  "RSES Original"
)


# ==============================================================================
# Table 4: Examples of Correlated-Factor Model with Cross-Loadings
#            and Correlated Errors
#
# NOTE: The two-factor cross-loading models in Table 4 do not converge when
# fitted to the RSES data because the RSES scale is essentially unidimensional
# in this sample (all items load on a single dominant factor). To demonstrate
# that these models are identified and estimable, we generate a simulated
# dataset (N = 500) from a population that has the two-factor cross-loading
# structure shown in Table 4.
# ==============================================================================

set.seed(42)
N <- 500
phi <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
F_sim <- mvrnorm(N, mu = c(0, 0), Sigma = phi)

# Loadings used to generate the simulated data:
#   x1 loads only on F1; x2 loads only on F2;
#   x3-x8 cross-load on both factors.
lambda1_sim <- c(1.0, 0.0, 0.7, 0.6, 0.7, 0.8, 0.7, 0.5)
lambda2_sim <- c(0.0, 1.0, 0.3, 0.4, 0.3, 0.2, 0.3, 0.5)

X_sim <- matrix(0, N, 8)
for (i in 1:8) {
  X_sim[, i] <- lambda1_sim[i] * F_sim[, 1] +
                lambda2_sim[i] * F_sim[, 2] +
                rnorm(N, 0, sqrt(0.4))
}
sim_data <- as.data.frame(X_sim)
names(sim_data) <- paste0("RSES_", 1:8)

# -----------------------------------------------------------------------------
# Table 4 Model (a): Over-identified (lavaan df = 13)
# x1 (RSES_1) loads exclusively on F1; x2 (RSES_2) loads exclusively on F2;
# all remaining items cross-load on both factors (Theorem 6).
#
# NOTE: lavaan returns df = 13 (8 indicators -> 36 moments - 23 parameters).
# The manuscript prints df = 15, which does not match this structure.
# -----------------------------------------------------------------------------
model_4a <- '
F1 =~ 1*RSES_1 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F2 =~ 1*RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ F2
'
results_list[["4a"]] <- fit_and_summarize(
  model_4a, sim_data,
  "Table 4 Model (a): Over-identified (Cross-Loading Rule, p >= 3m)",
  "Simulated"
)

# -----------------------------------------------------------------------------
# Table 4 Model (b): Over-identified (lavaan df = 13)
# EFA-model initial solution: factor correlation fixed to 0, and an
# additional cross-loading of x2 (RSES_2) on F1 is added.
#
# NOTE: lavaan returns df = 13. The manuscript prints df = 15.
# -----------------------------------------------------------------------------
model_4b <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F2 =~ 1*RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ 0*F2
'
results_list[["4b"]] <- fit_and_summarize(
  model_4b, sim_data,
  "Table 4 Model (b): Over-identified (EFA initial solution)",
  "Simulated"
)

# -----------------------------------------------------------------------------
# Table 4 Model (c): Over-identified (lavaan df = 1)
# Adds correlated residuals on top of Model (a) according to Theorem 8.
# With m = 2, the three base sets are:
#   x1 = {RSES_1, RSES_2}, x2 = {RSES_3, RSES_4}, x3 = {RSES_5, RSES_6}.
# Correlated residuals occur within each set but not across sets.
# The additional indicators RSES_7 and RSES_8 have correlated residuals
# with all indicators except those in one fully identified set.
#
# NOTE: lavaan returns df = 1 (36 moments - 35 parameters). The manuscript
# prints df = 3.
# -----------------------------------------------------------------------------
model_4c <- '
F1 =~ 1*RSES_1 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F2 =~ 1*RSES_2 + RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ F2
# correlated residuals within the three base sets
RSES_1 ~~ RSES_2
RSES_3 ~~ RSES_4
RSES_5 ~~ RSES_6
# RSES_7 correlated with all items except those in set x1 = {RSES_1, RSES_2}
RSES_7 ~~ RSES_3 + RSES_4 + RSES_5 + RSES_6 + RSES_8
# RSES_8 correlated with all items except those in set x3 = {RSES_5, RSES_6}
RSES_8 ~~ RSES_1 + RSES_2 + RSES_3 + RSES_4 + RSES_7
'
results_list[["4c"]] <- fit_and_summarize(
  model_4c, sim_data,
  "Table 4 Model (c): Over-identified (Cross-Loading-Plus-Correlated-Error Rule)",
  "Simulated"
)


# ==============================================================================
# Table 5: Examples of CTCU Models with Multiple Factors
# (RSES data; items 1-9)
#
# NOTE: The models converge on the RSES data, but the degrees of freedom
# reported by lavaan do not match those printed in the paper. With the number
# of indicators shown in the diagrams, the lavaan df values are the
# mathematically correct ones. The df values in the paper appear inconsistent
# with the stated number of free parameters and the number of indicators.
# ==============================================================================

# -----------------------------------------------------------------------------
# Table 5 Model (a): Over-identified (lavaan df = 15; paper reports df = 24)
# Three-factor CTCU model with three method effects.
# Each method effect is a set of three items, one per factor:
#   method 1: {RSES_1, RSES_4, RSES_7}
#   method 2: {RSES_2, RSES_5, RSES_8}
#   method 3: {RSES_3, RSES_6, RSES_9}
# Correlated residuals occur within each method effect but not across
# method effects (Theorem 8 with m = 3).
# -----------------------------------------------------------------------------
model_5a <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3
F2 =~ 1*RSES_4 + RSES_5 + RSES_6
F3 =~ 1*RSES_7 + RSES_8 + RSES_9
F1 ~~ F2 + F3
F2 ~~ F3
# method effect 1
RSES_1 ~~ RSES_4 + RSES_7
RSES_4 ~~ RSES_7
# method effect 2
RSES_2 ~~ RSES_5 + RSES_8
RSES_5 ~~ RSES_8
# method effect 3
RSES_3 ~~ RSES_6 + RSES_9
RSES_6 ~~ RSES_9
'
results_list[["5a"]] <- fit_and_summarize(
  model_5a, rses_data,
  "Table 5 Model (a): Over-identified (3 factors, 3 method effects)",
  "RSES Original"
)

# -----------------------------------------------------------------------------
# Table 5 Model (b): Over-identified (lavaan df = 11; paper reports df = 15)
# Two-factor CTCU model with two method effects plus one extra method-effect
# set needed for Theorem 8.
# Base sets for Condition 1 (no cross-set correlated residuals):
#   x1 = {RSES_1, RSES_2}, x2 = {RSES_3, RSES_4}, x3 = {RSES_5, RSES_6}.
# Extra set x4 = {RSES_7, RSES_8}: items have correlated residuals with
# each other and with items in set x2, but not with items in x1 or x3.
#
# NOTE ON q/df: lavaan reports 25 free parameters and df = 11, whereas the
# paper reports q = 23 and df = 15. The paper's counts appear inconsistent
# with the model structure shown in the diagram.
# -----------------------------------------------------------------------------
model_5b <- '
F1 =~ 1*RSES_1 + RSES_2 + RSES_3 + RSES_4
F2 =~ 1*RSES_5 + RSES_6 + RSES_7 + RSES_8
F1 ~~ F2
# within-set correlated residuals for x1, x2, x3, x4
RSES_1 ~~ RSES_2
RSES_3 ~~ RSES_4
RSES_5 ~~ RSES_6
RSES_7 ~~ RSES_8
# extra set x4 correlated with set x2 only
RSES_7 ~~ RSES_3 + RSES_4
RSES_8 ~~ RSES_3 + RSES_4
'
results_list[["5b"]] <- fit_and_summarize(
  model_5b, rses_data,
  "Table 5 Model (b): Over-identified (2 factors, method effects)",
  "RSES Original"
)


# ==============================================================================
# Final summary table
# ==============================================================================
cat("\n\n============================================================\n")
cat("SUMMARY OF ALL MODELS\n")
cat("============================================================\n\n")
results_df <- do.call(rbind, results_list)
print(results_df, row.names = FALSE)

# Also save the summary as a clean CSV file
write.csv(results_df, "model_summary.csv", row.names = FALSE)
cat("\nSummary table also saved to: model_summary.csv\n")
