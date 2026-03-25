# Measurement invariance testing for the PBAT two-factor model
# Compares Spain and Sweden across configural, metric, and scalar models
# Results saved to ./results/invariance_results.RData and ./results/invariance_fit.csv

library(lavaan)
library(semTools)

# ------------------------------------------------------------------------------
# Short labels for the 14 items used in the two-factor model
# ------------------------------------------------------------------------------
pos_items <- c("Motivation+", "Health+", "Attention+", "Social+",
               "Affect+", "Challenge+", "Cognition+")
neg_items <- c("Motivation-", "Health-", "Attention-", "Social-",
               "Affect-", "Challenge-", "Cognition-")
pbat_14   <- c(pos_items, neg_items)

# lavaan-safe column names: replace + with p, - with n
pbat_14_safe <- gsub("-", "n", pbat_14, fixed = TRUE)
pbat_14_safe <- gsub("+", "p", pbat_14_safe, fixed = TRUE)

# CFA model using lavaan-safe names
model <- "
  Positive =~ Motivationp + Healthp + Attentionp + Socialp +
              Affectp + Challengep + Cognitionp
  Negative =~ Motivationn + Healthn + Attentionn + Socialn +
              Affectn + Challengen + Cognitionn
"

# ------------------------------------------------------------------------------
# Load and prepare Spanish data
# ------------------------------------------------------------------------------
df_es <- read.csv("./data/data_es.csv")
df_es <- na.omit(df_es)

# Long-label index → short-label mapping (confirmed in 02_Descriptives.R + 03_Boruta_es.R)
short_labels_18 <- c(
  "Motivation+", "Health+", "Attention+", "Social+", "Affect+", "Challenge+", "Cognition+",
  "Motivation-", "Health-", "Attention-", "Social-", "Affect-", "Challenge-", "Cognition-",
  "Variation-",  "Variation+", "Retention-", "Retention+"
)

# PBAT_order maps CSV columns 4:21 → long-label indices → short labels
PBAT_order <- c(16, 11, 5, 17, 13, 2, 14, 3, 8, 18, 6, 15, 7, 10, 4, 1, 9, 12)
colnames(df_es)[4:21] <- short_labels_18[PBAT_order]

df_es <- df_es[, pbat_14]
colnames(df_es) <- pbat_14_safe
df_es$country <- "Spain"

# ------------------------------------------------------------------------------
# Load and prepare Swedish data
# ------------------------------------------------------------------------------
df_sw <- read.csv("./data/data_sw.csv")
df_sw <- na.omit(df_sw)

# Swedish CSV columns 3:20 are already in the order below (from 03b_Boruta_sw.R)
sw_col_labels <- c(
  "Variation+", "Social-",   "Affect+",    "Retention-", "Challenge-", "Health+",
  "Cognition-", "Attention+","Motivation-","Retention+", "Challenge+", "Variation-",
  "Cognition+", "Attention-","Social+",    "Motivation+","Health-",    "Affect-"
)
colnames(df_sw)[3:20] <- sw_col_labels

df_sw <- df_sw[, pbat_14]
colnames(df_sw) <- pbat_14_safe
df_sw$country <- "Sweden"

# ------------------------------------------------------------------------------
# Combine datasets
# ------------------------------------------------------------------------------
df_combined <- rbind(df_es, df_sw)
df_combined$country <- factor(df_combined$country)

cat("Sample sizes:\n")
print(table(df_combined$country))

# ------------------------------------------------------------------------------
# Step 1 — Single-group CFA per country
# ------------------------------------------------------------------------------
fit_es <- cfa(model, data = df_combined[df_combined$country == "Spain", ],
              estimator = "MLR")
fit_sw <- cfa(model, data = df_combined[df_combined$country == "Sweden", ],
              estimator = "MLR")

cat("\n--- Step 1: Single-group fit (Spain) ---\n")
print(fitMeasures(fit_es, c("cfi.robust", "rmsea.robust", "srmr")))

cat("\n--- Step 1: Single-group fit (Sweden) ---\n")
print(fitMeasures(fit_sw, c("cfi.robust", "rmsea.robust", "srmr")))

# ------------------------------------------------------------------------------
# Steps 2–4 — Multi-group invariance sequence
# ------------------------------------------------------------------------------
fit_configural <- cfa(model, data = df_combined, group = "country",
                      estimator = "MLR")

fit_metric     <- cfa(model, data = df_combined, group = "country",
                      group.equal = "loadings",
                      estimator = "MLR")

fit_scalar     <- cfa(model, data = df_combined, group = "country",
                      group.equal = c("loadings", "intercepts"),
                      estimator = "MLR")

# ------------------------------------------------------------------------------
# Step 5 — Partial scalar invariance
# Free intercepts identified by lavTestScore() on the scalar model, in
# descending order of X2: Socialp (39.5), Challengep (12.4), Healthp (10.7),
# Socialn (10.3), Cognitionn (9.8). Stopping criterion: dCFI vs metric <= .010.
# ------------------------------------------------------------------------------
partial_free <- c("Socialp", "Challengep", "Healthp", "Socialn", "Cognitionn")

fit_partial  <- cfa(model, data = df_combined, group = "country",
                    group.equal   = c("loadings", "intercepts"),
                    group.partial = paste0(partial_free, "~1"),
                    estimator = "MLR")

cat("\n--- Step 5: Partial scalar fit ---\n")
print(fitMeasures(fit_partial, c("cfi.robust", "rmsea.robust", "srmr")))

# ------------------------------------------------------------------------------
# Collect fit indices
# ------------------------------------------------------------------------------
fit_indices <- function(fit, label) {
  fm <- fitMeasures(fit, c("cfi.robust", "rmsea.robust", "rmsea.ci.lower.robust",
                            "rmsea.ci.upper.robust", "srmr", "chisq.scaled",
                            "df.scaled", "pvalue.scaled"))
  data.frame(
    Model     = label,
    CFI       = round(fm["cfi.robust"], 3),
    RMSEA     = round(fm["rmsea.robust"], 3),
    RMSEA_low = round(fm["rmsea.ci.lower.robust"], 3),
    RMSEA_up  = round(fm["rmsea.ci.upper.robust"], 3),
    SRMR      = round(fm["srmr"], 3),
    Chi2      = round(fm["chisq.scaled"], 3),
    df        = fm["df.scaled"],
    p         = round(fm["pvalue.scaled"], 3),
    row.names = NULL
  )
}

fit_table <- rbind(
  fit_indices(fit_configural, "Configural"),
  fit_indices(fit_metric,     "Metric"),
  fit_indices(fit_scalar,     "Scalar"),
  fit_indices(fit_partial,    "Partial Scalar")
)

# Add delta columns (each model vs. its predecessor)
fit_table$dCFI <- c(
  NA,
  round(fit_table$CFI[2] - fit_table$CFI[1], 3),
  round(fit_table$CFI[3] - fit_table$CFI[2], 3),
  round(fit_table$CFI[4] - fit_table$CFI[2], 3)
)
fit_table$dRMSEA <- c(
  NA,
  round(fit_table$RMSEA[2] - fit_table$RMSEA[1], 3),
  round(fit_table$RMSEA[3] - fit_table$RMSEA[2], 3),
  round(fit_table$RMSEA[4] - fit_table$RMSEA[2], 3)
)

cat("\n--- Invariance fit table ---\n")
print(fit_table)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------
write.csv(fit_table,
          "./results/invariance_fit.csv",
          row.names = FALSE)

save(fit_configural, fit_metric, fit_scalar, fit_partial,
     fit_es, fit_sw, fit_table,
     file = "./results/invariance_results.RData")

cat("\nResults saved to ./results/invariance_fit.csv and ./results/invariance_results.RData\n")
