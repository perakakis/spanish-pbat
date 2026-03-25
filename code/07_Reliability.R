# Reliability indices for the Spanish PBAT dataset
# Cronbach's alpha and McDonald's omega where appropriate
# Results saved to ./results/reliability_results.csv

library(psych)

df <- read.csv("./data/data_es.csv")
df <- na.omit(df)

# ------------------------------------------------------------------------------
# PBAT item mapping: CSV cols 4:21 → short labels via PBAT_order
# ------------------------------------------------------------------------------
short_labels_18 <- c(
  "Motivation+", "Health+", "Attention+", "Social+", "Affect+",
  "Challenge+", "Cognition+",
  "Motivation-", "Health-", "Attention-", "Social-", "Affect-",
  "Challenge-", "Cognition-",
  "Variation-", "Variation+", "Retention-", "Retention+"
)
PBAT_order <- c(16, 11, 5, 17, 13, 2, 14, 3, 8, 18, 6, 15, 7, 10, 4, 1, 9, 12)
colnames(df)[4:21] <- short_labels_18[PBAT_order]

# NEED mapping: NEED_order <- c(1,2,5,6,3,4) applied to NEED_labels
# Result: NEED1=Autonomy Sat, NEED2=Autonomy Frus,
#         NEED3=Competence Sat, NEED4=Competence Frus,
#         NEED5=Connection Sat, NEED6=Connection Frus
colnames(df)[31:36] <- c(
  "Autonomy.Sat", "Autonomy.Frus",
  "Competence.Sat", "Competence.Frus",
  "Connection.Sat", "Connection.Frus"
)

# ------------------------------------------------------------------------------
# Helper: extract alpha and omega from psych output into a tidy row
# ------------------------------------------------------------------------------
get_reliability <- function(items_df, scale_name, run_omega = FALSE) {
  a <- alpha(items_df, warnings = FALSE)
  alpha_val <- round(a$total$raw_alpha, 3)

  if (run_omega) {
    suppressWarnings(
      om <- omega(items_df, nfactors = 1, plot = FALSE)
    )
    omega_val <- round(om$omega.tot, 3)
  } else {
    omega_val <- NA
  }

  data.frame(
    Scale   = scale_name,
    n_items = ncol(items_df),
    alpha   = alpha_val,
    omega   = omega_val,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Vitality (VITAL1, VITAL2, VITAL3) — alpha + omega
# ------------------------------------------------------------------------------
cat("--- Vitality ---\n")
vital_df <- df[, c("VITAL1", "VITAL2", "VITAL3")]
print(alpha(vital_df, warnings = FALSE)$total)
suppressWarnings(print(omega(vital_df, nfactors = 1, plot = FALSE)))

# ------------------------------------------------------------------------------
# 2. Need Satisfaction (NEED1, NEED3, NEED5) — alpha only
# ------------------------------------------------------------------------------
cat("\n--- Need Satisfaction ---\n")
sat_df <- df[, c("Autonomy.Sat", "Competence.Sat", "Connection.Sat")]
print(alpha(sat_df, warnings = FALSE)$total)

# ------------------------------------------------------------------------------
# 3. Need Frustration (NEED2, NEED4, NEED6) — alpha only
# ------------------------------------------------------------------------------
cat("\n--- Need Frustration ---\n")
frus_df <- df[, c("Autonomy.Frus", "Competence.Frus", "Connection.Frus")]
print(alpha(frus_df, warnings = FALSE)$total)

# ------------------------------------------------------------------------------
# 4. PBAT Positive items (7 items) — alpha only
# ------------------------------------------------------------------------------
cat("\n--- PBAT Positive ---\n")
pos_items <- c("Motivation+", "Health+", "Attention+", "Social+",
               "Affect+", "Challenge+", "Cognition+")
pos_df <- df[, pos_items]
print(alpha(pos_df, warnings = FALSE)$total)

# ------------------------------------------------------------------------------
# 5. PBAT Negative items (7 items) — alpha only
# ------------------------------------------------------------------------------
cat("\n--- PBAT Negative ---\n")
neg_items <- c("Motivation-", "Health-", "Attention-", "Social-",
               "Affect-", "Challenge-", "Cognition-")
neg_df <- df[, neg_items]
print(alpha(neg_df, warnings = FALSE)$total)

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------
results <- rbind(
  get_reliability(vital_df, "Vitality",          run_omega = TRUE),
  get_reliability(sat_df,   "Need Satisfaction",  run_omega = FALSE),
  get_reliability(frus_df,  "Need Frustration",   run_omega = FALSE),
  get_reliability(pos_df,   "PBAT Positive",      run_omega = FALSE),
  get_reliability(neg_df,   "PBAT Negative",      run_omega = FALSE)
)

cat("\n--- Reliability summary ---\n")
print(results)

write.csv(results, "./results/reliability_results.csv", row.names = FALSE)
cat("\nResults saved to ./results/reliability_results.csv\n")
