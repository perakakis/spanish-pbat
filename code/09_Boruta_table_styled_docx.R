# Generate styled (heatmap-shaded) Word table of Boruta results
# Mirrors the visual design of figures/Boruta_table.pdf:
#   lower rank (more important) = darker grey, higher rank = lighter grey, NA = white
# Output: tables/Boruta_styled_table.docx

library(flextable)
library(officer)

# ---- Load data ---------------------------------------------------------------
d <- read.csv("./tables/Boruta_results.csv", stringsAsFactors = FALSE,
              na.strings = c("", "NA"))

outcome_cols <- c("Sad", "Anxious", "Stressed", "Angry", "NoSupp", "VITAL", "HEALTH")

# ---- Colour helper: rank -> hex grey ----------------------------------------
# Matches ggplot scale_fill_gradient(low="grey", high="white"):
#   rank 1 (most important) -> "#808080" (grey), max rank -> "#FFFFFF" (white)
rank_to_hex <- function(rank_vals) {
  all_ranks <- unlist(d[, outcome_cols], use.names = FALSE)
  rmin <- min(all_ranks, na.rm = TRUE)   # should be 1
  rmax <- max(all_ranks, na.rm = TRUE)
  
  sapply(rank_vals, function(r) {
    if (is.na(r)) return("#FFFFFF")
    # proportion: rank 1 -> 0 (dark grey), rank_max -> 1 (white)
    p <- (r - rmin) / (rmax - rmin)
    # interpolate between #808080 (128,128,128) and #FFFFFF (255,255,255)
    v <- round(128 + p * (255 - 128))
    sprintf("#%02X%02X%02X", v, v, v)
  })
}

# Pre-compute all fill colours
fill_matrix <- matrix(NA_character_, nrow = nrow(d), ncol = length(outcome_cols))
for (j in seq_along(outcome_cols)) {
  fill_matrix[, j] <- rank_to_hex(d[[outcome_cols[j]]])
}

# ---- Display values: NA -> em-dash ------------------------------------------
d_display <- d
for (col in outcome_cols) {
  d_display[[col]] <- ifelse(is.na(d[[col]]), "\u2014",
                             as.character(d[[col]]))
}

# ---- Build flextable ---------------------------------------------------------
ft <- flextable(d_display)

# Column labels
ft <- set_header_labels(ft,
  Predictor = "PBAT Item",
  Sad       = "Sad",
  Anxious   = "Anxious",
  Stressed  = "Stressed",
  Angry     = "Angry",
  NoSupp    = "No Support",
  VITAL     = "Vitality",
  HEALTH    = "Health"
)

# Spanner header
ft <- add_header_row(ft,
  values    = c("", "Emotional Symptoms (STOP-D)", "Well-being"),
  colwidths = c(1, 5, 2)
)

# ---- Apply per-cell shading --------------------------------------------------
for (i in seq_len(nrow(d))) {
  for (j in seq_along(outcome_cols)) {
    ft <- bg(ft, i = i, j = j + 1,   # j+1 because column 1 is Predictor
             bg = fill_matrix[i, j], part = "body")
  }
}

# ---- Typography & layout -----------------------------------------------------
ft <- font(ft, fontname = "Calibri", part = "all")
ft <- fontsize(ft, size = 10, part = "all")
ft <- bold(ft, bold = TRUE, part = "header")

# Borders
ft <- border_remove(ft)
outer_border <- fp_border(color = "#444444", width = 1.5)
inner_border <- fp_border(color = "#BBBBBB", width = 0.5)

ft <- border_outer(ft, border = outer_border, part = "all")
ft <- border_inner_h(ft, border = inner_border, part = "body")
ft <- border_inner_v(ft, border = inner_border, part = "body")
ft <- hline_bottom(ft, border = fp_border(color = "#444444", width = 1.5), part = "header")

# Header background: dark slate to match the heatmap aesthetic
ft <- bg(ft, bg = "#4A4A4A", part = "header")
ft <- color(ft, color = "white", part = "header")

# Alignment
ft <- align(ft, align = "center", part = "header")
ft <- align(ft, j = 1, align = "left",   part = "body")
ft <- align(ft, j = 2:8, align = "center", part = "body")

# Column widths
ft <- width(ft, j = 1, width = 3.0)
ft <- width(ft, j = 2:8, width = 0.78)

# Row height & padding
ft <- height_all(ft, height = 0.28, part = "body")
ft <- padding(ft, padding.top = 4, padding.bottom = 4,
              padding.left = 4, padding.right = 4, part = "all")

# ---- Legend note -------------------------------------------------------------
note_text <- paste0(
  "Note. Cell shading reflects Boruta importance rank: darker grey = lower rank ",
  "(stronger predictor), white = item not selected (\u2014). ",
  "STOP-D = Symptom Tracking of Psychological Disorders."
)

# ---- Assemble Word document --------------------------------------------------
doc <- read_docx()

# Page margins (narrow)
doc <- body_set_default_section(doc,
  prop_section(
    page_margins = page_mar(top = 1, bottom = 1, left = 1, right = 1)
  )
)

doc <- body_add_par(doc, "Table X", style = "Normal")
doc <- body_add_par(doc,
  "Boruta Feature Importance Rankings for PBAT Items Across Outcome Variables",
  style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, note_text, style = "Normal")

out_path <- "./tables/Boruta_styled_table.docx"
print(doc, target = out_path)
message("Saved: ", out_path)
