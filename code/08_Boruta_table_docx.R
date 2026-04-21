# Generate APA-style Word table of Boruta results
# Requires: flextable, officer

library(flextable)
library(officer)

# ---- Load data ---------------------------------------------------------------
d <- read.csv("./tables/Boruta_results.csv", stringsAsFactors = FALSE)

# Replace NA with em-dash for display
d[is.na(d)] <- "\u2014"

# ---- Column headers (APA-friendly) ------------------------------------------
col_labels <- c(
  Predictor  = "PBAT Item",
  Sad        = "Sad",
  Anxious    = "Anxious",
  Stressed   = "Stressed",
  Angry      = "Angry",
  NoSupp     = "No Support",
  VITAL      = "Vitality",
  HEALTH     = "Health"
)

# ---- Build flextable ---------------------------------------------------------
ft <- flextable(d)

# Apply column labels
ft <- set_header_labels(ft, values = col_labels)

# Add a spanner header: two groups of outcome columns
ft <- add_header_row(
  ft,
  values = c("", "Emotional Symptoms (STOP-D)", "Well-being"),
  colwidths = c(1, 5, 2)
)

# ---- APA formatting ----------------------------------------------------------
# Font
ft <- font(ft, fontname = "Times New Roman", part = "all")
ft <- fontsize(ft, size = 11, part = "all")

# Borders: APA uses top/bottom lines on header rows and a bottom line on the table
ft <- border_remove(ft)

thick <- fp_border(color = "black", width = 1.5)
thin  <- fp_border(color = "black", width = 0.5)

ft <- hline_top(ft,    border = thick, part = "header")
ft <- hline(ft,        border = thin,  part = "header", i = 1)   # under spanner row
ft <- hline_bottom(ft, border = thick, part = "header")
ft <- hline_bottom(ft, border = thick, part = "body")

# Alignment
ft <- align(ft, align = "center", part = "header")
ft <- align(ft, j = 1, align = "left",   part = "body")
ft <- align(ft, j = 2:8, align = "center", part = "body")

# Column widths (inches)
ft <- width(ft, j = 1, width = 3.2)
ft <- width(ft, j = 2:8, width = 0.75)

# Padding
ft <- padding(ft, padding.top = 3, padding.bottom = 3, part = "all")

# Bold header text
ft <- bold(ft, bold = TRUE, part = "header")

# ---- Caption / Note ----------------------------------------------------------
caption_text <- "Table X"
title_text   <- "Boruta Feature Importance Rankings for PBAT Items Across Outcome Variables"
note_text    <- "Note. Rankings indicate the importance rank of each PBAT item as a predictor of each outcome variable identified by the Boruta algorithm. Items are grouped by subscale: Positive Activation (rows 1\u20137) and Negative Avoidance (rows 8\u201314), Variation (rows 15\u201316), and Retention (rows 17\u201318). \u2014 = item not selected as an important predictor. STOP-D = Symptom Tracking of Psychological Disorders."

# ---- Assemble Word document --------------------------------------------------
doc <- read_docx()

doc <- body_add_par(doc, caption_text, style = "Normal")
doc <- body_add_par(doc, title_text,   style = "Normal")
doc <- body_add_flextable(doc, ft)
doc <- body_add_par(doc, note_text,    style = "Normal")

# Italic "Note." prefix
doc <- body_replace_all_text(doc, "Note.", "Note.", warn = FALSE)  # placeholder; styled below

# ---- Save --------------------------------------------------------------------
out_path <- "./tables/Boruta_APA_table.docx"
print(doc, target = out_path)
message("Saved: ", out_path)
