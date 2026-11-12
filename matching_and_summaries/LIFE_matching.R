# LIFE resting EEG inventory, ID matching and repeated-row diagnostics
# Run: Rscript /data/p_02865/code/matching_and_summaries/LIFE_matching.R
# Requires: readxl. Distinct, nonempty IDmatch rows are preserved (not merged).

library(readxl)

life <- "/data/p_02865/data/LIFE"
out <- file.path(life, "derivatives/matching_and_summaries")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
cell_values <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- NA_character_
  x
}
clean_id <- function(x) {
  x <- toupper(cell_values(x))
  x[x == "NA" & !is.na(x)] <- NA_character_
  x
}
strip_trailing <- function(x) {
  x <- sub("[_-]+$", "", x)
  x <- sub("[_-]2$", "", x)
  sub("[_-]+$", "", x)
}
nonempty <- function(paths) {
  info <- file.info(paths)
  !is.na(info$size) & info$size > 0
}
header_link_ok <- function(path, field) {
  lines <- tryCatch(readLines(path, n = 100L, warn = FALSE),
                    error = function(e) character())
  entry <- grep(paste0("^\\s*", field, "\\s*="), lines,
                value = TRUE, perl = TRUE, ignore.case = TRUE)
  if (!length(entry)) return(FALSE)
  target <- gsub("\\\\", "/", trimws(sub("^[^=]*=", "", entry[1])))
  if (!grepl("^(/|[A-Za-z]:/)", target)) target <- file.path(dirname(path), target)
  isTRUE(nonempty(target))
}
classify <- function(x) {
  x <- cell_values(x)
  n <- length(unique(x[!is.na(x)]))
  if (n == 0) "all_missing" else if (n > 1) "conflict" else
    if (anyNA(x)) "complementary" else "same"
}
show_values <- function(x) {
  x <- unique(cell_values(x))
  paste(x[!is.na(x)], collapse = " | ")
}

# 1. Identify recordings in all three EEG folders.
eeg_dirs <- c(
  ruhe = file.path(life, "raw/EEG/ruhe"),
  esf_data = file.path(life, "raw/EEG/esf_data_2014-08-21/EEG/Ruhe"),
  fehlende_Daten = file.path(life, "raw/EEG/fehlende_Daten")
)
if (any(!dir.exists(eeg_dirs))) stop("EEG folders missing: ",
                                    paste(eeg_dirs[!dir.exists(eeg_dirs)], collapse = ", "))
get_files <- function(pattern) {
  do.call(rbind, lapply(names(eeg_dirs), function(src) {
    paths <- list.files(eeg_dirs[[src]], pattern = pattern,
                        full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    data.frame(source = rep(src, length(paths)), path = paths)
  }))
}
all_vhdr <- get_files("\\.vhdr$")
if (!nrow(all_vhdr)) stop("No VHDR recordings found.")
aux <- get_files("\\.(eeg|vmrk)$")
# An orphan .eeg/.vmrk has no same-stem .vhdr in the same folder.
strip_ext <- function(x) sub("(?i)\\.(vhdr|vmrk|eeg)$", "", x, perl = TRUE)
orphan <- aux[!strip_ext(aux$path) %in% strip_ext(all_vhdr$path), , drop = FALSE]

rec <- all_vhdr
rec$vhdr_path <- rec$path
rec$vhdr_file <- basename(rec$path)
stem <- sub("(?i)\\.vhdr$", "", rec$vhdr_file, perl = TRUE)
corrected_stem <- stem
# Handle filename formatting errors by pattern, without participant-specific IDs.
# Some files omit "_Ruhe", contain a typo ("_Ruhen"/"_Ruihe"),
# or omit the standard "LI" prefix.
corrected_stem <- sub("(?i)_(Ruhen|Ruihe)$", "_Ruhe", stem, perl = TRUE)
missing_ruhe <- grepl("(?i)^LI[0-9]+[A-Z]?$", corrected_stem, perl = TRUE)
corrected_stem[missing_ruhe] <- paste0(corrected_stem[missing_ruhe], "_Ruhe")
known_typo <- stem != corrected_stem
no_li <- grepl("(?i)^[0-9]+[A-Z]?_Ruhe(?:_[0-9]+)?$", corrected_stem, perl = TRUE)
recognized <- grepl("(?i)^LI[A-Z0-9_-]+_Ruhe(?:_[0-9]+)?$",
                    corrected_stem, perl = TRUE) | no_li
unrecognized <- rec[!recognized, , drop = FALSE]
rec <- rec[recognized, , drop = FALSE]
corrected_stem <- corrected_stem[recognized]
known_typo <- known_typo[recognized]
no_li <- no_li[recognized]
rec$filename_rule <- ifelse(known_typo, "corrected_typo",
                            ifelse(no_li, "no_LI_prefix", "standard"))
rec$raw_life_id <- clean_id(sub("(?i)_Ruhe(?:_[0-9]+)?$", "", corrected_stem, perl = TRUE))
rec$suffix <- ifelse(grepl("(?i)_Ruhe_[0-9]+$", corrected_stem, perl = TRUE),
                     sub("(?i)^.*_Ruhe_([0-9]+)$", "\\1", corrected_stem, perl = TRUE),
                     "none")
rec$eeg_path <- sub("(?i)\\.vhdr$", ".eeg", rec$vhdr_path, perl = TRUE)
rec$vmrk_path <- sub("(?i)\\.vhdr$", ".vmrk", rec$vhdr_path, perl = TRUE)
rec$complete <- nonempty(rec$vhdr_path) & nonempty(rec$eeg_path) & nonempty(rec$vmrk_path)
rec$header_links_ok <- vapply(rec$vhdr_path, function(p) {
  header_link_ok(p, "DataFile") && header_link_ok(p, "MarkerFile")
}, logical(1))
rec$ready <- rec$complete & rec$header_links_ok

# 2. Read IDmatch, remove empty data rows and exact duplicate rows.
idmatch <- read_excel(file.path(life, "tabular/IDmatch.xlsx"),
                      sheet = "Zusammenführen1", .name_repair = "minimal")
if (!"Column1" %in% names(idmatch)) stop("IDmatch lacks Column1.")
ids_all <- clean_id(idmatch$Column1)
fields <- grep("^Table3[._]", names(idmatch), value = TRUE)
if (!length(fields)) stop("No Table3 fields found.")
find_field <- function(suffix) {
  hits <- grep(paste0("^Table3[._]", suffix, "$"), names(idmatch), value = TRUE)
  if (length(hits) != 1L) stop("Expected one Table3 field named ", suffix)
  hits
}
key_cols <- vapply(c("TEILNEHMER_SIC", "EEG_GRUPPE", "TEILNEHMER_GESCHLECHT_EEG"),
                   find_field, character(1))
# All three key fields must be missing. To avoid losing clinical data, also
# require that ALL Table3 fields are empty (only Column1/Column2 may remain).
all_missing <- function(cols) Reduce(`&`, lapply(idmatch[cols], function(x) is.na(cell_values(x))))
empty_keys <- all_missing(key_cols)
empty_data <- empty_keys & all_missing(fields)
partial_keys <- empty_keys & !empty_data  # Keep and report these, if present.

# Compare all columns after treating blank strings as missing, within each ID.
duplicate_row <- rep(FALSE, nrow(idmatch))
duplicate_of <- rep(NA_integer_, nrow(idmatch))
remaining <- which(!empty_data & !is.na(ids_all))
groups <- split(remaining, ids_all[remaining])
for (idx in groups) {
  if (length(idx) < 2L) next
  comparison <- as.data.frame(lapply(idmatch[idx, , drop = FALSE], cell_values),
                              check.names = FALSE)
  comparison[["Column1"]] <- clean_id(comparison[["Column1"]])
  is_copy <- duplicated(comparison)
  if (!any(is_copy)) next
  duplicate_row[idx[is_copy]] <- TRUE
  for (j in which(is_copy)) {
    earlier <- which(vapply(seq_len(j - 1L), function(k)
      identical(lapply(comparison, `[[`, k), lapply(comparison, `[[`, j)),
      logical(1)))[1]
    duplicate_of[idx[j]] <- idx[earlier] + 1L  # Excel row, header included.
  }
}
keep_idx <- which(!empty_data & !duplicate_row)
clean_data <- idmatch[keep_idx, , drop = FALSE]
ids <- ids_all[keep_idx]

removals <- data.frame(
  excel_row = c(which(empty_data), which(duplicate_row)) + 1L,
  life_id = ids_all[c(which(empty_data), which(duplicate_row))],
  reason = c(rep("empty_Table3_data", sum(empty_data)),
             rep("identical_duplicate", sum(duplicate_row))),
  duplicate_of_excel_row = c(rep(NA_integer_, sum(empty_data)),
                             duplicate_of[duplicate_row])
)
removals <- removals[order(removals$excel_row), , drop = FALSE]
save_csv(removals, "IDmatch_row_removals.csv")

# 3. Match EEG participant IDs to the cleaned IDmatch table.
valid_ids <- unique(ids[!is.na(ids)])
raw_ids <- unique(rec$raw_life_id)
lookup <- data.frame(raw_life_id = raw_ids, life_id = raw_ids,
                     match_method = "unmatched", candidate_ids = "")
for (i in seq_len(nrow(lookup))) {
  raw <- lookup$raw_life_id[i]
  if (raw %in% valid_ids) {
    lookup$match_method[i] <- "exact"
    next
  }
  base <- strip_trailing(raw)
  alternative <- if (startsWith(base, "LI")) sub("^LI", "", base) else paste0("LI", base)
  candidates <- valid_ids[strip_trailing(valid_ids) %in% c(base, alternative)]
  lookup$candidate_ids[i] <- paste(candidates, collapse = " | ")
  if (length(candidates) == 1L) {
    lookup$life_id[i] <- candidates
    lookup$match_method[i] <- if (strip_trailing(candidates) == base)
      "suffix_cleanup" else "LI_prefix_cleanup"
  } else if (length(candidates) > 1L) {
    lookup$match_method[i] <- "ambiguous_review"
  }
}
rec$life_id <- lookup$life_id[match(rec$raw_life_id, lookup$raw_life_id)]
rec$match_method <- lookup$match_method[match(rec$raw_life_id, lookup$raw_life_id)]

# Select one recording per person: complete/valid first, then Ruhe_2,
# then the main ruhe folder. All recordings remain in the inventory.
ord <- order(rec$life_id, !rec$ready, rec$suffix != "2",
             rec$source != "ruhe", rec$vhdr_path)
rec$selected <- FALSE
rec$selected[ord[!duplicated(rec$life_id[ord])]] <- TRUE
rec$path <- NULL
save_csv(rec, "EEG_recordings_inventory.csv")

part <- rec[rec$selected, c("life_id", "vhdr_file", "vhdr_path", "source", "suffix", "ready")]
names(part)[2:6] <- c("selected_vhdr_file", "selected_vhdr_path",
                       "selected_source", "selected_suffix", "selected_ready")
part <- part[order(part$life_id), , drop = FALSE]
part$n_recordings <- as.integer(table(rec$life_id)[part$life_id])
part$n_ready_recordings <- as.integer(tapply(rec$ready, rec$life_id, sum)[part$life_id])
part$has_ruhe_2 <- part$life_id %in% rec$life_id[rec$suffix == "2"]
part$n_sources <- as.integer(tapply(rec$source, rec$life_id,
                                   function(x) length(unique(x)))[part$life_id])
part$n_idmatch_rows_original <- as.integer(table(ids_all[!is.na(ids_all)])[part$life_id])
part$n_idmatch_rows_original[is.na(part$n_idmatch_rows_original)] <- 0L
part$n_idmatch_rows <- as.integer(table(ids[!is.na(ids)])[part$life_id])
part$n_idmatch_rows[is.na(part$n_idmatch_rows)] <- 0L
part$matched <- part$n_idmatch_rows > 0L & part$n_ready_recordings > 0L
save_csv(part, "EEG_participants_inventory.csv")

# Keep ALL distinct, nonempty rows (including multiple EEG groups per participant).
matched_ids <- part$life_id[part$matched]
matched_idx <- which(!is.na(ids) & ids %in% matched_ids)
matched_data <- clean_data[matched_idx, , drop = FALSE]
save_csv(matched_data, "IDmatch_EEG.csv")
save_csv(data.frame(field_name = sub("^Table3[._]", "", fields),
                    original_column = fields), "IDmatch_variables.csv")

# 4. Compare remaining repeated rows: complementary vs conflicting values.
dup_ids <- part$life_id[part$matched & part$n_idmatch_rows > 1L]
date_col <- find_field("EEG_DATUM")
group_col <- find_field("EEG_GRUPPE")
summaries <- vector("list", length(dup_ids))
decisions <- list()
k <- 0L
for (i in seq_along(dup_ids)) {
  id <- dup_ids[i]
  idx <- which(ids == id)  # which() avoids artificial NA-indexed rows.
  rows <- clean_data[idx, , drop = FALSE]
  excel_rows <- keep_idx[idx] + 1L
  states <- vapply(rows, classify, character(1))
  summaries[[i]] <- data.frame(
    life_id = id,
    n_rows_before = part$n_idmatch_rows_original[match(id, part$life_id)],
    n_rows_retained = nrow(rows),
    source_excel_rows = paste(excel_rows, collapse = " | "),
    n_complementary_columns = sum(states == "complementary"),
    n_conflicting_columns = sum(states == "conflict"),
    eeg_date_status = classify(rows[[date_col]]),
    eeg_dates = show_values(rows[[date_col]]),
    eeg_group_status = classify(rows[[group_col]]),
    eeg_groups = show_values(rows[[group_col]]),
    n_eeg_recordings = part$n_recordings[match(id, part$life_id)],
    eeg_files = paste(rec$vhdr_file[rec$life_id == id], collapse = " | "),
    different_columns = paste(names(states)[states %in% c("complementary", "conflict")],
                              collapse = "; ")
  )
  for (j in which(states %in% c("complementary", "conflict"))) {
    x <- cell_values(rows[[j]])
    k <- k + 1L
    decisions[[k]] <- data.frame(
      life_id = id, column = names(rows)[j], decision = states[j],
      source_excel_rows = paste(excel_rows, collapse = " | "),
      values_by_excel_row = paste(paste0(excel_rows, ": ",
                                       ifelse(is.na(x), "<missing>", x)), collapse = " | ")
    )
  }
}
dup_summary <- if (length(summaries)) do.call(rbind, summaries) else
  data.frame(life_id = character(), n_rows_retained = integer(),
             eeg_date_status = character(), eeg_group_status = character())
row_decisions <- if (length(decisions)) do.call(rbind, decisions) else
  data.frame(life_id = character(), column = character(), decision = character(),
             source_excel_rows = character(), values_by_excel_row = character())
save_csv(dup_summary, "IDmatch_duplicates.csv")
save_csv(row_decisions, "IDmatch_row_decisions.csv")

# 5. Small audit files for imperfect filename/ID matches and orphan files.
review <- lookup[lookup$match_method != "exact", , drop = FALSE]
review$recording_files <- vapply(review$raw_life_id, function(id)
  paste(rec$vhdr_file[rec$raw_life_id == id], collapse = " | "), character(1))
save_csv(review, "EEG_ID_matching_review.csv")
issues <- data.frame(issue = character(), source = character(), filename = character(), path = character())
add_issue <- function(issue, x, path_col) {
  data.frame(issue = rep(issue, nrow(x)), source = x$source,
             filename = basename(x[[path_col]]), path = x[[path_col]])
}
issues <- rbind(issues,
                add_issue("unrecognized_vhdr", unrecognized, "vhdr_path"),
                add_issue("incomplete_or_invalid_vhdr", rec[!rec$ready, , drop = FALSE], "vhdr_path"),
                add_issue("orphan_eeg_or_vmrk", orphan, "path"))
save_csv(issues, "EEG_file_issues.csv")

# 6. Concise console and text report.
report <- capture.output({
  cat("========== EEG ==========" , "\n")
  cat("VHDR files found:", nrow(all_vhdr), "| recognized:", nrow(rec),
      "| unrecognized:", nrow(unrecognized), "\n")
  cat("Recognized EEG recordings by folder:\n"); print(table(rec$source))
  cat("Filename exceptions recognized:", sum(rec$filename_rule != "standard"), "\n")
  if (any(rec$filename_rule != "standard"))
    print(rec$vhdr_file[rec$filename_rule != "standard"])
  cat("Complete file triplets:", sum(rec$complete), "/", nrow(rec), "\n")
  cat("Complete triplets with valid VHDR links:", sum(rec$ready), "/", nrow(rec), "\n")
  cat("Orphan .eeg/.vmrk files:", nrow(orphan), "\n")
  if (nrow(orphan)) print(orphan$path)
  cat("Unique EEG participants:", nrow(part), "\n")
  cat("Participants with >1 recording:", sum(part$n_recordings > 1),
      "| with _Ruhe_2:", sum(part$has_ruhe_2), "\n")

  cat("\n========== MATCHING ==========" , "\n")
  cat("Original IDmatch rows:", nrow(idmatch), "\n")
  cat("Empty Table3 rows removed:", sum(empty_data), "\n")
  cat("Rows with three missing key fields but other Table3 data (kept):",
      sum(partial_keys), "\n")
  cat("Exact duplicate rows removed:", sum(duplicate_row), "\n")
  cat("IDmatch rows remaining:", nrow(clean_data), "\n")
  cat("Matched unique EEG participants:", sum(part$matched), "\n")
  cat("Matched IDmatch rows retained:", nrow(matched_data), "\n")
  cat("Participants with only empty IDmatch rows:",
      sum(part$n_idmatch_rows_original > 0 & part$n_idmatch_rows == 0), "\n")
  cat("EEG participants without usable IDmatch data:", sum(part$n_idmatch_rows == 0), "\n")
  if (any(part$n_idmatch_rows == 0)) print(part$life_id[part$n_idmatch_rows == 0])
  cat("IDs resolved by suffix/LI-prefix correction:",
      sum(lookup$match_method %in% c("suffix_cleanup", "LI_prefix_cleanup")), "\n")
  cat("Ambiguous matching IDs:", sum(lookup$match_method == "ambiguous_review"), "\n")
  cat("IDmatch rows per matched participant:\n")
  print(table(part$n_idmatch_rows[part$matched]))

  cat("\n========== REPEATED IDmatch ROWS ==========" , "\n")
  cat("Matched participants with repeated rows BEFORE cleanup:",
      sum(part$n_idmatch_rows_original > 1 & part$n_ready_recordings > 0), "\n")
  cat("Matched participants with repeated rows AFTER cleanup:", length(dup_ids), "\n")
  cat("Participants resolved to a single row by cleanup:",
      sum(part$n_idmatch_rows_original > 1 & part$n_idmatch_rows == 1 &
            part$n_ready_recordings > 0), "\n")
  cat("Remaining repeated-row counts:\n")
  print(table(dup_summary$n_rows_retained))
  cat("Participants with complementary columns:",
      sum(dup_summary$n_complementary_columns > 0), "\n")
  cat("Participants with conflicting nonmissing columns:",
      sum(dup_summary$n_conflicting_columns > 0), "\n")
  cat("EEG_DATUM status:\n"); print(table(dup_summary$eeg_date_status))
  cat("EEG_GRUPPE status:\n"); print(table(dup_summary$eeg_group_status))
  cat("EEG_GRUPPE values/patterns:\n")
  if (nrow(dup_summary)) print(head(sort(table(dup_summary$eeg_groups), decreasing = TRUE), 12))
  cat("Multiple EEG recordings vs repeated IDmatch rows:\n")
  print(table(multiple_eeg = part$n_recordings[part$matched] > 1,
              repeated_idmatch = part$n_idmatch_rows[part$matched] > 1))
  cat("Complementary fields logged:", sum(row_decisions$decision == "complementary"),
      "| conflicting fields logged:", sum(row_decisions$decision == "conflict"), "\n")
  cat("Table3 variables:", length(fields), "\n")
  cat("Distinct rows remain separate; no data from different rows were merged.\n")
})
cat(paste(report, collapse = "\n"), "\n")
writeLines(report, file.path(out, "LIFE_matching_summary.txt"))
cat("\nSaved outputs to:", out, "\n")
