# Run AFTER LIFE_generate_data_manual_review.R and manual review:
#   Rscript LIFE_resolve_matching.R [matching_output_dir] [generator_script]
# Review files (created once; never overwritten by the generator):
#   matching_and_summaries/manual/data_unresolved_rows_manual.csv
#   matching_and_summaries/manual/data_match_issues_manual.csv
# Enter 1 in select_primary/select for ONE candidate per case.
# For a case with no acceptable candidate, enter keep_missing in resolution.
# This script validates reviews, reruns the original data assembly with the
# selections applied, and only then creates derivatives/data.csv.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L)
  stop("Usage: Rscript LIFE_resolve_matching.R [matching_output_dir] [generator_script]")
life <- "/data/p_02865/data/LIFE"
out <- if (length(args) >= 1L) args[[1L]] else
  file.path(life, "derivatives/matching_and_summaries")

self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(self)) dirname(normalizePath(sub("^--file=", "", self[[1L]]))) else getwd()
resolver_script <- if (length(self)) normalizePath(sub("^--file=", "", self[[1L]])) else NA_character_
generator <- if (length(args) == 2L) args[[2L]] else
  file.path(script_dir, "04_LIFE_generate_data_manual_review.R")
if (!file.exists(generator)) stop("Generator script not found: ", generator)
if (!dir.exists(out)) stop("Matching output folder not found: ", out)

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing required review file: ", path)
  x <- read.csv(path, colClasses = "character", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = "", fileEncoding = "UTF-8-BOM")
  x[] <- lapply(x, function(col) {
    col <- trimws(as.character(col))
    col[is.na(col)] <- ""
    col
  })
  x
}

# Spreadsheet UTF-7 exports may encode underscores as +AF8- and hyphens as +AC0-.
# Decode only review CSVs; leave all other source-data readers unchanged.
read_review_csv <- function(path) {
  x <- read_csv(path)
  decode_review <- function(v) {
    v <- gsub("+AF8-", "_", v, fixed = TRUE)
    gsub("+AC0-", "-", v, fixed = TRUE)
  }
  names(x) <- decode_review(names(x))
  if (anyDuplicated(names(x)))
    stop("Duplicate column names after decoding review CSV: ", path)
  x[] <- lapply(x, decode_review)
  x
}

# Some spreadsheet programs encode CSV *quotes* as UTF-7 +ACI-, in
# addition to encoding underscores (+AF8-) and hyphens (+AC0-). The quotes
# must be decoded BEFORE read.csv(): otherwise embedded commas shift columns
# and a candidate value can be mistaken for source_path.
#
# The generator also reads the active manual files in a separate R process.
# Normalize each affected file once (with an exact-byte backup) so both
# scripts parse the same decisions. Do not change any selection values.
normalize_manual_utf7 <- function(path) {
  if (!file.exists(path)) stop("Missing required review file: ", path)
  bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  contents <- rawToChar(bytes)
  if (!grepl("\\+(?:ACI|AF8|AC0)-", contents, perl = TRUE))
    return(invisible(FALSE))

  decoded <- suppressWarnings(iconv(contents, from = "UTF-7", to = "UTF-8",
                                    sub = NA_character_))
  if (is.na(decoded))
    stop("Could not decode UTF-7 CSV: ", path,
         ". Preserve the file and export it as CSV UTF-8.")

  original_md5 <- unname(tools::md5sum(path))
  backup <- paste0(path, ".utf7_original_", original_md5, ".bak")
  if (!file.exists(backup) && !file.copy(path, backup, overwrite = FALSE))
    stop("Could not back up original manual CSV: ", path)
  if (!identical(unname(tools::md5sum(backup)), original_md5))
    stop("Existing UTF-7 backup does not match original manual CSV: ", backup)

  temporary <- tempfile(pattern = ".normalized_review_", tmpdir = dirname(path),
                        fileext = ".csv")
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  writeBin(charToRaw(enc2utf8(decoded)), temporary)
  parsed <- tryCatch(read_csv(temporary), error = function(e)
    stop("Decoded manual CSV cannot be parsed: ", path, ": ", conditionMessage(e)))
  required <- if (grepl("data_unresolved_rows_manual[.]csv$", path))
    c("life_id", "input_csv_row", "select_primary", "source_path", "source_md5")
  else c("source", "life_id", "issue", "field", "source_row",
         "candidate_life_id", "select", "resolution", "source_path", "source_md5")
  if (!all(required %in% names(parsed)))
    stop("Decoded manual CSV is missing columns: ",
         paste(setdiff(required, names(parsed)), collapse = ", "),
         ". Original file and backup remain available: ", path)

  # Rename on the same filesystem is atomic on Linux; retain the original
  # bytes in the .bak file even if the replacement fails.
  if (!file.rename(temporary, path))
    stop("Could not replace manual CSV with normalized UTF-8 version: ", path,
         ". Original backup: ", backup)
  cat("Normalized UTF-7 manual CSV to UTF-8 (backup: ", backup, ")\n", sep = "")
  invisible(TRUE)
}
manual_dir <- file.path(out, "manual")
normalize_manual_utf7(file.path(manual_dir, "data_unresolved_rows_manual.csv"))
normalize_manual_utf7(file.path(manual_dir, "data_match_issues_manual.csv"))

primary_template <- read_review_csv(file.path(out, "data_unresolved_rows_template.csv"))
issue_template <- read_review_csv(file.path(out, "data_match_issues_template.csv"))
manual_dir <- file.path(out, "manual")
primary <- read_review_csv(file.path(manual_dir, "data_unresolved_rows_manual.csv"))
issues <- read_review_csv(file.path(manual_dir, "data_match_issues_manual.csv"))

check_template <- function(reference, edited, editable, key_cols, label,
                           display_only = character()) {
  if (!identical(names(reference), names(edited)))
    stop(label, ": manual CSV columns differ from the generated template")
  immutable <- setdiff(names(reference), editable)
  if (!all(key_cols %in% immutable)) stop("Missing template key columns: ", label)
  if (!all(display_only %in% immutable))
    stop("Unknown display-only columns in template: ", label)
  make_key <- function(x) do.call(paste, c(x[key_cols], sep = "\035"))
  key_ref <- make_key(reference)
  key_edit <- make_key(edited)
  if (anyDuplicated(key_ref) || anyDuplicated(key_edit))
    stop(label, ": duplicate candidate keys detected")
  if (!setequal(key_ref, key_edit))
    stop(label, ": candidate rows changed. Review the refreshed template and ",
         "reconcile the persistent manual file before continuing")
  manual_order <- match(key_ref, key_edit)
  edited <- edited[manual_order, , drop = FALSE]
  rownames(edited) <- NULL
  rownames(reference) <- NULL
  for (field in immutable) {
    if (identical(reference[[field]], edited[[field]])) next
    # Compare the actual text values and report the first mismatching CSV row.
    a <- reference[[field]]
    b <- edited[[field]]
    different <- which(!( (is.na(a) & is.na(b)) |
      (!is.na(a) & !is.na(b) & a == b) ))
    if (!length(different)) next
    i <- different[[1L]]
    # The first line is the header, so data-frame row i is CSV line i+1.
    where <- paste0("template CSV line ", i + 1L,
                    ", manual CSV line ", manual_order[i] + 1L)
    key_description <- paste(vapply(key_cols, function(col)
      paste0(col, "=", edited[[col]][i]), character(1L)), collapse = ", ")
    detail <- paste0("; ", where, "; candidate: ", key_description,
                     "; template value: ", encodeString(a[i], quote = '"'),
                     "; manual value: ", encodeString(b[i], quote = '"'))
    if (field %in% display_only) {
      warning(label, ": display-only column differs: ", field,
              " (", length(different), " row(s))", detail,
              ". Source-row and fingerprint checks remain enabled.",
              call. = FALSE)
    } else {
      stop(label, ": read-only column changed: ", field, detail,
           ". Restore the template value; only change review columns")
    }
  }
  edited
}

primary <- check_template(primary_template, primary,
  editable = c("select_primary", "review_note"),
  key_cols = c("life_id", "input_csv_row"), label = "Primary EEG/IDmatch review")
issues <- check_template(issue_template, issues,
  editable = c("select", "resolution", "review_note"),
  key_cols = c("source", "life_id", "issue", "field", "source_row", "candidate_life_id"),
  label = "External-source review", display_only = "candidate_value")

check_hashes <- function(x) {
  if (!nrow(x)) return(invisible(NULL))
  if (!all(c("source_path", "source_md5") %in% names(x)))
    stop("Missing source fingerprints in manual template")
  for (path in unique(x$source_path)) {
    expected <- unique(x$source_md5[x$source_path == path])
    if (!length(path) || !nzchar(path) || length(expected) != 1L || !nzchar(expected))
      stop("Invalid source fingerprint for ", path)
    if (!file.exists(path)) stop("Reviewed source no longer exists: ", path)
    actual <- unname(tools::md5sum(path))
    if (is.na(actual) || actual != expected)
      stop("Reviewed source file changed since template generation: ", path,
           ". Re-run the generator and reconcile your choices")
  }
}
check_hashes(primary)
check_hashes(issues)

# The review may also depend on the participant inventory and selected mapping,
# even when the conflicting measurement originates from a different source file.
check_context <- function(x) {
  if (!nrow(x)) return(invisible(NULL))
  checks <- list(mapping_md5 = file.path(out, "IDmatch_variables_manual.csv"),
                 matching_input_md5 = file.path(out, "IDmatch_EEG.csv"))
  for (column in intersect(names(checks), names(x))) {
    path <- checks[[column]]
    if (!file.exists(path)) stop("Matching context disappeared: ", path)
    expected <- unique(x[[column]])
    if (length(expected) != 1L ||
        unname(tools::md5sum(path)) != expected)
      stop("Matching inputs changed since review: ", path,
           ". Refresh the templates and reconcile your choices")
  }
}
check_context(primary)
check_context(issues)

if (any(!primary$select_primary %in% c("", "0", "1")))
  stop("select_primary must be blank, 0 or 1")
if (any(!issues$select %in% c("", "0", "1")))
  stop("select must be blank, 0 or 1")
if (any(!issues$resolution %in% c("", "keep_missing")))
  stop("resolution must be blank or keep_missing")

# Exactly one chosen input row per participant that remained duplicated.
for (id in unique(primary$life_id)) {
  z <- primary[primary$life_id == id, , drop = FALSE]
  if (sum(z$select_primary == "1") != 1L)
    stop("Select exactly one primary IDmatch row (1) for ", id)
}

# Source problems can be resolved by selecting a candidate, or explicitly
# acknowledging that all corresponding measurements should remain missing.
case_key <- function(x) paste(x$source, x$life_id, x$issue, x$field,
                             sep = "\035")
for (key in unique(case_key(issues))) {
  z <- issues[case_key(issues) == key, , drop = FALSE]
  chosen <- which(z$select == "1")
  mark_missing <- any(z$resolution == "keep_missing")
  label <- paste(z$source[1L], z$life_id[1L], z$issue[1L], z$field[1L], sep = " / ")
  if (length(chosen) > 1L)
    stop("Choose at most one candidate per issue: ", label)
  if (length(chosen) == 1L) {
    if (mark_missing) stop("Cannot both select a candidate and keep_missing: ", label)
    if (!nzchar(z$source_row[chosen]) && !nzchar(z$candidate_life_id[chosen]))
      stop("The chosen issue row has no candidate: ", label)
  } else if (!mark_missing) {
    stop("Unreviewed issue: ", label, ". Select 1 or mark keep_missing")
  }
}

# A single GAD-7 questionnaire source row cannot belong to two participants.
selected_gad <- issues[issues$source == "GAD7" & issues$select == "1", , drop = FALSE]
if (anyDuplicated(selected_gad$source_row))
  stop("A GAD-7 source row has been selected for more than one participant")

# Review summary is written only after the validation and generator succeed.
picks <- which(primary$select_primary == "1")
primary_audit <- data.frame(
  category = rep("primary_row", length(picks)),
  source = rep("IDmatch_EEG.csv", length(picks)), life_id = primary$life_id[picks],
  issue = rep("unresolved_participant_row", length(picks)),
  field = rep("", length(picks)), source_row = primary$input_csv_row[picks],
  candidate_life_id = rep("", length(picks)),
  decision = rep("selected", length(picks)),
  note = primary$review_note[picks], stringsAsFactors = FALSE)
issue_audit <- issues[issues$select == "1" | issues$resolution == "keep_missing", ,
                      drop = FALSE]
if (nrow(issue_audit)) {
  issue_audit <- issue_audit[!duplicated(case_key(issue_audit)), , drop = FALSE]
  issue_audit <- data.frame(category = "external_source", source = issue_audit$source,
    life_id = issue_audit$life_id, issue = issue_audit$issue,
    field = issue_audit$field, source_row = issue_audit$source_row,
    candidate_life_id = issue_audit$candidate_life_id,
    decision = ifelse(issue_audit$select == "1", "selected", "keep_missing"),
    note = issue_audit$review_note, stringsAsFactors = FALSE)
} else issue_audit <- primary_audit[FALSE, , drop = FALSE]
audit <- rbind(primary_audit, issue_audit)

cat("Manual review validated: ", length(unique(primary$life_id)),
    " unresolved participant(s), ", length(unique(case_key(issues))),
    " external issue(s).\n", sep = "")
# Keep the active manual filenames stable, but archive a numbered, immutable
# snapshot of each successful resolution. The latest numbered snapshot is NOT
# used as input automatically; only manual/*.csv above are read.
history_root <- file.path(manual_dir, "history")
dir.create(history_root, recursive = TRUE, showWarnings = FALSE)
existing_runs <- basename(list.dirs(history_root, recursive = FALSE, full.names = TRUE))
existing_runs <- existing_runs[grepl("^run_[0-9]+$", existing_runs)]
run_numbers <- suppressWarnings(as.integer(sub("^run_", "", existing_runs)))
run_numbers <- run_numbers[!is.na(run_numbers)]
next_run <- if (length(run_numbers)) max(run_numbers) + 1L else 1L
run_dir <- file.path(history_root, sprintf("run_%03d", next_run))
if (file.exists(run_dir)) stop("History run already exists: ", run_dir)
staging_dir <- tempfile(pattern = ".resolution_in_progress_", tmpdir = history_root)
if (!dir.create(staging_dir)) stop("Cannot create resolution staging directory")
copy_snapshot <- function(source, destination, required = TRUE) {
  if (!file.exists(source)) {
    if (required) stop("Cannot archive missing file: ", source)
    return(invisible(FALSE))
  }
  if (!file.copy(source, destination, overwrite = FALSE))
    stop("Cannot archive: ", source, " -> ", destination)
  invisible(TRUE)
}
copy_snapshot(file.path(manual_dir, "data_unresolved_rows_manual.csv"), staging_dir)
copy_snapshot(file.path(manual_dir, "data_match_issues_manual.csv"), staging_dir)
previous_data <- file.path(life, "derivatives", "data.csv")
had_previous_data <- file.exists(previous_data)
if (had_previous_data)
  copy_snapshot(previous_data, file.path(staging_dir, "data_before_resolution.csv"))
restore_previous_data <- function() {
  if (had_previous_data) {
    if (!file.copy(file.path(staging_dir, "data_before_resolution.csv"),
                   previous_data, overwrite = TRUE))
      warning("Could not restore previous data.csv; backup is in ", staging_dir)
  } else if (file.exists(previous_data)) {
    if (!file.remove(previous_data))
      warning("Could not remove incomplete data.csv: ", previous_data)
  }
}
cat("Recomputing the full data table with reviewed selections...\n")
status <- system2(file.path(R.home("bin"), "Rscript"),
  args = c(shQuote(normalizePath(generator)), shQuote(normalizePath(out))),
  env = "LIFE_MATCHING_RESOLVE=1", stdout = "", stderr = "")
if (!identical(as.integer(status), 0L)) {
  restore_previous_data()
  stop("Data regeneration failed; previous data.csv restored when available. ",
       "Check the preceding messages. Manual snapshots remain at ", staging_dir)
}

final_path <- file.path(life, "derivatives", "data.csv")
final_data <- tryCatch(read_csv(final_path), error = function(e) e)
if (inherits(final_data, "error") || anyDuplicated(final_data$life_id)) {
  restore_previous_data()
  stop("Final data could not be validated; previous data.csv restored when available. ",
       "Manual snapshots remain at ", staging_dir)
}

write.csv(audit, file.path(out, "data_manual_resolution_audit.csv"),
          na = "", row.names = FALSE)

# Avoid archiving repeated runs when decisions, generated outputs and code are
# identical to the most recent successful resolution. Changes create a new run.
generator_hash <- unname(tools::md5sum(generator))
resolver_hash <- if (!is.na(resolver_script))
  unname(tools::md5sum(resolver_script)) else ""
latest_run <- if (length(run_numbers))
  file.path(history_root, sprintf("run_%03d", max(run_numbers))) else NULL
if (!is.null(latest_run)) {
  current_files <- c(
    data_unresolved_rows_manual.csv = file.path(manual_dir, "data_unresolved_rows_manual.csv"),
    data_match_issues_manual.csv = file.path(manual_dir, "data_match_issues_manual.csv"),
    data.csv = final_path,
    data_manual_resolution_audit.csv = file.path(out, "data_manual_resolution_audit.csv"),
    data_final_row_provenance.csv = file.path(out, "data_final_row_provenance.csv"),
    data_final_match_issues.csv = file.path(out, "data_final_match_issues.csv"),
    data_final_summary.txt = file.path(out, "data_final_summary.txt"))
  identical_file <- function(current, archived) {
    file.exists(current) && file.exists(archived) &&
      identical(unname(tools::md5sum(current)), unname(tools::md5sum(archived)))
  }
  saved_manifest <- file.path(latest_run, "manifest.txt")
  saved_lines <- if (file.exists(saved_manifest)) readLines(saved_manifest) else character()
  same_code <- paste0("generator_md5=", generator_hash) %in% saved_lines &&
    nzchar(resolver_hash) && paste0("resolver_md5=", resolver_hash) %in% saved_lines
  same_outputs <- all(vapply(names(current_files), function(name)
    identical_file(current_files[[name]], file.path(latest_run, name)), logical(1L)))
  if (same_code && same_outputs) {
    unlink(staging_dir, recursive = TRUE)
    cat("Final matched data: ", final_path, "\n", sep = "")
    cat("Manual decisions audit: ",
        file.path(out, "data_manual_resolution_audit.csv"), "\n", sep = "")
    cat("No new snapshot: decisions, code and outputs match ",
        basename(latest_run), "\n", sep = "")
    quit(save = "no", status = 0L)
  }
}
copy_snapshot(final_path, file.path(staging_dir, "data.csv"))
copy_snapshot(file.path(out, "data_manual_resolution_audit.csv"), staging_dir)
copy_snapshot(file.path(out, "data_final_row_provenance.csv"), staging_dir)
copy_snapshot(file.path(out, "data_final_match_issues.csv"), staging_dir)
copy_snapshot(file.path(out, "data_final_summary.txt"), staging_dir)
writeLines(c(paste0("resolution_run=", basename(run_dir)),
             paste0("resolved_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
             paste0("generator_md5=", generator_hash),
             paste0("resolver_md5=", resolver_hash),
             paste0("final_data_md5=", unname(tools::md5sum(final_path)))),
           file.path(staging_dir, "manifest.txt"))
if (!file.rename(staging_dir, run_dir))
  stop("Resolution succeeded but archive could not be finalized: ", staging_dir)
cat("Final matched data: ", final_path, "\n", sep = "")
cat("Manual decisions audit: ",
    file.path(out, "data_manual_resolution_audit.csv"), "\n", sep = "")
cat("Versioned resolution snapshot: ", run_dir, "\n", sep = "")
