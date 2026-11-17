# LIFE resting EEG inventory, ID matching and repeated-row diagnostics
# Run: Rscript /data/p_02865/code/matching_and_summaries/LIFE_matching.R
# Requires: readxl, openxlsx. Distinct, nonempty IDmatch rows are preserved (not merged).

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
# Normalize recording numbers appearing before "_Ruhe"
# (e.g., "_2_Ruhe" -> "_Ruhe_2").
corrected_stem <- sub(
  "(?i)_([0-9]+)_Ruhe$",
  "_Ruhe_\\1",
  corrected_stem,
  perl = TRUE
)
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
fields <- grep("^Table3[._]", names(idmatch), value = TRUE)
if (!length(fields)) stop("No Table3 fields found.")
find_field <- function(suffix) {
  hits <- grep(paste0("^Table3[._]", suffix, "$"), names(idmatch), value = TRUE)
  if (length(hits) != 1L) stop("Expected one Table3 field named ", suffix)
  hits
}
key_cols <- vapply(c("TEILNEHMER_SIC", "EEG_GRUPPE", "TEILNEHMER_GESCHLECHT_EEG"),
                   find_field, character(1))
# 2a. Recover ONLY EEG participants with no populated Table3 data.
# The original IDmatch.xlsx remains untouched. This block does not merge
# existing populated participant rows or alter any downstream matching logic.
fill_missing_idmatch <- function(idmatch, rec, fields, life, out) {
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("IDmatch recovery requires openxlsx: install.packages('openxlsx')")

  source_file <- file.path(life, "tabular/IDmatch.xlsx")
  recovered_file <- file.path(life, "tabular/IDmatch_recovered.xlsx")
  behavior_dir <- file.path(life, "tabular/Behavior/data")
  if (!dir.exists(behavior_dir)) stop("Behavior directory missing: ", behavior_dir)
  if (normalizePath(source_file, mustWork = TRUE) ==
      normalizePath(recovered_file, mustWork = FALSE))
    stop("Recovery output may not overwrite IDmatch.xlsx")

  # Work with the same text representation used by the CSV matching pipeline.
  # Only NEW cells are written to the recovered workbook; preexisting Excel
  # cells, types and other worksheets are not rewritten.
  idmatch[] <- lapply(idmatch, as.character)
  origin_ids <- clean_id(idmatch$Column1)
  source_fields <- sub("^Table3[._]", "", fields)
  normalize_header <- function(x) {
    x <- toupper(trimws(as.character(x)))
    sub("^TABLE3[._]", "", x)
  }
  canonical <- function(x) {
    x <- clean_id(x)
    x <- sub("_(RUHE|PILOT)(?:_[0-9]+)?$", "", x, perl = TRUE)
    x <- strip_trailing(x)
    sub("^LI", "", x)
  }
  unique_values <- function(x) unique(x[!is.na(x) & nzchar(x)])
  blank_table <- function(columns) as.data.frame(
    setNames(replicate(length(columns), character(), simplify = FALSE), columns),
    stringsAsFactors = FALSE, check.names = FALSE)
  audit <- function(x, name) save_csv(x, name)

  # Resolve EEG names against *original* Column1 values using exact match,
  # then a unique LI-prefix / recording-suffix correction.
  eeg_ids <- unique(clean_id(rec$raw_life_id))
  eeg_ids <- eeg_ids[!is.na(eeg_ids)]
  original_unique <- unique_values(origin_ids)
  # Resolve on unique canonical keys in a vectorized lookup.
  original_canonical <- canonical(original_unique)
  unique_original <- !duplicated(original_canonical) &
    !duplicated(original_canonical, fromLast = TRUE)
  canonical_to_original <- setNames(original_unique[unique_original],
                                    original_canonical[unique_original])
  resolved <- eeg_ids
  missing_exact <- !eeg_ids %in% original_unique
  alternative <- unname(canonical_to_original[canonical(eeg_ids[missing_exact])])
  resolved[which(missing_exact)] <- ifelse(is.na(alternative),
                                           eeg_ids[missing_exact], alternative)
  ambiguous_raw <- missing_exact & canonical(eeg_ids) %in%
    original_canonical[!unique_original]
  resolved[ambiguous_raw] <- NA_character_
  if (anyNA(resolved))
    warning("Ambiguous EEG-to-IDmatch IDs excluded from recovery: ",
            paste(eeg_ids[is.na(resolved)], collapse = " | "))
  all_owners <- unique_values(resolved)
  original_is_empty <- Reduce(`&`, lapply(idmatch[fields],
    function(x) is.na(cell_values(x))))
  target_ids <- all_owners[vapply(all_owners, function(id) {
    r <- which(origin_ids == id)
    !length(r) || all(original_is_empty[r])
  }, logical(1))]

  summary_cols <- c("life_id", "status", "original_excel_rows",
                    "recovered_excel_row", "new_participant_row",
                    "source_rows", "source_files", "fields_filled",
                    "conflicting_fields")
  values_cols <- c("life_id", "field", "value", "file", "sheet",
                   "excel_row", "excel_column")
  conflict_cols <- c("life_id", "field", "values", "source_locations")
  error_cols <- c("file", "sheet", "error")
  unmapped_cols <- c("life_id", "file", "sheet", "excel_column",
                     "source_header", "source_value")

  if (!length(target_ids)) {
    audit(blank_table(summary_cols), "IDmatch_recovery_summary.csv")
    audit(blank_table(values_cols), "IDmatch_recovery_sources.csv")
    audit(blank_table(conflict_cols), "IDmatch_recovery_conflicts.csv")
    audit(blank_table(error_cols), "IDmatch_recovery_errors.csv")
    audit(blank_table(unmapped_cols), "IDmatch_recovery_unmapped.csv")
    message("IDmatch recovery: no EEG participants need recovery.")
    return(idmatch)
  }
  cat("IDmatch recovery candidates:", length(target_ids), "\n")

  # Build identifier aliases for ALL EEG participants, not only targets.
  # This prevents a SIC / LI variant used by two people from being accepted.
  keys <- data.frame(owner = character(), key = character(), source = character())
  add_key <- function(owner, value, source) {
    owner <- clean_id(owner)
    value <- clean_id(value)
    ok <- !is.na(owner) & !is.na(value)
    if (any(ok)) keys <<- rbind(keys, data.frame(
      owner = owner[ok], key = value[ok], source = rep(source, sum(ok))))
  }
  add_variants <- function(owner, value, source) {
    add_key(owner, value, source)
    base <- canonical(value)
    add_key(owner, base, paste0(source, ":without_LI"))
    add_key(owner, ifelse(is.na(base), NA_character_, paste0("LI", base)),
            paste0(source, ":with_LI"))
  }
  add_variants(resolved[!is.na(resolved)], eeg_ids[!is.na(resolved)],
               "EEG_filename")
  add_variants(all_owners, all_owners, "Column1")

  # Column2 and TEILNEHMER_SIC can contain aliases that appear in Behavior.
  # Also reserve aliases of people outside the EEG subset to detect collisions.
  owner_canonical <- canonical(all_owners)
  unique_owners <- !duplicated(owner_canonical) &
    !duplicated(owner_canonical, fromLast = TRUE)
  canonical_to_owner <- setNames(all_owners[unique_owners],
                                 owner_canonical[unique_owners])
  resolve_owner <- function(id) {
    id <- clean_id(id)
    result <- id
    valid <- !is.na(id)
    exact <- match(id[valid], all_owners)
    mapped <- unname(canonical_to_owner[canonical(id[valid])])
    result[valid] <- ifelse(!is.na(exact), all_owners[pmax(exact, 1L)],
                            ifelse(is.na(mapped), id[valid], mapped))
    result
  }
  full_owner <- resolve_owner(origin_ids)
  sic_name <- fields[normalize_header(source_fields) == "TEILNEHMER_SIC"]
  if ("Column2" %in% names(idmatch))
    add_key(full_owner, idmatch$Column2, "Zusammenführen1:Column2")
  if (length(sic_name) == 1L)
    add_key(full_owner, idmatch[[sic_name]], "Zusammenführen1:SIC")

  # Aliases from the IDmatch sheet (the same mapping used by generate_data).
  fallback <- list()
  if ("IDmatch" %in% readxl::excel_sheets(source_file)) {
    alt_headers <- names(readxl::read_excel(source_file, sheet = "IDmatch",
                         n_max = 0, .name_repair = "minimal"))
    if (length(alt_headers) >= 2L) {
      alt_types <- rep("skip", length(alt_headers))
      alt_types[1:2] <- "text"
      alternate <- readxl::read_excel(source_file, sheet = "IDmatch",
        col_types = alt_types, .name_repair = "minimal")
      alias_life <- clean_id(alternate[[1L]])
      alias_value <- clean_id(alternate[[2L]])
      alias_owner <- resolve_owner(alias_life)
      add_key(alias_owner, alias_value, "IDmatch:column2")
      eligible <- !is.na(alias_value) & alias_owner %in% target_ids
      fallback <- split(alias_value[eligible], alias_owner[eligible])
    }
  }
  keys <- keys[!duplicated(keys[c("owner", "key")]), , drop = FALSE]
  owners_by_key <- split(keys$owner, keys$key)
  ambiguous <- names(owners_by_key)[vapply(owners_by_key,
    function(x) length(unique(x)) > 1L, logical(1))]
  keys <- keys[!keys$key %in% ambiguous, , drop = FALSE]
  keys <- keys[keys$owner %in% target_ids, , drop = FALSE]
  if (!nrow(keys)) stop("No unambiguous search keys for missing EEG participants")

  # Scan workbooks once, and accept ONLY exact ID-cell matches. For each
  # worksheet, find its header row by matching raw variable names to Table3.*
  # names. The column-to-column mapping is unique and never based on position.
  workbooks <- sort(list.files(behavior_dir, pattern = "\\.xlsx$",
                               full.names = TRUE, recursive = TRUE,
                               ignore.case = TRUE))
  workbooks <- workbooks[!grepl("_Datajoin", basename(workbooks),
                               ignore.case = TRUE) &
                           !startsWith(basename(workbooks), "~$")]
  if (!length(workbooks)) stop("No eligible Behavior XLSX files")
  field_lookup <- normalize_header(source_fields)
  if (anyDuplicated(field_lookup))
    stop("Ambiguous Table3 variable names after header normalization")
  candidates <- list()
  unmapped <- list()
  errors <- list()
  source_hit <- setNames(integer(length(target_ids)), target_ids)
  source_files <- setNames(vector("list", length(target_ids)), target_ids)
  add_error <- function(path, sheet, err) {
    errors[[length(errors) + 1L]] <<- data.frame(
      file = path, sheet = sheet, error = err)
  }
  for (path in workbooks) {
    relative <- substring(path, nchar(behavior_dir) + 2L)
    sheets <- tryCatch(readxl::excel_sheets(path), error = function(e) {
      add_error(relative, "", conditionMessage(e)); NULL
    })
    if (is.null(sheets)) next
    for (sheet in sheets) {
      data <- tryCatch(readxl::read_excel(path, sheet = sheet,
        col_names = FALSE, col_types = "text", .name_repair = "minimal",
        range = readxl::cell_limits(c(1L, 1L), c(NA_integer_, NA_integer_))),
        error = function(e) {
          add_error(relative, sheet, conditionMessage(e)); NULL
        })
      if (is.null(data) || !nrow(data) || !ncol(data)) next
      mat <- as.matrix(data)
      # Source files normally have headers in the first row; also allow
      # a few introductory rows and reject sheets without recognizable fields.
      header_scores <- vapply(seq_len(min(12L, nrow(mat))), function(r)
        sum(normalize_header(mat[r, ]) %in% field_lookup, na.rm = TRUE), integer(1))
      if (!length(header_scores) || max(header_scores) == 0L) next
      header_row <- which.max(header_scores)
      source_names <- normalize_header(mat[header_row, ])
      mapped <- match(source_names, field_lookup)
      # Duplicate variable headers in one sheet cannot be mapped safely.
      duplicated_header <- !is.na(mapped) &
        (duplicated(mapped) | duplicated(mapped, fromLast = TRUE))
      mapped[duplicated_header] <- NA_integer_
      field_columns <- which(!is.na(mapped))
      if (!length(field_columns)) next

      # Find ID keys in ANY column, as in data_search.R. Header rows do not
      # count as participant data; a source row must resolve to one person.
      nr <- nrow(mat)
      flat_keys <- clean_id(as.vector(mat))
      positions <- which(flat_keys %in% keys$key)
      if (!length(positions)) next
      excel_rows <- (positions - 1L) %% nr + 1L
      owners <- keys$owner[match(flat_keys[positions], keys$key)]
      valid <- excel_rows > header_row
      excel_rows <- excel_rows[valid]
      owners <- owners[valid]
      if (!length(excel_rows)) next
      for (r in unique(excel_rows)) {
        owners_here <- unique(owners[excel_rows == r])
        if (length(owners_here) != 1L) {
          add_error(relative, sheet, paste0("Ambiguous owners in Excel row ", r))
          next
        }
        id <- owners_here[1L]
        source_hit[id] <- source_hit[id] + 1L
        source_files[[id]] <- c(source_files[[id]], relative)
        for (c in which(is.na(mapped))) {
          value <- cell_values(mat[r, c])
          if (is.na(value)) next
          header <- cell_values(mat[header_row, c])
          if (is.na(header)) next
          unmapped[[length(unmapped) + 1L]] <- data.frame(
            life_id = id, file = relative, sheet = sheet,
            excel_column = c, source_header = header,
            source_value = value, stringsAsFactors = FALSE)
        }
        for (c in field_columns) {
          value <- cell_values(mat[r, c])
          if (is.na(value)) next
          candidates[[length(candidates) + 1L]] <- data.frame(
            life_id = id, field = fields[mapped[c]], value = value,
            file = relative, sheet = sheet, excel_row = r,
            excel_column = c, stringsAsFactors = FALSE)
        }
      }
    }
  }
  sources <- if (length(candidates)) do.call(rbind, candidates) else
    blank_table(values_cols)
  problems <- if (length(errors)) do.call(rbind, errors) else
    blank_table(error_cols)
  audit(sources, "IDmatch_recovery_sources.csv")
  audit(problems, "IDmatch_recovery_errors.csv")
  audit(if (length(unmapped)) do.call(rbind, unmapped) else
          blank_table(unmapped_cols), "IDmatch_recovery_unmapped.csv")

  # Fill each source field only if its nonempty source values agree.
  # If multiple source rows disagree, keep the field blank and audit it.
  conflicts <- list()
  updates <- list()
  summary_rows <- list()
  if (nrow(sources)) {
    by_field <- split(seq_len(nrow(sources)),
                      paste(sources$life_id, sources$field, sep = "\037"))
  } else by_field <- list()
  resolved_fields <- list()
  for (group in by_field) {
    values <- unique_values(sources$value[group])
    id <- sources$life_id[group[1L]]
    field <- sources$field[group[1L]]
    if (length(values) == 1L) {
      resolved_fields[[length(resolved_fields) + 1L]] <- data.frame(
        life_id = id, field = field, value = values)
    } else if (length(values) > 1L) {
      origins <- paste(unique(paste0(sources$file[group], " [",
        sources$sheet[group], "] row ", sources$excel_row[group])),
        collapse = " | ")
      conflicts[[length(conflicts) + 1L]] <- data.frame(
        life_id = id, field = field, values = paste(values, collapse = " | "),
        source_locations = origins)
    }
  }
  resolved_df <- if (length(resolved_fields)) do.call(rbind, resolved_fields) else
    data.frame(life_id = character(), field = character(), value = character())
  conflicts_df <- if (length(conflicts)) do.call(rbind, conflicts) else
    blank_table(conflict_cols)
  audit(conflicts_df, "IDmatch_recovery_conflicts.csv")

  # Populate the FIRST empty existing row; append a new row only when the
  # person is entirely absent from Zusammenführen1. Other empty duplicate
  # rows are left for the existing duplicate/empty-row cleanup.
  for (id in target_ids) {
    existing <- which(origin_ids == id)
    is_new <- !length(existing)
    if (is_new) {
      idmatch[nrow(idmatch) + 1L, ] <- NA_character_
      row <- nrow(idmatch)
      idmatch$Column1[row] <- id
      updates[[length(updates) + 1L]] <- data.frame(
        row = row, column = match("Column1", names(idmatch)), value = id)
    } else row <- existing[1L]

    # Supplement missing Column2 only from a unique IDmatch-sheet alias.
    # Do not choose arbitrarily if the sheet contains multiple aliases.
    aliases <- unique_values(fallback[[id]])
    if ("Column2" %in% names(idmatch) && is.na(cell_values(idmatch$Column2[row])) &&
        length(aliases) == 1L && !aliases %in% ambiguous) {
      idmatch$Column2[row] <- aliases
      updates[[length(updates) + 1L]] <- data.frame(
        row = row, column = match("Column2", names(idmatch)), value = aliases)
    }

    recovered <- resolved_df[resolved_df$life_id == id, , drop = FALSE]
    inserted <- 0L
    for (j in seq_len(nrow(recovered))) {
      field <- recovered$field[j]
      if (!is.na(cell_values(idmatch[[field]][row]))) next
      value <- recovered$value[j]
      idmatch[[field]][row] <- value
      updates[[length(updates) + 1L]] <- data.frame(
        row = row, column = match(field, names(idmatch)), value = value)
      inserted <- inserted + 1L
    }
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      life_id = id,
      status = if (inserted) "recovered" else "no_matching_columns_with_data",
      original_excel_rows = paste(existing + 1L, collapse = " | "),
      recovered_excel_row = row + 1L,
      new_participant_row = is_new,
      source_rows = source_hit[[id]],
      source_files = paste(unique(source_files[[id]]), collapse = " | "),
      fields_filled = inserted,
      conflicting_fields = sum(conflicts_df$life_id == id),
      stringsAsFactors = FALSE)
  }
  recovery_summary <- do.call(rbind, summary_rows)
  audit(recovery_summary, "IDmatch_recovery_summary.csv")

  # Save a NEW workbook with all original sheets. Write ONLY newly recovered
  # cells; leave original cells and all other participant rows untouched.
  # Group adjacent columns to reduce Excel writer overhead.
  wb <- openxlsx::loadWorkbook(source_file)
  changed <- if (length(updates)) do.call(rbind, updates) else
    data.frame(row = integer(), column = integer(), value = character())
  if (nrow(changed)) {
    changed <- changed[order(changed$row, changed$column), , drop = FALSE]
    for (r in unique(changed$row)) {
      piece <- changed[changed$row == r, , drop = FALSE]
      split_at <- cumsum(c(TRUE, diff(piece$column) != 1L))
      for (run in split(piece, split_at))
        openxlsx::writeData(wb, sheet = "Zusammenführen1",
          x = matrix(run$value, nrow = 1L), startRow = r + 1L,
          startCol = run$column[1L], colNames = FALSE, rowNames = FALSE)
    }
  }
  openxlsx::saveWorkbook(wb, recovered_file, overwrite = TRUE)
  cat("IDmatch recovery:", sum(recovery_summary$status == "recovered"), "/",
      length(target_ids), "participants received Table3 data;",
      sum(recovery_summary$new_participant_row), "rows appended;",
      nrow(conflicts_df), "conflicting fields left blank.\n")
  cat("Recovered workbook:", recovered_file, "\n")
  if (nrow(problems)) warning("Behavior workbook read/ID errors: ", nrow(problems),
                               "; see IDmatch_recovery_errors.csv")
  if (any(recovery_summary$status != "recovered"))
    warning("Some EEG IDs still lack Table3 data; see IDmatch_recovery_summary.csv")
  idmatch
}

idmatch <- fill_missing_idmatch(idmatch, rec, fields, life, out)
ids_all <- clean_id(idmatch$Column1)

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
