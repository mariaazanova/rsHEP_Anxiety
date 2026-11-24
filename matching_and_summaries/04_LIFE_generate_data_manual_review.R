# Run after LIFE_matching.R to prepare manual review:
#   Rscript LIFE_generate_data_manual_review.R
# After marking the generated files in matching_and_summaries/manual/, run:
#   Rscript LIFE_resolve_matching.R
# General participant-level data ONLY; pipeline HEP/ECG and preprocessing
# statuses are joined AFTER the common exclusions by LIFE_pipeline_matching.R.
# Requires readxl:
#   install.packages("readxl")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript 04_LIFE_generate_data_manual_review.R [matching_output_dir]")
life <- "/data/p_02865/data/LIFE"
out <- if (length(args)) args[[1L]] else file.path(life, "derivatives/matching_and_summaries")
if (!dir.exists(out)) stop("Missing matching output directory: ", out)

read_csv <- function(path) read.csv(path, colClasses = "character",
  check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM",
  na.strings = "")
# Spreadsheet UTF-7 exports may encode underscores as +AF8-.
# Decode only review CSVs; leave all other source-data readers unchanged.
read_review_csv <- function(path) {
  x <- read_csv(path)
  decode_review <- function(v) gsub("+AF8-", "_", v, fixed = TRUE)
  names(x) <- decode_review(names(x))
  if (anyDuplicated(names(x)))
    stop("Duplicate column names after decoding review CSV: ", path)
  x[] <- lapply(x, decode_review)
  x
}

optional_csv <- function(name) {
  path <- file.path(out, name)
  if (file.exists(path)) read_csv(path) else NULL
}
write_csv <- function(x, name) write.csv(x, file.path(out, name),
                                          na = "", row.names = FALSE)
normalize <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
normalize_id <- function(x) {
  x <- toupper(normalize(x))
  x[!is.na(x) & x == "NA"] <- NA_character_
  x
}
show <- function(x) if (length(x) == 1L && !is.na(x)) as.character(x) else ""
joined <- function(x, sep = " | ") {
  x <- unique(as.character(x[!is.na(x)]))
  paste(x[nzchar(x)], collapse = sep)
}
empty_table <- function(cols) {
  z <- setNames(rep(list(character()), length(cols)), cols)
  as.data.frame(z, check.names = FALSE, stringsAsFactors = FALSE)
}
bind_events <- function(lst, cols) {
  if (!length(lst)) return(empty_table(cols))
  result <- do.call(rbind, lst)
  rownames(result) <- NULL
  result[, cols, drop = FALSE]
}
conflicting <- function(a, b, cols) {
  cols <- intersect(cols, intersect(names(a), names(b)))
  cols[!is.na(a[cols]) & !is.na(b[cols]) & a[cols] != b[cols]]
}
changes <- function(a, b, cols) {
  cols <- intersect(cols, intersect(names(a), names(b)))
  cols[(is.na(a[cols]) != is.na(b[cols])) |
       (!is.na(a[cols]) & !is.na(b[cols]) & a[cols] != b[cols])]
}
contains <- function(x, pattern) grepl(pattern, x, perl = TRUE, ignore.case = TRUE)
yes <- function(x) !is.na(x) & toupper(trimws(as.character(x))) %in%
  c("TRUE", "T", "YES", "Y", "1")

# ---- Manual review: preparation pass versus validated resolution pass --------
resolution_mode <- identical(Sys.getenv("LIFE_MATCHING_RESOLVE", ""), "1")
manual_dir <- file.path(out, "manual")
manual_primary_path <- file.path(manual_dir, "data_unresolved_rows_manual.csv")
manual_issues_path <- file.path(manual_dir, "data_match_issues_manual.csv")
primary_review <- issue_review <- NULL
consumed_decisions <- character()
if (resolution_mode) {
  if (!file.exists(manual_primary_path) || !file.exists(manual_issues_path))
    stop("Missing manual decisions in ", manual_dir)
  primary_review <- read_review_csv(manual_primary_path)
  issue_review <- read_review_csv(manual_issues_path)
}
review_choice <- function(source, id, issue, field = "") {
  if (!resolution_mode || !nrow(issue_review)) return(NULL)
  blank <- function(x) ifelse(is.na(x), "", x)
  ii <- which(blank(issue_review$source) == source &
    blank(issue_review$life_id) == id &
    blank(issue_review$issue) == issue &
    blank(issue_review$field) == field & yes(issue_review$select))
  if (!length(ii)) return(NULL)
  if (length(ii) != 1L) stop("Multiple manual choices for ", source, "/", id, "/", issue)
  consumed_decisions <<- unique(c(consumed_decisions,
    paste(source, id, issue, field, sep = "\035")))
  issue_review[ii, , drop = FALSE]
}
review_source_row <- function(source, id, issue, field = "") {
  choice <- review_choice(source, id, issue, field)
  if (is.null(choice)) return(NA_integer_)
  row <- suppressWarnings(as.integer(choice$source_row))
  if (is.na(row) || row < 2L) stop("Invalid reviewed source row for ", id)
  row - 1L  # source row numbers include the header, just like input_csv_row
}

# ---- 1. Mapping: keep only requested variables, plus required core fields --
raw <- optional_csv("IDmatch_EEG.csv")
manual <- optional_csv("IDmatch_variables_manual.csv")
if (is.null(raw) || is.null(manual))
  stop("Missing IDmatch_EEG.csv or idmatch_variables_manual.csv in ", out)
# Some spreadsheet programs save underscores as UTF-7 +AF8- strings.
decode <- function(x) gsub("+AF8-", "_", x, fixed = TRUE)
names(manual) <- decode(names(manual))
manual[] <- lapply(manual, decode)
required_map <- c("original_column", "keep", "new_name")
if (!all(required_map %in% names(manual)))
  stop("Mapping lacks: ", joined(setdiff(required_map, names(manual))))
manual[] <- lapply(manual, normalize)
manual$keep[is.na(manual$keep)] <- "0"
if (any(!manual$keep %in% c("0", "1"))) stop("keep must be 0 or 1")
if (anyDuplicated(names(raw))) stop("Duplicate columns in IDmatch_EEG.csv")
identifier_sources <- c("Column1", "Column2", "Table3.TEILNEHMER_SIC")
if (!all(identifier_sources %in% names(raw)))
  stop("Missing one or more required identifier columns: ",
       joined(setdiff(identifier_sources, names(raw))))
forced <- c("Table3.TEILNEHMER_SIC", "Table3.EEG_GRUPPE", "Table3.EEG_DATUM")
selection <- manual[manual$keep == "1" |
                      (!is.na(manual$original_column) &
                         manual$original_column %in% forced), , drop = FALSE]
if (anyNA(selection$original_column) || anyNA(selection$new_name))
  stop("Selected variables must have both original_column and new_name")
if (anyDuplicated(selection$original_column) || anyDuplicated(selection$new_name))
  stop("Selected original_column and new_name must each be unique")
if (any(selection$original_column %in% c("Column1", "Column2")))
  stop("Remove Column1/Column2 from keep=1 mapping; they are always retained")
if (any(selection$new_name %in% c("life_id", "Column2")))
  stop("life_id and Column2 are reserved output names")
if (!all(forced %in% selection$original_column))
  stop("SIC, EEG_GRUPPE and EEG_DATUM must exist in manual mapping")
missing_sources <- setdiff(selection$original_column, names(raw))
if (length(missing_sources))
  stop("Variables not in IDmatch_EEG.csv: ", joined(missing_sources))
name_for <- function(original) {
  idx <- match(original, selection$original_column)
  if (is.na(idx)) NA_character_ else selection$new_name[idx]
}
sic <- name_for("Table3.TEILNEHMER_SIC")
grp <- name_for("Table3.EEG_GRUPPE")
dt <- name_for("Table3.EEG_DATUM")
age <- name_for("Table3.age_EEG")
if (any(grepl(paste0("^(?:", grp, "|", dt, ")[2-9][0-9]*$"),
              selection$new_name)))
  stop("Names such as group2/date_eeg2 are reserved for alternative sessions")

values <- data.frame(life_id = normalize_id(raw$Column1),
  Column2 = normalize(raw$Column2), raw[selection$original_column],
  check.names = FALSE, stringsAsFactors = FALSE)
names(values)[-(1:2)] <- selection$new_name
values[] <- lapply(values, normalize)
if (anyNA(values$life_id)) stop("Blank participant life_id: cannot reconcile")
if (!nrow(values)) stop("IDmatch_EEG.csv has no rows")
values <- values[, c("life_id", "Column2", sic,
                     setdiff(names(values), c("life_id", "Column2", sic))), drop = FALSE]
orig <- setNames(selection$original_column, selection$new_name)

# Core EEG evidence: actual session date/age and EEG_RUHE/Befund fields.
# EEG_SUBSTANZ variables are acquisition-related but do not by themselves
# prove which EEG recording a tabular row belongs to.
core_eeg <- unique(c(dt, age[!is.na(age)], names(orig)[
  contains(orig, "^Table3[.]EEG_RUHE[_]|^Table3[.]BEFUND_EEG_")]))
core_eeg <- intersect(core_eeg, names(values))
core_eeg <- setdiff(core_eeg, grp)
# All selected EEG-associated measures must follow the primary EEG row.
eeg_session <- unique(c(grp, dt, age[!is.na(age)], names(orig)[
  contains(orig, "^Table3[.]EEG_|^Table3[.]BEFUND_EEG_")]))
eeg_session <- intersect(eeg_session, names(values))
# VIGALL must come ONLY from the primary tabular row, even for close matches.
vigall <- names(orig)[contains(orig, "VIGALL") |
                        contains(names(orig), "^vigall")]
vigall <- intersect(vigall, names(values))
eeg_session <- setdiff(eeg_session, vigall)
identifiers <- c("life_id", "Column2", sic)

# Columns empty in the entire input may be dropped from the final CSV.
all_empty <- setdiff(names(values)[vapply(values, function(x)
  all(is.na(x)), logical(1))], identifiers)


# ---- 2. Prior matching diagnostics and EEG recording paths (NO dates) ------
inventory <- optional_csv("EEG_recordings_inventory.csv")
participants <- optional_csv("EEG_participants_inventory.csv")
prior_duplicates <- optional_csv("IDmatch_duplicates.csv")
prior_issues <- optional_csv("EEG_file_issues.csv")
if (is.null(inventory)) {
  warning("EEG_recordings_inventory.csv missing; recording paths unavailable")
  inventory <- data.frame(life_id = character(), vhdr_file = character(),
                          vhdr_path = character(), selected = character(),
                          ready = character(), suffix = character())
}
if (!all(c("life_id", "vhdr_file", "vhdr_path") %in% names(inventory)))
  stop("EEG_recordings_inventory.csv must have life_id, vhdr_file, vhdr_path")
for (field in c("selected", "ready", "suffix"))
  if (!field %in% names(inventory))
    inventory[[field]] <- rep(NA_character_, nrow(inventory))
inventory$life_id <- normalize_id(inventory$life_id)
inventory$is_2 <- contains(inventory$vhdr_file, "_2[.]vhdr$")
inventory$is_selected <- yes(inventory$selected)
inventory$is_ready <- yes(inventory$ready)

# Previously rejected _2 filenames are included in the *review table* only.
# They are not silently promoted to a validated, selected recording.
extra_files <- data.frame(life_id = character(), vhdr_file = character(),
  vhdr_path = character(), selected = character(), ready = character(),
  suffix = character(), is_2 = logical(), is_selected = logical(),
  is_ready = logical(), stringsAsFactors = FALSE)
if (!is.null(prior_issues) &&
    all(c("issue", "filename", "path") %in% names(prior_issues))) {
  unresolved_files <- prior_issues[prior_issues$issue == "unrecognized_vhdr" &
    contains(prior_issues$filename, "_2[.]vhdr$"), , drop = FALSE]
  for (i in seq_len(nrow(unresolved_files))) {
    f <- unresolved_files$filename[i]
    # e.g. LI12345678_Pilot_2.vhdr -> LI12345678
    proposed <- normalize_id(sub("(?i)_([^_]+)_2[.]vhdr$", "", f, perl = TRUE))
    if (!(proposed %in% values$life_id)) next
    if (any(inventory$life_id == proposed & inventory$vhdr_file == f)) next
    extra_files[nrow(extra_files) + 1L, ] <- list(
      proposed, f, unresolved_files$path[i], "FALSE", "FALSE", "2",
      TRUE, FALSE, FALSE)
  }
}
recordings <- rbind(inventory[, names(extra_files), drop = FALSE], extra_files)


# Manual resolution is only needed for otherwise ambiguous rows.
# No association between _2 and a particular tabular row is inferred from the
# filename alone; the user can explicitly specify that row here.
# Older date-derived overrides are deliberately ignored to preserve unresolved rows.
# Set TRUE only when data_primary_row_overrides.csv has been MANUALLY checked.
use_manual_overrides <- resolution_mode
overrides <- if (use_manual_overrides)
  primary_review[yes(primary_review$select_primary), c("life_id", "input_csv_row"),
                 drop = FALSE] else NULL
if (!is.null(overrides)) {
  needed <- c("life_id", "input_csv_row")
  if (!all(needed %in% names(overrides)))
    stop("data_primary_row_overrides.csv must have life_id,input_csv_row")
  overrides$life_id <- normalize_id(overrides$life_id)
  overrides$input_csv_row <- suppressWarnings(as.integer(overrides$input_csv_row))
  if (anyNA(overrides$input_csv_row) || anyNA(overrides$life_id) ||
      anyDuplicated(overrides$life_id))
    stop("Override IDs and row numbers must be present and unique")
  if (any(overrides$input_csv_row < 2L | overrides$input_csv_row > nrow(values) + 1L))
    stop("Invalid input_csv_row in override; first data row is CSV row 2")
  if (any(values$life_id[overrides$input_csv_row - 1L] != overrides$life_id))
    stop("An override input_csv_row does not belong to its life_id")
}
override_row <- function(id) {
  if (is.null(overrides)) return(NA_integer_)
  i <- match(id, overrides$life_id)
  if (is.na(i)) NA_integer_ else overrides$input_csv_row[i]
}

# ---- 3. Row reconciliation -------------------------------------------------
cell_actions <- integer()
record <- function(i) list(id = values$life_id[i],
  vals = unlist(values[i, -1L, drop = FALSE], use.names = TRUE),
  primary = as.integer(i + 1L), sources = as.integer(i + 1L),
  alts = list(), reason = "single_input_row")
# Only count decisions; do not save columns/values in verbose audit tables.
log_cell <- function(main, other, col, before, supplied, action, rule) {
  if (!action %in% names(cell_actions)) cell_actions[action] <<- 1L
  else cell_actions[action] <<- cell_actions[action] + 1L
}
source_ids <- function(r) joined(sort(r$sources))

has_core <- function(r) any(!is.na(r$vals[core_eeg]))
core_count <- function(r) sum(!is.na(r$vals[core_eeg]))
choose_best <- function(rows, id) {
  ov <- override_row(id)
  if (!is.na(ov)) {
    w <- which(vapply(rows, function(r) ov %in% r$sources, logical(1)))
    if (length(w) == 1L) return(w)
    if (length(w) > 1L) stop("Reviewed override appears in several records for ", id)
    # Another incompatible identifier partition is handled separately below.
  }
  # Only used when every row is already mutually compatible outside
  # group/date; tie-breaking here cannot discard a conflicting EEG/VIGALL value.
  available_vig <- vigall[vapply(vigall, function(col)
    any(vapply(rows, function(r) !is.na(r$vals[[col]]), logical(1L))), logical(1L))]
  joint <- function(r) {
    has_core(r) && !is.na(r$vals[[dt]]) &&
      (is.na(age) || !is.na(r$vals[[age]])) &&
      length(available_vig) > 0L && all(!is.na(r$vals[available_vig]))
  }
  scores <- t(vapply(rows, function(r) c(
    joint = as.integer(joint(r)),
    has_eeg = as.integer(has_core(r)),
    date_age = as.integer(!is.na(r$vals[[dt]]) &&
                            (is.na(age) || !is.na(r$vals[[age]]))),
    vig = sum(!is.na(r$vals[vigall])),
    core = core_count(r),
    all = sum(!is.na(r$vals)),
    input = -r$primary), numeric(7L)))
  order(-scores[, "joint"], -scores[, "has_eeg"],
        -scores[, "date_age"], -scores[, "vig"],
        -scores[, "core"], -scores[, "all"], -scores[, "input"])[1L]
}

# Choose a primary EEG/VIGALL row only when it has a *clear* advantage.
# No _2 filename or recording metadata is used to select a tabular row.
select_primary <- function(rows, id) {
  ov <- override_row(id)
  if (!is.na(ov) && any(vapply(rows, function(r) ov %in% r$sources, logical(1L))))
    return(list(index = choose_best(rows, id),
                reason = "reviewed_primary_override"))
  has_eeg <- vapply(rows, has_core, logical(1L))
  if (sum(has_eeg) == 1L)
    return(list(index = which(has_eeg), reason = "only_row_with_core_eeg_information"))

  available_vig <- vigall[vapply(vigall, function(col)
    any(vapply(rows, function(r) !is.na(r$vals[[col]]), logical(1L))), logical(1L))]
  full_joint <- vapply(rows, function(r)
    has_core(r) && !is.na(r$vals[[dt]]) &&
      (is.na(age) || !is.na(r$vals[[age]])) &&
      length(available_vig) > 0L && all(!is.na(r$vals[available_vig])), logical(1L))
  if (sum(full_joint) == 1L)
    return(list(index = which(full_joint),
                reason = "only_row_with_complete_eeg_and_vigall"))

  # If one record covers *every* nonmissing EEG/VIGALL field found in the
  # other candidates (and at least one more), it wins. Conflicting VALUES in
  # fields covered by both are logged, never copied from the secondary row.
  signal <- unique(c(core_eeg, vigall))
  cover <- lapply(rows, function(r) !is.na(r$vals[signal]))
  dominated <- vapply(seq_along(rows), function(i) {
    if (sum(has_eeg) && !has_eeg[i]) return(FALSE)
    all(vapply(seq_along(rows), function(j)
      all(!cover[[j]] | cover[[i]]), logical(1L))) &&
      any(vapply(seq_along(rows), function(j)
        any(cover[[i]] & !cover[[j]]), logical(1L)))
  }, logical(1L))
  if (sum(dominated) == 1L)
    return(list(index = which(dominated),
                reason = "dominant_eeg_vigall_completeness"))

  if (!any(has_eeg)) {
    with_vig <- vapply(rows, function(r)
      any(!is.na(r$vals[vigall])), logical(1L))
    if (sum(with_vig) == 1L)
      return(list(index = which(with_vig),
                  reason = "only_row_with_vigall_no_eeg"))
    return(list(index = NA_integer_, reason = "no_clear_eeg_vigall_primary"))
  }
  list(index = NA_integer_, reason = "competing_eeg_vigall_rows_no_clear_primary")
}

# Promote all alternate date/group pairs to separate columns: group2/date_eeg2,
# group3/date_eeg3, ... . Keep GROUP and DATE from a single same source row.
add_alt <- function(main, other) {
  g <- other$vals[[grp]]
  d <- other$vals[[dt]]
  # Only a genuine nonmissing disagreement needs an alternative visit pair.
  # A second row with missing group/date is NOT another recorded session.
  if (!length(conflicting(main$vals, other$vals, c(grp, dt)))) return(main)
  dup <- vapply(main$alts, function(a)
    identical(g, a$group) && identical(d, a$date), logical(1))
  if (any(dup)) return(main)
  main$alts[[length(main$alts) + 1L]] <-
    list(group = g, date = d, source_rows = source_ids(other))
  log_cell(main, other, paste0(grp, "/", dt),
           paste(show(main$vals[[grp]]), show(main$vals[[dt]]), sep = " / "),
           paste(show(g), show(d), sep = " / "),
           "saved_alternative_group_date", "preserve_secondary_visit")
  main
}

# Absorb selected rows without creating synthetic EEG/VIGALL combinations.
# Different nonmissing values are retained from primary and logged.
absorb <- function(main, other, rule, keep_alt = TRUE) {
  if (length(conflicting(main$vals, other$vals, c("Column2", sic))))
    stop("Internal error: merging different identifiers for ", main$id)
  if (keep_alt) main <- add_alt(main, other)
  different_visits <- length(conflicting(main$vals, other$vals, c(grp, dt, age))) > 0L
  for (col in names(main$vals)) {
    a <- main$vals[[col]]
    b <- other$vals[[col]]
    if (is.na(a) && is.na(b)) next
    if (!is.na(a) && !is.na(b)) {
      if (a != b) log_cell(main, other, col, a, b,
                           "kept_primary_conflict", rule)
      next
    }
    if (!is.na(a)) next
    if (col %in% vigall) {
      log_cell(main, other, col, a, b, "skipped_secondary_vigall", rule)
    } else if (col %in% eeg_session && different_visits) {
      log_cell(main, other, col, a, b, "skipped_secondary_eeg_session", rule)
    } else if (col %in% c(grp, dt) && length(main$alts)) {
      # Do not mix alternate visit date/group into the primary paired fields.
      log_cell(main, other, col, a, b, "kept_primary_visit_pair", rule)
    } else {
      main$vals[[col]] <- b
      log_cell(main, other, col, a, b, "filled_from_secondary", rule)
    }
  }
  main$sources <- c(main$sources, other$sources)
  main
}

# Distinct SICs (or conflicting Column2) cannot automatically share a row.
# If there are several known values, missing-identifier rows remain separate.
partition_on <- function(rows, col) {
  found <- unique(vapply(rows, function(r) show(r$vals[[col]]), character(1)))
  found <- found[nzchar(found)]
  if (length(found) < 2L) return(list(rows))
  keys <- vapply(rows, function(r) {
    a <- show(r$vals[[col]])
    if (nzchar(a)) paste0("known:", a) else "MISSING:REVIEW"
  }, character(1))
  unname(split(rows, factor(keys, levels = unique(keys))))
}

# All rows of a partition can be merged via (a) ONLY when no pair has
# conflicting nonmissing cells outside EEG group/date. Missing vs observed
# is complementary, not a conflict. Also covers identical duplicates.
can_merge_group_date <- function(rows) {
  if (length(rows) < 2L) return(TRUE)
  cols <- setdiff(names(rows[[1L]]$vals), c(grp, dt))
  for (i in seq_len(length(rows) - 1L))
    for (j in seq.int(i + 1L, length(rows)))
      if (length(conflicting(rows[[i]]$vals, rows[[j]]$vals, cols)))
        return(FALSE)
  TRUE
}

ids <- unique(values$life_id)
final <- list()
for (id in ids) {
  source_records <- lapply(which(values$life_id == id), record)
  sic_parts <- partition_on(source_records, sic)
  for (sp in sic_parts) for (part in partition_on(sp, "Column2")) {
    if (length(part) == 1L) {
      final[[length(final) + 1L]] <- part[[1L]]
      next
    }
    # Rule A: no conflicts other than group/date, even when some cells are NA.
    # Pick the EEG/VIGALL-complete row; NEVER backfill VIGALL from another row.
    if (can_merge_group_date(part)) {
      best <- choose_best(part, id)
      main <- part[[best]]
      main$reason <- if (is.na(override_row(id)) ||
        !any(vapply(part, function(r) override_row(id) %in% r$sources, logical(1L))))
        "only_group_date_disagreement_or_complementary" else "reviewed_primary_override"
      for (j in setdiff(seq_along(part), best)) {
        other <- part[[j]]
        rule <- if (!length(changes(main$vals, other$vals, names(main$vals))))
          "exact_duplicate" else "group_date_or_complementary"
        main <- absorb(main, other, rule)
      }
      final[[length(final) + 1L]] <- main
      next
    }
    # Rule B: prefer a row with EEG and VIGALL acquired together.
    # If no uniquely superior row exists, retain the competing records.
    decision <- select_primary(part, id)
    preferred <- decision$index
    reason <- decision$reason
    if (is.na(preferred)) {
      for (r in part) {
        r$reason <- paste0("UNRESOLVED:", reason)
        final[[length(final) + 1L]] <- r
      }
      next
    }
    main <- part[[preferred]]
    main$reason <- reason
    for (j in setdiff(seq_along(part), preferred))
      main <- absorb(main, part[[j]], reason)
    final[[length(final) + 1L]] <- main
  }
}

# Exactly once: each input CSV data row must appear in the final provenance.
final <- final[order(vapply(final, function(r) min(r$sources), integer(1)))]
coverage <- unlist(lapply(final, `[[`, "sources"), use.names = FALSE)
if (!identical(sort(as.integer(coverage)), seq_len(nrow(values)) + 1L))
  stop("AUDIT ERROR: an original input row was lost or used more than once")

# After auditing every original row, a manual primary choice may discard OTHER
# identifier partitions. Never absorb a row with conflicting SIC/Column2.
manual_dropped <- list()
manual_dropped_n <- 0L
if (resolution_mode && !is.null(overrides) && nrow(overrides)) {
  for (id in overrides$life_id) {
    positions <- which(vapply(final, function(r) r$id == id, logical(1L)))
    if (length(positions) < 2L) next
    chosen <- positions[vapply(final[positions], function(r)
      override_row(id) %in% r$sources, logical(1L))]
    if (length(chosen) != 1L) stop("Manual primary row cannot uniquely resolve ", id)
    discarded <- setdiff(positions, chosen)
    for (j in discarded) {
      manual_dropped_n <- manual_dropped_n + length(final[[j]]$sources)
      manual_dropped[[length(manual_dropped) + 1L]] <-
        data.frame(life_id = id, selected_input_csv_row = override_row(id),
          discarded_input_csv_rows = source_ids(final[[j]]),
          reason = "different_identifier_partition_not_merged", stringsAsFactors = FALSE)
    }
    final <- final[-discarded]
  }
}

# ---- 4. Assemble selected rows (no output yet) -----------------------------
output_cols <- setdiff(names(values), all_empty)
if (!all(identifiers %in% output_cols)) stop("Identifiers must never be dropped")
max_alt <- max(vapply(final, function(r) length(r$alts), integer(1L)))
alt_cols <- if (max_alt) unlist(lapply(seq_len(max_alt), function(k)
  c(paste0(grp, k + 1L), paste0(dt, k + 1L)))) else character()
if (any(alt_cols %in% output_cols)) stop("Alternative names collide with existing fields")
result_data <- as.data.frame(setNames(rep(list(rep(NA_character_, length(final))),
                                length(c(output_cols, alt_cols))),
                                c(output_cols, alt_cols)), check.names = FALSE)
for (i in seq_along(final)) {
  r <- final[[i]]
  result_data$life_id[i] <- r$id
  for (col in setdiff(output_cols, "life_id"))
    result_data[[col]][i] <- r$vals[[col]]
  if (length(r$alts)) for (j in seq_along(r$alts)) {
    result_data[[paste0(grp, j + 1L)]][i] <- r$alts[[j]]$group
    result_data[[paste0(dt, j + 1L)]][i] <- r$alts[[j]]$date
  }
}
# Remove fully empty alternates, too.
if (length(alt_cols)) {
  alt_empty <- alt_cols[vapply(result_data[alt_cols], function(x)
    all(is.na(x)), logical(1L))]
  result_data <- result_data[, setdiff(names(result_data), alt_empty), drop = FALSE]
}

# Remove the entire requested EP-EW interval in the manually defined order.
# Include any variables newly selected between these two endpoints in future maps.
start <- which(manual$new_name == "coronary_disease_angina_pectoris_2")
end <- which(manual$new_name == "diabetes_question_2")
if (length(start) != 1L || length(end) != 1L || end < start)
  stop("Cannot identify EP-EW boundaries in idmatch_variables_manual.csv")
block_sources <- manual$original_column[seq.int(start, end)]
block <- selection$new_name[selection$original_column %in% block_sources]
result_data <- result_data[, !names(result_data) %in% block, drop = FALSE]

# non-unique describes the FINAL table, not repetitions in the initial IDmatch.
# This is a deliberately non-syntactic CSV header; use [[...]] in R.
id_counts <- table(result_data$life_id)
result_data[["non-unique"]] <- as.integer(id_counts[result_data$life_id] > 1L)
result_data <- result_data[, c("life_id", "Column2", sic, "non-unique",
  setdiff(names(result_data), c("life_id", "Column2", sic, "non-unique"))), drop = FALSE]

# ---- 5. Safe matching of participant IDs across external data sources -----
# Exact matching first. For unmatched IDs, allow LIFE's documented LI-prefix /
# trailing suffix corrections, but NEVER guess when a fallback is ambiguous.
main_ids <- unique(result_data$life_id)
canon_id <- function(x) {
  x <- normalize_id(x)
  x <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", x, perl = TRUE)
  x <- sub("[_-]+2$", "", x, perl = TRUE)
  x <- sub("[_-]+$", "", x, perl = TRUE)
  sub("^LI", "", x)
}
canonical_main <- canon_id(main_ids)
issues <- list()
add_issue <- function(source, participant, type, detail) {
  issues[[length(issues) + 1L]] <<- data.frame(
    source = source, life_id = participant, issue = type,
    details = detail, stringsAsFactors = FALSE)
}
resolve_ids <- function(x, source) {
  raw_id <- normalize_id(x)
  ans <- match(raw_id, main_ids)
  missing <- which(!is.na(raw_id) & is.na(ans))
  for (i in missing) {
    candidates <- which(!is.na(canonical_main) &
                        canonical_main == canon_id(raw_id[i]))
    if (length(candidates) == 1L) ans[i] <- candidates
    if (length(candidates) > 1L) {
      selected <- review_choice(source, raw_id[i], "ambiguous_id_alias")
      target <- if (is.null(selected)) NA_integer_ else
        match(selected$candidate_life_id, main_ids)
      if (!is.na(target) && target %in% candidates) ans[i] <- target else
        add_issue(source, raw_id[i], "ambiguous_id_alias",
                  paste(main_ids[candidates], collapse = " | "))
    }
  }
  ans
}
add_column <- function(name, indexed_values) {
  if (name %in% names(result_data)) stop("New column already exists: ", name)
  result_data[[name]] <<- indexed_values[match(result_data$life_id, main_ids)]
}
coverage <- list()
add_coverage <- function(source, column, matched, nonmissing, duplicates, note = "") {
  coverage[[length(coverage) + 1L]] <<- data.frame(
    source = source, column = column, participants_found = matched,
    participants_with_value = nonmissing, source_duplicate_ids = duplicates,
    note = note, stringsAsFactors = FALSE)
}

# ---- 6. CER-S: read the correctly separated XLSX, never repeated headers ---
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Please install readxl: install.packages('readxl')")
ecga_path <- file.path(life, "tabular/AMPS_CER-S/ECGAna_allSubj.xlsx")
if (!file.exists(ecga_path)) stop("Missing CER-S Excel table: ", ecga_path)
ecga <- as.data.frame(readxl::read_excel(ecga_path, sheet = 1,
  col_types = "text", .name_repair = "minimal"),
  check.names = FALSE, stringsAsFactors = FALSE)
names(ecga) <- trimws(names(ecga))
ecga[] <- lapply(ecga, normalize)
beats <- c("N Beats", "V Beats", "S Beats")
if (!all(beats %in% names(ecga)))
  stop("CER-S Excel missing beat columns: ", paste(setdiff(beats, names(ecga)), collapse = ", "))

# An occasional embedded header ("N Beats", "V Beats", "S Beats") or
# 'Excluded' status is NOT a beat count. Accept only nonnegative integers.
clean_beat <- function(x) {
  x <- normalize(x)
  numeric_cell <- !is.na(x) & grepl("^[0-9]+(?:\\.0+)?$", x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  out[numeric_cell] <- sub("\\.0+$", "", x[numeric_cell])
  out
}
invalid_beats <- setNames(vapply(beats, function(col)
  sum(!is.na(ecga[[col]]) & is.na(clean_beat(ecga[[col]]))), integer(1L)), beats)
for (col in beats) ecga[[col]] <- clean_beat(ecga[[col]])
# Automatically locate the filename column; preserve its original header.
file_col <- names(ecga)[vapply(ecga, function(x) {
  x <- basename(normalize(x))
  any(!is.na(x) & grepl("(?i)^ECG_.+[.]hea[.]csv$", x, perl = TRUE))
}, logical(1L))]
if (length(file_col) != 1L)
  stop("Expected exactly one CER-S filename column with ECG_LI..._Ruhe.hea.csv; ",
       "candidates: ", paste(file_col, collapse = ", "))
ecga_files <- basename(normalize(ecga[[file_col]]))
ecga_ids <- sub("(?i)[.]hea[.]csv$", "", ecga_files, perl = TRUE)
ecga_ids <- sub("(?i)^ECG_", "", ecga_ids, perl = TRUE)
ecga_ids <- sub("(?i)_(Ruhe|Pilot)(?:_[0-9]+)?$", "", ecga_ids, perl = TRUE)
valid_ecga <- !is.na(ecga_files) &
  grepl("(?i)^ECG_.+[.]hea[.]csv$", ecga_files, perl = TRUE)
ecga_keys <- rep(NA_integer_, nrow(ecga))
ecga_keys[valid_ecga] <- resolve_ids(ecga_ids[valid_ecga], "ECGAna_allSubj.xlsx")
matched_ecga <- unique(ecga_keys[!is.na(ecga_keys)])
ecga_repeated <- sum(table(ecga_keys[!is.na(ecga_keys)]) > 1L)
beat_values <- as.data.frame(setNames(lapply(beats, function(x)
  rep(NA_character_, length(main_ids))), beats), check.names = FALSE)
for (k in matched_ecga) {
  inds <- which(ecga_keys == k)
  # Skip embedded header / excluded-only records when a numeric record exists.
  usable <- inds[rowSums(!is.na(ecga[inds, beats, drop = FALSE])) > 0L]
  if (!length(usable)) {
    add_issue("ECGAna_allSubj.xlsx", main_ids[k], "no_numeric_beat_counts",
              joined(ecga_files[inds]))
    next
  }
  inds <- usable
  # When several files exist, prefer _Ruhe_2, then _Ruhe; do not use arbitrary
  # source row ordering if equally preferred files report different results.
  is2 <- grepl("(?i)_Ruhe_2[.]hea[.]csv$", ecga_files[inds], perl = TRUE)
  isruhe <- grepl("(?i)_Ruhe(?:_[0-9]+)?[.]hea[.]csv$",
                  ecga_files[inds], perl = TRUE)
  priority <- as.integer(is2) * 10L + as.integer(isruhe)
  chosen <- inds[priority == max(priority)]
  if (length(chosen) > 1L) {
    v <- ecga[chosen, beats, drop = FALSE]
    v[] <- lapply(v, normalize)
    if (any(vapply(v, function(x)
      length(unique(x[!is.na(x)])) > 1L, logical(1L)))) {
      manual_index <- review_source_row("ECGAna_allSubj.xlsx", main_ids[k],
        "conflicting_equal_priority_ecg_files")
      if (!is.na(manual_index) && manual_index %in% chosen) {
        chosen <- manual_index
      } else {
        add_issue("ECGAna_allSubj.xlsx", main_ids[k], "conflicting_equal_priority_ecg_files",
                  paste(basename(ecga_files[chosen]), collapse = " | "))
        next
      }
    }
    # Equal-priority duplicates with compatible values: choose the most complete.
    # A manually chosen conflicting row is already unique and must be retained.
    if (length(chosen) > 1L) {
      completeness <- rowSums(!is.na(v))
      chosen <- chosen[which.max(completeness)]
    }
  }
  for (col in beats) beat_values[[col]][k] <- normalize(ecga[[col]][chosen])
}
for (j in seq_along(beats)) {
  new <- c("ecg_n_beats", "ecg_v_beats", "ecg_s_beats")[j]
  vals <- beat_values[[beats[j]]]
  add_column(new, vals)
  add_coverage("ECGAna_allSubj.xlsx", new, length(matched_ecga),
               sum(!is.na(vals)), ecga_repeated,
               "_Ruhe_2 preferred for duplicated ECG filenames")
}

# ---- 7. GAD-7: match questionnaire SIC to Column2/sic, then IDmatch alias --
if (!requireNamespace("readxl", quietly = TRUE))
  stop("Please install readxl: install.packages('readxl')")
gad_path <- file.path(life, "tabular/Behavior/data/PV0491_T00019_NODUP.xlsx")
if (!file.exists(gad_path)) stop("Missing GAD-7 workbook: ", gad_path)
# Read identifiers as TEXT to preserve leading zeroes and exact ID spelling.
gad <- as.data.frame(readxl::read_excel(gad_path, sheet = 1,
  col_types = "text", .name_repair = "minimal"),
  check.names = FALSE, stringsAsFactors = FALSE)
names(gad) <- trimws(names(gad))
gad_items <- paste0("GAD7_GAD7_", seq_len(7L))
gad_columns <- c(gad_items, "GAD7_DATUM")
if (!all(c("GAD7_SIC", gad_columns) %in% names(gad)))
  stop("GAD-7 workbook missing columns: ",
       paste(setdiff(c("GAD7_SIC", gad_columns), names(gad)), collapse = ", "))
gad[] <- lapply(gad, normalize)

# Some Excel date cells imported as text are actually Excel serial numbers,
# e.g. 40604.5060416667 (the fractional part is time of day). Export as
# YYYY-MM-DD using the conventional Excel 1900 date system.
format_gad_date <- function(x, date_field = "GAD7_DATUM") {
  x <- normalize(x)
  output <- x
  serial <- !is.na(x) & grepl("^[0-9]+(?:\\.[0-9]+)?$", x, perl = TRUE)
  if (any(serial)) {
    ii <- which(serial)
    n <- suppressWarnings(as.numeric(x[ii]))
    valid <- is.finite(n) & n >= 20000 & n <= 80000
    output[ii[valid]] <- format(as.Date(floor(n[valid]),
                                   origin = "1899-12-30"), "%Y-%m-%d")
    output[ii[!valid]] <- NA_character_
  }
  other <- which(!is.na(x) & !serial)
  for (fmt in c("%Y-%m-%d", "%d.%m.%Y", "%d/%m/%Y", "%Y/%m/%d")) {
    if (!length(other)) break
    # The first ten characters also cover ISO datetime strings.
    trial <- suppressWarnings(as.Date(substr(x[other], 1L, 10L), format = fmt))
    good <- which(!is.na(trial))
    if (length(good)) output[other[good]] <- format(trial[good], "%Y-%m-%d")
    other <- other[-good]
  }
  if (length(other)) warning("Unrecognized ", date_field, " value(s): ",
                              paste(unique(x[other]), collapse = " | "))
  output[other] <- NA_character_
  output
}
# The GAD-7 and MEDANAMEEG workbooks use the Excel 1904 date system.
# Their serialized dates need the 1462-day offset after 1900-system parsing.
shift_date_1462 <- function(x) {
  as.character(as.Date(x, format = "%Y-%m-%d") + 1462L)
}
gad$GAD7_DATUM <- shift_date_1462(format_gad_date(gad$GAD7_DATUM))
gad_key <- normalize_id(gad$GAD7_SIC)

# Include identifiers from all compatible absorbed rows, but NEVER reuse an
# identifier belonging to an incompatible partition manually discarded above.
# Never use life_id itself as a GAD7_SIC alias.
retained_rows <- unique(unlist(lapply(final, `[[`, "sources"), use.names = FALSE)) - 1L
alias_values <- values[retained_rows, , drop = FALSE]
direct_aliases <- lapply(main_ids, function(id) {
  x <- alias_values[alias_values$life_id == id, c("Column2", sic), drop = FALSE]
  unique(na.omit(normalize_id(unlist(x, use.names = FALSE))))
})

# Direct matches to either Column2 or SIC have priority over any fallback.
gad_hits <- lapply(direct_aliases, function(keys)
  which(!is.na(gad_key) & gad_key %in% keys))
gad_method <- ifelse(lengths(gad_hits) > 0L, "direct_Column2_or_sic", "unmatched")

# The IDmatch sheet may contain an additional alias in its second column.
# Use it only for LIFE IDs without a direct questionnaire match.
fallback_aliases <- replicate(length(main_ids), character(), simplify = FALSE)
need_fallback <- which(lengths(gad_hits) == 0L)
if (length(need_fallback)) {
  idmatch_path <- file.path(life, "tabular/IDmatch.xlsx")
  if (file.exists(idmatch_path) && "IDmatch" %in% readxl::excel_sheets(idmatch_path)) {
    alt <- as.data.frame(readxl::read_excel(idmatch_path,
      sheet = "IDmatch", col_types = "text", .name_repair = "minimal"),
      check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(alt) < 2L) stop("IDmatch.xlsx / IDmatch must have at least two columns")
    alias_life <- normalize_id(alt[[1L]])
    alias_value <- normalize_id(alt[[2L]])
    # Earlier LIFE_matching used exact IDs, then unique LI-prefix/suffix cleanup.
    # Apply the same restricted fallback only when its result is unambiguous.
    pos <- match(alias_life, main_ids)
    for (j in which(!is.na(alias_life) & is.na(pos))) {
      hits <- which(!is.na(canonical_main) &
                    canonical_main == canon_id(alias_life[j]))
      if (length(hits) == 1L) pos[j] <- hits
    }
    for (j in which(!is.na(pos) & !is.na(alias_value))) {
      k <- pos[j]
      fallback_aliases[[k]] <- unique(c(fallback_aliases[[k]], alias_value[j]))
    }
    for (k in need_fallback) {
      hits <- which(!is.na(gad_key) & gad_key %in% fallback_aliases[[k]])
      if (length(hits)) {
        gad_hits[[k]] <- hits
        gad_method[k] <- "IDmatch_sheet_alias"
      }
    }
  } else {
    warning("GAD-7 fallback unavailable: missing IDmatch.xlsx or 'IDmatch' sheet")
  }
}

# A questionnaire row must not be claimed by multiple LIFE participants.
claimed <- unlist(gad_hits, use.names = FALSE)
claim_counts <- tabulate(claimed, nbins = nrow(gad))
gad_selected <- rep(NA_integer_, length(main_ids))
for (k in seq_along(main_ids)) {
  hit <- unique(gad_hits[[k]])
  if (!length(hit)) next
  gad_problem <- if (any(claim_counts[hit] > 1L))
    "shared_questionnaire_identifier" else if (length(hit) > 1L)
    "multiple_questionnaire_rows" else NA_character_
  if (!is.na(gad_problem)) {
    manual_index <- review_source_row("GAD7", main_ids[k], gad_problem)
    if (!is.na(manual_index) && manual_index %in% hit) {
      gad_selected[k] <- manual_index
    } else if (gad_problem == "shared_questionnaire_identifier") {
      add_issue("GAD7", main_ids[k], gad_problem,
                "At least one GAD7_SIC matches multiple LIFE participants")
    } else {
      add_issue("GAD7", main_ids[k], gad_problem,
                paste("Matching Excel rows:", paste(hit + 1L, collapse = " | ")))
    }
    next
  }
  gad_selected[k] <- hit
}
if (anyDuplicated(gad_selected[!is.na(gad_selected)]))
  stop("Manual GAD-7 decisions assign one questionnaire row to multiple participants")

for (col in gad_columns) {
  matched_values <- rep(NA_character_, length(main_ids))
  has_match <- which(!is.na(gad_selected))
  if (length(has_match)) matched_values[has_match] <- gad[[col]][gad_selected[has_match]]
  output_col <- if (col == "GAD7_DATUM") "date_gad7" else col
  add_column(output_col, matched_values)  # Same value on EVERY final row of a LIFE ID.
  add_coverage("PV0491_T00019_NODUP.xlsx", output_col, length(has_match),
               sum(!is.na(matched_values)), sum(table(gad_key[!is.na(gad_key)]) > 1L),
               "GAD7_SIC via Column2/sic; IDmatch alias only if unmatched")
}
gad_direct_n <- sum(!is.na(gad_selected) & gad_method == "direct_Column2_or_sic")
gad_fallback_n <- sum(!is.na(gad_selected) & gad_method == "IDmatch_sheet_alias")
gad_ambiguous_n <- sum(lengths(gad_hits) > 0L & is.na(gad_selected))

# ---- 8. MEDANAMEEG: neurological history, keyed by LIFE participant ID ---
# No fuzzy matching beyond the same exact/unique canonical-ID logic used above.
med_file <- file.path(life, "tabular/Behavior/data/PV0426_T00399_NODUP.xlsx")
if (!file.exists(med_file)) stop("Missing MEDANAMEEG workbook: ", med_file)
med <- as.data.frame(readxl::read_excel(med_file, sheet = 1,
  col_types = "text", .name_repair = "minimal"),
  check.names = FALSE, stringsAsFactors = FALSE)
names(med) <- trimws(names(med))
if (anyDuplicated(names(med))) stop("Duplicate MEDANAMEEG workbook headers")
med_sources <- c("MEDANAMEEG_DATUM", "MEDANAMEEG_MAE_IK",
  "MEDANAMEEG_MAE_SHT", "MEDANAMEEG_MAE_SF", "MEDANAMEEG_MAE_ME",
  "MEDANAMEEG_MAE_OP")
med_labels <- c("date_medanameeg", "Intrakranielle_Blutung",
  "schaedel_hirn_trauma", "schaedelfraktur", "meningitis", "head_brain_op")
if (!all(c("MEDANAMEEG_SIC", med_sources) %in% names(med)))
  stop("MEDANAMEEG workbook missing: ",
       paste(setdiff(c("MEDANAMEEG_SIC", med_sources), names(med)), collapse = ", "))
med[] <- lapply(med, normalize)
med$MEDANAMEEG_DATUM <- shift_date_1462(format_gad_date(
  med$MEDANAMEEG_DATUM, "MEDANAMEEG_DATUM"))
med_keys <- resolve_ids(med$MEDANAMEEG_SIC, basename(med_file))
med_matched <- unique(med_keys[!is.na(med_keys)])
med_repeated <- sum(table(med_keys[!is.na(med_keys)]) > 1L)
for (j in seq_along(med_sources)) {
  source_field <- med_sources[j]
  target_field <- med_labels[j]
  vals <- rep(NA_character_, length(main_ids))
  for (k in med_matched) {
    indices <- which(med_keys == k)
    nonmissing <- unique(med[[source_field]][indices])
    nonmissing <- nonmissing[!is.na(nonmissing)]
    if (length(nonmissing) == 1L) vals[k] <- nonmissing
    if (length(nonmissing) > 1L) {
      issue_type <- paste0("conflicting_", source_field, "_duplicates")
      manual_index <- review_source_row(basename(med_file), main_ids[k],
                                        issue_type, source_field)
      if (!is.na(manual_index) && manual_index %in% indices &&
          !is.na(med[[source_field]][manual_index]))
        vals[k] <- med[[source_field]][manual_index]
      else add_issue(basename(med_file), main_ids[k], issue_type,
                     paste(length(indices), "input rows; value left missing"))
    }
  }
  add_column(target_field, vals)
  add_coverage(basename(med_file), target_field, length(med_matched),
               sum(!is.na(vals)), med_repeated)
}

# ---- 8a. CNS medication exclusion: participant-level LIFE ID -------------
medications_path <- file.path(life,
  "derivatives/matching_and_summaries/medications/LIFE_medications_by_ID.csv")
if (!file.exists(medications_path))
  stop("Missing medication summary: ", medications_path)
medications <- read_csv(medications_path)
names(medications) <- trimws(names(medications))
if (!ncol(medications) || names(medications)[1L] != "life_id" ||
    !"cns_exclude_definite" %in% names(medications))
  stop("Medication table must start with life_id and contain cns_exclude_definite")
medication_keys <- resolve_ids(medications[[1L]], basename(medications_path))
medication_matched <- unique(medication_keys[!is.na(medication_keys)])
medication_repeated <- sum(table(medication_keys[!is.na(medication_keys)]) > 1L)
medication_values <- rep(NA_character_, length(main_ids))
medication_source <- normalize(medications$cns_exclude_definite)
for (k in medication_matched) {
  matches <- which(medication_keys == k)
  distinct <- unique(na.omit(medication_source[matches]))
  if (length(distinct) == 1L) medication_values[k] <- distinct
  if (length(distinct) > 1L) {
    manual_index <- review_source_row(basename(medications_path), main_ids[k],
                                       "conflicting_cns_exclude_definite", "cns_exclude_definite")
    if (!is.na(manual_index) && manual_index %in% matches &&
        !is.na(medication_source[manual_index]))
      medication_values[k] <- medication_source[manual_index]
    else add_issue(basename(medications_path), main_ids[k],
      "conflicting_cns_exclude_definite", "Conflicting duplicate rows; left missing")
  }
}
add_column("medications_cns_exclude_definite", medication_values)
add_coverage(basename(medications_path), "medications_cns_exclude_definite",
  length(medication_matched), sum(!is.na(medication_values)), medication_repeated)

# ---- 8b. Echo LVOT stroke volume: match pseudonym to Column2 or SIC ------
echo_path <- file.path(life, "tabular/Behavior/data/PV491_T01092_ECHO_LVOT.csv")
if (!file.exists(echo_path)) stop("Missing ECHO LVOT file: ", echo_path)
echo_fields <- c("ECHO_READ_STUDYDATETIME", "ECHO_READ_AGE", "ECHO_READ_LVOT_SV")
# The first column has a damaged header (e.g., '���pseudonym'). Identify it
# by POSITION, while validating all three named measurement columns. Try
# common delimiters/encodings rather than silently reading one giant column.
read_echo <- function(path) {
  for (enc in c("UTF-8-BOM", "latin1", "UTF-16LE")) {
    for (sep in c(",", ";", "\t")) {
      candidate <- suppressWarnings(tryCatch(read.table(path,
        header = TRUE, sep = sep, quote = '"', comment.char = "",
        colClasses = "character", check.names = FALSE, stringsAsFactors = FALSE,
        na.strings = "", fileEncoding = enc), error = function(e) NULL))
      if (is.null(candidate) || ncol(candidate) < 4L) next
      # Leave the garbled first header untouched; only named fields matter.
      names(candidate)[-1L] <- trimws(names(candidate)[-1L])
      if (all(echo_fields %in% names(candidate))) return(candidate)
    }
  }
  stop("Could not parse ECHO CSV into pseudonym and required columns: ",
       paste(echo_fields, collapse = ", "))
}
echo <- read_echo(echo_path)
if (anyDuplicated(names(echo)[-1L])) stop("Duplicate ECHO measurement headers")
echo[] <- lapply(echo, normalize)
echo_pseudonym <- normalize_id(echo[[1L]])
echo$ECHO_READ_STUDYDATETIME <- format_gad_date(
  echo$ECHO_READ_STUDYDATETIME, "ECHO_READ_STUDYDATETIME")

# Same pseudonym may appear in several original IDmatch rows of one LIFE
# participant. It must not be shared by DIFFERENT LIFE participants.
# direct_aliases were built above from all original Column2/SIC values.
alias_values <- unlist(direct_aliases, use.names = FALSE)
alias_owners <- rep(seq_along(main_ids), lengths(direct_aliases))
alias_map <- split(alias_owners, alias_values)
echo_keys <- rep(NA_integer_, nrow(echo))
for (i in which(!is.na(echo_pseudonym))) {
  candidates <- unique(alias_map[[echo_pseudonym[i]]])
  if (length(candidates) == 1L) echo_keys[i] <- candidates
  if (length(candidates) > 1L) {
    selected <- review_choice(basename(echo_path), echo_pseudonym[i],
                              "shared_echo_pseudonym")
    target <- if (is.null(selected)) NA_integer_ else
      match(selected$candidate_life_id, main_ids)
    if (!is.na(target) && target %in% candidates) echo_keys[i] <- target else
      add_issue(basename(echo_path), echo_pseudonym[i],
        "shared_echo_pseudonym", paste(main_ids[candidates], collapse = " | "))
  }
}
echo_matched <- unique(echo_keys[!is.na(echo_keys)])
echo_repeated <- sum(table(echo_keys[!is.na(echo_keys)]) > 1L)
next_echo_name <- function(base, existing) {
  if (!(base %in% existing)) return(base)
  n <- 2L
  while (paste0(base, "_", n) %in% existing) n <- n + 1L
  paste0(base, "_", n)
}
echo_labels <- vapply(
  c("date_echo", "age_echo", "echo_lvot_sv"),
  next_echo_name, character(1L),
  existing = names(result_data)
)
for (j in seq_along(echo_fields)) {
  source_col <- echo_fields[j]
  output_col <- echo_labels[j]
  vals <- rep(NA_character_, length(main_ids))
  for (k in echo_matched) {
    matches <- which(echo_keys == k)
    distinct <- unique(na.omit(echo[[source_col]][matches]))
    if (length(distinct) == 1L) vals[k] <- distinct
    if (length(distinct) > 1L) {
      issue_type <- paste0("conflicting_", source_col, "_duplicates")
      manual_index <- review_source_row(basename(echo_path), main_ids[k],
                                        issue_type, source_col)
      if (!is.na(manual_index) && manual_index %in% matches &&
          !is.na(echo[[source_col]][manual_index]))
        vals[k] <- echo[[source_col]][manual_index]
      else add_issue(basename(echo_path), main_ids[k], issue_type,
        paste(length(matches), "ECHO rows; value left missing"))
    }
  }
  add_column(output_col, vals)
  add_coverage(basename(echo_path), output_col, length(echo_matched),
               sum(!is.na(vals)), echo_repeated,
               "ECHO pseudonym matched to original Column2/SIC aliases")
}

# ---- 9. HR/HRV: import all measurements by participant ID ---------------
# Keep the names from the source file; only prefix when a name already exists.
hrv_path <- file.path(life, "derivatives/HR_HRV_all_participants.csv")
if (!file.exists(hrv_path)) stop("Missing HR/HRV file: ", hrv_path)
hrv <- read_csv(hrv_path)
names(hrv) <- trimws(names(hrv))
if (anyDuplicated(names(hrv))) stop("Duplicate column names in HR/HRV file")
hrv_id_candidates <- c("subject_id", "life_id", "participant_id", "SubjectID", "subject", "Subject", "ID", "id")
hrv_id_col <- intersect(hrv_id_candidates, names(hrv))
if (!length(hrv_id_col)) {
  warning("No recognized HR/HRV ID header; assuming the FIRST column is the ID: ",
          names(hrv)[1L])
  hrv_id_col <- names(hrv)[1L]
}
hrv_id_col <- hrv_id_col[1L]
hrv_keys <- resolve_ids(hrv[[hrv_id_col]], "HR_HRV_all_participants.csv")
hrv_matched <- unique(hrv_keys[!is.na(hrv_keys)])
if (!length(hrv_matched))
  warning("No HR/HRV identifiers matched LIFE participant IDs. Check ID format and source.")
hrv_repeated <- sum(table(hrv_keys[!is.na(hrv_keys)]) > 1L)
hrv_fields <- setdiff(names(hrv), hrv_id_col)
if (!length(hrv_fields)) warning("HR/HRV table has no measurement columns")
hrv_names <- setNames(character(length(hrv_fields)), hrv_fields)
for (src in hrv_fields) {
  label <- if (src %in% names(result_data)) paste0("HR_HRV_", src) else src
  if (label %in% names(result_data) || label %in% hrv_names)
    stop("Ambiguous HR/HRV output column: ", label)
  hrv_names[src] <- label
  vals <- rep(NA_character_, length(main_ids))
  original <- normalize(hrv[[src]])
  for (k in hrv_matched) {
    matches <- which(hrv_keys == k)
    nonmissing <- unique(original[matches][!is.na(original[matches])])
    if (length(nonmissing) == 1L) vals[k] <- nonmissing
    if (length(nonmissing) > 1L) {
      manual_index <- review_source_row("HR_HRV_all_participants.csv", main_ids[k],
        "conflicting_duplicate_measurements", src)
      if (!is.na(manual_index) && manual_index %in% matches &&
          !is.na(original[manual_index]))
        vals[k] <- original[manual_index]
      else add_issue("HR_HRV_all_participants.csv", main_ids[k],
        "conflicting_duplicate_measurements",
        paste(src, ": ", length(nonmissing), " different values; left missing"))
    }
  }
  add_column(label, vals)
  add_coverage("HR_HRV_all_participants.csv", label, length(hrv_matched),
               sum(!is.na(vals)), hrv_repeated)
}

# ---- 10. Derived participant-level values ---------------------------------
# Scores are calculated for each output row; external measurements already
# repeat identically across unresolved duplicate rows of the same LIFE ID.
# All questionnaire items must be coded 0,1,2,3. Unexpected values become NA.
item_numeric <- function(x) {
  x <- normalize(x)
  n <- suppressWarnings(as.numeric(gsub(",", ".", x, fixed = TRUE)))
  n[!is.finite(n) | is.na(n) | !(n %in% 0:3)] <- NA_real_
  n
}
positive_number <- function(x) {
  n <- suppressWarnings(as.numeric(gsub(",", ".", normalize(x), fixed = TRUE)))
  n[!is.finite(n)] <- NA_real_
  n
}
nonmissing_row_sum <- function(x) {
  out <- rowSums(x, na.rm = TRUE)
  out[rowSums(!is.na(x)) == 0L] <- NA_real_
  out
}

# VIGALL drowsiness: B+C > 900 one-second segments (strictly >15 min).
# VIGALL comes from the same original row selected for EEG; never backfilled.
if (!all(c("vigall_b", "vigall_c") %in% names(result_data)))
  stop("Expected vigall_b and vigall_c in selected IDmatch variables")
vig_b <- positive_number(result_data$vigall_b)
vig_c <- positive_number(result_data$vigall_c)
result_data$drowsy_15min <- ifelse(!is.na(vig_b) & !is.na(vig_c),
  as.integer(vig_b + vig_c > 900), NA_integer_)

# GAD-7: requested available-item sum, with a separate preregistration validity
# flag. A partially observed score is NOT valid for the registered models.
if (!all(gad_items %in% names(result_data)))
  stop("GAD-7 items missing after questionnaire matching")
gad_mat <- as.data.frame(lapply(result_data[gad_items], item_numeric))
gad_valid_count <- rowSums(!is.na(gad_mat))
result_data$all_GAD7 <- as.integer(gad_valid_count == 7L)
result_data$GAD7 <- nonmissing_row_sum(gad_mat)
result_data$gad7_group <- ifelse(result_data$all_GAD7 == 1L &
                                  result_data$GAD7 <= 4, "low",
                           ifelse(result_data$all_GAD7 == 1L &
                                    result_data$GAD7 >= 10, "high", NA_character_))

# IDS-SR30: score ONE value from each alternative appetite and weight pair,
# using max(item 11, item 12), max(item 13, item 14).
# LIFE-specific mapping verified by the data owner: IDS_I09A is scored item 9
# (severity of mood variation, 0-3). IDS_I09B is NOT included in IDS total.
# Valid IDS-SR total requires all 28 scored items (range 0-84).
ids_regular <- sprintf("IDS_I%02d", c(1:10, 15:30))  # 26 ordinary items
ids_slots <- c(ids_regular, "IDS_appetite", "IDS_weight") # 28 scored slots
ids_mat <- matrix(NA_real_, nrow = nrow(result_data), ncol = length(ids_slots),
                  dimnames = list(NULL, ids_slots))
for (field in intersect(setdiff(ids_regular, "IDS_I09"), names(result_data)))
  ids_mat[, field] <- item_numeric(result_data[[field]])
# The input field IDS_I09A represents the standard scored item 9 in LIFE.
if ("IDS_I09A" %in% names(result_data)) {
  ids_mat[, "IDS_I09"] <- item_numeric(result_data[["IDS_I09A"]])
} else {
  warning("IDS_I09A (LIFE item 9) is missing: IDS scores may be incomplete.")
}
ids_pair <- function(a, b) {
  a <- if (a %in% names(result_data)) item_numeric(result_data[[a]]) else
    rep(NA_real_, nrow(result_data))
  b <- if (b %in% names(result_data)) item_numeric(result_data[[b]]) else
    rep(NA_real_, nrow(result_data))
  out <- pmax(a, b, na.rm = TRUE)
  out[is.infinite(out)] <- NA_real_
  out
}
ids_mat[, "IDS_appetite"] <- ids_pair("IDS_I11", "IDS_I12")
ids_mat[, "IDS_weight"] <- ids_pair("IDS_I13", "IDS_I14")
result_data$IDS_n_items <- rowSums(!is.na(ids_mat))
result_data$all_IDS <- as.integer(result_data$IDS_n_items == 28L)
result_data$IDS <- nonmissing_row_sum(ids_mat)

# Extrasystoles: the CER-S table leaves S/V blank when none were detected.
# Keep the imported ecg_n_beats/ecg_s_beats/ecg_v_beats unchanged: missing
# S/V means zero ONLY for deriving the count, not for the source columns.
required_beats <- c("ecg_s_beats", "ecg_v_beats", "ecg_n_beats")
if (!all(required_beats %in% names(result_data)))
  stop("Cannot calculate extrasystoles; beat columns are missing")
s <- positive_number(result_data$ecg_s_beats)
v <- positive_number(result_data$ecg_v_beats)
n <- positive_number(result_data$ecg_n_beats)
# Any recorded numeric N/S/V value establishes a usable CER-S measurement;
# if all three are absent (e.g., Excluded), do NOT turn absence into zero.
beat_information <- !is.na(n) | !is.na(s) | !is.na(v)
result_data$extrasystolic_N <- ifelse(
  beat_information, rowSums(cbind(s, v), na.rm = TRUE), NA_real_)
result_data$extrasystolic_proportion <- ifelse(
  !is.na(n) & n > 0 & !is.na(result_data$extrasystolic_N),
  result_data$extrasystolic_N / n, NA_real_)

# Compute the across-participant reference on unique measured LIFE IDs only.
# This INCLUDES genuine zero counts when N is available and S/V are blank,
# and EXCLUDES participants whose entire CER-S record is unmeasured.
extra_unique <- !duplicated(result_data$life_id)
extra_reference <- result_data$extrasystolic_N[extra_unique]
extra_reference <- extra_reference[!is.na(extra_reference)]
extra_mean <- if (length(extra_reference)) mean(extra_reference) else NA_real_
extra_sd <- if (length(extra_reference) >= 2L) sd(extra_reference) else NA_real_
result_data$extrasystolic_N_SD <- if (is.finite(extra_sd) && extra_sd > 0)
  (result_data$extrasystolic_N - extra_mean) / extra_sd else
  rep(NA_real_, nrow(result_data))
# Preregistered one-sided exclusion: only more than 4 SD ABOVE the mean.
result_data$extrasystolic_exclusion <- ifelse(
  is.na(result_data$extrasystolic_N_SD), NA_integer_,
  as.integer(result_data$extrasystolic_N_SD > 4)
)

# ---- 11. Save ONLY useful compact outputs ------------------------------
# Do not silently ignore a manual choice if the candidates changed as a
# consequence of selecting a different participant/identifier partition.
if (resolution_mode && nrow(issue_review)) {
  selected <- issue_review[yes(issue_review$select), , drop = FALSE]
  chosen_keys <- paste(selected$source, selected$life_id,
                       selected$issue, ifelse(is.na(selected$field), "", selected$field),
                       sep = "\035")
  unused <- setdiff(chosen_keys, consumed_decisions)
  if (length(unused))
    stop("One or more selected external-source decisions were not used; ",
         "review source candidates after the primary-row choice")
}
# Only explicitly acknowledged missing/unresolvable source problems may remain.
# Manual selections that unexpectedly produce NEW conflicts fail before data.csv.
if (resolution_mode && length(issues)) {
  blank <- function(x) ifelse(is.na(x), "", x)
  issue_field <- function(source, problem, details) {
    if (identical(source, "HR_HRV_all_participants.csv") &&
        identical(problem, "conflicting_duplicate_measurements"))
      return(trimws(sub(":.*$", "", details)))
    if (identical(problem, "conflicting_cns_exclude_definite"))
      return("cns_exclude_definite")
    if (grepl("^conflicting_.*_duplicates$", problem))
      return(sub("_duplicates$", "", sub("^conflicting_", "", problem)))
    ""
  }
  for (entry in issues) {
    source <- entry$source[[1L]]
    id <- entry$life_id[[1L]]
    problem <- entry$issue[[1L]]
    field <- issue_field(source, problem, entry$details[[1L]])
    reviewed <- which(blank(issue_review$source) == source &
      blank(issue_review$life_id) == id & blank(issue_review$issue) == problem &
      blank(issue_review$field) == field &
      blank(issue_review$resolution) == "keep_missing")
    if (!length(reviewed)) stop("New or unreviewed matching issue after manual selection: ",
                                source, " / ", id, " / ", problem, " / ", field)
  }
}
if (anyDuplicated(names(result_data))) stop("Duplicate output column names")
if (nrow(result_data) != length(final)) stop("Internal error: row count changed")
if (!identical(result_data[["non-unique"]],
               as.integer(table(result_data$life_id)[result_data$life_id] > 1L)))
  stop("Internal error: non-unique flag does not match final IDs")
# Preparation writes an intermediate; only validated resolution writes data.csv.
if (resolution_mode && anyDuplicated(result_data$life_id))
  stop("Unresolved duplicate LIFE participants remain after manual primary selection")
data_path <- if (resolution_mode) file.path(life, "derivatives", "data.csv") else
  file.path(out, "data_pre_review.csv")
write.csv(result_data, data_path, na = "", row.names = FALSE)

# Row provenance: original source row numbers + choices, without all cell values.
provenance <- do.call(rbind, lapply(seq_along(final), function(i) {
  r <- final[[i]]
  data.frame(data_csv_row = i + 1L, life_id = r$id,
    non_unique = result_data[["non-unique"]][i],
    primary_input_csv_row = r$primary, all_input_csv_rows = source_ids(r),
    n_input_rows_merged = length(r$sources), selection_rule = r$reason,
    n_alternative_eeg_dates = length(r$alts), stringsAsFactors = FALSE)
}))
write_csv(provenance, if (resolution_mode) "data_final_row_provenance.csv" else
  "data_row_provenance.csv")

file_info <- function(id) {
  f <- recordings[!is.na(recordings$life_id) & recordings$life_id == id, , drop = FALSE]
  list(paths = joined(f$vhdr_path),
       selected = joined(f$vhdr_path[f$is_selected]))
}
unresolved <- list()
for (id in unique(result_data$life_id[result_data[["non-unique"]] == 1L])) {
  positions <- which(result_data$life_id == id)
  f <- file_info(id)
  unresolved[[length(unresolved) + 1L]] <- data.frame(
    life_id = id, n_rows = length(positions),
    data_csv_rows = joined(positions + 1L),
    primary_input_csv_rows = joined(provenance$primary_input_csv_row[positions]),
    eeg_dates = joined(paste0(positions + 1L, ": ",
                       ifelse(is.na(result_data[[dt]][positions]), "<missing>",
                              result_data[[dt]][positions]))),
    eeg_groups = joined(paste0(positions + 1L, ": ",
                        ifelse(is.na(result_data[[grp]][positions]), "<missing>",
                               result_data[[grp]][positions]))),
    recording_paths = f$paths, selected_recording_path = f$selected,
    selection_rules = joined(provenance$selection_rule[positions]),
    stringsAsFactors = FALSE)
}
unresolved_cols <- c("life_id", "n_rows", "data_csv_rows", "primary_input_csv_rows",
  "eeg_dates", "eeg_groups", "recording_paths", "selected_recording_path",
  "selection_rules")
write_csv(bind_events(unresolved, unresolved_cols),
  if (resolution_mode) "data_final_unresolved_rows.csv" else "data_unresolved_rows.csv")
# Dates and candidate SICs/Excel row numbers help manual matching, without
# dumping full source rows into the audit. All values reflect original matches.
issue_cols <- c("source", "life_id", "issue", "details")
issue_table <- bind_events(issues, issue_cols)
issue_table$eeg_date <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_date <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_sic <- rep(NA_character_, nrow(issue_table))
issue_table$gad7_excel_rows <- rep(NA_character_, nrow(issue_table))
for (i in seq_len(nrow(issue_table))) {
  id <- issue_table$life_id[i]
  k <- match(id, main_ids)
  if (!is.na(k)) {
    # Include dates across unresolved LIFE rows, plus saved alternative visits.
    positions <- which(result_data$life_id == id)
    date_fields <- intersect(unique(c(dt, grep(paste0("^", dt, "[2-9][0-9]*$"),
                                       names(result_data), value = TRUE))),
                             names(result_data))
    issue_table$eeg_date[i] <- joined(unlist(result_data[positions, date_fields,
                                                  drop = FALSE], use.names = FALSE))
    candidates <- unique(gad_hits[[k]])
  } else if (issue_table$source[i] == "GAD7") {
    # For a questionnaire-side identifier that could not resolve to LIFE.
    candidates <- which(!is.na(gad_key) & gad_key == normalize_id(id))
  } else {
    candidates <- integer()
  }
  if (length(candidates)) {
    issue_table$gad7_date[i] <- joined(gad$GAD7_DATUM[candidates])
    issue_table$gad7_sic[i] <- joined(gad_key[candidates])
    issue_table$gad7_excel_rows[i] <- joined(candidates + 1L)
  }
}
write_csv(issue_table, if (resolution_mode) "data_final_match_issues.csv" else
  "data_match_issues.csv")

# ---- 11a. Candidate-level manual review templates (never overwrite edits) ---
if (!resolution_mode) {
  dir.create(manual_dir, recursive = TRUE, showWarnings = FALSE)
  input_path <- file.path(out, "IDmatch_EEG.csv")
  primary_ids <- unique(result_data$life_id[result_data[["non-unique"]] == 1L])
  original_idx <- which(values$life_id %in% primary_ids)
  primary_paths <- vapply(values$life_id[original_idx], function(id)
    file_info(id)$paths, character(1L))
  primary_template <- data.frame(
    life_id = values$life_id[original_idx],
    input_csv_row = as.character(original_idx + 1L),
    Column2 = values$Column2[original_idx],
    sic = values[[sic]][original_idx],
    eeg_group = values[[grp]][original_idx],
    eeg_date = values[[dt]][original_idx],
    core_eeg_fields = as.character(rowSums(!is.na(values[original_idx, core_eeg,
                                                     drop = FALSE]))),
    vigall_fields = as.character(rowSums(!is.na(values[original_idx, vigall,
                                                   drop = FALSE]))),
    recording_paths = unname(primary_paths),
    source_path = rep(input_path, length(original_idx)),
    source_md5 = rep(unname(tools::md5sum(input_path)), length(original_idx)),
    mapping_md5 = rep(unname(tools::md5sum(file.path(out,
      "IDmatch_variables_manual.csv"))), length(original_idx)),
    select_primary = "", review_note = "", stringsAsFactors = FALSE)
  # source_path/hash are repeated on every candidate to detect stale decisions.
  write_csv(primary_template, "data_unresolved_rows_template.csv")
  if (!file.exists(manual_primary_path))
    write.csv(primary_template, manual_primary_path, na = "", row.names = FALSE)

  issue_paths <- c(
    "ECGAna_allSubj.xlsx" = ecga_path, "GAD7" = gad_path,
    setNames(med_file, basename(med_file)),
    setNames(medications_path, basename(medications_path)),
    setNames(echo_path, basename(echo_path)),
    "HR_HRV_all_participants.csv" = hrv_path)
  issue_hashes <- tools::md5sum(unname(issue_paths))
  matching_hash <- unname(tools::md5sum(input_path))
  mapping_hash <- unname(tools::md5sum(file.path(out, "IDmatch_variables_manual.csv")))
  review_candidates <- list()
  add_candidates <- function(source, id, issue_type, field = "", rows = integer(),
                             owner_ids = character(), source_ids = character(),
                             preview = character()) {
    size <- max(length(rows), length(owner_ids), 1L)
    if (length(rows) && length(rows) != size) stop("Bad candidate row lengths")
    if (length(owner_ids) && length(owner_ids) != size) stop("Bad owner ID lengths")
    path <- unname(issue_paths[source])
    if (length(path) != 1L || is.na(path) || !file.exists(path))
      stop("Missing source file for review issue: ", source)
    expand <- function(x) if (length(x)) as.character(x) else rep("", size)
    review_candidates[[length(review_candidates) + 1L]] <<- data.frame(
      source = source, life_id = id, issue = issue_type, field = field,
      source_row = if (length(rows)) as.character(rows + 1L) else rep("", size),
      candidate_life_id = expand(owner_ids),
      candidate_identifier = expand(source_ids),
      candidate_value = expand(preview), source_path = path,
      source_md5 = unname(issue_hashes[path]),
      matching_input_md5 = matching_hash, mapping_md5 = mapping_hash,
      select = "", resolution = if (issue_type == "no_numeric_beat_counts")
        "keep_missing" else "", review_note = "", stringsAsFactors = FALSE)
  }
  for (i in seq_len(nrow(issue_table))) {
    source <- issue_table$source[i]
    id <- issue_table$life_id[i]
    problem <- issue_table$issue[i]
    k <- match(id, main_ids)
    if (problem == "ambiguous_id_alias") {
      possible <- main_ids[!is.na(canonical_main) & canonical_main == canon_id(id)]
      add_candidates(source, id, problem, owner_ids = possible,
                     source_ids = rep(id, length(possible)), preview = possible)
    } else if (problem == "shared_echo_pseudonym") {
      possible <- main_ids[unique(alias_map[[id]])]
      add_candidates(source, id, problem, owner_ids = possible,
                     source_ids = rep(id, length(possible)), preview = possible)
    } else if (source == "GAD7") {
      inds <- if (!is.na(k)) unique(gad_hits[[k]]) else integer()
      add_candidates(source, id, problem, rows = inds,
                     source_ids = gad_key[inds],
                     preview = if (length(inds)) paste0("date=", gad$GAD7_DATUM[inds])
                       else character())
    } else if (source == "ECGAna_allSubj.xlsx") {
      inds <- if (!is.na(k)) which(ecga_keys == k) else integer()
      if (problem == "conflicting_equal_priority_ecg_files") {
        inds <- inds[rowSums(!is.na(ecga[inds, beats, drop = FALSE])) > 0L]
        if (length(inds)) {
          is2 <- grepl("(?i)_Ruhe_2[.]hea[.]csv$", ecga_files[inds], perl = TRUE)
          isruhe <- grepl("(?i)_Ruhe(?:_[0-9]+)?[.]hea[.]csv$",
                          ecga_files[inds], perl = TRUE)
          p <- as.integer(is2) * 10L + as.integer(isruhe)
          inds <- inds[p == max(p)]
        }
      } else inds <- integer()  # no_numeric_beat_counts: no selectable value
      add_candidates(source, id, problem, rows = inds,
        source_ids = ecga_files[inds],
        preview = if (length(inds)) paste0("N=", ecga[[beats[1L]]][inds], "; V=",
                         ecga[[beats[2L]]][inds], "; S=", ecga[[beats[3L]]][inds])
          else character())
    } else if (source == basename(med_file)) {
      field <- sub("_duplicates$", "", sub("^conflicting_", "", problem))
      inds <- if (!is.na(k)) which(med_keys == k & !is.na(med[[field]])) else integer()
      add_candidates(source, id, problem, field, inds,
                     source_ids = med$MEDANAMEEG_SIC[inds], preview = med[[field]][inds])
    } else if (source == basename(medications_path)) {
      field <- "cns_exclude_definite"
      inds <- if (!is.na(k)) which(medication_keys == k &
                                   !is.na(medication_source)) else integer()
      add_candidates(source, id, problem, field, inds,
                     source_ids = medications[[1L]][inds],
                     preview = medication_source[inds])
    } else if (source == basename(echo_path)) {
      field <- sub("_duplicates$", "", sub("^conflicting_", "", problem))
      inds <- if (!is.na(k)) which(echo_keys == k & !is.na(echo[[field]])) else integer()
      add_candidates(source, id, problem, field, inds,
                     source_ids = echo_pseudonym[inds], preview = echo[[field]][inds])
    } else if (source == "HR_HRV_all_participants.csv") {
      field <- trimws(sub(":.*$", "", issue_table$details[i]))
      inds <- if (!is.na(k)) which(hrv_keys == k & !is.na(hrv[[field]])) else integer()
      add_candidates(source, id, problem, field, inds,
                     source_ids = hrv[[hrv_id_col]][inds], preview = hrv[[field]][inds])
    } else stop("Unrecognized matching issue for manual template: ", source, "/", problem)
  }
  issue_template_cols <- c("source", "life_id", "issue", "field", "source_row",
    "candidate_life_id", "candidate_identifier", "candidate_value",
    "source_path", "source_md5", "matching_input_md5", "mapping_md5",
    "select", "resolution", "review_note")
  issue_template <- unique(bind_events(review_candidates, issue_template_cols))
  write_csv(issue_template, "data_match_issues_template.csv")
  if (!file.exists(manual_issues_path))
    write.csv(issue_template, manual_issues_path, na = "", row.names = FALSE)
}

summary_coverage <- bind_events(coverage, c("source", "column", "participants_found",
                                            "participants_with_value", "source_duplicate_ids", "note"))
count_cells <- function(action) if (action %in% names(cell_actions))
  unname(cell_actions[[action]]) else 0L
n_original_repeated <- sum(table(values$life_id) > 1L)
n_final_repeated <- sum(table(result_data$life_id) > 1L)
report <- capture.output({
  cat("========== LIFE analysis data preparation ==========\n")
  cat("Matched IDmatch input rows:         ", nrow(values), "\n")
  cat("Unique participants:                ", length(main_ids), "\n")
  cat("Participants with repeated inputs: ", n_original_repeated, "\n")
  cat("Rows in ", basename(data_path), ": ", nrow(result_data), "\n", sep = "")
  cat("Participants with one final row:    ", sum(table(result_data$life_id) == 1L), "\n")
  cat("Participants still non-unique:      ", n_final_repeated, "\n")
  cat("Original rows absorbed:            ",
      nrow(values) - nrow(result_data) - manual_dropped_n, "\n")
  if (resolution_mode)
    cat("Original rows deliberately discarded:", manual_dropped_n, "\n")
  cat("Complementary cells filled:        ", count_cells("filled_from_secondary"), "\n")
  cat("Primary conflicting cells kept:    ", count_cells("kept_primary_conflict"), "\n")
  cat("Secondary VIGALL cells NOT used:   ", count_cells("skipped_secondary_vigall"), "\n")
  cat("Other session EEG cells NOT used:  ", count_cells("skipped_secondary_eeg_session"), "\n")
  cat("Alternative group/date pairs saved:", count_cells("saved_alternative_group_date"), "\n")
  cat("Fully empty selected columns dropped:", length(all_empty), "\n")
  cat("EP-EW columns removed:               ", paste(block, collapse = ", "), "\n")
  cat("Final output columns:               ", ncol(result_data), "\n")
  cat("MEDANAMEEG matched participants:      ", length(med_matched), "\n")
  cat("MEDANAMEEG source duplicate IDs:     ", med_repeated, "\n")
  cat("CNS medication matched participants:  ", length(medication_matched), "\n")
  cat("CNS medication flags with values:     ", sum(!is.na(medication_values)), "\n")
  cat("ECHO LVOT matched participants:       ", length(echo_matched), "\n")
  cat("ECHO source duplicate IDs:           ", echo_repeated, "\n")
  cat("HR/HRV columns imported:             ", length(hrv_fields), "\n")
  cat("GAD-7 complete (unique participants):", sum(result_data$all_GAD7[!duplicated(result_data$life_id)] == 1L), "\n")
  cat("IDS complete (unique participants):  ", sum(result_data$all_IDS[!duplicated(result_data$life_id)] == 1L), "\n")
  cat("VIGALL B+C >900 (unique participants):", sum(result_data$drowsy_15min[extra_unique] == 1L, na.rm = TRUE), "\n")
  cat("VIGALL B/C missing (unique subjects): ", sum(is.na(result_data$drowsy_15min[extra_unique])), "\n")
  cat("CER-S without usable beat counts:    ", sum(!beat_information[extra_unique]), "\n")
  cat("CER-S with N but missing S and V (assumed zero): ",
      sum(!is.na(n[extra_unique]) & is.na(s[extra_unique]) & is.na(v[extra_unique])), "\n")
  cat("Extrasystole SD reference participants: ", length(extra_reference), "\n")
  cat("Observed extrasystolic_N = 0 (unique participants): ",
      sum(result_data$extrasystolic_N[extra_unique] == 0, na.rm = TRUE), "\n")
  cat("GAD-7 groups high/low (unique participants):\n")
  print(table(result_data$gad7_group[!duplicated(result_data$life_id)], useNA = "ifany"))
  cat("Extrasystole count mean (unique participants):", extra_mean, "\n")
  cat("Extrasystole count SD (unique participants):  ", extra_sd, "\n")
  cat("Extrasystole >4-SD exclusions (unique participants):",
      sum(result_data$extrasystolic_exclusion[extra_unique] == 1L, na.rm = TRUE), "\n")
  cat("Pipeline HEP/ECG/preprocessing imported: NO (deferred until post-exclusion)\n")
  cat("Prior EEG recording inventory loaded:", !is.null(inventory), "\n")
  if (!is.null(participants)) cat("Earlier LIFE matching participants: ", nrow(participants), "\n")
  if (!is.null(prior_duplicates)) cat("Earlier LIFE repeated-row cases:    ", nrow(prior_duplicates), "\n")
  cat("\n---------- External data coverage ----------\n")
  print(summary_coverage[, c("column", "participants_found", "participants_with_value",
                             "source_duplicate_ids"), drop = FALSE], row.names = FALSE)
  cat("\nCER-S filenames interpreted from column:", file_col, "\n")
  cat("CER-S invalid/non-numeric count cells ignored:",
      paste(paste(names(invalid_beats), invalid_beats, sep = "="), collapse = "; "), "\n")
  cat("Final rows flagged non-unique:", sum(result_data[["non-unique"]] == 1L), "\n")
  cat("GAD-7 matched directly:              ", gad_direct_n, "\n")
  cat("GAD-7 matched via IDmatch fallback:  ", gad_fallback_n, "\n")
  cat("GAD-7 ambiguous cases left missing:  ", gad_ambiguous_n, "\n")
  cat("GAD-7 participants without a match:  ", sum(lengths(gad_hits) == 0L), "\n")
  cat("External matching issues needing review:", length(issues), "\n")
  if (length(issues)) {
    print(table(vapply(issues, function(x) x$source[1L], character(1L))))
  }
  cat("\nNon-unique participants and EEG dates: see data_unresolved_rows.csv\n")
  cat("EEG/VIGALL primary row is preserved; VIGALL is NEVER backfilled.\n")
  cat("External source values are replicated across all surviving rows of an ID.\n")
})
writeLines(report, file.path(out,
  if (resolution_mode) "data_final_summary.txt" else "data_summary.txt"))
if (resolution_mode) {
  dropped_cols <- c("life_id", "selected_input_csv_row",
                    "discarded_input_csv_rows", "reason")
  write_csv(bind_events(manual_dropped, dropped_cols), "data_manual_dropped_rows.csv")
}

cat(paste(report, collapse = "\n"), "\n")
cat("\nSaved analysis data to: ", data_path, "\n", sep = "")
if (!resolution_mode){
  cat("Saved diagnostics and editable manual review files to: ", out, "\n", sep = "")
}else{
  cat("Saved final matching provenance, issues, and audit to: ", out, "\n", sep = "")
}
