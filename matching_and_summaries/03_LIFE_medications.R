#!/usr/bin/env Rscript
#
# Usage:
#   Rscript LIFE_medications.R [IDmatch.xlsx] [output_dir] [ATC_reference]
#
# Required R package: readxl. Optional: httr + jsonlite for online lookups.
# Outputs in the output directory:
#   LIFE_medications_by_ID.csv       keys + medication-derived columns only
#   medications_atc_lookup.csv       distinct ATC labels, matches, CNS status
#   medications_cns_decisions.csv    persistent, editable manual overrides
#   medications_cns_yes_qc.csv       all YES codes, labels, categories, reasons
#   medications_atc_unmatched.csv   unmatched/invalid ATC codes
#   medications_atc_summary.txt     quantitative matching and exclusion QC
#   medications_cns_narrative.txt   exclusions by category + actual drug labels

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L) {
  stop('Usage: Rscript LIFE_medications.R [IDmatch.xlsx] [output_dir] [ATC_reference]')
}
life <- '/data/p_02865/data/LIFE'
input <- if (length(args) >= 1L) args[[1L]] else file.path(life, 'tabular/IDmatch_recovered.xlsx')
out_dir <- if (length(args) >= 2L) args[[2L]] else
  file.path(life, 'derivatives/matching_and_summaries/medications')
ref_path <- if (length(args) >= 3L) args[[3L]] else ''
if (!nzchar(ref_path)) {
  default_ref <- file.path(life, 'tabular/atc_reference.csv')
  if (file.exists(default_ref)) ref_path <- default_ref
}
USE_ONLINE <- tolower(Sys.getenv('LIFE_ATC_ONLINE', '1')) %in% c('1', 'true', 'yes')
if (USE_ONLINE && (!requireNamespace('httr', quietly = TRUE) ||
                   !requireNamespace('jsonlite', quietly = TRUE))) {
  warning('httr/jsonlite not installed: online lookup disabled; local and cache used')
  USE_ONLINE <- FALSE
}
API_BASE <- 'https://tx.fhir.org/r4/CodeSystem/$lookup'
API_SYSTEM <- 'http://www.whocc.no/atc'

if (!requireNamespace('readxl', quietly = TRUE)) {
  stop("Install readxl: install.packages('readxl')")
}
if (!file.exists(input)) stop('IDmatch workbook not found: ', input)
if (nzchar(ref_path) && !file.exists(ref_path)) stop('ATC reference missing: ', ref_path)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out_dir)) stop('Could not create output directory: ', out_dir)

p <- function(name) file.path(out_dir, name)
read_csv <- function(path) read.csv(path, colClasses = 'character',
  check.names = FALSE, stringsAsFactors = FALSE,
  fileEncoding = 'UTF-8-BOM', na.strings = '')
write_csv <- function(x, path) write.csv(x, path, na = '',
  row.names = FALSE, fileEncoding = 'UTF-8')
normalize <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x)] <- NA_character_
  x
}
normalize_id <- function(x) {
  x <- toupper(normalize(x))
  x[!is.na(x) & x == 'NA'] <- NA_character_
  x
}
concat <- function(x, sep = ' | ') {
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (length(x)) paste(x, collapse = sep) else NA_character_
}
pct <- function(num, den) if (den > 0) sprintf('%.1f%%', 100 * num / den) else 'n/a'

# ATC hierarchy: A, A01, A01A, A01AA, A01AA01.
is_atc <- function(x) grepl('^[A-Z]([0-9]{2}([A-Z]([A-Z]([0-9]{2})?)?)?)?$', x)
parse_atc <- function(cell) {
  cell <- normalize(cell)
  if (is.na(cell) || identical(cell, '#')) return(character())
  parts <- unlist(strsplit(cell, '[#;,[:space:]]+', perl = TRUE), use.names = FALSE)
  parts <- unique(toupper(trimws(parts)))
  parts[nzchar(parts)]
}

# -------- 1. Read the same raw worksheet and identifiers as LIFE matching -----
raw <- as.data.frame(readxl::read_excel(input, sheet = 'Zusammenführen1',
  .name_repair = 'minimal'), check.names = FALSE,
  stringsAsFactors = FALSE)
if (anyDuplicated(names(raw))) stop('Duplicate headers in IDmatch sheet')
if (!'Column1' %in% names(raw)) stop('Missing Column1 (LIFE identifier)')

field <- function(suffix, required = TRUE) {
  hits <- grep(paste0('^Table3[._]', suffix, '$'), names(raw), value = TRUE)
  if (length(hits) > 1L) stop('Ambiguous Table3 field: ', suffix)
  if (!length(hits)) {
    if (required) stop('Required Table3 field not found: ', suffix)
    return(NA_character_)
  }
  hits
}
med_col <- field('ADULT_MEDA_H_ATC')
sic_col <- field('TEILNEHMER_SIC')
grp_col <- field('EEG_GRUPPE')
date_col <- field('EEG_DATUM')
# Matching code removes entirely empty Table3 rows and exact duplicates.
fields <- grep('^Table3[._]', names(raw), value = TRUE)
nonempty_table3 <- rep(FALSE, nrow(raw))
for (nm in fields) nonempty_table3 <- nonempty_table3 | !is.na(normalize(raw[[nm]]))
life_id_all <- normalize_id(raw[['Column1']])
valid <- !is.na(life_id_all) & nonempty_table3
# Duplicate check uses ALL columns, not only medications: do not drop visits
# with other differing clinical/EEG fields.
kept_idx <- which(valid)
if (length(kept_idx)) {
  cmp <- as.data.frame(lapply(raw[kept_idx, , drop = FALSE], normalize),
                       check.names = FALSE, stringsAsFactors = FALSE)
  cmp[['Column1']] <- life_id_all[kept_idx]
  duplicated_full <- duplicated(cmp)
  kept_idx <- kept_idx[!duplicated_full]
} else {
  duplicated_full <- logical()
}
raw <- raw[kept_idx, , drop = FALSE]
life_ids <- life_id_all[kept_idx]
source_rows <- kept_idx + 1L    # Excel rows include the first header row
row_atc <- lapply(raw[[med_col]], parse_atc)
all_codes <- sort(unique(unlist(row_atc, use.names = FALSE)))
valid_codes <- all_codes[is_atc(all_codes)]

# -------- 2. Local official ATC reference + previous cache + optional API -----
# A simple reference with atc_code,medication_name is sufficient.
read_ref <- function(path) {
  if (!nzchar(path)) return(setNames(character(), character()))
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c('xlsx', 'xls')) {
    tab <- as.data.frame(readxl::read_excel(path, col_names = FALSE,
      col_types = 'text', .name_repair = 'minimal'), check.names = FALSE)
  } else if (ext %in% c('csv', 'tsv', 'txt')) {
    # Try semicolon as well as comma; official German exports often use ';'.
    first <- readLines(path, n = 1L, warn = FALSE, encoding = 'UTF-8')
    sep <- if (ext == 'tsv') '\t' else
      if (length(first) && lengths(regmatches(first, gregexpr(';', first, fixed = TRUE))) >
          lengths(regmatches(first, gregexpr(',', first, fixed = TRUE)))) ';' else ','
    tab <- read.table(path, header = TRUE, sep = sep, quote = '"',
      colClasses = 'character', check.names = FALSE, stringsAsFactors = FALSE,
      fileEncoding = 'UTF-8-BOM', na.strings = '', comment.char = '', fill = TRUE)
  } else stop('Reference must be CSV, TSV, XLS or XLSX: ', path)

  header_key <- function(x) tolower(gsub('[^a-z0-9]', '',
    iconv(as.character(x), to = 'ASCII//TRANSLIT', sub = '')))
  select_cols <- function(hdr) {
    hdr <- header_key(hdr)
    i <- which(hdr %in% c('atc', 'atccode', 'code', 'atcschluessel'))
    j <- which(hdr %in% c('medicationname', 'atcname', 'atcbedeutung',
      'atcbezeichnung', 'wirkstoff', 'bezeichnung', 'display', 'name', 'bedeutung'))
    if (!length(i) || !length(j)) return(NULL)
    j <- setdiff(j, i[[1L]])
    if (!length(j)) return(NULL)
    c(i[[1L]], j[[1L]])
  }
  cols <- select_cols(names(tab))
  if (is.null(cols)) {
    for (i in seq_len(min(80L, nrow(tab)))) {
      hdr <- unlist(tab[i, , drop = TRUE], use.names = FALSE)
      found <- select_cols(hdr)
      if (length(found) == 2L) {
        cols <- found
        tab <- tab[-seq_len(i), , drop = FALSE]
        break
      }
    }
  }
  if (is.null(cols)) stop('Cannot identify ATC + name reference columns; use: atc_code,medication_name')
  codes <- toupper(normalize(tab[[cols[[1L]]]]))
  labels <- normalize(tab[[cols[[2L]]]])
  ok <- !is.na(codes) & !is.na(labels) & is_atc(codes)
  codes <- codes[ok]
  labels <- labels[ok]
  conflicting <- unique(codes)[vapply(unique(codes), function(x)
    length(unique(labels[codes == x])) > 1L, logical(1L))]
  if (length(conflicting)) stop('Conflicting names for ATC codes in reference: ',
                                paste(head(conflicting, 20L), collapse = ', '))
  keep <- !duplicated(codes)
  setNames(labels[keep], codes[keep])
}
ref <- read_ref(ref_path)
lookup_path <- p('medications_atc_lookup.csv')
cache <- setNames(character(), character())
cache_version <- setNames(character(), character())
if (file.exists(lookup_path)) {
  old <- read_csv(lookup_path)
  if (all(c('atc_code', 'medication_name', 'match_status') %in% names(old))) {
    good <- !is.na(old$medication_name) &
      old$match_status %in% c('matched_local', 'matched_cache', 'matched_api') &
      !duplicated(old$atc_code)
    cache <- setNames(old$medication_name[good], old$atc_code[good])
    if ('source_version' %in% names(old))
      cache_version <- setNames(old$source_version[good], old$atc_code[good])
  }
}

fhir_lookup <- function(code) {
  if (!requireNamespace('httr', quietly = TRUE) ||
      !requireNamespace('jsonlite', quietly = TRUE)) {
    stop("Install optional API packages: install.packages(c('httr','jsonlite'))")
  }
  last_error <- 'API request failed'
  for (attempt in 1:2) {
    response <- tryCatch(httr::GET(API_BASE,
      query = list(system = API_SYSTEM, code = code, `_format` = 'json'),
      httr::add_headers(Accept = 'application/fhir+json'), httr::timeout(12)),
      error = function(e) e)
    if (inherits(response, 'error')) {
      last_error <- conditionMessage(response)
    } else if (httr::status_code(response) != 200L) {
      last_error <- paste0('HTTP ', httr::status_code(response))
      if (httr::status_code(response) %in% c(400L, 404L, 422L)) break
    } else {
      body <- tryCatch(jsonlite::fromJSON(httr::content(response, as = 'text',
        encoding = 'UTF-8'), simplifyVector = FALSE), error = function(e) e)
      if (!inherits(body, 'error') && identical(body$resourceType, 'Parameters')) {
        pars <- body$parameter
        get_par <- function(name) {
          hit <- Filter(function(z) identical(z$name, name), pars)
          if (!length(hit) || is.null(hit[[1L]]$valueString)) return(NA_character_)
          as.character(hit[[1L]]$valueString)
        }
        label <- get_par('display')
        if (!is.na(label) && nzchar(label)) return(list(name = label,
          version = get_par('version'), status = 'matched_api', error = ''))
        last_error <- 'API returned no display name'
      } else last_error <- 'Unexpected FHIR API response'
    }
    if (attempt < 2L) Sys.sleep(0.4)
  }
  # An HTTP 404/400 may mean only this code is missing; a connection failure
  # indicates the entire service is unavailable.
  list(name = NA_character_, version = NA_character_, status = 'unmatched_api',
       error = last_error)
}

lookup <- data.frame(atc_code = all_codes,
  medication_name = rep(NA_character_, length(all_codes)),
  match_status = rep(NA_character_, length(all_codes)),
  source_version = rep(NA_character_, length(all_codes)),
  lookup_error = rep(NA_character_, length(all_codes)),
  stringsAsFactors = FALSE)
api_available <- USE_ONLINE
api_errors <- 0L
for (i in seq_len(nrow(lookup))) {
  code <- lookup$atc_code[[i]]
  if (code == 'VITMIN') {
    lookup$medication_name[[i]] <- 'Vitamins (unspecified)'
    lookup$match_status[[i]] <- 'special_placeholder'
  } else if (!is_atc(code)) {
    lookup$match_status[[i]] <- 'invalid_code'
  } else if (code %in% names(ref)) {
    lookup$medication_name[[i]] <- ref[[code]]
    lookup$match_status[[i]] <- 'matched_local'
    lookup$source_version[[i]] <- basename(ref_path)
  } else if (code %in% names(cache)) {
    lookup$medication_name[[i]] <- cache[[code]]
    lookup$match_status[[i]] <- 'matched_cache'
    if (code %in% names(cache_version))
      lookup$source_version[[i]] <- cache_version[[code]]
  } else if (api_available) {
    result <- fhir_lookup(code)
    lookup$medication_name[[i]] <- result$name
    lookup$match_status[[i]] <- result$status
    lookup$source_version[[i]] <- result$version
    lookup$lookup_error[[i]] <- result$error
    # Prevent large batches of failed requests on an offline compute node.
    if (grepl('^(API request failed|Could not resolve|Timeout|Failed to connect)',
              result$error, ignore.case = TRUE)) {
      api_errors <- api_errors + 1L
      if (api_errors >= 1L) api_available <- FALSE
    }
    if (api_available) Sys.sleep(0.12)
  } else {
    lookup$match_status[[i]] <- if (USE_ONLINE) 'api_unavailable' else 'not_in_reference'
  }
}

# -------- 3. LIBERAL SCREENING: only definite CNS exclusions ---------------
# Preregistration section 2.2.4:
# "Participants would be excluded if they were taking medications that affect
# the central nervous system, e.g., opioids, antidepressants, anti-epileptic
# drugs, sedatives, anxiolytics, and anti-Parkinsonian drugs."
#
# The aim is to KEEP participants unless an ATC code gives convincing evidence
# of a CNS-active drug in those categories or a very closely related class.
# No broad group-N screening, no inference from possible CNS side effects,
# and no default 'review' classification for every plausible drug.
# Exclusions are based on ATC drug classification, not verified exposure timing.
#
# N05A (antipsychotics), N06BA (ADHD stimulants), and CNS-directed dementia
# medications also have unequivocal central indications despite not being
# spelled out as examples in the preregistration's non-exhaustive 'e.g.' list.
# These explicit additions can be removed from the rules if a narrower
# definition is desired, without changing parsing or output.
CNS_YES_RULES <- data.frame(
  prefix = c('N02A', 'N03A', 'N04', 'N05A', 'N05B', 'N05C',
             'N06A', 'N06CA', 'N06BA', 'N06DA', 'N06DX01',
             'N07BC', 'N01AH',
             'R05DA01', 'R05DA03', 'R05DA04', 'R05DA05', 'R05DA06',
             'R05DA08', 'R05DA10', 'R05DA12', 'R05DA13', 'R05DA14'),
  reason = c(
    'Opioid analgesic (preregistration: opioids)',
    'Antiepileptic drug (preregistration)',
    'Anti-Parkinsonian drug (preregistration)',
    'Antipsychotic / unequivocally CNS-directed drug',
    'Anxiolytic (preregistration)',
    'Hypnotic or sedative (preregistration)',
    'Antidepressant (preregistration)',
    'Antidepressant/psycholeptic combination',
    'Centrally acting ADHD medication / psychostimulant',
    'CNS-directed anti-dementia cholinesterase inhibitor',
    'Memantine / CNS-directed anti-dementia medication',
    'Opioid substitution treatment (preregistration: opioids)',
    'Opioid anesthetic (preregistration: opioids)',
    rep('Opioid antitussive containing a definite opioid ingredient', 10L)
  ), stringsAsFactors = FALSE)

# Each rule has one explicit reporting category for participant-level QC.
# Use a consistent category even where several ATC branches represent a class.
CNS_YES_RULES$category <- c(
  'Opioids', 'Antiepileptics', 'Anti-Parkinsonian drugs',
  'Antipsychotics', 'Anxiolytics', 'Hypnotics and sedatives',
  'Antidepressants', 'Antidepressants',
  'ADHD psychostimulants', 'Anti-dementia drugs', 'Anti-dementia drugs',
  'Opioids', 'Opioids', rep('Opioids', 10L)
)
stopifnot(nrow(CNS_YES_RULES) == length(CNS_YES_RULES$category))

# For a LIBERAL exclusion strategy, borderline botanical treatments and
# melatonin products are deliberately not auto-excluded. This is a screening
# convention, not a claim of zero CNS effects. These exception prefixes must
# be applied BEFORE the broad N05B/N05C/N06A rule.
# German ATC codes include e.g. N05BP03 lavender oil, N05CP01 valerian,
# N06AP01 St. John's wort. WHO equivalents: N05BX05, N05CM09, N06AX25.
# N05CH01 melatonin and N05CX combination products are also outside automatic
# exclusion; other N05CH receptor agonists remain YES as definite hypnotics.
# Only well-identified R05DA opioid ingredients (not the entire antitussive
# class) are YES, while non-opioid cough suppressants are not.
CNS_NO_EXACT <- c('N05BX05', 'N05CM09', 'N06AX25', 'N05CH01')
CNS_NO_PREFIX <- c('N05BP', 'N05BH', 'N05CP', 'N06AP', 'N05CX')

screen_code <- function(code, label, status) {
  if (code == 'VITMIN')
    return(c('no', 'Unspecified vitamins: not a definite CNS exclusion'))
  if (!is_atc(code))
    return(c('unknown', 'Malformed / non-ATC code: NOT auto-excluded'))
  if (is.na(label) || !status %in% c('matched_local', 'matched_cache', 'matched_api'))
    return(c('unknown', 'ATC code could not be identified: NOT auto-excluded'))

  if (code %in% CNS_NO_EXACT || any(startsWith(code, CNS_NO_PREFIX))) {
    return(c('no', 'Liberal screening exception: botanical/melatonin preparation; NOT auto-excluded'))
  }

  hit <- which(startsWith(code, CNS_YES_RULES$prefix))
  if (length(hit)) {
    best <- hit[which.max(nchar(CNS_YES_RULES$prefix[hit]))]
    return(c('yes', CNS_YES_RULES$reason[[best]]))
  }

  # Crucial: all other identified drugs are 'no' for screening purposes,
  # even if they could have central side effects. No speculative exclusions.
  c('no', 'No definite CNS exclusion under liberal preregistration rules')
}
screen <- lapply(seq_len(nrow(lookup)), function(i)
  screen_code(lookup$atc_code[[i]], lookup$medication_name[[i]],
              lookup$match_status[[i]]))
lookup$suggested_cns <- vapply(screen, `[[`, character(1L), 1L)
lookup$cns_reason <- vapply(screen, `[[`, character(1L), 2L)

# Review/override file is persistent: manual values may be yes, no, or review.
review_path <- p('medications_cns_decisions.csv')
if (file.exists(review_path)) {
  decisions <- read_csv(review_path)
  if (!all(c('atc_code', 'final_cns', 'reviewer_note') %in% names(decisions)))
    stop('Existing medications_cns_decisions.csv has invalid structure')
  if (anyDuplicated(decisions$atc_code)) stop('Duplicate ATC codes in review CSV')
} else decisions <- data.frame(atc_code = character(),
  final_cns = character(), reviewer_note = character(), stringsAsFactors = FALSE)
new <- setdiff(lookup$atc_code, decisions$atc_code)
if (length(new)) decisions <- rbind(decisions, data.frame(atc_code = new,
  final_cns = NA_character_, reviewer_note = NA_character_, stringsAsFactors = FALSE))
decisions$final_cns <- tolower(normalize(decisions$final_cns))
if (any(!is.na(decisions$final_cns) & !decisions$final_cns %in%
        c('yes', 'no', 'review', 'unknown')))
  stop('Invalid final_cns (use yes/no/review/unknown or blank)')
match_idx <- match(decisions$atc_code, lookup$atc_code)
decisions$medication_name <- lookup$medication_name[match_idx]
decisions$suggested_cns <- lookup$suggested_cns[match_idx]
decisions$cns_reason <- lookup$cns_reason[match_idx]
decisions <- decisions[, c('atc_code', 'medication_name', 'suggested_cns',
  'cns_reason', 'final_cns', 'reviewer_note')]
write_csv(decisions, review_path)
manual <- decisions$final_cns[match(lookup$atc_code, decisions$atc_code)]
lookup$cns_exclusion <- ifelse(!is.na(manual), manual, lookup$suggested_cns)
lookup$cns_source <- ifelse(!is.na(manual), 'manual_override', 'liberal_definite_only_rule')
manual_note <- decisions$reviewer_note[match(lookup$atc_code, decisions$atc_code)]
lookup$cns_final_reason <- ifelse(!is.na(manual) & !is.na(manual_note) &
  nzchar(manual_note), manual_note, lookup$cns_reason)
# A YES override outside automatic YES groups needs a truthful reason, even
# if the reviewer has not entered an explanatory note yet.
manual_yes_without_reason <- !is.na(manual) & manual == 'yes' &
  lookup$suggested_cns != 'yes' & (is.na(manual_note) | !nzchar(manual_note))
lookup$cns_final_reason[manual_yes_without_reason] <-
  'Manual YES override outside automatic criteria; reviewer rationale missing'

# Category labels are for explanatory reporting, not extra exclusion rules.
# If a reviewer manually marks a code YES outside the predefined categories,
# keep it visible as an explicitly manual, uncategorized inclusion.
get_final_category <- function(code, final_cns) {
  if (is.na(final_cns) || final_cns != 'yes') return(NA_character_)
  if (code %in% CNS_NO_EXACT || any(startsWith(code, CNS_NO_PREFIX)))
    return('Other (manual YES override)')
  hit <- which(startsWith(code, CNS_YES_RULES$prefix))
  if (!length(hit)) return('Other (manual YES override)')
  CNS_YES_RULES$category[[hit[which.max(nchar(CNS_YES_RULES$prefix[hit]))]]]
}
lookup$cns_category <- vapply(seq_len(nrow(lookup)), function(i)
  get_final_category(lookup$atc_code[[i]], lookup$cns_exclusion[[i]]),
  character(1L))

# -------- 4. Participant/visit-row ATC summaries (NO large data CSV) --------
idx <- setNames(seq_len(nrow(lookup)), lookup$atc_code)
row_names <- function(tokens) {
  if (!length(tokens)) return(NA_character_)
  names <- lookup$medication_name[unname(idx[tokens])]
  names[is.na(names)] <- paste0('[unmatched: ', tokens[is.na(names)], ']')
  concat(names)
}
row_status <- function(tokens) {
  if (!length(tokens)) return(NA_character_)  # '#' is missing, not a negative
  statuses <- lookup$cns_exclusion[unname(idx[tokens])]
  if ('yes' %in% statuses) return('yes')
  if ('review' %in% statuses) return('review')  # only if manually requested
  if ('unknown' %in% statuses || anyNA(statuses)) return('unknown')
  'no'
}
row_yes <- function(tokens) {
  if (!length(tokens)) return(NA_character_)
  concat(tokens[lookup$cns_exclusion[unname(idx[tokens])] == 'yes'])
}

# -------- 5. Identifier-keyed medication columns (NO full data CSV) ---------
# Keep the existing join keys and one row per retained IDmatch workbook row.
values_or_na <- function(col) if (is.na(col)) rep(NA_character_, nrow(raw)) else normalize(raw[[col]])
row_class <- vapply(row_atc, row_status, character(1L))
output <- data.frame(
  life_id = life_ids,
  Column2 = if ('Column2' %in% names(raw)) normalize(raw[['Column2']]) else rep(NA_character_, nrow(raw)),
  sic = values_or_na(sic_col),
  group = values_or_na(grp_col),
  date_eeg = values_or_na(date_col),
  idmatch_excel_row = source_rows,
  medications_atc = normalize(raw[[med_col]]),
  medication_names = vapply(row_atc, row_names, character(1L)),
  cns_medication_exclusion = row_class,
  # 'review', 'unknown', and missing do not create an automatic exclusion.
  cns_exclude_definite = as.integer(row_class %in% 'yes'),
  cns_yes_atc_codes = vapply(row_atc, row_yes, character(1L)),
  stringsAsFactors = FALSE)
write_csv(output, p('LIFE_medications_by_ID.csv'))

# -------- 6. QC counts: codes, workbook rows AND unique participants --------
lookup$n_idmatch_rows <- vapply(lookup$atc_code, function(code)
  sum(vapply(row_atc, function(x) code %in% x, logical(1L))), integer(1L))
lookup$n_participants <- vapply(lookup$atc_code, function(code) {
  present <- vapply(row_atc, function(x) code %in% x, logical(1L))
  length(unique(life_ids[present]))
}, integer(1L))
lookup <- lookup[order(-lookup$n_idmatch_rows, lookup$atc_code), , drop = FALSE]
write_csv(lookup, lookup_path)
write_csv(lookup[lookup$match_status %in% c('invalid_code', 'unmatched_api',
  'api_unavailable', 'not_in_reference'), , drop = FALSE],
  p('medications_atc_unmatched.csv'))
yes_codes <- lookup[lookup$cns_exclusion == 'yes',
  c('atc_code', 'medication_name', 'cns_category', 'cns_final_reason',
    'cns_source', 'match_status', 'n_idmatch_rows', 'n_participants'), drop = FALSE]
write_csv(yes_codes, p('medications_cns_yes_qc.csv'))

# Build distinct life_id x YES ATC code pairs to avoid counting the same
# participant more than once per class or medicine, even if IDmatch has
# several different visit rows for that person.
rows_with_yes <- which(output$cns_exclude_definite == 1L)
if (length(rows_with_yes)) {
  yes_pairs <- do.call(rbind, lapply(rows_with_yes, function(i) {
    codes <- row_atc[[i]]
    codes <- codes[lookup$cns_exclusion[match(codes, lookup$atc_code)] == 'yes']
    if (!length(codes)) return(NULL)
    data.frame(life_id = rep(life_ids[[i]], length(codes)),
               atc_code = codes, stringsAsFactors = FALSE)
  }))
  yes_pairs <- unique(yes_pairs)
} else yes_pairs <- NULL
if (is.null(yes_pairs)) yes_pairs <- data.frame(
  life_id = character(), atc_code = character(), stringsAsFactors = FALSE)
yes_pairs$cns_category <- lookup$cns_category[match(yes_pairs$atc_code, lookup$atc_code)]

total_participants <- length(unique(life_ids))
excluded_participants <- length(unique(yes_pairs$life_id))
matched_statuses <- c('matched_local', 'matched_cache', 'matched_api')
observed <- lookup[lookup$atc_code %in% valid_codes, , drop = FALSE]
matched <- observed$match_status %in% matched_statuses
matched_unique <- sum(matched)
all_mentions <- sum(observed$n_idmatch_rows)
matched_mentions <- sum(observed$n_idmatch_rows[matched])
matched_yes_unique <- sum(observed$cns_exclusion == 'yes' & matched)
yes_mentions <- sum(observed$n_idmatch_rows[observed$cns_exclusion == 'yes' & matched])
row_has_yes <- output$cns_exclude_definite == 1L

summary_lines <- capture.output({
  cat('========== LIFE ATC medications / liberal CNS screening ==========\n')
  cat('Input workbook:             ', input, '\n', sep = '')
  cat('Sheet (same as LIFE matching): Zusammenführen1\n')
  cat('ATC field:                  ', med_col, '\n', sep = '')
  cat('ATC reference:              ', if (nzchar(ref_path)) ref_path else '(cache / API)', '\n', sep = '')
  cat('Online fallback:            ', USE_ONLINE, '\n', sep = '')
  cat('Original IDmatch rows:      ', nrow(raw) + sum(!valid) + sum(duplicated_full), '\n', sep = '')
  cat('Kept IDmatch rows:          ', nrow(output), '\n', sep = '')
  cat('Distinct LIFE participants:', total_participants, '\n')
  cat('ATC rows with codes:        ', sum(lengths(row_atc) > 0L), '\n', sep = '')
  cat('ATC rows without codes:     ', sum(lengths(row_atc) == 0L), '\n', sep = '')
  cat('Unique valid ATC codes:     ', nrow(observed), '\n', sep = '')
  cat('Successfully matched codes: ', matched_unique, ' / ', nrow(observed),
      ' (', pct(matched_unique, nrow(observed)), ')\n', sep = '')
  cat('Valid ATC code mentions:   ', all_mentions, '\n', sep = '')
  cat('Matched code mentions:     ', matched_mentions, ' / ', all_mentions,
      ' (', pct(matched_mentions, all_mentions), ')\n', sep = '')
  cat('Definite YES codes (matched): ', matched_yes_unique, ' / ', matched_unique,
      ' (', pct(matched_yes_unique, matched_unique), ')\n', sep = '')
  cat('YES code mentions:         ', yes_mentions, ' / ', matched_mentions,
      ' matched (', pct(yes_mentions, matched_mentions), ')\n', sep = '')
  cat('Total YES codes:           ', sum(observed$cns_exclusion == 'yes'), '\n', sep = '')
  cat('Manual CNS overrides:      ', sum(lookup$cns_source == 'manual_override'), '\n', sep = '')
  cat('IDmatch rows flagged:      ', sum(row_has_yes), '\n', sep = '')
  cat('Unique participants flagged:', excluded_participants, ' / ', total_participants,
      ' (', pct(excluded_participants, total_participants), ')\n', sep = '')
  cat('Special VITMIN tokens:      ', sum(vapply(row_atc, function(x) 'VITMIN' %in% x, logical(1L))), ' rows\n', sep = '')
  cat('Unmatched valid ATC codes: ', sum(!matched), '\n', sep = '')
  cat('Invalid ATC tokens:        ', sum(lookup$match_status == 'invalid_code'), '\n', sep = '')
  cat('\nCNS status by retained IDmatch row:\n')
  print(table(ifelse(is.na(output$cns_medication_exclusion), 'missing',
                     output$cns_medication_exclusion)))
  cat('\nOnly unambiguous CNS medication classes are automatic YES exclusions.\n')
  cat('Unknown and missing medication codes are NOT automatically exclusions.\n')
  cat('Manual decisions persist in medications_cns_decisions.csv.\n')
  cat('Category breakdown and medication labels: medications_cns_narrative.txt\n')
})
writeLines(summary_lines, p('medications_atc_summary.txt'), useBytes = TRUE)
cat(paste(summary_lines, collapse = '\n'), '\n', sep = '')

# -------- 7. Narrative report: per-category unique participant counts ------
# A participant taking antidepressants and anxiolytics appears in both groups,
# but is counted ONCE in the grand total. Ingredient names are used as
# returned by the ATC lookup; an ATC class-level entry may have a class label.
category_order <- c('Opioids', 'Antidepressants', 'Antiepileptics',
  'Hypnotics and sedatives', 'Anxiolytics', 'Anti-Parkinsonian drugs',
  'Antipsychotics', 'ADHD psychostimulants', 'Anti-dementia drugs',
  'Other (manual YES override)')
categories <- unique(yes_pairs$cns_category)
categories <- c(category_order[category_order %in% categories],
                sort(setdiff(categories, category_order)))

category_descriptions <- c(
  'Opioids' = 'opioid medications',
  'Antidepressants' = 'antidepressants',
  'Antiepileptics' = 'antiepileptic drugs',
  'Hypnotics and sedatives' = 'hypnotics or sedatives',
  'Anxiolytics' = 'anxiolytics',
  'Anti-Parkinsonian drugs' = 'anti-Parkinsonian drugs',
  'Antipsychotics' = 'antipsychotic drugs',
  'ADHD psychostimulants' = 'centrally acting ADHD psychostimulants',
  'Anti-dementia drugs' = 'CNS-directed anti-dementia drugs',
  'Other (manual YES override)' = 'other manually confirmed medications'
)

narrative <- capture.output({
  cat('LIFE MEDICATION EXCLUSIONS - DESCRIPTIVE NARRATIVE\n')
  cat('=================================================\n\n')
  cat(sprintf('We flagged %d of %d unique participants (%s) for exclusion based on\n',
              excluded_participants, total_participants,
              pct(excluded_participants, total_participants)))
  cat('definite CNS-relevant medications recorded.\n')
  cat('These are medication-screening exclusions, not a record that participants\n')
  cat('have already been dropped from the subsequent analysis.\n\n')
  cat('Category counts represent unique participants in that category.\n')
  cat('A participant may contribute to multiple categories, while the overall\n')
  cat('total counts each LIFE ID only once.\n\n')
  if (!length(categories)) cat('No definite YES medication categories found.\n')
  for (category in categories) {
    these <- yes_pairs[yes_pairs$cns_category == category, , drop = FALSE]
    n_people <- length(unique(these$life_id))
    action <- if (category %in% names(category_descriptions))
      category_descriptions[[category]] else category
    cat(sprintf('%d participants were flagged due to %s.\n', n_people, action))
    cat('The specific ATC medication labels in this category were:\n')
    codes <- unique(these$atc_code)
    counts <- vapply(codes, function(code)
      length(unique(these$life_id[these$atc_code == code])), integer(1L))
    codes <- codes[order(-counts, codes)]
    for (code in codes) {
      item <- lookup[match(code, lookup$atc_code), , drop = FALSE]
      label <- item$medication_name[[1L]]
      if (is.na(label) || !nzchar(label)) label <- paste0('[unmatched: ', code, ']')
      n <- length(unique(these$life_id[these$atc_code == code]))
      cat(sprintf('  - %s (ATC %s): %d participant%s\n',
                  label, code, n, if (n == 1L) '' else 's'))
    }
    cat('\n')
  }
  cat('Interpretation note: ATC labels describe substances or therapeutic\n')
  cat('classes, not necessarily commercial brand names. Group-level ATC\n')
  cat('codes may have category labels rather than individual ingredients.\n')
  cat('The underlying medication inventory does not establish exposure dates.\n')
  cat('Only definite YES drugs are counted; other medications remain included.\n')
  cat('Review the full YES list and reasons in medications_cns_yes_qc.csv.\n')
})
writeLines(narrative, p('medications_cns_narrative.txt'), useBytes = TRUE)
cat('Narrative report: ', p('medications_cns_narrative.txt'), '\n', sep = '')
cat('Saved outputs in: ', out_dir, '\n', sep = '')
