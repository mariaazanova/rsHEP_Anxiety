#!/usr/bin/env Rscript
# 06b_LIFE_match_add_controls_after_exclusion.R
# Run AFTER 06_exclusion_common_only_updated.R and BEFORE 07_LIFE_pipeline_matching_clean.R.
#
# 1. Reproduce the conservative matching procedure of 09 against PV0696_datajoin.
# 2. Keep only high-confidence reciprocal one-to-one confirmed_unique SIC matches.
# 3. Join *additional* control variables from the separate PV0696_*_NODUP.xlsx files:
#      echo_lvedv_mod_bp, echo_lvesv_mod_bp, lvedv_lvesv, height, diabetes_final.
# 4. Preserve all existing cells and rows. Write a pre-join snapshot, then add the
#    five columns to data_after_exclusion.csv, which is read by script 07.
# 5. Produce an audit comparing valid/nonempty measurements in the complete LIFE
#    source files with coverage in the final matched post-exclusion participants.
#
# If 06 left several rows for a life_id, construct a consensus for matching:
# conflicting fields become unavailable for matching. Do NOT confirm a LIFE ID
# with ANY conflicting nonmissing duplicate-row fields (conservative safety gate).
# Those IDs remain present in the final data with empty newly-added controls.
#
# Diabetes definition: ADA_DS (binary 0/1), T01404/ADA2021; never DS_COMB3.
# Height 996, 998, 999 are documented invalid codes and treated as missing.
# Stroke-volume difference requires both MOD_BP volumes observed.
#
# Usage: Rscript 06b_LIFE_match_add_controls_after_exclusion.R [life_directory]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 1L) stop("Usage: Rscript 06b_LIFE_match_add_controls_after_exclusion.R [LIFE_directory]")
life      <- if (length(args)) args[[1L]] else "/data/p_02865/data/LIFE"
deriv     <- file.path(life, "derivatives")
base_out  <- file.path(deriv, "matching_and_summaries")
a_path    <- file.path(deriv, "data_after_exclusion.csv")
map_path  <- file.path(base_out, "IDmatch_variables_manual.csv")
dj_path   <- file.path(life, "tabular/Behavior/PV0696_FILES/data/PV0696_datajoin.xlsx")
prev_path <- file.path(base_out, "datajoin_diagnostic", "datajoin_participant_candidates.csv")
out       <- file.path(base_out, "datajoin_after_exclusion")
for (p in c(a_path, map_path, dj_path)) if (!file.exists(p)) stop("Missing input: ", p)
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl: install.packages('readxl')")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- Settings ----------------------------------------------------------------
# Dates: data_after_exclusion.csv holds ISO dates on the TRUE calendar (readxl applied the
# 1904 flag of IDmatch). Datajoin serials are read with the 1900 origin and shifted; the
# shift (-1462 / 0 / +1462) is chosen from the data, per date column, with a global fallback.
A_SERIAL_ORIGIN  <- "1904-01-01"   # only used if a raw serial number is found in data_after_exclusion
DJ_SERIAL_ORIGIN <- "1899-12-30"
SHIFTS           <- c(-1462L, 0L, 1462L)
DATE_MIN_OVERLAP <- 5L
DATE_WINNER_RATIO <- 1.5
# Candidate generation (values shared by too many rows are not used to PROPOSE pairs)
MAX_KEY_ROWS        <- 30L        # max rows per value in either file
MAX_PAIRS_PER_FIELD <- 300000     # skip a field if it would propose more pairs
MIN_ANCHOR_HITS     <- 2L         # a pair must be proposed by >= this many columns
MAX_TOTAL_PAIRS     <- 3000000    # keep the pairs with most anchor hits if above this
RARE_MAX            <- 5L         # value counts as "specific" if in <= 5 rows in both files
# Decision thresholds (sweep table shows sensitivity)
MIN_AGREE   <- 15L                # agreeing columns
MAX_CONF    <- 0.15               # conflicting / jointly observed columns
MIN_MARGIN  <- 10                 # score gap to the 2nd best candidate (both directions)
UNRELIABLE_CONFLICT <- 0.30       # column dropped in pass 2 if > 30% of confirmed pairs disagree
UNRELIABLE_MIN_N    <- 30L
# Never used for matching (derived / pipeline columns of data_after_exclusion)
DERIVED_NAMES <- c("life_id", "Column2", "non-unique", "common_exclude", "overall_exclude",
  "exclusion_reasons", "gad7_group", "GAD7", "IDS", "IDS_n_items", "all_GAD7", "all_IDS",
  "drowsy_15min", "extrasystolic_N", "extrasystolic_proportion", "extrasystolic_N_SD",
  "extrasystolic_exclusion", "ecg_n_beats", "ecg_v_beats", "ecg_s_beats",
  "medications_cns_exclude_definite")
# Columns created in 04_LIFE_generate_data_manual_review.R from other workbooks, renamed there
ALIASES <- data.frame(
  A_col = c("date_gad7", "date_medanameeg", "Intrakranielle_Blutung", "schaedel_hirn_trauma",
            "schaedelfraktur", "meningitis", "head_brain_op", "date_echo", "age_echo", "echo_lvot_sv"),
  B_col = c("GAD7_DATUM", "MEDANAMEEG_DATUM", "MEDANAMEEG_MAE_IK", "MEDANAMEEG_MAE_SHT",
            "MEDANAMEEG_MAE_SF", "MEDANAMEEG_MAE_ME", "MEDANAMEEG_MAE_OP", "ECHO_READ_STUDYDATETIME",
            "ECHO_READ_AGE", "ECHO_READ_LVOT_SV"), stringsAsFactors = FALSE)

# ---- Helpers -----------------------------------------------------------------
norm <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | !nzchar(x) | toupper(x) %in% c("NA", "<NA>", "NULL")] <- NA_character_
  x
}
decode <- function(x) gsub("+AF8-", "_", x, fixed = TRUE)
write_csv <- function(x, f) write.csv(x, file.path(out, f), row.names = FALSE, na = "", fileEncoding = "UTF-8")
pct <- function(x, n) sprintf("%d (%.1f%%)", x, if (n > 0) 100 * x / n else 0)
norm_value <- function(x) {
  x <- toupper(gsub("\\s+", " ", norm(x)))
  num <- !is.na(x) & grepl("^[-+]?[0-9]*[.,]?[0-9]+$", x)
  if (any(num)) {
    n <- suppressWarnings(as.numeric(sub(",", ".", x[num])))
    x[num] <- ifelse(is.finite(n), as.character(signif(n, 10)), x[num])
  }
  x
}
parse_dates <- function(x, origin) {
  x <- norm(x); res <- rep(NA_character_, length(x))
  is_serial <- !is.na(x) & grepl("^[0-9]+(\\.[0-9]+)?$", x)
  s <- which(is_serial)
  if (length(s)) {
    num <- suppressWarnings(as.numeric(x[s]))
    ok <- is.finite(num) & num >= 1 & num <= 100000
    res[s[ok]] <- as.character(as.Date(floor(num[ok]), origin = origin))
  }
  rest <- which(!is.na(x) & is.na(res) & !is_serial)
  for (fmt in c("%Y-%m-%d", "%Y/%m/%d", "%d.%m.%Y", "%d/%m/%Y", "%d-%m-%Y")) {
    if (!length(rest)) break
    d <- suppressWarnings(as.Date(substr(x[rest], 1L, 10L), format = fmt))
    good <- !is.na(d)
    if (any(good)) res[rest[good]] <- as.character(d[good])
    rest <- rest[!good]
  }
  ymd <- which(!is.na(x) & is.na(res) & grepl("^[12][0-9]{7}$", x))
  if (length(ymd)) res[ymd] <- as.character(suppressWarnings(as.Date(x[ymd], format = "%Y%m%d")))
  res[!is.na(res) & (res < "1900-01-01" | res > "2100-12-31")] <- NA_character_
  res
}
shift_dates <- function(x, days) {
  if (days == 0L) return(x)
  i <- which(!is.na(x)); r <- x
  if (length(i)) r[i] <- as.character(as.Date(x[i]) + as.integer(days))
  r
}
overlap <- function(a, b) {
  ta <- table(a[!is.na(a)]); tb <- table(b[!is.na(b)])
  k <- intersect(names(ta), names(tb))
  if (!length(k)) 0L else as.integer(sum(pmin(as.integer(ta[k]), as.integer(tb[k]))))
}

# ---- 1. Read inputs ----------------------------------------------------------
data_full <- read.csv(a_path, colClasses = "character", check.names = FALSE,
                      stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", na.strings = "")
names(data_full) <- trimws(names(data_full))
if (anyDuplicated(names(data_full))) stop("Duplicate headers in data_after_exclusion.csv")
if (!"life_id" %in% names(data_full)) stop("data_after_exclusion.csv needs life_id")
NEW_COLS <- c("echo_lvedv_mod_bp", "echo_lvesv_mod_bp", "lvedv_lvesv",
              "height", "diabetes_final")
if (any(NEW_COLS %in% names(data_full))) {
  stop("New controls already exist in data_after_exclusion.csv: ",
       paste(intersect(NEW_COLS, names(data_full)), collapse = ", "),
       ". Re-run script 06 to regenerate a clean post-exclusion dataset.")
}
full_ids <- toupper(norm(data_full$life_id))
if (anyNA(full_ids)) stop("Missing or blank life_id in post-exclusion data")
if (!nrow(data_full)) stop("data_after_exclusion.csv has no rows")
# Script 06 sometimes retains unresolved multiple rows per participant. A is an
# analysis-only one-row-per-ID matching view. We do NOT change these data rows.
uid <- unique(full_ids)
original_first_rows <- match(uid, full_ids) + 1L
id_rows <- split(seq_along(full_ids), factor(full_ids, levels = uid))
A <- data_full[match(uid, full_ids), , drop = FALSE]
A[] <- lapply(A, norm)
A$life_id <- uid
repeated <- which(lengths(id_rows) > 1L)
duplicate_conflicts <- list()
for (k in repeated) {
  ix <- id_rows[[k]]
  for (column in setdiff(names(A), "life_id")) {
    vals <- unique(norm(data_full[[column]][ix]))
    vals <- vals[!is.na(vals)]
    if (length(vals) > 1L) {
      A[[column]][k] <- NA_character_  # do not match on a conflicting field
      duplicate_conflicts[[length(duplicate_conflicts) + 1L]] <- data.frame(
        life_id = uid[k], column = column, values = paste(vals, collapse = " | "),
        original_csv_rows = paste(ix + 1L, collapse = " | "), stringsAsFactors = FALSE)
    } else if (length(vals) == 1L) {
      A[[column]][k] <- vals
    } else {
      A[[column]][k] <- NA_character_
    }
  }
}
dup_qc <- if (length(duplicate_conflicts)) do.call(rbind, duplicate_conflicts) else
  data.frame(life_id = character(), column = character(), values = character(),
             original_csv_rows = character(), stringsAsFactors = FALSE)
conflicting_duplicate_ids <- unique(dup_qc$life_id)
write_csv(dup_qc, "datajoin_after_exclusion_duplicate_input_conflicts.csv")
nA <- nrow(A)  # unique post-exclusion participants, not raw CSV rows
dj <- as.data.frame(readxl::read_excel(dj_path, sheet = 1, col_types = "text", .name_repair = "minimal"),
                    check.names = FALSE, stringsAsFactors = FALSE)
names(dj) <- trimws(names(dj))
if (anyDuplicated(names(dj))) stop("Duplicate headers in datajoin")
if (!"SIC" %in% names(dj)) stop("Datajoin needs SIC")
dj[] <- lapply(dj, norm); dj$SIC <- toupper(dj$SIC)
if (anyNA(dj$SIC) || anyDuplicated(dj$SIC))
  stop("Datajoin SIC values must be nonempty and unique for safe 1-to-1 matching")
nB <- nrow(dj)

manual <- read.csv(map_path, colClasses = "character", check.names = FALSE, stringsAsFactors = FALSE,
                   fileEncoding = "UTF-8-BOM", na.strings = "")
names(manual) <- decode(names(manual)); manual[] <- lapply(manual, decode); manual[] <- lapply(manual, norm)
if (!all(c("original_column", "keep", "new_name") %in% names(manual))) stop("Mapping lacks original_column/keep/new_name")
manual$keep[is.na(manual$keep)] <- "0"
forced <- c("Table3.TEILNEHMER_SIC", "Table3.EEG_GRUPPE", "Table3.EEG_DATUM")  # as in script 04
sel <- manual[!is.na(manual$new_name) & !is.na(manual$original_column) &
              (manual$keep == "1" | manual$original_column %in% forced), , drop = FALSE]

# ---- 2. Column map: data_after_exclusion name  <->  datajoin name -----------------
empty_map <- function() data.frame(A_col = character(), B_col = character(), source = character(),
                                   stringsAsFactors = FALSE)
m1 <- data.frame(A_col = sel$new_name, B_col = sub("^Table3[._]", "", sel$original_column),
                 source = rep("manual_mapping", nrow(sel)), stringsAsFactors = FALSE)
m1 <- m1[m1$A_col %in% names(A) & m1$B_col %in% names(dj), , drop = FALSE]
m1 <- m1[!duplicated(m1$A_col) & !duplicated(m1$B_col), , drop = FALSE]
al <- ALIASES[ALIASES$A_col %in% names(A) & ALIASES$B_col %in% names(dj) &
              !ALIASES$A_col %in% m1$A_col & !ALIASES$B_col %in% m1$B_col, , drop = FALSE]
m2 <- if (nrow(al)) data.frame(A_col = al$A_col, B_col = al$B_col, source = "alias_table",
                               stringsAsFactors = FALSE) else empty_map()
same <- setdiff(intersect(names(A), names(dj)), c(m1$A_col, m2$A_col, m1$B_col, m2$B_col, "life_id", "Column2"))
m3 <- if (length(same)) data.frame(A_col = same, B_col = same, source = "same_name",
                                   stringsAsFactors = FALSE) else empty_map()
cmap <- rbind(m1, m2, m3)
cmap$is_id_like <- grepl("(^|_)(SIC|ID|PSEUDONYM)($|_)", cmap$B_col, ignore.case = TRUE)
cmap$is_derived <- cmap$A_col %in% DERIVED_NAMES
cmap$nonmissing_A <- vapply(cmap$A_col, function(c) sum(!is.na(A[[c]])), integer(1))
cmap$nonmissing_B <- vapply(cmap$B_col, function(c) sum(!is.na(dj[[c]])), integer(1))
cmap$is_date <- grepl("DATUM|STUDYDATETIME", cmap$B_col, ignore.case = TRUE)
cmap$skip_reason <- ifelse(cmap$is_id_like, "identifier_like",
                    ifelse(cmap$is_derived, "derived_column",
                    ifelse(cmap$nonmissing_A == 0L, "empty_in_data_after_exclusion",
                    ifelse(cmap$nonmissing_B == 0L, "empty_in_datajoin", ""))))
cmap$use <- cmap$skip_reason == ""
um <- cmap[cmap$use, , drop = FALSE]; K <- nrow(um)
if (K < 5L) stop("Fewer than 5 usable mapped columns (", K, "): check IDmatch_variables_manual.csv")

# Columns of data_after_exclusion without counterpart, and datajoin columns not in data_after_exclusion
unm <- setdiff(names(A), cmap$A_col)
unm_tab <- data.frame(column = unm, nonmissing = vapply(unm, function(c) sum(!is.na(A[[c]])), integer(1)),
  reason = ifelse(unm %in% DERIVED_NAMES, "derived/pipeline column",
           ifelse(unm %in% sel$new_name & !sub("^Table3[._]", "", sel$original_column[match(unm, sel$new_name)]) %in% names(dj),
                  "IDmatch-only column (original not in datajoin)",
                  "no counterpart (other source / not in datajoin)")), stringsAsFactors = FALSE)
dj_new <- setdiff(names(dj), c(cmap$B_col, "SIC"))
dj_new_tab <- data.frame(datajoin_column = dj_new,
  n_nonmissing = vapply(dj_new, function(c) sum(!is.na(dj[[c]])), integer(1)),
  stringsAsFactors = FALSE)
write_csv(cmap, "datajoin_after_exclusion_column_map.csv")
write_csv(unm_tab, "datajoin_after_exclusion_unmapped_columns.csv")
write_csv(dj_new_tab, "datajoin_columns_not_in_data_after_exclusion.csv")

# ---- 3. Dates: choose the datajoin shift from the data ------------------------------
dkeys <- which(um$is_date)
a_iso <- b_raw <- vector("list", K)
ov <- matrix(0L, nrow = length(dkeys), ncol = length(SHIFTS), dimnames = list(um$A_col[dkeys], SHIFTS))
for (r in seq_along(dkeys)) {
  k <- dkeys[r]
  a_iso[[k]] <- parse_dates(A[[um$A_col[k]]], A_SERIAL_ORIGIN)
  b_raw[[k]] <- parse_dates(dj[[um$B_col[k]]], DJ_SERIAL_ORIGIN)
  for (s in seq_along(SHIFTS)) ov[r, s] <- overlap(a_iso[[k]], shift_dates(b_raw[[k]], SHIFTS[s]))
}
global_shift <- 1462L
if (length(dkeys) && max(colSums(ov)) > 0) global_shift <- SHIFTS[which.max(colSums(ov))]
col_shift <- rep(global_shift, K); shift_why <- rep("global", K)
for (r in seq_along(dkeys)) {
  o <- sort(ov[r, ], decreasing = TRUE)
  lead <- as.integer(names(o)[1L])
  if (o[1L] >= DATE_MIN_OVERLAP && o[1L] > o[2L] && o[1L] >= DATE_WINNER_RATIO * max(1, o[2L])) {
    col_shift[dkeys[r]] <- lead; shift_why[dkeys[r]] <- "column_winner"
  }
}
shift_tab <- data.frame(column = um$A_col[dkeys], datajoin_column = um$B_col[dkeys],
  overlap_minus1462 = ov[, 1L], overlap_0 = ov[, 2L], overlap_plus1462 = ov[, 3L],
  applied_shift = col_shift[dkeys], decided_by = shift_why[dkeys], row.names = NULL)
write_csv(shift_tab, "datajoin_after_exclusion_date_shifts.csv")

# ---- 4. Normalised values -> integer codes ---------------------------------------
CA <- CB <- FA <- FB <- vector("list", K)
for (k in seq_len(K)) {
  if (um$is_date[k]) {
    a <- a_iso[[k]]; b <- shift_dates(b_raw[[k]], col_shift[k])
  } else {
    a <- norm_value(A[[um$A_col[k]]]); b <- norm_value(dj[[um$B_col[k]]])
  }
  lev <- unique(c(a[!is.na(a)], b[!is.na(b)]))
  CA[[k]] <- match(a, lev); CB[[k]] <- match(b, lev)
  FA[[k]] <- tabulate(CA[[k]][!is.na(CA[[k]])], length(lev))
  FB[[k]] <- tabulate(CB[[k]][!is.na(CB[[k]])], length(lev))
}

# ---- 5. Candidate pairs from selective values (counted per pair) ------------------
acc_key <- numeric(); acc_n <- numeric(); buf <- list()
collapse <- function() {
  if (!length(buf)) return(invisible())
  allk <- c(acc_key, unlist(buf, use.names = FALSE))
  alln <- c(acc_n, rep(1, length(allk) - length(acc_key)))
  o <- order(allk); allk <- allk[o]; alln <- alln[o]
  r <- rle(allk); ends <- cumsum(r$lengths); cs <- cumsum(alln)
  acc_key <<- r$values
  acc_n <<- cs[ends] - c(0, cs[ends[-length(ends)]])
  buf <<- list()
}
anchor_log <- vector("list", K)
for (k in seq_len(K)) {
  ta <- FA[[k]]; tb <- FB[[k]]
  okc <- which(ta > 0 & tb > 0 & ta <= MAX_KEY_ROWS & tb <= MAX_KEY_ROWS)
  nprop <- if (length(okc)) sum(as.numeric(ta[okc]) * tb[okc]) else 0
  why <- if (!length(okc)) "no_usable_shared_value" else if (nprop > MAX_PAIRS_PER_FIELD) "too_many_pairs" else ""
  anchor_log[[k]] <- data.frame(column = um$A_col[k], usable_values = length(okc),
                                proposed_pairs = nprop, skipped = why, stringsAsFactors = FALSE)
  if (nzchar(why)) next
  ia <- which(CA[[k]] %in% okc); ib <- which(CB[[k]] %in% okc)
  m <- merge(data.frame(i = ia, code = CA[[k]][ia]), data.frame(j = ib, code = CB[[k]][ib]), by = "code")
  buf[[length(buf) + 1L]] <- (m$i - 1) * nB + m$j
  if (length(buf) >= 40L) collapse()
}
collapse()
anchor_tab <- do.call(rbind, anchor_log)
write_csv(anchor_tab, "datajoin_after_exclusion_anchor_qc.csv")
if (!length(acc_key)) stop("No candidate pairs could be proposed. See datajoin_after_exclusion_anchor_qc.csv")
rows_any_pair <- length(unique((acc_key - 1) %/% nB + 1))
keep <- acc_n >= MIN_ANCHOR_HITS
if (!any(keep)) stop("No pair is proposed by >= ", MIN_ANCHOR_HITS, " columns; lower MIN_ANCHOR_HITS")
ak <- acc_key[keep]; an <- acc_n[keep]
if (length(ak) > MAX_TOTAL_PAIRS) { o <- order(-an)[seq_len(MAX_TOTAL_PAIRS)]; ak <- ak[o]; an <- an[o] }
PI <- as.integer((ak - 1) %/% nB + 1); PJ <- as.integer((ak - 1) %% nB + 1)
rows_with_pair <- length(unique(PI))

# ---- 6. Scoring and ranking ------------------------------------------------------
score_pairs <- function(cols) {
  n <- length(PI)
  agree <- conf <- dag <- dcf <- spec <- integer(n)
  for (k in cols) {
    a <- CA[[k]][PI]; b <- CB[[k]][PJ]
    both <- !is.na(a) & !is.na(b)
    eq <- both & (a == b); ne <- both & (a != b)
    agree <- agree + eq; conf <- conf + ne
    if (um$is_date[k]) { dag <- dag + eq; dcf <- dcf + ne }
    rare <- eq & (FA[[k]][a] <= RARE_MAX) & (FB[[k]][a] <= RARE_MAX)
    spec <- spec + (rare %in% TRUE)
  }
  nb <- agree + conf
  data.frame(i = PI, j = PJ, hits = an, agree = agree, conf = conf, n_both = nb,
             conf_frac = ifelse(nb > 0, conf / nb, 1), date_agree = dag, date_conf = dcf,
             specific = spec, score = agree + 2 * spec - 5 * conf)
}
rank_best <- function(P, by) {
  other <- if (by == "i") "j" else "i"
  o <- order(P[[by]], -P$score, -P$agree); Q <- P[o, , drop = FALSE]
  g <- Q[[by]]; first <- !duplicated(g)
  same_next <- c(g[-1L] == g[-length(g)], FALSE)
  nxt_score <- c(Q$score[-1L], NA_real_); nxt_other <- c(Q[[other]][-1L], NA)
  b <- Q[first, , drop = FALSE]
  has2 <- same_next[first]
  b$second_score <- ifelse(has2, nxt_score[first], NA_real_)
  b$second_other <- ifelse(has2, nxt_other[first], NA)
  b$margin <- b$score - pmax(ifelse(is.na(b$second_score), 0, b$second_score), 0)
  b
}
status_fun <- function(b, min_agree, max_conf, min_margin) {
  ok <- b$agree >= min_agree & b$conf_frac <= max_conf
  ifelse(!ok, "weak_or_conflicting",
  ifelse(b$margin < min_margin, "ambiguous_close_candidates",
  ifelse(!b$recip, "datajoin_row_fits_other_participant_better",
  ifelse(b$margin_B < min_margin, "ambiguous_close_candidates", "confirmed_unique"))))
}
run_matching <- function(cols) {
  P <- score_pairs(cols)
  bA <- rank_best(P, "i"); bB <- rank_best(P, "j")
  mB <- match(bA$j, bB$j)
  bA$recip <- bB$i[mB] == bA$i
  bA$margin_B <- bB$margin[mB]
  bA$status <- status_fun(bA, MIN_AGREE, MAX_CONF, MIN_MARGIN)
  list(P = P, bA = bA, bB = bB)
}
col_agreement <- function(pi, pj, cols) {
  do.call(rbind, lapply(cols, function(k) {
    a <- CA[[k]][pi]; b <- CB[[k]][pj]
    both <- !is.na(a) & !is.na(b)
    data.frame(column = um$A_col[k], datajoin_column = um$B_col[k], source = um$source[k],
      n_both = sum(both), n_agree = sum(both & a == b), n_conflict = sum(both & a != b),
      conflict_rate = ifelse(sum(both) > 0, round(sum(both & a != b) / sum(both), 3), NA),
      A_missing_datajoin_present = sum(is.na(a) & !is.na(b)),
      A_present_datajoin_missing = sum(!is.na(a) & is.na(b)), stringsAsFactors = FALSE)
  }))
}
all_cols <- seq_len(K)
r1 <- run_matching(all_cols)
conf1 <- r1$bA[r1$bA$status == "confirmed_unique" &
                !(A$life_id[r1$bA$i] %in% conflicting_duplicate_ids), ]
ca1 <- col_agreement(conf1$i, conf1$j, all_cols)
unrel <- ca1$n_both >= UNRELIABLE_MIN_N & !is.na(ca1$conflict_rate) & ca1$conflict_rate > UNRELIABLE_CONFLICT
r2 <- if (any(unrel)) run_matching(all_cols[!unrel]) else r1
bA <- r2$bA; bB <- r2$bB
# Conservative: if source rows for a participant disagree, do not silently
# confirm that participant from a consensus of potentially different records.
if (length(conflicting_duplicate_ids))
  bA$status[A$life_id[bA$i] %in% conflicting_duplicate_ids] <- "duplicate_input_conflict"
final_cols <- if (any(unrel)) all_cols[!unrel] else all_cols
conf <- bA[bA$status == "confirmed_unique", ]
ca <- col_agreement(conf$i, conf$j, all_cols)
ca$used_in_final_scoring <- all_cols %in% final_cols
write_csv(ca, "datajoin_after_exclusion_column_agreement.csv")

# ---- 7. Result tables ---------------------------------------------------------
mi <- match(seq_len(nA), bA$i)
res <- data.frame(life_id = A$life_id, data_after_exclusion_row = original_first_rows,
  status = ifelse(is.na(mi), "no_candidate", bA$status[mi]),
  SIC = dj$SIC[bA$j[mi]], datajoin_excel_row = bA$j[mi] + 1L,
  n_agree = bA$agree[mi], n_conflict = bA$conf[mi], conflict_frac = round(bA$conf_frac[mi], 3),
  date_agree = bA$date_agree[mi], date_conflict = bA$date_conf[mi], specific_agree = bA$specific[mi],
  score = bA$score[mi], second_best_score = bA$second_score[mi], margin = bA$margin[mi],
  second_best_SIC = dj$SIC[bA$second_other[mi]], datajoin_row_best_for_this_participant = bA$recip[mi],
  margin_from_datajoin_side = bA$margin_B[mi], stringsAsFactors = FALSE)
write_csv(res, "datajoin_after_exclusion_participant_matches.csv")
cw <- res[res$status == "confirmed_unique", c("life_id", "SIC", "datajoin_excel_row", "n_agree",
                                             "n_conflict", "conflict_frac", "margin")]
stopifnot(!anyNA(cw$SIC), !anyNA(cw$life_id),
          !anyDuplicated(cw$SIC), !anyDuplicated(cw$life_id))
write_csv(cw, "datajoin_after_exclusion_crosswalk_confirmed.csv")

sw <- expand.grid(min_agree = c(5L, 15L, 30L, 60L), max_conflict_frac = c(0.05, 0.10, 0.15, 0.25),
                  min_margin = c(5, 10, 20))
sw$confirmed <- vapply(seq_len(nrow(sw)), function(r)
  sum(status_fun(bA, sw$min_agree[r], sw$max_conflict_frac[r], sw$min_margin[r]) == "confirmed_unique" &
      !(A$life_id[bA$i] %in% conflicting_duplicate_ids)), 1L)
sw$pct_of_rows <- round(100 * sw$confirmed / nA, 1)
write_csv(sw, "datajoin_after_exclusion_threshold_sweep.csv")

# datajoin side
b_conf <- unique(conf$j); b_cand <- unique(bB$j)
side <- data.frame(SIC = dj$SIC, datajoin_excel_row = seq_len(nB) + 1L,
  status = ifelse(seq_len(nB) %in% b_conf, "matched_confirmed",
           ifelse(seq_len(nB) %in% b_cand, "has_candidate_not_confirmed", "no_candidate_in_data_after_exclusion")),
  stringsAsFactors = FALSE)
write_csv(side, "datajoin_after_exclusion_SIC_side.csv")

# comparison with the earlier IDmatch-based diagnostic (08)
prev_lines <- "Earlier diagnostic file not found (skipped)."
if (file.exists(prev_path)) {
  pv <- read.csv(prev_path, colClasses = "character", check.names = FALSE, na.strings = "")
  if (all(c("life_id", "datajoin_sic", "candidate_status") %in% names(pv))) {
    pv <- pv[pv$candidate_status == "strong_unique_candidate_review", ]
    pv$life_id <- toupper(pv$life_id); pv$datajoin_sic <- toupper(pv$datajoin_sic)
    pv <- pv[!duplicated(pv$life_id), ]
    in_A <- pv$life_id %in% A$life_id
    mm <- match(pv$life_id, res$life_id)
    new_sic <- ifelse(in_A, res$SIC[mm], NA); new_st <- ifelse(in_A, res$status[mm], NA)
    prev_lines <- c(
      paste("Earlier strong candidates (IDmatch-based):", nrow(pv)),
      paste("  removed by exclusions (not in data_after_exclusion):", sum(!in_A)),
      paste("  still present:", sum(in_A)),
      paste("    confirmed now with the SAME SIC:", sum(in_A & new_st == "confirmed_unique" & new_sic == pv$datajoin_sic, na.rm = TRUE)),
      paste("    confirmed now with a DIFFERENT SIC (check!):", sum(in_A & new_st == "confirmed_unique" & new_sic != pv$datajoin_sic, na.rm = TRUE)),
      paste("    not confirmed now:", sum(in_A & new_st != "confirmed_unique", na.rm = TRUE)),
      paste("  confirmed now but not 'strong' earlier:",
            sum(res$status == "confirmed_unique" & !res$life_id %in% pv$life_id)))
    chk <- data.frame(life_id = pv$life_id, earlier_SIC = pv$datajoin_sic, now_status = new_st, now_SIC = new_sic)
    write_csv(chk[in_A & (is.na(chk$now_SIC) | chk$now_SIC != chk$earlier_SIC), ],
              "datajoin_after_exclusion_vs_earlier_differences.csv")
  }
}

# ---- 8. Text report ------------------------------------------------------------
st_tab <- table(factor(res$status, levels = c("confirmed_unique", "ambiguous_close_candidates",
  "datajoin_row_fits_other_participant_better", "weak_or_conflicting", "duplicate_input_conflict", "no_candidate")))
report <- capture.output({
  cat("========== LIFE: data_after_exclusion  vs  PV0696_datajoin (diagnostic) ==========\n")
  cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
  cat("---------- 1. INPUT ----------\n")
  cat("data_after_exclusion:", a_path, "\n  original CSV rows:", nrow(data_full),
      "| matching rows (unique participants):", nA, "| unique life_id:", length(unique(A$life_id)),
      "| life_id with >1 original row:", length(repeated),
      "| conflicting duplicate participant IDs (not confirmed):", length(conflicting_duplicate_ids),
      "| columns:", ncol(A), "\n")
  cat("datajoin:", dj_path, "\n  rows:", nB, "| unique SIC:", length(unique(dj$SIC)), "| columns:", ncol(dj), "\n")
  cat("Mapping file:", map_path, "(selected keep=1/forced rows:", nrow(sel), ")\n\n")
  cat("---------- 2. COLUMN CORRESPONDENCE ----------\n")
  cat("Mapped pairs:", nrow(cmap), "| by manual mapping:", sum(cmap$source == "manual_mapping"),
      "| alias table:", sum(cmap$source == "alias_table"), "| same name:", sum(cmap$source == "same_name"), "\n")
  cat("Used for matching:", K, "| skipped:", sum(!cmap$use), "\n")
  print(table(skip_reason = cmap$skip_reason[!cmap$use]))
  cat("Date columns used:", length(dkeys), "\n")
  cat("data_after_exclusion columns without counterpart:", nrow(unm_tab), "\n")
  print(table(unm_tab$reason))
  cat("datajoin columns NOT in data_after_exclusion (= candidates to add later):", length(dj_new), "\n\n")
  cat("---------- 3. DATE SHIFT OF DATAJOIN (chosen from data) ----------\n")
  cat("Global shift:", global_shift, "days (sum of overlaps per shift:",
      paste(SHIFTS, colSums(ov), sep = " -> ", collapse = " | "), ")\n")
  print(shift_tab, row.names = FALSE)
  if (length(dkeys) && max(colSums(ov)) == 0) cat("!! No date overlap at any shift. Check date parsing in both files.\n")
  cat("\n---------- 4. CANDIDATE PAIRS ----------\n")
  cat("Fields skipped (too many pairs / no usable values):", sum(nzchar(anchor_tab$skipped)), "\n")
  cat("Participants with any proposed pair (>=1 column):", pct(rows_any_pair, nA), "\n")
  cat("Pairs proposed by >=", MIN_ANCHOR_HITS, "columns:", length(PI),
      "| participants covered:", pct(rows_with_pair, nA), "\n\n")
  cat("---------- 5. COLUMN RELIABILITY ----------\n")
  cat("Pass 1 confirmed:", nrow(conf1), "\n")
  cat("Columns dropped in pass 2 (>", UNRELIABLE_CONFLICT * 100, "% disagreement among confirmed pairs):",
      sum(unrel), "\n", sep = "")
  if (any(unrel)) print(ca1[unrel, c("column", "datajoin_column", "n_both", "conflict_rate")], row.names = FALSE)
  cat("Columns with the most disagreement among final confirmed pairs (top 15; n_both >= 30):\n")
  top <- ca[ca$n_both >= 30, ]; top <- top[order(-top$conflict_rate), ]
  print(head(top[, c("column", "datajoin_column", "n_both", "n_conflict", "conflict_rate")], 15), row.names = FALSE)
  cat("\n---------- 6. RESULT (default: min_agree", MIN_AGREE, "| max_conflict", MAX_CONF,
      "| min_margin", MIN_MARGIN, ") ----------\n")
  print(st_tab)
  cat("\n>>> ROWS OF data_after_exclusion THAT CAN BE MATCHED TO THE DATAJOIN:", pct(nrow(cw), nA), "\n")
  cat("    (one row per life_id; crosswalk: datajoin_after_exclusion_crosswalk_confirmed.csv)\n")
  cat("Confirmed pairs: median agreeing columns:",
      if (nrow(cw)) median(cw$n_agree) else NA_real_,
      "| median conflicts:", if (nrow(cw)) median(cw$n_conflict) else NA_real_,
      "| max conflict fraction:", if (nrow(cw)) max(cw$conflict_frac) else NA_real_, "\n")
  cat("Cells that the datajoin would ADD in shared columns (A missing, datajoin present), sum over columns:",
      sum(ca$A_missing_datajoin_present), "\n")
  cat("Cells where both are filled but differ, sum over columns:", sum(ca$n_conflict), "\n\n")
  cat("Threshold sweep (confirmed participants):\n")
  print(sw[sw$min_margin == MIN_MARGIN, ], row.names = FALSE)
  cat("(full grid incl. other margins in datajoin_after_exclusion_threshold_sweep.csv)\n\n")
  cat("---------- 7. DATAJOIN SIDE ----------\n")
  print(table(side$status))
  cat("Datajoin rows without a data_after_exclusion participant are mostly new people or excluded participants.\n\n")
  cat("---------- 8. COMPARISON WITH EARLIER IDmatch-BASED DIAGNOSTIC ----------\n")
  cat(paste(prev_lines, collapse = "\n"), "\n\n")
  cat("---------- 9. FILES ----------\n")
  cat(paste0("  ", list.files(out), collapse = "\n"), "\n")
  cat("\nMatching diagnostic complete; control variables will be joined below.\n")
})
writeLines(report, file.path(out, "datajoin_after_exclusion_summary.txt"), useBytes = TRUE)
cat(paste(report, collapse = "\n"), "\n")
cat("\nSaved outputs to:", out, "\n")


# ---- 9. Join five independent controls on confirmed SICs only ----------------
# This script regenerates the confirmed crosswalk above, so a stale crosswalk
# from a previous exclusion run cannot be used accidentally.
source_dir <- file.path(life, "tabular", "Behavior", "PV0696_FILES", "data")
echo_path   <- file.path(source_dir, "PV0696_T01092_NODUP.xlsx")
anthro_path <- file.path(source_dir, "PV0696_T00047_NODUP.xlsx")
ada_path    <- file.path(source_dir, "PV0696_T01404_NODUP.xlsx")
for (path in c(echo_path, anthro_path, ada_path)) {
  if (!file.exists(path)) stop("Missing control source workbook: ", path,
                               "\nNo post-exclusion dataset has been modified.")
}

# Explicit header resolution, not a fixed Excel column number. NODUP workbooks
# may expose either the long descriptive name or the instrument field code.
find_one <- function(df, candidates, description) {
  h <- toupper(trimws(names(df)))
  for (candidate in toupper(candidates)) {
    which_col <- which(h == candidate)
    if (length(which_col) == 1L) return(names(df)[which_col])
    if (length(which_col) > 1L) stop("Ambiguous header: ", description)
  }
  stop("Column not found for ", description, ". Expected one of: ",
       paste(candidates, collapse = ", "), "\nWorkbook headers include: ",
       paste(head(names(df), 20L), collapse = ", "))
}
read_source <- function(path, key_candidates, source_name) {
  x <- as.data.frame(readxl::read_excel(
    path, sheet = 1L, col_types = "text", .name_repair = "minimal"),
    check.names = FALSE, stringsAsFactors = FALSE)
  names(x) <- trimws(names(x))
  if (anyDuplicated(names(x))) stop("Duplicate headers in ", path)
  kc <- find_one(x, key_candidates, paste(source_name, "SIC"))
  keys <- toupper(norm(x[[kc]]))
  dup <- unique(keys[!is.na(keys) & duplicated(keys)])
  if (length(dup))
    stop("Source claims NODUP but has duplicate participant SICs in ", path,
         " (first values: ", paste(head(dup, 10L), collapse = ", "),
         "). Aborting instead of choosing an arbitrary record.")
  list(data = x, SIC = keys, key_column = kc, path = path, name = source_name)
}

E <- read_source(echo_path, c("SIC", "ECHO_READ_SIC", "T01092_F0002"), "echo")
H <- read_source(anthro_path, c("SIC", "ANTHRO_SIC", "T00047_F0004"), "anthropometry")
D <- read_source(ada_path, c("SIC", "ADA_SIC", "T01404_F0003"), "ADA diabetes")

col_ed <- find_one(E$data, c("ECHO_READ_LVEDV_MOD_BP", "T01092_F0375"), "echo LVEDV MOD_BP")
col_es <- find_one(E$data, c("ECHO_READ_LVESV_MOD_BP", "T01092_F0381"), "echo LVESV MOD_BP")
col_h  <- find_one(H$data, c("ANTHRO_F0010", "T00047_F0010"), "height")
col_d  <- find_one(D$data, c("ADA_DS", "T01404_F0013"), "diabetes status")

# Parsed values have a fixed missingness definition. We never replace missing
# values with 0, including no matched SIC or only one of the two echo volumes.
as_num <- function(x) {
  y <- suppressWarnings(as.numeric(gsub(",", ".", norm(x), fixed = TRUE)))
  y[!is.finite(y)] <- NA_real_
  y
}
raw_ed <- E$data[[col_ed]]
raw_es <- E$data[[col_es]]
raw_h  <- H$data[[col_h]]
raw_d  <- D$data[[col_d]]
v_ed <- as_num(raw_ed)
v_es <- as_num(raw_es)
v_h  <- as_num(raw_h)
v_d  <- norm(raw_d)
v_h[v_h %in% c(996, 998, 999)] <- NA_real_
# The supplied T01404 codebook gives only 0 and 1 for ADA_DS.
# Unrecognized diabetes codes must not be silently reinterpreted as 'no'.
invalid_diabetes <- !is.na(v_d) & !(v_d %in% c("0", "1"))
if (any(invalid_diabetes)) {
  stop("Unexpected values in ADA_DS: ",
       paste(unique(v_d[invalid_diabetes]), collapse = ", "),
       ". Verify the ADA workbook before joining. No dataset was modified.")
}
v_d[!(v_d %in% c("0", "1"))] <- NA_character_
v_sv <- ifelse(!is.na(v_ed) & !is.na(v_es), v_ed - v_es, NA_real_)

# m_cw maps every original post-exclusion row (even repeated life_id) to one
# confirmed SIC. Non-confirmed participants always receive missing controls.
m_cw <- match(full_ids, cw$life_id)
confirmed_sic <- cw$SIC[m_cw]
ix_e <- match(confirmed_sic, E$SIC)
ix_h <- match(confirmed_sic, H$SIC)
ix_d <- match(confirmed_sic, D$SIC)

# For diagnostics, analyze unique participants; don't let unresolved duplicate
# rows artificially inflate coverage numbers.
match_uid <- match(uid, cw$life_id)
uid_sic <- cw$SIC[match_uid]
ui_e <- match(uid_sic, E$SIC)
ui_h <- match(uid_sic, H$SIC)
ui_d <- match(uid_sic, D$SIC)

fmt_num <- function(x) {
  out <- as.character(x)
  out[is.na(x)] <- NA_character_
  out
}
added <- data.frame(
  echo_lvedv_mod_bp = fmt_num(v_ed[ix_e]),
  echo_lvesv_mod_bp = fmt_num(v_es[ix_e]),
  lvedv_lvesv = fmt_num(v_sv[ix_e]),
  height = fmt_num(v_h[ix_h]),
  diabetes_final = v_d[ix_d],
  stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(names(added), NEW_COLS)) stop("Internal error: new column order")
if (nrow(added) != nrow(data_full)) stop("Internal error: join changed row count")

# Valid-coverage table: raw nonempty cells, valid source measurements, and
# availability in the confirmed post-exclusion sample. The last count uses
# unique IDs, not rows, and includes *no* unconfirmed matches.
row_coverage <- function(output_col, instrument, path, source_col, raw, valid,
                         source_match, final_participant_values, comment = "") {
  data.frame(
    output_column = output_col,
    source = instrument,
    source_file = basename(path),
    source_column = source_col,
    source_rows_total = length(raw),
    source_raw_nonempty = sum(!is.na(norm(raw))),
    source_valid_nonmissing = sum(!is.na(valid)),
    confirmed_participants_with_source_row = sum(!is.na(source_match)),
    confirmed_participants_with_valid_value = sum(!is.na(final_participant_values)),
    postexclusion_participants_total = length(uid),
    postexclusion_confirmed_participants = nrow(cw),
    note = comment, stringsAsFactors = FALSE)
}

coverage <- rbind(
  row_coverage("echo_lvedv_mod_bp", "T01092", echo_path, col_ed, raw_ed,
               v_ed, ui_e, v_ed[ui_e]),
  row_coverage("echo_lvesv_mod_bp", "T01092", echo_path, col_es, raw_es,
               v_es, ui_e, v_es[ui_e]),
  row_coverage("lvedv_lvesv", "T01092", echo_path,
               paste0(col_ed, " - ", col_es),
               ifelse(!is.na(norm(raw_ed)) & !is.na(norm(raw_es)), "both", NA_character_),
               v_sv, ui_e, v_sv[ui_e],
               "Nonmissing iff both source echo volumes are valid and present"),
  row_coverage("height", "T00047", anthro_path, col_h, raw_h,
               v_h, ui_h, v_h[ui_h], "996/998/999 are missing by instrument documentation"),
  row_coverage("diabetes_final", "T01404", ada_path, col_d, raw_d,
               v_d, ui_d, v_d[ui_d], "ADA_DS: 0=no diabetes; 1=diabetes; NA=unknown")
)
write_csv(coverage, "datajoin_controls_coverage.csv")

# Additional audit checks, without imposing undocumented clinical exclusions.
extra_qc <- data.frame(
  check = c("echo_valid_both", "echo_nonpositive_diff", "height_996", "height_998",
            "height_999", "height_non_numeric_raw", "height_outside_50_to_250_cm",
            "diabetes_yes_initial", "diabetes_no_initial", "diabetes_yes_matched",
            "diabetes_no_matched", "unmatched_after_crosswalk", "repeated_input_ids",
            "repeated_input_ids_with_conflicts"),
  count = c(sum(!is.na(v_sv)), sum(!is.na(v_sv) & v_sv <= 0),
            sum(as_num(raw_h) == 996, na.rm = TRUE),
            sum(as_num(raw_h) == 998, na.rm = TRUE),
            sum(as_num(raw_h) == 999, na.rm = TRUE),
            sum(!is.na(norm(raw_h)) & is.na(as_num(raw_h))),
            sum(!is.na(v_h) & (v_h < 50 | v_h > 250)),
            sum(v_d == "1", na.rm = TRUE), sum(v_d == "0", na.rm = TRUE),
            sum(v_d[ui_d] == "1", na.rm = TRUE), sum(v_d[ui_d] == "0", na.rm = TRUE),
            length(uid) - nrow(cw), length(repeated), length(conflicting_duplicate_ids)),
  stringsAsFactors = FALSE)
write_csv(extra_qc, "datajoin_controls_additional_qc.csv")

# Append only. Existing field names and VALUES are kept as is.
result <- data_full
for (field in NEW_COLS) result[[field]] <- added[[field]]
if (nrow(result) != nrow(data_full) ||
    !identical(names(result)[seq_along(names(data_full))], names(data_full)) ||
    !identical(result[, names(data_full), drop = FALSE], data_full))
  stop("Internal error: an existing value, header, or row changed")

# ---- 10. Human-readable report and safe export -------------------------------
report_controls <- capture.output({
  cat("========== LIFE: high-confidence matching + post-exclusion controls =========\n")
  cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n")
  cat("Matching model: original script 09 (MIN_AGREE=", MIN_AGREE,
      ", MAX_CONF=", MAX_CONF, ", MIN_MARGIN=", MIN_MARGIN,
      "); reciprocal best and 1:1 only\n", sep = "")
  cat("Input post-exclusion rows:", nrow(data_full), "\n")
  cat("Unique post-exclusion life_id:", length(uid), "\n")
  cat("Confirmed unique life_id <-> datajoin SIC:", nrow(cw), "\n")
  cat("Unconfirmed (never assigned external values):", length(uid) - nrow(cw), "\n")
  cat("Multiple input rows for one life_id:", length(repeated), "\n")
  cat("IDs with conflicting input rows, withheld from matching:",
      length(conflicting_duplicate_ids), "\n\n")
  cat("Source columns and missingness (COUNTS ARE UNIQUE PARTICIPANTS ON JOIN):\n")
  print(coverage[, c("output_column", "source_column", "source_rows_total",
       "source_raw_nonempty", "source_valid_nonmissing",
       "confirmed_participants_with_source_row",
       "confirmed_participants_with_valid_value")], row.names = FALSE)
  cat("\nNODUP source key columns: echo=", E$key_column,
      "; anthro=", H$key_column, "; diabetes=", D$key_column, "\n", sep = "")
  cat("SICs found in T01092 (of confirmed):", sum(!is.na(ui_e)), "\n")
  cat("SICs found in T00047 (of confirmed):", sum(!is.na(ui_h)), "\n")
  cat("SICs found in T01404 (of confirmed):", sum(!is.na(ui_d)), "\n")
  cat("\nDiabetes coding: ADA_DS 0=absent, 1=present, missing=unknown.\n")
  cat("Height coding: 996/998/999=missing (source instrument codes).\n")
  cat("Stroke volume: LVEDV_MOD_BP - LVESV_MOD_BP (mL), both required.\n")
  cat("\nAdditional checks:\n")
  print(extra_qc, row.names = FALSE)
  cat("\nExisting input data are unchanged. Exactly 5 columns appended.\n")
  cat("Matching crosswalk:", file.path(out, "datajoin_after_exclusion_crosswalk_confirmed.csv"), "\n")
  cat("Coverage CSV:", file.path(out, "datajoin_controls_coverage.csv"), "\n")
  cat("Updated pipeline input:", a_path, "\n")
})
# Write the summary after the enriched CSV is saved successfully.

# Backup the original post-exclusion file before any update; then write the new
# CSV to a temporary path in the same filesystem and replace the final path.
backup <- file.path(deriv, "data_after_exclusion_before_controls.csv")
if (!file.copy(a_path, backup, overwrite = TRUE))
  stop("Unable to create pre-controls snapshot: ", backup)
tmp <- tempfile(pattern = "data_after_exclusion_controls_", tmpdir = deriv, fileext = ".csv")
write.csv(result, tmp, row.names = FALSE, na = "", fileEncoding = "UTF-8")
check <- read.csv(tmp, colClasses = "character", check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", na.strings = "")
if (nrow(check) != nrow(result) || !identical(names(check), names(result)) ||
    !identical(check[, names(data_full), drop = FALSE], data_full))
  stop("Output verification failed; original dataset preserved at: ", backup)
if (!file.rename(tmp, a_path)) stop("Could not atomically update ", a_path,
                                     "; original is backed up at: ", backup)
writeLines(report_controls, file.path(out, "datajoin_controls_join_summary.txt"))
cat(paste(report_controls, collapse = "\n"), "\n")
cat("\nSaved pre-join backup:", backup, "\n")
cat("Run 07_LIFE_pipeline_matching_clean.R next (it reads the enriched data_after_exclusion.csv).\n")
