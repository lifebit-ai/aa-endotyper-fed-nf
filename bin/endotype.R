
#!/usr/bin/env Rscript

library(tidyverse)
library(janitor)
library(tidyHeatmap)
library(ComplexHeatmap)

# paths <- list.files("mounted-data-readonly", full.names = TRUE)

# linked_file <- read_csv(paths[str_detect(paths, "demo_1737")])
# meas_file <- read_csv(paths[str_detect(paths, "aa_1737")])

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if we have the required arguments
if (length(args) < 2) {
  stop("Usage: Rscript endotype.R <linked_file.csv> <meas_file.csv>")
}

# Read the input files
linked_file <- read_csv(args[1])
meas_file <- read_csv(args[2])

comb_file <- linked_file %>% 
  inner_join(meas_file, by = "accession")

## accessory functions

aa_age_to_group <- function(age) {
  grp <- dplyr::case_when(
    is.na(age)    ~ "Missing",
    age < 8       ~ "8-12",
    age <= 12     ~ "8-12",
    age <= 17     ~ "13-17",
    TRUE          ~ "18-30"
  )
  factor(grp, levels = c("8-12","13-17","18-30","Missing"))
}

aa_first_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit)) hit[[1]] else NULL
}

# ───────────────────────────────────────────────────────────────────────────────
# 2) Canonicalize + (optionally) restrict visits + longify to Property/Value
#    Works for both formats:
#      - long:   analyte_name_col + analyte_value_col present (e.g., SDY524)
#      - wide:   analytes are columns (e.g., SDY1737)
# ───────────────────────────────────────────────────────────────────────────────
aa_prepare <- function(
    df,
    accession_col,
    sex_col,
    age_group_col = NULL,         # if not given, compute from age_years_col
    age_years_col = NULL,
    analyte_name_col = NULL,      # (long) e.g., "AutoAB"
    analyte_value_col = NULL,     # (long) e.g., "Value"
    analyte_cols = NULL          # (wide) e.g., c("GAD65","IA_2ic","MIAA","Zn_T8")
    #visit_cols = c("Visit_Num","Visitnum","Visit"),
    #restrict_visits = NULL        # numeric or -1 (for textual baseline)
) {
  stopifnot(accession_col %in% names(df), sex_col %in% names(df))
  
  out <- df %>%
    dplyr::rename(
      Accession = !!rlang::sym(accession_col),
      Sex       = !!rlang::sym(sex_col)
    )
  
  # ---- Age grouping (canonical labels)
  AGE_LEVELS <- c("0-7","8-12","13-17","18-30",">30","Missing")
  
  # replace aa_age_to_group() inside aa_prepare() with:
  aa_age_to_group <- function(x) {
    xx <- trimws(as.character(x))
    xx[xx == ""] <- NA_character_
    # map pre-labeled buckets safely (no "" key)
    map <- c(
      "0-7" = "0-7", "<8" = "0-7", "< 8" = "0-7",
      "8-12" = "8-12",
      "13-17" = "13-17",
      "18-30" = "18-30",
      ">30" = ">30", "> 30" = ">30", "30+" = ">30", "31+" = ">30",
      "Missing" = "Missing", "missing" = "Missing"
    )
    lbl <- unname(map[xx])  # unmatched → NA
    # numeric fallback for anything not matched above
    num  <- suppressWarnings(as.numeric(xx))
    fill <- is.na(lbl) & !is.na(num)
    lbl[fill & num < 8]               <- "0-7"
    lbl[fill & num >= 8  & num <= 12] <- "8-12"
    lbl[fill & num >= 13 & num <= 17] <- "13-17"
    lbl[fill & num >= 18 & num <= 30] <- "18-30"
    lbl[fill & num > 30]              <- ">30"
    lbl[is.na(lbl)] <- "Missing"
    
    factor(lbl, levels = c("0-7","8-12","13-17","18-30",">30","Missing"))
  }
  
  if (!is.null(age_group_col) && age_group_col %in% names(out)) {
    out <- out %>% dplyr::mutate(Age_Group = aa_age_to_group(.data[[age_group_col]]))
  } else if (!is.null(age_years_col) && age_years_col %in% names(out)) {
    out <- out %>% dplyr::mutate(Age_Group = aa_age_to_group(.data[[age_years_col]]))
  } else {
    out <- out %>% dplyr::mutate(Age_Group = factor("Missing", levels = AGE_LEVELS))
  }
  
  # ---- Visit filter (numeric or textual baseline)
  # vcol <- intersect(visit_cols, names(out))[1]
  # if (!is.na(vcol) && length(vcol) && !is.null(restrict_visits)) {
  #   vraw <- out[[vcol]]
  #   vnum <- suppressWarnings(as.numeric(as.character(vraw)))
  #   
  #   if (!all(is.na(vnum))) {
  #     # numeric/labelled visits
  #     keep <- !is.na(vnum) & vnum %in% restrict_visits
  #     out  <- out[keep, , drop = FALSE]
  #     out[[vcol]] <- vnum[keep]   # <-- align length with filtered rows
  #   } else {
  #     # textual visits: allow only baseline via -1
  #     if (all(restrict_visits == -1)) {
  #       vtxt <- tolower(trimws(as.character(vraw)))
  #       base_syn <- c("baseline","screening","day -1","day-1","v-1","-1")
  #       keep <- vtxt %in% base_syn
  #       out  <- out[keep, , drop = FALSE]
  #     } else {
  #       stop("Textual visit column cannot match non-baseline visits.")
  #     }
  #   }
  # }
  
  
  # ---- Long input
  if (!is.null(analyte_name_col) && !is.null(analyte_value_col) &&
      analyte_name_col %in% names(out) && analyte_value_col %in% names(out)) {
    
    tidy <- out %>%
      dplyr::transmute(
        Accession, Sex, Age_Group,
        Property = .data[[analyte_name_col]],
        Value    = suppressWarnings(as.numeric(.data[[analyte_value_col]]))
      )
    
  } else {
    # ---- Wide input → long
    id_like <- c("Accession","Sex","Age_Group", age_years_col)
    
    # (A) If analyte_cols explicitly provided, use them and force numeric
    if (!is.null(analyte_cols)) {
      present <- intersect(analyte_cols, names(out))
      if (!length(present)) stop("None of the specified analyte_cols are present.")
      out[present] <- lapply(out[present], function(x) suppressWarnings(as.numeric(x)))
      value_cols <- present
    } else {
      # (B) Coerce common antibody names if present
      common <- intersect(c("gad65","ia_2ic","miaa","zn_t8"), names(out))
      if (length(common)) {
        out[common] <- lapply(out[common], function(x) suppressWarnings(as.numeric(x)))
      }
      num_cols   <- names(out)[vapply(out, is.numeric, logical(1))]
      id_present <- intersect(id_like, names(out))
      value_cols <- setdiff(num_cols, id_present)
    }
    
    if (!length(value_cols))
      stop('No numeric analyte columns to melt in wide→long step. Supply analyte_cols=c("gad65","ia_2ic","miaa","zn_t8").')
    
    tidy <- out %>%
      dplyr::select(Accession, Sex, Age_Group, dplyr::all_of(value_cols)) %>%
      tidyr::pivot_longer(dplyr::all_of(value_cols), names_to = "Property", values_to = "Value")
  }
  
  if (!nrow(tidy)) stop("No rows left after visit filtering/preparation.")
  
  tidy %>%
    dplyr::group_by(Accession, Sex, Age_Group, Property) %>%
    dplyr::summarise(Value = suppressWarnings(mean(Value, na.rm = TRUE)), .groups = "drop")
}


# ───────────────────────────────────────────────────────────────────────────────
# 3) Build matrix (+ simple NA fill for clustering only)
# ───────────────────────────────────────────────────────────────────────────────
aa_build_matrix <- function(tidy_df) {
  wide <- tidy_df %>%
    pivot_wider(names_from = Property, values_from = Value)
  
  row_annot <- wide %>% select(Accession, Sex, Age_Group)
  mat <- wide %>%
    select(-Sex, -Age_Group) %>%
    tibble::column_to_rownames("Accession") %>%
    as.matrix()
  
  list(mat = mat, annot = row_annot)
}


aa_fill_for_clustering <- function(mat) {
  if (!nrow(mat) || !ncol(mat)) return(mat)
  # column medians
  cm <- apply(mat, 2, function(x) {
    m <- suppressWarnings(stats::median(x, na.rm = TRUE))
    if (is.infinite(m) || is.na(m)) 0 else m
  })
  # fill NAs
  for (j in seq_len(ncol(mat))) {
    na_idx <- is.na(mat[, j])
    if (any(na_idx)) mat[na_idx, j] <- cm[j]
  }
  mat
}

# ───────────────────────────────────────────────────────────────────────────────
# 4) Two small plotters: ward.D2 and k-means (kmeans done outside to avoid NA issues)
# ───────────────────────────────────────────────────────────────────────────────
aa_plot_ward <- function(tidy_df, scale_mode = "column",
                         dist_method = "minkowski", minkowski_p = 2) {
  dist_fun <- function(x) stats::dist(x, method = dist_method, p = minkowski_p)
  
  hm <- tidy_df %>%
    tidyHeatmap::heatmap(
      .row    = Accession,
      .column = Property,
      .value  = Value,
      scale   = scale_mode,
      clustering_distance_rows    = dist_fun,
      clustering_method_rows      = "ward.D2",
      clustering_distance_columns = dist_fun,
      clustering_method_columns   = "ward.D2"
    )
  
  if ("Sex" %in% names(tidy_df))       hm <- hm %>% annotation_tile(Sex)
  if ("Age_Group" %in% names(tidy_df)) hm <- hm %>% annotation_tile(Age_Group)
  
  hm
}


aa_plot_kmeans <- function(tidy_df, km_rows = 4, km_cols = 2, scale_mode = "column") {
  mm <- aa_build_matrix(tidy_df)
  if (!nrow(mm$mat)||!ncol(mm$mat)) stop("Empty matrix for k-means heatmap.")
  
  mat_fill <- aa_fill_for_clustering(mm$mat)
  
  # safe k
  kr <- max(1, min(km_rows, nrow(mat_fill)))
  kc <- max(1, min(km_cols, ncol(mat_fill)))
  
  # robust kmeans (fallback to 1 cluster if it errors)
  kmeans_safe <- function(x, centers) {
    out <- try(stats::kmeans(x, centers = centers, iter.max = 50), silent = TRUE)
    if (inherits(out, "try-error")) list(cluster = rep(1L, nrow(x))) else out
  }
  
  set.seed(1)
  rkm <- kmeans_safe(mat_fill,  kr)
  ckm <- kmeans_safe(t(mat_fill), kc)
  
  row_lab <- paste0("R", as.integer(rkm$cluster))
  col_lab <- paste0("C", as.integer(ckm$cluster))
  
  # named split vectors aligned to the matrix order
  row_split_vec <- factor(row_lab, levels = unique(row_lab))
  names(row_split_vec) <- rownames(mat_fill)
  row_split_vec <- row_split_vec[rownames(mm$mat)]
  col_split_vec <- factor(col_lab, levels = unique(col_lab))
  names(col_split_vec) <- colnames(mat_fill)
  col_split_vec <- col_split_vec[colnames(mm$mat)]
  
  hm <- tidy_df %>%
    tidyHeatmap::heatmap(
      .row         = Accession,
      .column      = Property,
      .value       = Value,
      scale        = scale_mode,
      row_split    = row_split_vec,
      column_split = col_split_vec,
      row_order    = rownames(mm$mat),
      column_order = colnames(mm$mat)
    )
  
  if ("Sex" %in% names(tidy_df))       hm <- hm %>% tidyHeatmap::annotation_tile(Sex)
  if ("Age_Group" %in% names(tidy_df)) hm <- hm %>% tidyHeatmap::annotation_tile(Age_Group)
  
  hm
}

# ───────────────────────────────────────────────────────────────────────────────
# 5) Tiny wrapper
# ───────────────────────────────────────────────────────────────────────────────
aa_heatmaps <- function(tidy_df,
                        scale_mode = "column",
                        dist_method = "minkowski", minkowski_p = 2,
                        k_rows = 4, k_cols = 2) {
  list(
    hm_ward   = aa_plot_ward(tidy_df, scale_mode, dist_method, minkowski_p),
    hm_kmeans = aa_plot_kmeans(tidy_df, k_rows, k_cols, scale_mode)
  )
}


# ───────────────────────────────────────────────────────────────────────────────
# 6) Running Analyses
# ───────────────────────────────────────────────────────────────────────────────

tidied_data <- aa_prepare(df = comb_file, 
                          accession_col = "accession", 
                          sex_col = "sex", 
                          age_years_col = "age", 
                          analyte_cols = c("gad65", "ia_2ic", "miaa", "zn_t8"))

heatmaps <- aa_heatmaps(tidy_df = tidied_data, 
                        scale_mode = "column")

save_pdf(heatmaps$hm_ward, filename = "sdy_ward.pdf")
save_pdf(heatmaps$hm_kmeans, filename = "sdy_kmeans.pdf")