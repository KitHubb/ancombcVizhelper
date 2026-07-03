.as_ancombc2_df <- function(x, object_name) {
  if (base::is.null(x)) {
    base::stop("`out$", object_name, "` is missing. Check the ANCOM-BC2 execution options.")
  }

  x <- base::as.data.frame(x, check.names = FALSE)

  if (!"taxon" %in% base::names(x)) {
    if (base::is.null(base::rownames(x))) {
      base::stop("`out$",
                 object_name,
                 "` does not contain a taxon column or row names.")
    }

    x$taxon <- base::rownames(x)
  }

  tibble::as_tibble(x)
}


.select_lfc_cols <- function(x,
                             prefix = NULL,
                             require_diff = FALSE) {
  lfc_cols <- base::names(x)[base::startsWith(base::names(x), "lfc_")]

  if (!base::is.null(prefix) && base::nzchar(prefix)) {
    suffix <- base::sub("^lfc_", "", lfc_cols)
    lfc_cols <- lfc_cols[base::startsWith(suffix, prefix)]
  }

  if (require_diff) {
    suffix <- base::sub("^lfc_", "", lfc_cols)

    lfc_cols <- lfc_cols[base::paste0("diff_", suffix) %in% base::names(x)]
  }

  if (base::length(lfc_cols) == 0) {
    base::stop("No usable lfc_ columns were found for `prefix = ", prefix, "`.")
  }

  lfc_cols
}


.extract_structural_zeros <- function(out) {
  zero_ind <- out$zero_ind

  empty_result <- function() {
    tibble::tibble(
      taxon = base::character(),
      structural_zero = base::logical(),
      structural_zero_groups = base::character()
    )
  }

  if (base::is.null(zero_ind) || base::length(zero_ind) == 0) {
    return(empty_result())
  }

  zero_ind <- base::as.data.frame(zero_ind, check.names = FALSE)

  taxa <- if ("taxon" %in% base::names(zero_ind)) {
    base::as.character(zero_ind$taxon)
  } else {
    base::rownames(zero_ind)
  }

  zero_cols <- base::grep("^structural_zero", base::names(zero_ind), value = TRUE)

  if (base::is.null(taxa) ||
      base::length(zero_cols) == 0) {
    return(empty_result())
  }

  zero_mat <- base::sapply(zero_cols, function(col) {
    base::as.logical(zero_ind[[col]])
  })

  if (base::is.null(base::dim(zero_mat))) {
    zero_mat <- base::matrix(zero_mat,
                             ncol = 1,
                             dimnames = base::list(NULL, zero_cols))
  }

  has_structural_zero <- base::rowSums(zero_mat, na.rm = TRUE) > 0

  if (!base::any(has_structural_zero)) {
    return(empty_result())
  }

  structural_groups <- base::apply(zero_mat[has_structural_zero, , drop = FALSE], 1, function(x) {
    base::paste(base::colnames(zero_mat)[base::which(x)], collapse = "; ")
  })

  tibble::tibble(
    taxon = taxa[has_structural_zero],
    structural_zero = TRUE,
    structural_zero_groups = structural_groups
  )
}


.prepare_ancombc2_lfc_data <- function(out,
                                       result = base::c("res", "res_pair", "res_dunn", "res_global", "res_trend"),
                                       prefix = NULL,
                                       sensitivity = base::c("keep", "robust_only"),
                                       show_all = TRUE) {
  result <- base::match.arg(result)
  sensitivity <- base::match.arg(sensitivity)

  if (!base::is.list(out) || base::is.null(out$res)) {
    base::stop(
      "`out` must be the full result object returned by `ANCOMBC::ancombc2()`. ",
      "Individual result tables such as `out$res` are not accepted."
    )
  }

  if (result %in% base::c("res_global", "res_trend") &&
      (base::is.null(prefix) || !base::nzchar(prefix))) {
    base::stop("`result = '",
               result,
               "'` requires `prefix` to identify the group-related LFCs to display.")
  }

  result_tbl <- .as_ancombc2_df(out[[result]], result)

  structural_zero <- .extract_structural_zeros(out)

  if (result %in% base::c("res", "res_pair", "res_dunn")) {
    lfc_cols <- .select_lfc_cols(x = result_tbl,
                                 prefix = prefix,
                                 require_diff = TRUE)

    suffix <- base::sub("^lfc_", "", lfc_cols)
    diff_cols <- base::paste0("diff_", suffix)
    ss_cols <- base::paste0("passed_ss_", suffix)

    missing_cols <- base::setdiff(base::c(diff_cols, ss_cols), base::names(result_tbl))

    if (base::length(missing_cols) > 0) {
      base::stop(
        "The following ANCOM-BC2 columns are missing: ",
        base::paste(missing_cols, collapse = ", "),
        ". Check whether `pseudo_sens = TRUE` was used."
      )
    }

    df_lfc <- result_tbl |>
      dplyr::select(dplyr::all_of(base::c("taxon", lfc_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(lfc_cols),
        names_to = "comparison",
        values_to = "lfc"
      ) |>
      dplyr::mutate(comparison = base::sub("^lfc_", "", .data$comparison))

    df_diff <- result_tbl |>
      dplyr::select(dplyr::all_of(base::c("taxon", diff_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(diff_cols),
        names_to = "comparison",
        values_to = "significant"
      ) |>
      dplyr::mutate(
        comparison = base::sub("^diff_", "", .data$comparison),
        significant = base::as.logical(.data$significant)
      )

    df_ss <- result_tbl |>
      dplyr::select(dplyr::all_of(base::c("taxon", ss_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(ss_cols),
        names_to = "comparison",
        values_to = "passed_ss"
      ) |>
      dplyr::mutate(
        comparison = base::sub("^passed_ss_", "", .data$comparison),
        passed_ss = base::as.logical(.data$passed_ss)
      )

    df_final <- df_lfc |>
      dplyr::left_join(df_diff, by = base::c("taxon", "comparison")) |>
      dplyr::left_join(df_ss, by = base::c("taxon", "comparison"))

    comparison_order <- suffix

    inference_note <- base::switch(result,
                                   res =
                                     "Primary ANCOM-BC2 result: FDR is controlled within each coefficient, but mdFDR is not controlled across multiple group comparisons.",
                                   res_pair =
                                     "Pairwise directional ANCOM-BC2 result: diff_* reflects mdFDR-controlled pairwise inference.",
                                   res_dunn =
                                     "Dunnett-type ANCOM-BC2 result: diff_* reflects mdFDR-controlled comparisons against the reference group.")

    if (result == "res" && base::length(lfc_cols) > 1) {
      base::warning(
        "`result = 'res'` does not control mdFDR across multiple group comparisons. ",
        "Use `res_pair` or `res_dunn` for directional comparisons across multiple groups.",
        call. = FALSE
      )
    }
  }

  if (result == "res_global") {
    lfc_source <- .as_ancombc2_df(out$res, "res")

    lfc_cols <- .select_lfc_cols(x = lfc_source,
                                 prefix = prefix,
                                 require_diff = TRUE)

    required_cols <- base::c("diff_abn", "passed_ss")

    missing_cols <- base::setdiff(required_cols, base::names(result_tbl))

    if (base::length(missing_cols) > 0) {
      base::stop("The following required columns are missing from `res_global`: ",
                 base::paste(missing_cols, collapse = ", "))
    }

    df_lfc <- lfc_source |>
      dplyr::select(dplyr::all_of(base::c("taxon", lfc_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(lfc_cols),
        names_to = "comparison",
        values_to = "lfc"
      ) |>
      dplyr::mutate(comparison = base::sub("^lfc_", "", .data$comparison))

    df_status <- result_tbl |>
      dplyr::transmute(
        taxon = .data$taxon,
        significant = base::as.logical(.data$diff_abn),
        passed_ss = base::as.logical(.data$passed_ss)
      )

    df_final <- df_lfc |>
      dplyr::left_join(df_status, by = "taxon")

    comparison_order <- base::sub("^lfc_", "", lfc_cols)

    inference_note <-
      "Global ANCOM-BC2 result: taxa are selected using the taxon-level global test; displayed LFCs are group coefficients from out$res."
  }

  if (result == "res_trend") {
    lfc_cols <- .select_lfc_cols(x = result_tbl,
                                 prefix = prefix,
                                 require_diff = FALSE)

    required_cols <- base::c("diff_abn", "passed_ss")

    missing_cols <- base::setdiff(required_cols, base::names(result_tbl))

    if (base::length(missing_cols) > 0) {
      base::stop("The following required columns are missing from `res_trend`: ",
                 base::paste(missing_cols, collapse = ", "))
    }

    df_lfc <- result_tbl |>
      dplyr::select(dplyr::all_of(base::c("taxon", lfc_cols))) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(lfc_cols),
        names_to = "comparison",
        values_to = "lfc"
      ) |>
      dplyr::mutate(comparison = base::sub("^lfc_", "", .data$comparison))

    df_status <- result_tbl |>
      dplyr::transmute(
        taxon = .data$taxon,
        significant = base::as.logical(.data$diff_abn),
        passed_ss = base::as.logical(.data$passed_ss)
      )

    df_final <- df_lfc |>
      dplyr::left_join(df_status, by = "taxon")

    comparison_order <- base::sub("^lfc_", "", lfc_cols)

    inference_note <-
      "Trend ANCOM-BC2 result: taxa are selected using the taxon-level trend test; displayed LFCs describe the fitted group coefficients."
  }

  if (base::nrow(structural_zero) > 0) {
    df_final <- df_final |>
      dplyr::anti_join(structural_zero |>
                         dplyr::select(dplyr::all_of("taxon")), by = "taxon")
  }

  df_final <- df_final[
    base::is.finite(df_final$lfc),
    ,
    drop = FALSE
  ]

  df_final$significant <- base::as.logical(df_final$significant)
  df_final$passed_ss <- base::as.logical(df_final$passed_ss)

  df_final$significant[
    base::is.na(df_final$significant)
  ] <- FALSE

  df_final$passed_ss[
    base::is.na(df_final$passed_ss)
  ] <- FALSE

  df_final$robust <- df_final$significant &
    df_final$passed_ss

  df_final$display <- if (
    base::identical(sensitivity, "robust_only")
  ) {
    df_final$robust
  } else {
    df_final$significant
  }

  df_final$plot_value <- base::ifelse(
    df_final$display,
    df_final$lfc,
    0
  )

  df_final$label <- base::ifelse(
    df_final$display,
    base::sprintf("%.2f", df_final$lfc),
    ""
  )

  if (!show_all) {

    taxa_to_keep <- base::unique(
      base::as.character(
        df_final$taxon[df_final$display %in% TRUE]
      )
    )

    df_final <- df_final[
      base::as.character(df_final$taxon) %in% taxa_to_keep,
      ,
      drop = FALSE
    ]
  }

  if (base::nrow(df_final) == 0) {
    reason <- if (sensitivity == "robust_only") {
      "No taxa were significant and passed the pseudo-count sensitivity analysis."
    } else {
      "No taxa were significant in the selected result."
    }

    base::stop("No taxa are available for display. ",
               reason,
               " Set `show_all = TRUE` to retain non-significant taxa in the output data.")
  }

  taxon_order <- base::unique(base::as.character(df_final$taxon))

  df_final <- df_final |>
    dplyr::mutate(
      comparison = base::factor(.data$comparison, levels = comparison_order),
      taxon = base::factor(.data$taxon, levels = base::rev(taxon_order))
    )

  base::list(
    data = df_final,
    structural_zero = structural_zero,
    result_type = result,
    inference_note = inference_note
  )
}
