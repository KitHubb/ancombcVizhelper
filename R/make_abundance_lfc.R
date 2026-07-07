#' Plot TSS abundance with ANCOM-BC2 log-fold changes
#'
#' Creates a composite figure containing group-wise average relative abundance
#' (TSS; mean +/- SE), ANCOM-BC2 log-fold changes with standard errors, and
#' the corresponding ANCOM-BC2 p and q values. The abundance panel is
#' descriptive only; no additional abundance-based hypothesis test is run.
#'
#' The function uses \code{.prepare_ancombc2_lfc_data()} to select taxa and
#' retrieve ANCOM-BC2 statistics. Consequently, \code{show_all} and
#' \code{sensitivity} follow the same selection rules as \code{make_barplots()}.
#' For pairwise/Dunn-type results, this is based on the helper's comparison-wise
#' differential-abundance result. Opacity of the LFC bars mirrors
#' \code{make_barplots()}: pseudo-count sensitivity-robust results are fully
#' opaque, whereas remaining results use \code{non_robust_alpha}.
#'
#' @param ps A \code{phyloseq} object containing the count table, taxonomy,
#'   and sample metadata.
#' @param out Full result object returned by \code{ANCOMBC::ancombc2()}.
#' @param result Character string specifying the ANCOM-BC2 result table to use.
#'   One of \code{"res"}, \code{"res_pair"}, \code{"res_dunn"},
#'   \code{"res_global"}, or \code{"res_trend"}.
#' @param prefix Optional model-coefficient prefix used to select LFC columns.
#'   For example, \code{"bmi"} selects coefficients beginning with
#'   \code{lfc_bmi}.
#' @param group Name of the grouping variable in \code{sample_data(ps)} used
#'   for the descriptive abundance panel.
#' @param tax_level Name of the taxonomy column in \code{tax_table(ps)} used
#'   for the abundance panel, such as \code{"Family"} or \code{"Genus"}.
#' @param comparison A single ANCOM-BC2 comparison/coefficient to display.
#'   When \code{groupnames = FALSE}, this is usually a group level such as
#'   \code{"lean"}; when \code{groupnames = TRUE}, it may include the model
#'   prefix, such as \code{"bmi_lean"}.
#' @param sensitivity Either \code{"keep"} to retain ANCOM-BC2 significant
#'   results, or \code{"robust_only"} to retain only results that also pass
#'   pseudo-count sensitivity analysis.
#' @param show_all Logical. Passed to \code{.prepare_ancombc2_lfc_data()}.
#'   If \code{FALSE}, retain taxa selected by the ANCOM-BC2 display rule; if
#'   \code{TRUE}, retain all available non-structural-zero taxa.
#' @param groupnames Logical. Passed to \code{.prepare_ancombc2_lfc_data()}.
#'   If \code{FALSE}, comparison labels omit the model prefix.
#' @param group_colors Optional named character vector assigning colors to the
#'   two \code{abundance_groups}. The reference group color is used for
#'   negative LFCs and the selected comparison group color is used for positive
#'   LFCs. Defaults to \code{c(reference = "blue", comparison = "red")} in
#'   the order of \code{abundance_groups}.
#' @param abundance_groups Character vector of length two specifying the groups
#'   shown in the abundance panel, in the order \code{c(reference, comparison)}.
#'   The second element must match \code{comparison}.
#' @param non_robust_alpha Numeric scalar between 0 and 1. Opacity used for LFC
#'   bars that are not pseudo-count sensitivity-robust. The default, \code{0.5},
#'   matches \code{make_barplots()}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{plot}{A patchwork object containing the abundance, LFC, p-value,
#'   and q-value panels.}
#'   \item{data}{The ANCOM-BC2 statistics used in the selected comparison.}
#'   \item{abundance_data}{Group-wise TSS abundance summaries used for the
#'   abundance panel.}
#'   \item{comparison}{The displayed ANCOM-BC2 comparison.}
#'   \item{reference_group}{The reference group inferred from
#'   \code{abundance_groups}.}
#'   \item{result_type}{The selected ANCOM-BC2 result table.}
#'   \item{inference_note}{Interpretive note provided by the ANCOM-BC2 helper.}
#' }
#'
#' @details The abundance panel is calculated after total-sum scaling within
#' each sample and displays mean relative abundance with standard error. Its
#' colors are matched to LFC direction: negative LFCs use the reference group
#' color and positive LFCs use the selected comparison group color.
#'
#' @examples
#' \dontrun{
#' fig <- make_abundance_lfc(
#'   ps = ps,
#'   out = output,
#'   result = "res_dunn",
#'   prefix = "bmi",
#'   group = "bmi",
#'   tax_level = "Family",
#'   comparison = "lean",
#'   sensitivity = "keep",
#'   show_all = FALSE,
#'   abundance_groups = c("overweight", "lean"),
#'   group_colors = c(overweight = "blue", lean = "red")
#' )
#'
#' fig$plot
#' }
#'
#' @family ANCOM-BC2 visualisation functions
#' @importFrom rlang .data .env
#' @export
make_abundance_lfc <- function(
    ps,
    out,
    result = base::c("res", "res_pair", "res_dunn", "res_global", "res_trend"),
    prefix = NULL,
    group,
    tax_level,
    comparison = NULL,
    sensitivity = base::c("keep", "robust_only"),
    show_all = FALSE,
    groupnames = FALSE,
    group_colors = NULL,
    abundance_groups = NULL,
    non_robust_alpha = 0.5
) {
  result <- base::match.arg(result)
  sensitivity <- base::match.arg(sensitivity)

  if (!base::is.numeric(non_robust_alpha) ||
      base::length(non_robust_alpha) != 1 ||
      base::is.na(non_robust_alpha) ||
      non_robust_alpha < 0 ||
      non_robust_alpha > 1) {
    base::stop("`non_robust_alpha` must be a single number between 0 and 1.")
  }

  prepared <- .prepare_ancombc2_lfc_data(
    out = out,
    result = result,
    prefix = prefix,
    sensitivity = sensitivity,
    show_all = show_all,
    groupnames = groupnames
  )

  available_comparisons <- base::unique(
    base::as.character(prepared$data$comparison)
  )

  if (base::is.null(comparison)) {
    if (base::length(available_comparisons) != 1) {
      base::stop(
        "Specify exactly one `comparison`. Available values: ",
        base::paste(available_comparisons, collapse = ", ")
      )
    }
    comparison <- available_comparisons
  }

  comparison_requested <- base::as.character(comparison)

  if (base::length(comparison_requested) != 1) {
    base::stop("`comparison` must contain exactly one ANCOM-BC2 coefficient.")
  }

  if (!comparison_requested %in% available_comparisons) {
    base::stop(
      "`comparison = ", comparison_requested, "` was not found. Available values: ",
      base::paste(available_comparisons, collapse = ", ")
    )
  }

  stat_df <- prepared$data |>
    dplyr::mutate(
      taxon = base::as.character(.data$taxon),
      comparison = base::as.character(.data$comparison),
      robust = base::as.logical(.data$robust)
    ) |>
    dplyr::filter(.data$comparison == .env$comparison_requested) |>
    dplyr::arrange(.data$lfc)

  if (base::nrow(stat_df) == 0) {
    base::stop("No ANCOM-BC2 data are available for the selected comparison.")
  }

  if (base::anyDuplicated(stat_df$taxon)) {
    base::stop(
      "The selected comparison contains duplicate taxon labels. ",
      "Check the ANCOM-BC2 result table before plotting."
    )
  }

  taxon_order <- stat_df$taxon

  sample_metadata <- data.frame(
    phyloseq::sample_data(ps)
  )

  if (!group %in% base::names(sample_metadata)) {
    base::stop("`group = ", group, "` is not present in sample_data(ps).")
  }

  sample_metadata$Sample <- base::rownames(sample_metadata)

  sample_df_all <- sample_metadata |>
    dplyr::transmute(
      Sample = .data$Sample,
      .group = base::as.character(.data[[group]])
    ) |>
    dplyr::filter(!base::is.na(.data$.group))

  available_groups <- base::unique(sample_df_all$.group)

  if (base::is.null(abundance_groups)) {
    abundance_groups <- available_groups
  } else {
    abundance_groups <- base::as.character(abundance_groups)

    if (base::anyDuplicated(abundance_groups)) {
      base::stop("`abundance_groups` must not contain duplicate group labels.")
    }

    missing_groups <- base::setdiff(abundance_groups, available_groups)

    if (base::length(missing_groups) > 0) {
      base::stop(
        "The following `abundance_groups` are not present in sample_data(ps): ",
        base::paste(missing_groups, collapse = ", "),
        ". Available values: ",
        base::paste(available_groups, collapse = ", ")
      )
    }
  }

  if (base::length(abundance_groups) != 2) {
    base::stop(
      "For a two-group abundance panel, `abundance_groups` must contain exactly ",
      "the reference group and the selected comparison group."
    )
  }

  if (!comparison_requested %in% abundance_groups) {
    base::stop(
      "`comparison = ", comparison_requested, "` must match one value in ",
      "`abundance_groups`. For example: abundance_groups = c(reference_group, ",
      comparison_requested, ")."
    )
  }

  reference_group <- base::setdiff(abundance_groups, comparison_requested)

  if (base::length(reference_group) != 1) {
    base::stop("Unable to identify a single reference group from `abundance_groups`.")
  }

  sample_df <- sample_df_all |>
    dplyr::filter(.data$.group %in% .env$abundance_groups) |>
    dplyr::mutate(
      .group = base::factor(.data$.group, levels = abundance_groups)
    )

  observed_groups <- base::unique(base::as.character(sample_df$.group))
  missing_after_filter <- base::setdiff(abundance_groups, observed_groups)

  if (base::length(missing_after_filter) > 0) {
    base::stop(
      "No samples remain for: ",
      base::paste(missing_after_filter, collapse = ", "),
      ". Check `abundance_groups` against sample_data(ps)[[group]]."
    )
  }

  if (base::is.null(group_colors)) {
    group_colors <- stats::setNames(
      base::c("blue", "red"),
      abundance_groups
    )
  } else {
    if (base::is.null(base::names(group_colors))) {
      base::stop(
        "`group_colors` must be a named vector. For example: ",
        "c(reference = 'blue', comparison = 'red')."
      )
    }

    missing_color_groups <- base::setdiff(abundance_groups, base::names(group_colors))

    if (base::length(missing_color_groups) > 0) {
      base::stop(
        "`group_colors` is missing colors for: ",
        base::paste(missing_color_groups, collapse = ", ")
      )
    }

    group_colors <- group_colors[abundance_groups]
  }

  lfc_colors <- base::c(
    negative = base::unname(group_colors[[reference_group]]),
    positive = base::unname(group_colors[[comparison_requested]])
  )

  stat_df <- stat_df |>
    dplyr::mutate(
      direction = base::ifelse(.data$lfc >= 0, "positive", "negative"),
      alpha_value = base::ifelse(.data$robust, 1, non_robust_alpha),
      p_label = .format_ancombc_p_value(.data$p),
      q_label = .format_ancombc_p_value(.data$q),
      taxon = base::factor(.data$taxon, levels = taxon_order),
      direction = base::factor(
        .data$direction,
        levels = base::c("negative", "positive")
      )
    )

  rel_ps <- phyloseq::transform_sample_counts(
    ps,
    function(x) {
      if (base::sum(x) == 0) x else x / base::sum(x)
    }
  )

  abundance_long <- phyloseq::psmelt(rel_ps) |>
    dplyr::transmute(
      Sample = base::as.character(.data$Sample),
      taxon = base::as.character(.data[[tax_level]]),
      relative_abundance = .data$Abundance
    ) |>
    dplyr::filter(.data$taxon %in% taxon_order) |>
    dplyr::group_by(.data$Sample, .data$taxon) |>
    dplyr::summarise(
      relative_abundance = base::sum(.data$relative_abundance),
      .groups = "drop"
    )

  abundance_summary <- tidyr::crossing(
    sample_df,
    taxon = taxon_order
  ) |>
    dplyr::left_join(
      abundance_long,
      by = base::c("Sample", "taxon")
    ) |>
    dplyr::mutate(
      relative_abundance = tidyr::replace_na(.data$relative_abundance, 0)
    ) |>
    dplyr::group_by(.data$taxon, .data$.group) |>
    dplyr::summarise(
      Mean = base::mean(.data$relative_abundance),
      SD = stats::sd(.data$relative_abundance),
      N = dplyr::n(),
      SE = .data$SD / base::sqrt(.data$N),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      taxon = base::factor(.data$taxon, levels = taxon_order),
      .group = base::factor(.data$.group, levels = abundance_groups)
    )

  abundance_plot <- ggplot2::ggplot(
    abundance_summary,
    ggplot2::aes(
      x = .data$taxon,
      y = .data$Mean,
      fill = .data$.group
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.65,
      color = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = base::pmax(.data$Mean - .data$SE, 0),
        ymax = .data$Mean + .data$SE
      ),
      position = ggplot2::position_dodge(width = 0.75),
      width = 0.18,
      linewidth = 0.4,
      color = "black",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      breaks = abundance_groups,
      drop = FALSE
    ) +
    ggplot2::scale_x_discrete(
      limits = taxon_order,
      drop = FALSE
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = base::c(0, 0.05))
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Average relative abundance (%)",
      fill = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 9),
      legend.key.size = grid::unit(0.35, "cm"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 2, 5.5, 5.5)
    )

  lfc_plot <- ggplot2::ggplot(
    stat_df,
    ggplot2::aes(
      x = .data$taxon,
      y = .data$lfc,
      fill = .data$direction,
      alpha = .data$alpha_value
    )
  ) +
    ggplot2::geom_col(
      width = 0.65,
      color = "black"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$lfc - .data$se,
        ymax = .data$lfc + .data$se
      ),
      width = 0.18,
      linewidth = 0.4,
      color = "black",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = lfc_colors,
      guide = "none"
    ) +
    ggplot2::scale_alpha_identity(guide = "none") +
    ggplot2::scale_x_discrete(
      limits = taxon_order,
      drop = FALSE
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Log fold change"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5.5, 2, 5.5, 2)
    )

  make_p_column <- function(data, label_col, panel_label) {
    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = 1,
        y = .data$taxon,
        label = .data[[label_col]]
      )
    ) +
      ggplot2::geom_text(size = 3) +
      ggplot2::scale_y_discrete(
        limits = taxon_order,
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        limits = base::c(0.5, 1.5),
        breaks = NULL
      ) +
      ggplot2::labs(x = panel_label, y = NULL) +
      ggplot2::theme_bw(base_size = 10) +
      ggplot2::theme(
        axis.text = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 10
        ),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5.5, 2, 5.5, 2)
      )
  }

  p_plot <- make_p_column(stat_df, "p_label", "p")
  q_plot <- make_p_column(stat_df, "q_label", "q")

  final_plot <- abundance_plot + lfc_plot + p_plot + q_plot +
    patchwork::plot_layout(
      widths = base::c(1.9, 1.1, 0.35, 0.35)
    )

  base::list(
    plot = final_plot,
    data = stat_df,
    abundance_data = abundance_summary,
    comparison = comparison_requested,
    reference_group = reference_group,
    result_type = result,
    inference_note = prepared$inference_note
  )
}

.format_ancombc_p_value <- function(x) {
  base::ifelse(
    base::is.na(x),
    "",
    base::ifelse(x < 0.001, "<0.001", base::sprintf("%.3f", x))
  )
}
