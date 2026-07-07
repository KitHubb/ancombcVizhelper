#' Create ANCOM-BC2 log-fold change bar plots
#'
#' Creates one bar plot per ANCOM-BC2 comparison. Positive log-fold changes
#' are shown in red and negative log-fold changes in blue. Error bars represent
#' plus or minus one standard error.
#'
#' @param out Full result object returned by `ANCOMBC::ancombc2()`.
#' @param result Result table to visualize.
#' @param prefix Prefix used to identify the target model coefficient.
#' @param groupnames Logical; include the group-variable prefix in comparison
#'   labels (e.g., `bmi_lean`). When `FALSE`, only group levels are shown
#'   (e.g., `lean`).
#' @param title Plot title.
#' @param sensitivity Whether to retain all significant results or only
#'   pseudo-count sensitivity-robust results.
#' @param show_all Logical; retain all taxa or only taxa with at least one
#'   displayed result.
#' @param group_order Taxon ordering method: `"none"`, `"mean"`, or a
#'   specific comparison name.
#' @param order Taxon order direction.
#'
#' @return A list with `data`, `structural_zero`, `result_type`,
#'   `inference_note`, and `plot`.
#' @importFrom rlang .data
#' @export
make_barplots <- function(out,
                          result = base::c("res", "res_pair", "res_dunn", "res_global", "res_trend"),
                          prefix = NULL,
                          title = "ANCOM-BC2 log fold changes",
                          sensitivity = base::c("keep", "robust_only"),
                          show_all = FALSE,
                          group_order = "mean",
                          order = base::c("asc", "desc"),
                          groupnames = FALSE) {
  result <- base::match.arg(result)
  sensitivity <- base::match.arg(sensitivity)
  order <- base::match.arg(order)

  prepared <- .prepare_ancombc2_lfc_data(
    out = out,
    result = result,
    prefix = prefix,
    sensitivity = sensitivity,
    show_all = show_all,
    groupnames = groupnames
  )

  df_final <- prepared$data |>
    dplyr::mutate(
      taxon = base::as.character(.data$taxon),
      comparison = base::as.character(.data$comparison)
    )

  if (sensitivity == "robust_only") {
    df_final <- df_final |>
      dplyr::mutate(bar_lfc = base::ifelse(.data$robust, .data$lfc, NA_real_))
  } else {
    df_final <- df_final |>
      dplyr::mutate(bar_lfc = .data$lfc)
  }

  df_final <- df_final |>
    dplyr::mutate(
      se_lower = .data$bar_lfc -
        .data$se,
      se_upper = .data$bar_lfc +
        .data$se,
      direction = base::ifelse(.data$bar_lfc >= 0, "positive", "negative"),
      alpha_value = base::ifelse(.data$robust, 1, 0.5)
    )

  displayed_lfc <- df_final |>
    dplyr::filter(base::is.finite(.data$bar_lfc))

  if (base::nrow(displayed_lfc) == 0) {
    base::stop("No log-fold changes are available for display. ",
               "There were no significant contrasts or no contrasts passed the sensitivity criterion.")
  }

  common_min <- base::min(displayed_lfc$se_lower, 0, na.rm = TRUE)

  common_max <- base::max(displayed_lfc$se_upper, 0, na.rm = TRUE)

  common_range <- common_max - common_min

  if (!base::is.finite(common_range) || common_range == 0) {
    common_range <- 1
  }

  padding <- common_range * 0.06

  common_limits <- base::c(common_min - padding, common_max + padding)

  groups <- base::levels(prepared$data$comparison)

  if (base::is.null(groups)) {
    groups <- base::unique(df_final$comparison)
  }

  if (group_order == "none") {
    taxon_order <- base::unique(df_final$taxon)

  } else if (group_order == "mean") {
    taxon_summary <- df_final |>
      dplyr::group_by(.data$taxon) |>
      dplyr::summarise(
        mean_lfc = base::mean(.data$bar_lfc, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(mean_lfc = base::ifelse(
        base::is.nan(.data$mean_lfc),
        NA_real_,
        .data$mean_lfc
      ))

    if (order == "asc") {
      taxon_order <- taxon_summary |>
        dplyr::arrange(base::is.na(.data$mean_lfc),
                       .data$mean_lfc) |>
        dplyr::pull(.data$taxon)
    } else {
      taxon_order <- taxon_summary |>
        dplyr::arrange(base::is.na(.data$mean_lfc),
                       dplyr::desc(.data$mean_lfc)) |>
        dplyr::pull(.data$taxon)
    }

  } else if (group_order %in% groups) {
    taxon_summary <- df_final |>
      dplyr::filter(.data$comparison == group_order) |>
      dplyr::select(dplyr::all_of(base::c("taxon", "bar_lfc")))

    if (order == "asc") {
      taxon_order <- taxon_summary |>
        dplyr::arrange(base::is.na(.data$bar_lfc),
                       .data$bar_lfc) |>
        dplyr::pull(.data$taxon)
    } else {
      taxon_order <- taxon_summary |>
        dplyr::arrange(base::is.na(.data$bar_lfc),
                       dplyr::desc(.data$bar_lfc)) |>
        dplyr::pull(.data$taxon)
    }

  } else {
    base::stop(
      "`group_order` must be 'none', 'mean', or one of the following comparisons: ",
      base::paste(groups, collapse = ", ")
    )
  }

  df_final <- df_final |>
    dplyr::mutate(
      taxon = base::factor(.data$taxon, levels = taxon_order),
      comparison = base::factor(.data$comparison, levels = groups),
      direction = base::factor(.data$direction, levels = base::c("negative", "positive"))
    )

  plots <- base::lapply(base::seq_along(groups), function(i) {
    g <- groups[i]

    panel_data <- df_final |>
      dplyr::filter(.data$comparison == g)

    p <- ggplot2::ggplot(
      panel_data,
      ggplot2::aes(
        x = .data$taxon,
        y = .data$bar_lfc,
        fill = .data$direction,
        alpha = .data$alpha_value
      )
    ) +
      ggplot2::geom_col(width = 0.7,
                        color = "black",
                        na.rm = TRUE) +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = .data$se_lower,
          ymax = .data$se_upper
        ),
        width = 0.18,
        color = "black",
        linewidth = 0.4,
        na.rm = TRUE
      ) +
      ggplot2::scale_fill_manual(values = base::c("negative" = "blue", "positive" = "red"),
                                 guide = "none") +
      ggplot2::scale_alpha_identity(guide = "none") +
      ggplot2::scale_y_continuous(limits = common_limits, expand = base::c(0, 0)) +
      ggplot2::coord_flip() +
      ggplot2::labs(title = g, x = NULL, y = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text.x = ggplot2::element_text(size = 10),
        plot.margin = ggplot2::margin(3, 6, 3, 3)
      )


    if (i == 1) {
      p <- p +
        ggplot2::scale_x_discrete(drop = FALSE)
    } else {
      p <- p +
        ggplot2::scale_x_discrete(drop = FALSE, labels = NULL) +
        ggplot2::theme(axis.ticks.y = ggplot2::element_blank())
    }

    p
  })

  axis_label <- patchwork::wrap_elements(
    full = grid::textGrob(
      "Log fold change",
      x = 0.5,
      y = 0.5,
      hjust = 0.5,
      vjust = 0.5,
      gp = grid::gpar(fontsize = 11)
    )
  )

  plot_row <- patchwork::wrap_plots(plots, nrow = 1) +
    patchwork::plot_annotation(title = title)

  final_plot <- (plot_row / axis_label) +
    patchwork::plot_layout(heights = base::c(1, 0.06))

  base::list(
    data = df_final,
    structural_zero = prepared$structural_zero,
    result_type = prepared$result_type,
    inference_note = prepared$inference_note,
    plot = final_plot
  )
}
