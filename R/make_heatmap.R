#' Create an ANCOM-BC2 log-fold change heatmap
#'
#' Creates a heatmap of ANCOM-BC2 log-fold changes. Non-significant cells
#' are displayed as white tiles. Structural-zero taxa are excluded because
#' their log-fold changes are not estimated by the standard model.
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
#'
#' @return A list with `data`, `structural_zero`, `result_type`,
#'   `inference_note`, and `plot`.
#' @importFrom rlang .data
#' @export
make_heatmap <- function(out,
                         result = base::c("res", "res_pair", "res_dunn", "res_global", "res_trend"),
                         prefix = NULL,
                         title = "ANCOM-BC2 log fold changes",
                         sensitivity = base::c("keep", "robust_only"),
                         show_all = TRUE,
                         groupnames = FALSE) {
  result <- base::match.arg(result)
  sensitivity <- base::match.arg(sensitivity)

  prepared <- .prepare_ancombc2_lfc_data(
    out = out,
    result = result,
    prefix = prefix,
    sensitivity = sensitivity,
    show_all = show_all,
    groupnames = groupnames
  )

  df_final <- prepared$data

  displayed_lfc <- df_final$lfc[df_final$display]

  max_abs_lfc <- if (base::length(displayed_lfc) == 0) {
    1
  } else {
    base::max(base::abs(displayed_lfc), na.rm = TRUE)
  }

  if (!base::is.finite(max_abs_lfc) || max_abs_lfc == 0) {
    max_abs_lfc <- 1
  }

  plot <- ggplot2::ggplot(
    df_final,
    ggplot2::aes(
      x = .data$comparison,
      y = .data$taxon,
      fill = .data$plot_value
    )
  ) +
    ggplot2::geom_tile(color = "black", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      limits = base::c(-max_abs_lfc, max_abs_lfc),
      name = "Log fold change"
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(df_final, .data$display),
      ggplot2::aes(
        label = .data$label,
        color = .data$robust
      ),
      fontface = "bold",
      size = 4
    ) +
    ggplot2::scale_color_manual(values = base::c("FALSE" = "black", "TRUE" = "aquamarine3"),
                                guide = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = title) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid = ggplot2::element_blank()
    )

  base::list(
    data = df_final,
    structural_zero = prepared$structural_zero,
    result_type = prepared$result_type,
    inference_note = prepared$inference_note,
    plot = plot
  )
}
