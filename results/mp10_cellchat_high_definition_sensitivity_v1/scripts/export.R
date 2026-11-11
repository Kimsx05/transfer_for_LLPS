# 导出一张独立科研图；KM 的风险表、GSEA 的多个轨道属于同一单图。
# 有图例：.pdf、.png、_legend_only.pdf、_panel_only.pdf。
# 无图例：.pdf、.png。不生成内容重复的 _full.pdf。

# 从已构建的 gtable 中查找实际显示的图例，而非仅检查 colour/fill 映射。
# 隐藏图例通常是 zeroGrob，必须排除；嵌套 patchwork 中的图例也会检查。
jin_find_legends <- function(plot_grob) {
  legends <- list()
  if (inherits(plot_grob, "gtable")) {
    for (index in seq_along(plot_grob$grobs)) {
      child <- plot_grob$grobs[[index]]
      if (grepl("^guide-box", plot_grob$layout$name[index]) &&
          inherits(child, "gtable") && length(child$grobs) > 0) {
        legends[[length(legends) + 1]] <- child
      } else {
        legends <- c(legends, jin_find_legends(child))
      }
    }
  }
  legends
}

# p：完整图对象，或 KM 模板返回的列表。
# path：不含扩展名的路径；width、height、panel_width 均以英寸计。
# 图例默认自动识别。legend_plot / panel_plot 仍可显式覆盖用于拆分的组件。
# 多个独立图例的外部组合图需明确提供拆分组件，不能静默遗漏。
# 返回不可见的命名路径向量。此函数只导出，不变更数据或统计结果。
jin_export <- function(
  p, path, width, height,
  legend_plot = NULL, panel_plot = NULL, panel_width = width,
  font = "Arial"
) {
  jin_setup(font)
  if (!capabilities("cairo")) stop("需要 Cairo PDF 支持，以嵌入所需字体。")

  # KM 可传整个返回列表，也可只传 $plot；两种调用都能找到顶部图例。
  if (is.list(p) && !inherits(p, "ggplot") &&
      all(c("plot", "legend_plot", "panel_plot") %in% names(p))) {
    if (is.null(legend_plot)) legend_plot <- p$legend_plot
    if (is.null(panel_plot)) panel_plot <- p$panel_plot
    p <- p$plot
  }
  export_parts <- attr(p, "jin_export_parts")
  if (!is.null(export_parts)) {
    if (is.null(legend_plot)) legend_plot <- export_parts$legend_plot
    if (is.null(panel_plot)) panel_plot <- export_parts$panel_plot
  }

  # 用 Cairo 测量文字，避免 PostScript 字体替换影响图例尺寸。
  measurement_file <- tempfile(fileext = ".pdf")
  grDevices::cairo_pdf(filename = measurement_file, family = font)
  measurement_device <- dev.cur()
  null_device_files <- character()
  previous_null_device <- cowplot::set_null_device(function(width, height) {
    temporary_pdf <- tempfile(fileext = ".pdf")
    null_device_files <<- c(null_device_files, temporary_pdf)
    grDevices::cairo_pdf(
      filename = temporary_pdf, width = width, height = height, family = font
    )
  })
  on.exit({
    cowplot::set_null_device(previous_null_device)
    dev.off(measurement_device)
    unlink(c(measurement_file, null_device_files))
  }, add = TRUE)

  # 先确认实际图例，再写文件；防止遗漏参数时只导出完整图而不拆分。
  build_grob <- function(plot) {
    if (inherits(plot, "patchwork")) {
      patchwork::patchworkGrob(plot)
    } else {
      ggplotGrob(plot)
    }
  }
  actual_legends <- jin_find_legends(build_grob(p))
  has_legend <- length(actual_legends) > 0
  legend_grob <- NULL
  if (has_legend) {
    if (is.null(legend_plot)) {
      if (length(actual_legends) != 1) {
        stop("图中有多个独立图例；请显式提供 legend_plot 和 panel_plot。")
      }
      legend_grob <- actual_legends[[1]]
    } else {
      supplied_legends <- jin_find_legends(build_grob(legend_plot))
      if (length(supplied_legends) != 1) {
        stop("legend_plot 必须包含一个完整、可提取的图例区域。")
      }
      legend_grob <- supplied_legends[[1]]
    }
    if (is.null(panel_plot)) {
      panel_plot <- if (inherits(p, "patchwork")) {
        p & theme(legend.position = "none")
      } else {
        p + theme(legend.position = "none")
      }
    }
    if (length(jin_find_legends(build_grob(panel_plot))) > 0) {
      stop("panel_plot 中仍含图例；请提供无图例主体。")
    }
  } else if (!is.null(legend_plot)) {
    stop("完整图实际没有图例，不应额外提供 legend_plot。")
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  save_pdf <- function(plot, file, width, height) {
    ggsave(
      file, plot, device = grDevices::cairo_pdf,
      width = width, height = height, bg = "white", limitsize = FALSE
    )
  }
  save_pdf(p, paste0(path, ".pdf"), width, height)
  ggsave(
    paste0(path, ".png"), p, device = ragg::agg_png,
    width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE
  )
  output_files <- c(pdf = paste0(path, ".pdf"), png = paste0(path, ".png"))

  if (has_legend) {
    # 沿用验收版的图例留白：宽高各增加 0.15 英寸，最小 0.5 英寸。
    legend_width <- convertWidth(sum(legend_grob$widths), "in", valueOnly = TRUE) + 0.15
    legend_height <- convertHeight(sum(legend_grob$heights), "in", valueOnly = TRUE) + 0.15
    save_pdf(
      legend_grob, paste0(path, "_legend_only.pdf"),
      max(legend_width, 0.5), max(legend_height, 0.5)
    )
    save_pdf(panel_plot, paste0(path, "_panel_only.pdf"), panel_width, height)
    output_files <- c(
      output_files,
      legend_only = paste0(path, "_legend_only.pdf"),
      panel_only = paste0(path, "_panel_only.pdf")
    )
  }
  invisible(output_files)
}
