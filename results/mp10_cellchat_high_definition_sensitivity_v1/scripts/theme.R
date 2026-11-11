# 公共主题、字体检查与绘图数据整理。
# 加载顺序：theme.R -> palettes.R -> export.R -> 所需模板的 plot.R。
# 本文件只处理绘图所需的格式，不筛选统计样本、不计算分数或检验结果。
suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(grid)
})

# 检查字体，并统一 geom_text / geom_label 的默认字体与字号。
# font：系统中已安装的字体族；默认必须是真实 Arial，而非替代字体。
# 返回：不可见的字体名；同时更新本次 R 会话的文字图层默认设置。
jin_setup <- function(font = "Arial") {
  if (!requireNamespace("systemfonts", quietly = TRUE)) {
    stop("请先安装 systemfonts，以检查所需字体是否可用。")
  }

  resolved_family <- systemfonts::font_info(family = font)$family[1]
  if (tolower(resolved_family) != tolower(font)) {
    stop("所需字体不可用：", font, "。请先安装该字体，不要使用静默字体替换。")
  }

  # ggplot2 的文字图层使用 mm，而 theme 的字号使用 pt；这里明确换算。
  update_geom_defaults("text", list(family = font, size = 12 / ggplot2::.pt))
  update_geom_defaults("label", list(family = font, size = 12 / ggplot2::.pt))
  invisible(font)
}

# 返回已验收的公共主题；各模板只覆盖自身特有的轴、边距或图例设置。
# base_size：以 pt 为单位的基础字号。常规线宽 0.4 mm 约等于 1.13 pt。
jin_theme <- function(font = "Arial", base_size = 12) {
  theme_classic(
    base_size = base_size,
    base_family = font,
    base_line_size = 0.4
  ) +
    theme(
      text = element_text(family = font, colour = "black"),
      axis.text = element_text(size = base_size, colour = "black"),
      axis.title = element_text(size = base_size),
      line = element_line(linewidth = 0.4),
      rect = element_rect(linewidth = 0.4),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption = element_blank(),
      strip.background = element_rect(
        fill = "#E5E5E5", colour = "black", linewidth = 0.4
      ),
      strip.text = element_text(size = base_size),
      plot.margin = margin(9, 12, 9, 9)
    )
}

# 将调用方列名映射为模板内部的统一列名，保留原来的行顺序及数值。
# columns 为命名字符向量，例如 c(x = "UMAP_1", y = "UMAP_2")：
# 向量的名称是内部列名，向量的值是输入数据中的列名。
# 返回独立的列子集，不通过引用修改调用方的 data.table。
jin_frame <- function(data, columns) {
  stopifnot(all(unname(columns) %in% names(data)))
  plot_data <- data.table::as.data.table(data)[, unname(columns), with = FALSE]
  data.table::setnames(plot_data, names(columns))
  plot_data
}

# 优先使用调用方指定的顺序；否则沿用因子水平或首次出现顺序。
# 不擅自对组名排序，也不静默丢弃缺失或未列出的组。
jin_order <- function(x, order = NULL) {
  if (is.null(order)) {
    order <- if (is.factor(x)) levels(x) else unique(as.character(x))
  }
  if (anyNA(x) || !all(as.character(x) %in% order)) {
    stop("存在缺失或未列出的分组；请提供已整理好的绘图数据和明确组序。")
  }
  order
}

# 只格式化已有 P 值：小于 0.001 时用科学计数法，否则保留三位小数。
# p = / P = 前缀由不同模板添加；本函数不执行统计检验。
jin_p_label <- function(p) {
  stopifnot(length(p) == 1, is.finite(p), p >= 0, p <= 1)
  if (p < 0.001) {
    formatC(p, format = "e", digits = 2)
  } else {
    formatC(p, format = "f", digits = 3)
  }
}

# 构建 UMAP 专用的 SCP blank 主题、箭头轴和等比例坐标范围。
# d 必须包含数值列 x、y。返回可通过 + 加入 ggplot 的组件列表。
# 保留组件叠加顺序，因为后加入的 theme 会覆盖前面的部分设置。
jin_blank_umap <- function(d, font = "Arial", xlab = "UMAP1", ylab = "UMAP2") {
  # 不再索引 SCP 返回值的内部层级。箭头使用公开 grid/ggplot2 接口，
  # 坐标、字号、线宽与 SCP 0.5.6 验收版完全一致；空白主题仍调用 SCP。
  coordinate_grob <- grobTree(gList(
    linesGrob(
      x = unit(c(0, 0.15), "npc"), y = unit(c(0, 0), "npc"),
      arrow = arrow(length = unit(0.02, "npc")),
      gp = gpar(lwd = 2, fontfamily = font)
    ),
    textGrob(
      label = xlab, x = unit(0, "npc"), y = unit(0, "npc"),
      vjust = 4 / 3, hjust = 0, gp = gpar(fontsize = 12, fontfamily = font)
    ),
    linesGrob(
      x = unit(c(0, 0), "npc"), y = unit(c(0, 0.15), "npc"),
      arrow = arrow(length = unit(0.02, "npc")),
      gp = gpar(lwd = 2, fontfamily = font)
    ),
    textGrob(
      label = ylab, x = unit(0, "npc"), y = unit(0, "npc"),
      vjust = -2 / 3, hjust = 0, rot = 90,
      gp = gpar(fontsize = 12, fontfamily = font)
    )
  ))
  blank_theme <- list(
    annotation_custom(coordinate_grob),
    SCP::theme_blank(add_coord = FALSE),
    coord_cartesian(clip = "off")
  )

  # 取两个方向中较大的跨度，并保留原版 8% 留白，使主体近似正方形。
  coordinate_span <- max(diff(range(d$x)), diff(range(d$y))) * 1.08
  list(
    blank_theme,
    jin_theme(font),
    SCP::theme_blank(add_coord = FALSE),
    theme(
      text = element_text(family = font),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.margin = margin(18, 15, 24, 24)
    ),
    coord_fixed(
      xlim = mean(range(d$x)) + c(-0.5, 0.5) * coordinate_span,
      ylim = mean(range(d$y)) + c(-0.5, 0.5) * coordinate_span,
      expand = FALSE,
      clip = "off"
    )
  )
}
