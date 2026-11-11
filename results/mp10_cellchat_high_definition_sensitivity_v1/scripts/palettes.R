# 最终验收版使用的颜色序列。
# 不在常量中固定项目组名；调用方通过 setNames(颜色, 组序) 关联新数据。
JIN_CLUSTER_COLORS <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
  "#FDBF6F", "#FF7F00", "#FB9A99"
)
JIN_STAGE_COLORS <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#FDBF6F", "#FB9A99")
JIN_KM_COLORS <- c("#D87869", "#1F78B4")
JIN_DEG_COLORS <- c(up = "#C8554E", down = "#427BA8", neutral = "#B8B8B8")
JIN_ORA_COLORS <- c("#FFF3B0", "#D84120")

# 为给定组序取得命名颜色向量。
# groups：按预期显示顺序排列的组名。
# colors：显式颜色表；无名称时按组序对应。省略时先查项目命名颜色表，再取 SCP 默认色。
# 返回值保留已有的颜色对应关系，不重新排序或重新设计配色。
jin_palette <- function(
  groups, colors = NULL, project_colors = getOption("jin.plotting.colors")
) {
  # 同一项目的所有分类图共用这一命名表；不按当前子集重新分配已有颜色。
  # 建议在开始绘图时设置 options(jin.plotting.colors = named_celltype_colors)。
  if (is.null(colors)) {
    if (!is.null(project_colors)) {
      if (is.null(names(project_colors)) || anyDuplicated(names(project_colors)) ||
          anyNA(project_colors) || any(!nzchar(names(project_colors)))) {
        stop("项目颜色表必须为名称唯一、无缺失的命名颜色向量。")
      }
      # 已定义的类型保留颜色；只有未定义类型才使用 SCP 的颜色。
      if (all(groups %in% names(project_colors))) {
        colors <- project_colors[groups]
      } else {
        all_groups <- unique(c(names(project_colors), as.character(groups)))
        colors <- SCP::palette_scp(factor(all_groups, levels = all_groups), type = "discrete")
        names(colors) <- all_groups
        colors[names(project_colors)] <- project_colors
        colors <- colors[groups]
      }
    } else {
      colors <- SCP::palette_scp(factor(groups, levels = groups), type = "discrete")
    }
  }

  if (is.null(names(colors))) {
    if (length(colors) != length(groups)) {
      stop("请为每个分组提供一个颜色。")
    }
    names(colors) <- groups
  }

  if (!all(groups %in% names(colors))) {
    stop("颜色表中缺少部分分组。")
  }
  colors
}
