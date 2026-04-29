#!/usr/bin/env Rscript
# ==============================================================================
# Script: 04_run_trajectory_analysis.R
# Authors: Alfredo Pherez-Farah (ORCID: 0000-0003-2213-3405); Willem de Koning (ORCID: 0000-0002-4594-8423)
# Purpose: Run Monocle3 trajectory analysis, pseudotime visualization, branch convergence analysis, and incoming-signature summaries.
# Inputs: Processed Seurat object with UMAP and final cluster annotations.
# Outputs: PNG trajectory plots and CSV summary tables.
# Dependencies: Seurat, monocle3, SeuratWrappers, igraph, Matrix, dplyr, tidyr, ggplot2, ggrepel, purrr, scales, patchwork.
# Environment: Main analysis environment (conda + renv).
# Notes:
# - Edit the configuration section before running.
# - This script uses the existing Seurat UMAP embedding for Monocle3 trajectory learning.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(igraph)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(purrr)
  library(scales)
  library(patchwork)
})

set.seed(1)

# ------------------------------------------------------------------------------
# 1. CONFIGURATION  ← EDIT THIS SECTION
# ------------------------------------------------------------------------------

input_file <- "path/to/combined_seurat_normalized.rds"
output_dir <- "path/to/trajectory_analysis"

cluster_column <- "final_clusters"
umap_reduction <- "umap"

num_dim <- 50
minimal_branch_len <- 14

# Trajectory groups
myeloid_clusters <- c("0", "2", "3", "4", "5", "13", "14", "17")
dc_clusters <- c("1", "7", "8", "9", "13")

# Root and convergence clusters
myeloid_root_clusters <- c("4")
dc_root_clusters <- c("1", "9")
dc_convergence_cluster <- "7"

# Convergence analysis
run_convergence_analysis <- TRUE
min_detection_rate <- 0.05
spearman_q_threshold <- 0.05
spearman_rho_threshold <- 0.2
near_progress_min <- 0.85
upstream_progress_max <- 0.20
top_trend_genes_to_plot <- 30

# Combined trajectory overlay
run_combined_overlay <- TRUE
myeloid_origin_clusters <- c(b1 = "2", b2 = "4")
dc_origin_clusters <- c(b1 = "1", b2 = "9")

# Incoming signature summary
run_incoming_signature_summary <- TRUE

# Density-weighted trajectory graph
run_density_graph <- TRUE
density_radius <- 0.25

# ------------------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ------------------------------------------------------------------------------

make_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

save_png <- function(plot, filename, width = 7, height = 5, dpi = 300) {
  ggsave(filename, plot = plot, width = width, height = height, dpi = dpi)
}

require_metadata <- function(object, columns) {
  missing_columns <- setdiff(columns, colnames(object@meta.data))
  if (length(missing_columns) > 0) {
    stop("Missing required metadata columns: ", paste(missing_columns, collapse = ", "))
  }
}

as_character_clusters <- function(object, cluster_column) {
  object[[cluster_column]] <- as.character(object[[cluster_column]][, 1])
  Idents(object) <- object[[cluster_column]][, 1]
  object
}

build_cds_from_seurat <- function(seurat_object, cluster_column, umap_reduction, num_dim) {
  cds <- as.cell_data_set(seurat_object)

  colData(cds)$seurat_clusters <- as.character(seurat_object@meta.data[[cluster_column]])

  seurat_umap <- Embeddings(seurat_object, reduction = umap_reduction)
  seurat_umap <- seurat_umap[colnames(cds), , drop = FALSE]
  reducedDims(cds)$UMAP <- seurat_umap

  cds <- preprocess_cds(cds, num_dim = num_dim)

  cds@clusters$UMAP <- list()
  cds@clusters$UMAP$clusters <- factor(colData(cds)$seurat_clusters)
  names(cds@clusters$UMAP$clusters) <- colnames(cds)

  cds@clusters$UMAP$partitions <- factor(rep(1, ncol(cds)))
  names(cds@clusters$UMAP$partitions) <- colnames(cds)

  cds
}

subset_cds_by_clusters <- function(cds, clusters) {
  cells_keep <- colnames(cds)[as.character(colData(cds)$seurat_clusters) %in% as.character(clusters)]
  cds_sub <- cds[, cells_keep]

  cds_sub@clusters$UMAP <- list()
  cds_sub@clusters$UMAP$clusters <- factor(colData(cds_sub)$seurat_clusters)
  names(cds_sub@clusters$UMAP$clusters) <- colnames(cds_sub)

  cds_sub@clusters$UMAP$partitions <- factor(rep(1, ncol(cds_sub)))
  names(cds_sub@clusters$UMAP$partitions) <- colnames(cds_sub)

  cds_sub
}

get_closest_vertex_df <- function(cds) {
  closest_vertex <- principal_graph_aux(cds)[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.data.frame(closest_vertex)
  colnames(closest_vertex) <- "closest_vertex"
  closest_vertex$cell_id <- rownames(closest_vertex)
  closest_vertex
}

get_majority_root_nodes <- function(cds, root_clusters) {
  closest_vertex <- get_closest_vertex_df(cds)

  root_nodes <- vapply(root_clusters, function(cluster_id) {
    cluster_cells <- colnames(cds)[as.character(colData(cds)$seurat_clusters) == as.character(cluster_id)]
    raw_node <- names(which.max(table(closest_vertex[cluster_cells, "closest_vertex"])))
    paste0("Y_", raw_node)
  }, character(1))

  names(root_nodes) <- as.character(root_clusters)
  root_nodes
}

get_node_coordinates <- function(cds) {
  dp_mst <- principal_graph_aux(cds)[["UMAP"]]$dp_mst

  node_coordinates <- if (is.list(dp_mst) && !is.null(dp_mst$Y)) {
    as.data.frame(t(dp_mst$Y))
  } else {
    as.data.frame(t(dp_mst))
  }

  node_coordinates <- tibble::rownames_to_column(node_coordinates, var = "name")
  colnames(node_coordinates)[2:3] <- c("x", "y")

  node_coordinates
}

learn_subset_trajectory <- function(
  cds,
  clusters,
  root_clusters,
  minimal_branch_len
) {
  cds_sub <- subset_cds_by_clusters(cds, clusters)

  cds_sub <- learn_graph(
    cds_sub,
    use_partition = FALSE,
    learn_graph_control = list(minimal_branch_len = minimal_branch_len)
  )

  root_nodes <- get_majority_root_nodes(cds_sub, root_clusters)

  cds_sub <- order_cells(
    cds_sub,
    reduction_method = "UMAP",
    root_pr_nodes = unname(root_nodes)
  )

  list(
    cds = cds_sub,
    root_nodes = root_nodes
  )
}

plot_root_nodes <- function(cds, plot, root_nodes) {
  node_coordinates <- get_node_coordinates(cds)

  highlight_nodes <- node_coordinates %>%
    filter(name %in% unname(root_nodes))

  plot +
    geom_point(
      data = highlight_nodes,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      color = "black",
      size = 4,
      shape = 21,
      fill = "red"
    )
}

get_graph_segments <- function(cds, label) {
  graph <- principal_graph(cds)[["UMAP"]]
  node_coordinates <- get_node_coordinates(cds)

  edges <- as.data.frame(igraph::as_edgelist(graph))
  colnames(edges) <- c("from", "to")

  edges %>%
    left_join(node_coordinates, by = c("from" = "name")) %>%
    rename(x = x, y = y) %>%
    left_join(node_coordinates, by = c("to" = "name"), suffix = c("", ".end")) %>%
    transmute(
      x = x,
      y = y,
      xend = x.end,
      yend = y.end,
      graph = label
    )
}

scale01 <- function(x) {
  (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE) + 1e-8)
}

# ------------------------------------------------------------------------------
# 3. LOAD DATA AND BUILD CDS
# ------------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

output_dir <- make_dir(output_dir)
trajectory_dir <- make_dir(file.path(output_dir, "monocle_trajectories"))
convergence_dir <- make_dir(file.path(output_dir, "convergence_analysis"))

seurat_object <- readRDS(input_file)
require_metadata(seurat_object, cluster_column)

if (!umap_reduction %in% Reductions(seurat_object)) {
  stop("UMAP reduction not found: ", umap_reduction)
}

seurat_object <- as_character_clusters(seurat_object, cluster_column)

cds <- build_cds_from_seurat(
  seurat_object = seurat_object,
  cluster_column = cluster_column,
  umap_reduction = umap_reduction,
  num_dim = num_dim
)

message("Loaded Seurat object and created Monocle3 CDS.")
message("Cells: ", ncol(cds))
message("Clusters: ", paste(sort(unique(colData(cds)$seurat_clusters)), collapse = ", "))

# ------------------------------------------------------------------------------
# 4. MYELOID TRAJECTORY
# ------------------------------------------------------------------------------

myeloid <- learn_subset_trajectory(
  cds = cds,
  clusters = myeloid_clusters,
  root_clusters = myeloid_root_clusters,
  minimal_branch_len = minimal_branch_len
)

p_myeloid_pseudotime <- plot_cells(
  myeloid$cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE
) +
  ggtitle("Myeloid trajectory pseudotime")

p_myeloid_clusters <- plot_cells(
  myeloid$cds,
  color_cells_by = "seurat_clusters",
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE
) +
  ggtitle("Myeloid trajectory by cluster")

p_myeloid_roots <- plot_root_nodes(
  cds = myeloid$cds,
  plot = p_myeloid_pseudotime,
  root_nodes = myeloid$root_nodes
) +
  ggtitle("Myeloid trajectory with root node")

save_png(
  p_myeloid_pseudotime,
  file.path(trajectory_dir, "myeloid_pseudotime.png"),
  width = 7,
  height = 6
)

save_png(
  p_myeloid_clusters,
  file.path(trajectory_dir, "myeloid_clusters.png"),
  width = 7,
  height = 6
)

save_png(
  p_myeloid_roots,
  file.path(trajectory_dir, "myeloid_pseudotime_with_root.png"),
  width = 7,
  height = 6
)

write.csv(
  data.frame(root_cluster = names(myeloid$root_nodes), root_node = unname(myeloid$root_nodes)),
  file.path(trajectory_dir, "myeloid_root_nodes.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 5. DC TRAJECTORY
# ------------------------------------------------------------------------------

dc <- learn_subset_trajectory(
  cds = cds,
  clusters = dc_clusters,
  root_clusters = dc_root_clusters,
  minimal_branch_len = minimal_branch_len
)

p_dc_pseudotime <- plot_cells(
  dc$cds,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE
) +
  ggtitle("DC trajectory pseudotime")

p_dc_clusters <- plot_cells(
  dc$cds,
  color_cells_by = "seurat_clusters",
  label_cell_groups = FALSE,
  show_trajectory_graph = TRUE
) +
  ggtitle("DC trajectory by cluster")

p_dc_roots <- plot_root_nodes(
  cds = dc$cds,
  plot = p_dc_pseudotime,
  root_nodes = dc$root_nodes
) +
  ggtitle("DC trajectory with root nodes")

save_png(
  p_dc_pseudotime,
  file.path(trajectory_dir, "dc_pseudotime.png"),
  width = 7,
  height = 6
)

save_png(
  p_dc_clusters,
  file.path(trajectory_dir, "dc_clusters.png"),
  width = 7,
  height = 6
)

save_png(
  p_dc_roots,
  file.path(trajectory_dir, "dc_pseudotime_with_roots.png"),
  width = 7,
  height = 6
)

write.csv(
  data.frame(root_cluster = names(dc$root_nodes), root_node = unname(dc$root_nodes)),
  file.path(trajectory_dir, "dc_root_nodes.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. CONVERGENCE ANALYSIS FOR DC TRAJECTORY
# ------------------------------------------------------------------------------

if (run_convergence_analysis) {
  message("Running convergence analysis for DC trajectory.")

  cds_dc <- dc$cds
  graph <- principal_graph(cds_dc)[["UMAP"]]
  closest_vertex <- get_closest_vertex_df(cds_dc)

  convergence_cells <- colnames(cds_dc)[
    as.character(colData(cds_dc)$seurat_clusters) == dc_convergence_cluster
  ]

  convergence_raw_node <- names(which.max(table(
    closest_vertex[convergence_cells, "closest_vertex"]
  )))
  convergence_node <- paste0("Y_", convergence_raw_node)

  root_nodes <- dc$root_nodes
  root1_node <- unname(root_nodes[1])
  root2_node <- unname(root_nodes[2])

  path1 <- shortest_paths(graph, from = root1_node, to = convergence_node, mode = "all")$vpath[[1]]
  path2 <- shortest_paths(graph, from = root2_node, to = convergence_node, mode = "all")$vpath[[1]]

  path1_nodes <- names(path1)
  path2_nodes <- names(path2)

  path1_only <- setdiff(path1_nodes, path2_nodes)
  path2_only <- setdiff(path2_nodes, path1_nodes)

  first_touch <- function(node) {
    path_to_convergence <- shortest_paths(
      graph,
      from = node,
      to = convergence_node,
      mode = "all"
    )$vpath[[1]]

    if (length(path_to_convergence) == 0) {
      return(NA_character_)
    }

    path_names <- names(path_to_convergence)

    branch1_hits <- which(path_names %in% path1_only)
    branch2_hits <- which(path_names %in% path2_only)

    if (length(branch1_hits) > 0 &&
        (length(branch2_hits) == 0 || min(branch1_hits) < min(branch2_hits))) {
      return("b1")
    }

    if (length(branch2_hits) > 0 &&
        (length(branch1_hits) == 0 || min(branch2_hits) < min(branch1_hits))) {
      return("b2")
    }

    distance_to_root1 <- as.numeric(distances(graph, v = node, to = root1_node, mode = "all"))
    distance_to_root2 <- as.numeric(distances(graph, v = node, to = root2_node, mode = "all"))

    if (is.finite(distance_to_root1) && is.finite(distance_to_root2)) {
      if (distance_to_root1 < distance_to_root2) return("b1")
      if (distance_to_root2 < distance_to_root1) return("b2")
    }

    NA_character_
  }

  node_names <- V(graph)$name
  branch_of_node <- setNames(
    vapply(node_names, first_touch, character(1)),
    node_names
  )

  cell_node_names <- paste0("Y_", closest_vertex$closest_vertex)
  names(cell_node_names) <- closest_vertex$cell_id

  b1_cells <- names(cell_node_names)[branch_of_node[cell_node_names] == "b1"]
  b2_cells <- names(cell_node_names)[branch_of_node[cell_node_names] == "b2"]

  neighborhood_nodes <- names(neighborhood(graph, order = 1, nodes = convergence_node, mode = "all")[[1]])
  neighborhood_raw <- sub("^Y_", "", neighborhood_nodes)

  b1_cells <- setdiff(
    b1_cells,
    closest_vertex$cell_id[closest_vertex$closest_vertex %in% neighborhood_raw]
  )

  b2_cells <- setdiff(
    b2_cells,
    closest_vertex$cell_id[closest_vertex$closest_vertex %in% neighborhood_raw]
  )

  message("Branch sizes after excluding immediate convergence neighborhood:")
  message("b1: ", length(b1_cells))
  message("b2: ", length(b2_cells))

  graph_distances <- distances(graph, v = V(graph), to = V(graph), mode = "all", weights = NA)
  distance_to_convergence <- graph_distances[, match(convergence_node, V(graph)$name)]
  names(distance_to_convergence) <- V(graph)$name

  cell_distance <- distance_to_convergence[cell_node_names]
  names(cell_distance) <- names(cell_node_names)

  b1_progress <- 1 - scale01(cell_distance[b1_cells])
  b2_progress <- 1 - scale01(cell_distance[b2_cells])

  mat <- counts(cds_dc)
  cpm <- t(t(mat) / Matrix::colSums(mat)) * 1e6
  logexpr <- log1p(cpm)

  detection_rate <- function(x) mean(x > 0, na.rm = TRUE)

  keep_genes <- rownames(logexpr)[
    (apply(logexpr[, b1_cells, drop = FALSE] > 0, 1, detection_rate) >= min_detection_rate) &
      (apply(logexpr[, b2_cells, drop = FALSE] > 0, 1, detection_rate) >= min_detection_rate)
  ]

  logexpr_filtered <- logexpr[keep_genes, , drop = FALSE]

  spearman_branch <- function(expr_mat, cells, progress_vec) {
    expr_sub <- expr_mat[, cells, drop = FALSE]
    progress <- progress_vec[cells]

    rhos <- numeric(nrow(expr_sub))
    p_values <- numeric(nrow(expr_sub))

    for (i in seq_len(nrow(expr_sub))) {
      expression <- as.numeric(expr_sub[i, ])

      if (sd(expression) == 0 || sd(progress) == 0) {
        rhos[i] <- 0
        p_values[i] <- 1
      } else {
        test <- suppressWarnings(
          cor.test(expression, progress, method = "spearman", exact = FALSE)
        )
        rhos[i] <- unname(test$estimate)
        p_values[i] <- test$p.value
      }
    }

    tibble(
      gene = rownames(expr_sub),
      rho = rhos,
      p = p_values,
      q = p.adjust(p_values, method = "BH")
    )
  }

  res1 <- spearman_branch(logexpr_filtered, b1_cells, b1_progress)
  res2 <- spearman_branch(logexpr_filtered, b2_cells, b2_progress)

  sig1 <- res1 %>%
    filter(q < spearman_q_threshold, rho > spearman_rho_threshold)

  sig2 <- res2 %>%
    filter(q < spearman_q_threshold, rho > spearman_rho_threshold)

  convergent <- inner_join(
    sig1 %>% select(gene, rho1 = rho, q1 = q),
    sig2 %>% select(gene, rho2 = rho, q2 = q),
    by = "gene"
  )

  near_convergence_cells <- closest_vertex$cell_id[
    closest_vertex$closest_vertex %in% neighborhood_raw
  ]

  b1_upstream <- names(sort(b1_progress, decreasing = FALSE))[
    seq_len(max(10, floor(upstream_progress_max * length(b1_progress))))
  ]

  b2_upstream <- names(sort(b2_progress, decreasing = FALSE))[
    seq_len(max(10, floor(upstream_progress_max * length(b2_progress))))
  ]

  peak_filter <- function(genes, logexpr_mat) {
    keep <- logical(length(genes))

    for (i in seq_along(genes)) {
      gene <- genes[i]
      expression <- as.numeric(logexpr_mat[gene, ])
      names(expression) <- colnames(logexpr_mat)

      near_mean <- mean(expression[intersect(near_convergence_cells, names(expression))], na.rm = TRUE)
      b1_mean <- mean(expression[intersect(b1_upstream, names(expression))], na.rm = TRUE)
      b2_mean <- mean(expression[intersect(b2_upstream, names(expression))], na.rm = TRUE)

      keep[i] <- is.finite(near_mean) && near_mean > b1_mean && near_mean > b2_mean
    }

    genes[keep]
  }

  convergent_genes <- peak_filter(convergent$gene, logexpr_filtered)

  convergent_ranked <- convergent %>%
    filter(gene %in% convergent_genes) %>%
    mutate(rho_sum = rho1 + rho2) %>%
    arrange(desc(rho_sum))

  write.csv(
    convergent_ranked,
    file.path(convergence_dir, "convergent_genes.csv"),
    row.names = FALSE
  )

  exclusive_b1 <- res1 %>%
    rename(rho1 = rho, p1 = p, q1 = q) %>%
    inner_join(
      res2 %>% rename(rho2 = rho, p2 = p, q2 = q),
      by = "gene"
    ) %>%
    filter(
      q1 < spearman_q_threshold,
      rho1 > spearman_rho_threshold,
      (rho2 <= 0 | q2 >= 0.20),
      (rho1 - rho2) > 0.20
    ) %>%
    arrange(desc(rho1 - rho2))

  write.csv(
    exclusive_b1,
    file.path(convergence_dir, "exclusive_branch1_genes.csv"),
    row.names = FALSE
  )

  mean_window <- function(gene, cells, progress, logexpr_mat) {
    progress_values <- progress[cells]
    expression <- as.numeric(logexpr_mat[gene, cells])

    near <- mean(expression[progress_values >= near_progress_min], na.rm = TRUE)
    upstream <- mean(expression[progress_values <= upstream_progress_max], na.rm = TRUE)

    c(
      near = near,
      upstream = upstream,
      delta = near - upstream
    )
  }

  delta_df <- lapply(rownames(logexpr_filtered), function(gene) {
    b1 <- mean_window(gene, b1_cells, b1_progress, logexpr_filtered)
    b2 <- mean_window(gene, b2_cells, b2_progress, logexpr_filtered)

    data.frame(
      gene = gene,
      delta_b1 = b1["delta"],
      delta_b2 = b2["delta"],
      near_b1 = b1["near"],
      near_b2 = b2["near"]
    )
  }) %>%
    bind_rows()

  exclusive_b1_shape <- delta_df %>%
    filter(
      is.finite(delta_b1),
      is.finite(delta_b2),
      delta_b1 > 0.30,
      delta_b2 < 0.10,
      (delta_b1 - delta_b2) > 0.25
    ) %>%
    arrange(desc(delta_b1 - delta_b2))

  write.csv(
    exclusive_b1_shape,
    file.path(convergence_dir, "exclusive_branch1_shape_genes.csv"),
    row.names = FALSE
  )

  exclusive_high_confidence <- inner_join(
    exclusive_b1 %>% select(gene, rho1, q1, rho2, q2),
    exclusive_b1_shape %>% select(gene, delta_b1, delta_b2),
    by = "gene"
  ) %>%
    arrange(desc((rho1 - rho2) + (delta_b1 - delta_b2)))

  write.csv(
    exclusive_high_confidence,
    file.path(convergence_dir, "exclusive_branch1_high_confidence_genes.csv"),
    row.names = FALSE
  )

  plot_gene_trend_pair <- function(gene, logexpr_mat = logexpr_filtered) {
    if (!gene %in% rownames(logexpr_mat)) return(NULL)

    df1 <- data.frame(
      branch = "b1",
      progress = b1_progress[b1_cells],
      expression = as.numeric(logexpr_mat[gene, b1_cells])
    )

    df2 <- data.frame(
      branch = "b2",
      progress = b2_progress[b2_cells],
      expression = as.numeric(logexpr_mat[gene, b2_cells])
    )

    trend_df <- bind_rows(df1, df2) %>%
      filter(is.finite(progress), is.finite(expression))

    rho1 <- res1$rho[match(gene, res1$gene)]
    rho2 <- res2$rho[match(gene, res2$gene)]

    ggplot(trend_df, aes(progress, expression)) +
      geom_point(alpha = 0.35, size = 0.6) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, span = 0.6) +
      facet_wrap(~branch, ncol = 2, scales = "fixed") +
      labs(
        title = paste0(gene, " — trend toward convergence"),
        subtitle = sprintf("Spearman rho: b1 = %.2f | b2 = %.2f", rho1, rho2),
        x = "Branch progress toward convergence",
        y = "log1p(CPM)"
      ) +
      theme_classic(base_size = 11)
  }

  trend_dir <- make_dir(file.path(convergence_dir, "gene_trends"))

  genes_to_plot <- head(convergent_ranked$gene, top_trend_genes_to_plot)

  for (i in seq_along(genes_to_plot)) {
    gene <- genes_to_plot[i]
    p <- plot_gene_trend_pair(gene)

    if (!is.null(p)) {
      save_png(
        p,
        file.path(trend_dir, sprintf("%02d_%s_trend_b1_b2.png", i, gene)),
        width = 9,
        height = 4
      )
    }
  }
}

# ------------------------------------------------------------------------------
# 7. COMBINED MYELOID + DC TRAJECTORY OVERLAY
# ------------------------------------------------------------------------------

if (run_combined_overlay) {
  combined_dir <- make_dir(file.path(output_dir, "combined_trajectory_overlay"))

  build_graph_for_clusters <- function(cds, target_clusters, origin_clusters, label) {
    trajectory <- learn_subset_trajectory(
      cds = cds,
      clusters = target_clusters,
      root_clusters = unname(origin_clusters),
      minimal_branch_len = minimal_branch_len
    )

    cds_sub <- trajectory$cds
    segments <- get_graph_segments(cds_sub, label)

    node_coordinates <- get_node_coordinates(cds_sub)

    root_points <- node_coordinates %>%
      filter(name %in% unname(trajectory$root_nodes)) %>%
      mutate(
        label = paste0(
          label,
          " origin cl ",
          names(origin_clusters)
        )
      )

    pseudotime_values <- monocle3::pseudotime(cds_sub)
    pseudotime_scaled <- scale01(pseudotime_values)

    pseudotime_df <- tibble(
      cell_id = names(pseudotime_scaled),
      pseudotime = as.numeric(pseudotime_scaled),
      subset = label
    )

    list(
      cds = cds_sub,
      segments = segments,
      root_points = root_points,
      pseudotime = pseudotime_df
    )
  }

  graph_myeloid <- build_graph_for_clusters(
    cds,
    myeloid_clusters,
    myeloid_origin_clusters,
    "Myeloid"
  )

  graph_dc <- build_graph_for_clusters(
    cds,
    dc_clusters,
    dc_origin_clusters,
    "DC"
  )

  umap_df <- as.data.frame(reducedDims(cds)$UMAP)
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  umap_df$cell_id <- colnames(cds)
  umap_df$cluster <- factor(colData(cds)$seurat_clusters)

  all_segments <- bind_rows(graph_myeloid$segments, graph_dc$segments)
  all_roots <- bind_rows(graph_myeloid$root_points, graph_dc$root_points)

  pseudotime_df <- bind_rows(graph_myeloid$pseudotime, graph_dc$pseudotime)
  pseudotime_map <- setNames(pseudotime_df$pseudotime, pseudotime_df$cell_id)
  umap_df$pseudotime <- unname(pseudotime_map[umap_df$cell_id])

  p_combined_clusters <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cluster), size = 0.25, alpha = 0.45) +
    geom_segment(
      data = subset(all_segments, graph == "Myeloid"),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "#1f77b4",
      alpha = 0.9
    ) +
    geom_segment(
      data = subset(all_segments, graph == "DC"),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "#d62728",
      alpha = 0.9
    ) +
    geom_point(
      data = subset(all_roots, grepl("^Myeloid", label)),
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.5,
      stroke = 0.6,
      fill = "#1f77b4",
      color = "black"
    ) +
    geom_point(
      data = subset(all_roots, grepl("^DC", label)),
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.5,
      stroke = 0.6,
      fill = "#d62728",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = all_roots,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 3,
      seed = 1
    ) +
    labs(
      title = "Combined trajectories with four origins",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 12)

  p_combined_pseudotime <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = pseudotime), size = 0.35, alpha = 0.7) +
    scale_color_viridis_c(
      option = "plasma",
      na.value = "grey85",
      name = "Scaled pseudotime"
    ) +
    geom_segment(
      data = subset(all_segments, graph == "Myeloid"),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "#1f77b4",
      alpha = 0.9
    ) +
    geom_segment(
      data = subset(all_segments, graph == "DC"),
      aes(x = x, y = y, xend = xend, yend = yend),
      inherit.aes = FALSE,
      linewidth = 0.6,
      color = "#d62728",
      alpha = 0.9
    ) +
    geom_point(
      data = all_roots,
      aes(x = x, y = y),
      inherit.aes = FALSE,
      shape = 21,
      size = 3.5,
      stroke = 0.6,
      fill = "white",
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = all_roots,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 3,
      seed = 1
    ) +
    labs(
      title = "Combined trajectories with pseudotime overlay",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 12)

  save_png(
    p_combined_clusters,
    file.path(combined_dir, "combined_four_origins_clusters.png"),
    width = 8,
    height = 6
  )

  save_png(
    p_combined_pseudotime,
    file.path(combined_dir, "combined_four_origins_pseudotime.png"),
    width = 8,
    height = 6
  )
}

# ------------------------------------------------------------------------------
# 8. INCOMING SIGNATURE SUMMARY
# ------------------------------------------------------------------------------

if (run_incoming_signature_summary && exists("convergent_ranked") && nrow(convergent_ranked) > 0) {
  signature_dir <- make_dir(file.path(output_dir, "incoming_signature_summary"))

  sig_genes <- convergent_ranked %>%
    filter(
      q1 < spearman_q_threshold,
      q2 < spearman_q_threshold,
      rho1 > spearman_rho_threshold,
      rho2 > spearman_rho_threshold
    ) %>%
    arrange(desc(rho1 + rho2)) %>%
    pull(gene) %>%
    unique()

  sig_genes <- intersect(sig_genes, rownames(logexpr_filtered))

  if (length(sig_genes) > 3) {
    t_grid <- seq(0, 1, length.out = 101)

    safe_loess <- function(x, y) {
      x <- as.numeric(x)
      y <- as.numeric(y)

      keep <- is.finite(x) & is.finite(y)
      x <- x[keep]
      y <- y[keep]

      if (length(unique(x)) < 6 || sum(y > 0) < 10) {
        return(rep(NA_real_, length(t_grid)))
      }

      tryCatch(
        predict(
          loess(y ~ x, span = 0.6, degree = 2, surface = "direct"),
          newdata = t_grid
        ),
        error = function(e) rep(NA_real_, length(t_grid))
      )
    }

    smooth_one_gene <- function(gene) {
      y1 <- as.numeric(logexpr_filtered[gene, b1_cells])
      x1 <- b1_progress[b1_cells]

      y2 <- as.numeric(logexpr_filtered[gene, b2_cells])
      x2 <- b2_progress[b2_cells]

      f1 <- safe_loess(x1, y1)
      f2 <- safe_loess(x2, y2)

      if (all(!is.finite(f1)) && all(!is.finite(f2))) {
        return(NULL)
      }

      all_fit <- c(f1, f2)
      fit_mean <- mean(all_fit, na.rm = TRUE)
      fit_sd <- sd(all_fit, na.rm = TRUE)
      if (fit_sd == 0) fit_sd <- 1

      tibble(
        gene = gene,
        t = rep(t_grid, 2),
        z = c((f1 - fit_mean) / fit_sd, (f2 - fit_mean) / fit_sd),
        branch = rep(c("b1", "b2"), each = length(t_grid))
      )
    }

    trend_df <- purrr::map(sig_genes, smooth_one_gene) %>%
      bind_rows() %>%
      filter(is.finite(z))

    summary_df <- trend_df %>%
      group_by(branch, t) %>%
      summarise(
        median = median(z, na.rm = TRUE),
        q25 = quantile(z, 0.25, na.rm = TRUE),
        q75 = quantile(z, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    p_signature <- ggplot() +
      annotate(
        "rect",
        xmin = 0.8,
        xmax = 1,
        ymin = -Inf,
        ymax = Inf,
        alpha = 0.06,
        fill = "grey50"
      ) +
      geom_line(
        data = trend_df,
        aes(t, z, group = interaction(gene, branch), color = branch),
        linewidth = 0.3,
        alpha = 0.12
      ) +
      geom_ribbon(
        data = summary_df,
        aes(t, ymin = q25, ymax = q75, fill = branch),
        alpha = 0.20
      ) +
      geom_line(
        data = summary_df,
        aes(t, median, color = branch),
        linewidth = 1.2
      ) +
      scale_color_manual(values = c(b1 = "#1f77b4", b2 = "#d62728"), guide = "none") +
      scale_fill_manual(values = c(b1 = "#1f77b4", b2 = "#d62728"), guide = "none") +
      labs(
        title = "Incoming signature definition",
        subtitle = "Thin lines = individual genes; thick line = branch median ± IQR",
        x = "Branch progress toward convergence node",
        y = "Standardized expression"
      ) +
      theme_classic(base_size = 11)

    seurat_object <- AddModuleScore(
      seurat_object,
      features = list(sig_genes),
      name = "IncomingSig"
    )

    score_name <- grep("^IncomingSig", colnames(seurat_object@meta.data), value = TRUE)[1]

    p_clusters <- DotPlot(
      seurat_object,
      features = score_name,
      group.by = cluster_column
    ) +
      coord_flip() +
      labs(
        title = "Incoming signature across clusters",
        x = "Cluster",
        y = NULL
      )

    fig <- p_signature / p_clusters + plot_layout(heights = c(1, 0.8))

    save_png(
      fig,
      file.path(signature_dir, "incoming_signature_summary.png"),
      width = 10,
      height = 9
    )

    save_png(
      p_signature,
      file.path(signature_dir, "incoming_signature_trends.png"),
      width = 8,
      height = 5
    )

    write.csv(
      data.frame(gene = sig_genes),
      file.path(signature_dir, "incoming_signature_genes.csv"),
      row.names = FALSE
    )
  } else {
    message("Skipping incoming signature summary; fewer than four signature genes passed filters.")
  }
}

# ------------------------------------------------------------------------------
# 9. DENSITY-WEIGHTED TRAJECTORY GRAPH
# ------------------------------------------------------------------------------

if (run_density_graph) {
  density_dir <- make_dir(file.path(output_dir, "density_weighted_graph"))

  plot_cells_with_density_graph <- function(
    cds,
    color_cells_by = "seurat_clusters",
    radius = 0.25
  ) {
    base_plot <- plot_cells(
      cds,
      color_cells_by = color_cells_by,
      label_cell_groups = TRUE,
      label_leaves = FALSE,
      label_branch_points = FALSE,
      label_roots = FALSE,
      show_trajectory_graph = FALSE
    )

    graph <- principal_graph(cds)[["UMAP"]]
    node_coordinates <- get_node_coordinates(cds)

    edges <- as.data.frame(igraph::as_edgelist(graph))
    colnames(edges) <- c("from", "to")

    edges <- edges %>%
      left_join(node_coordinates, by = c("from" = "name")) %>%
      rename(x1 = x, y1 = y) %>%
      left_join(node_coordinates, by = c("to" = "name")) %>%
      rename(x2 = x, y2 = y)

    cell_df <- as.data.frame(reducedDims(cds)$UMAP)
    colnames(cell_df)[1:2] <- c("x", "y")

    point_to_segment_distance <- function(px, py, x1, y1, x2, y2) {
      vx <- x2 - x1
      vy <- y2 - y1
      wx <- px - x1
      wy <- py - y1

      c1 <- vx * wx + vy * wy
      if (c1 <= 0) return(sqrt((px - x1)^2 + (py - y1)^2))

      c2 <- vx * vx + vy * vy
      if (c2 <= c1) return(sqrt((px - x2)^2 + (py - y2)^2))

      b <- c1 / c2
      bx <- x1 + b * vx
      by <- y1 + b * vy

      sqrt((px - bx)^2 + (py - by)^2)
    }

    edges$length <- sqrt((edges$x2 - edges$x1)^2 + (edges$y2 - edges$y1)^2)

    edges$density <- vapply(seq_len(nrow(edges)), function(i) {
      edge <- edges[i, ]

      distances_to_edge <- mapply(
        point_to_segment_distance,
        px = cell_df$x,
        py = cell_df$y,
        MoreArgs = list(
          x1 = edge$x1,
          y1 = edge$y1,
          x2 = edge$x2,
          y2 = edge$y2
        )
      )

      sum(distances_to_edge <= radius) / pmax(edge$length, 1e-8)
    }, numeric(1))

    edges$line_width <- scales::rescale(edges$density, to = c(0.4, 2.5))

    base_plot +
      geom_segment(
        data = edges,
        aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = line_width),
        inherit.aes = FALSE,
        color = "black",
        alpha = 0.9,
        lineend = "round"
      ) +
      scale_linewidth_identity()
  }

  p_density_dc <- plot_cells_with_density_graph(
    dc$cds,
    color_cells_by = "seurat_clusters",
    radius = density_radius
  ) +
    ggtitle("DC trajectory with density-weighted graph")

  save_png(
    p_density_dc,
    file.path(density_dir, "dc_density_weighted_graph.png"),
    width = 8,
    height = 6
  )
}

# ------------------------------------------------------------------------------
# 10. FINISH
# ------------------------------------------------------------------------------

message("Done. Outputs saved to: ", output_dir)