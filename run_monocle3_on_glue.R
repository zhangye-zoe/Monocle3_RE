#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(monocle3)
  library(Matrix)
  library(data.table)
  library(jsonlite)
  library(SingleCellExperiment)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)

parse_bool <- function(x) {
  if (is.null(x)) return(NULL)
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (x %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(sprintf("Cannot parse boolean value: %s", x))
}

parse_args <- function(args) {
  res <- list(
    input_dir = NULL,
    outdir = NULL,
    label_key = "cell_type",
    stage_order_json = NULL,
    root_labels = NULL,
    num_dim = 50,
    cluster_resolution = 1e-2,
    use_partition = TRUE,
    random_seed = 1234
  )
  
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1]] else NULL
    
    if (key == "--input-dir") res$input_dir <- val
    if (key == "--outdir") res$outdir <- val
    if (key == "--label-key") res$label_key <- val
    if (key == "--stage-order-json") res$stage_order_json <- val
    if (key == "--root-labels") res$root_labels <- val
    if (key == "--num-dim") res$num_dim <- as.integer(val)
    if (key == "--cluster-resolution") res$cluster_resolution <- as.numeric(val)
    if (key == "--use-partition") res$use_partition <- parse_bool(val)
    if (key == "--random-seed") res$random_seed <- as.integer(val)
    
    i <- i + 2
  }
  
  if (is.null(res$input_dir) || is.null(res$outdir)) {
    stop("Must provide --input-dir and --outdir")
  }
  res
}

# ---- helper: safely fetch closest vertex vector ----
get_closest_vertex_raw <- function(cds, reduction_method = "UMAP") {
  aux <- cds@principal_graph_aux[[reduction_method]]
  if (is.null(aux)) {
    stop(sprintf("cds@principal_graph_aux[['%s']] is NULL", reduction_method))
  }
  if (is.null(aux$pr_graph_cell_proj_closest_vertex)) {
    stop(sprintf("pr_graph_cell_proj_closest_vertex not found under principal_graph_aux[['%s']]", reduction_method))
  }
  
  closest_vertex <- aux$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), , drop = FALSE])
  closest_vertex <- as.character(closest_vertex[, 1])
  names(closest_vertex) <- colnames(cds)
  closest_vertex
}

# ---- helper: map closest vertex ids to actual principal node names ----
map_closest_vertex_to_pr_nodes <- function(cds, closest_vertex_raw, reduction_method = "UMAP") {
  g <- principal_graph(cds)[[reduction_method]]
  if (is.null(g)) {
    stop(sprintf("principal_graph(cds)[['%s']] is NULL", reduction_method))
  }
  
  valid_pr_nodes <- igraph::V(g)$name
  if (is.null(valid_pr_nodes) || length(valid_pr_nodes) == 0) {
    stop("No valid principal graph node names found.")
  }
  
  mapped <- closest_vertex_raw
  
  # Case 1: already principal node names like Y_1
  already_valid <- mapped %in% valid_pr_nodes
  
  # Case 2: numeric/internal ids -> map by vertex sequence order
  suppressWarnings(num_idx <- as.integer(mapped))
  can_map_num <- !is.na(num_idx) & num_idx >= 1 & num_idx <= length(valid_pr_nodes)
  
  mapped[can_map_num] <- valid_pr_nodes[num_idx[can_map_num]]
  
  # Re-check after numeric mapping
  final_valid <- mapped %in% valid_pr_nodes
  
  cat("Closest vertex raw: length =", length(closest_vertex_raw), "\n")
  cat("Closest vertex raw NA count:", sum(is.na(closest_vertex_raw)), "\n")
  cat("Already valid principal-node names:", sum(already_valid, na.rm = TRUE), "\n")
  cat("Mapped numeric/internal ids:", sum(can_map_num, na.rm = TRUE), "\n")
  cat("Final valid principal-node names:", sum(final_valid, na.rm = TRUE), "\n")
  
  if (sum(final_valid, na.rm = TRUE) == 0) {
    cat("First few raw closest vertices:\n")
    print(head(closest_vertex_raw))
    cat("First few valid principal graph nodes:\n")
    print(head(valid_pr_nodes))
    stop("Failed to map closest vertices to valid principal graph node names.")
  }
  
  mapped[!final_valid] <- NA_character_
  mapped
}

get_roots_from_labels <- function(cds, label_key, stage_order_map, explicit_root_labels = NULL, reduction_method = "UMAP") {
  if (!(label_key %in% colnames(colData(cds)))) {
    stop(sprintf("label_key '%s' not found in colData(cds)", label_key))
  }
  
  labels <- as.character(colData(cds)[[label_key]])
  cat("Unique labels in cds:\n")
  print(sort(unique(labels)))
  
  closest_vertex_raw <- get_closest_vertex_raw(cds, reduction_method = reduction_method)
  closest_vertex <- map_closest_vertex_to_pr_nodes(
    cds,
    closest_vertex_raw,
    reduction_method = reduction_method
  )
  
  cat("Mapped closest_vertex length:", length(closest_vertex), "\n")
  cat("Mapped closest_vertex NA count:", sum(is.na(closest_vertex)), "\n")
  
  if (!is.null(explicit_root_labels)) {
    root_label_vec <- strsplit(explicit_root_labels, ",")[[1]]
    root_label_vec <- trimws(root_label_vec)
    cat("Requested root labels:\n")
    print(root_label_vec)
    
    cell_ids <- colnames(cds)[labels %in% root_label_vec]
    cat("Matched root cells:", length(cell_ids), "\n")
    if (length(cell_ids) > 0) {
      cat("First few matched root cells:\n")
      print(head(cell_ids))
      cat("Their labels:\n")
      print(table(labels[labels %in% root_label_vec]))
    }
  } else {
    stage_scores <- as.numeric(stage_order_map[labels])
    cat("Stage score summary:\n")
    print(summary(stage_scores))
    
    earliest_stage <- min(stage_scores[is.finite(stage_scores)], na.rm = TRUE)
    cat("Earliest stage:", earliest_stage, "\n")
    
    cell_ids <- colnames(cds)[is.finite(stage_scores) & stage_scores == earliest_stage]
    cat("Matched earliest-stage cells:", length(cell_ids), "\n")
  }
  
  if (length(cell_ids) == 0) {
    stop("No root cells found from stage order / root labels.")
  }
  
  root_vertices <- closest_vertex[cell_ids]
  cat("Root vertex NA count after mapping:", sum(is.na(root_vertices)), "\n")
  cat("First few mapped root vertices:\n")
  print(head(root_vertices))
  
  root_vertices <- root_vertices[!is.na(root_vertices)]
  
  if (length(root_vertices) == 0) {
    stop("Matched root cells exist, but none could be mapped to valid principal graph nodes.")
  }
  
  tab <- sort(table(root_vertices), decreasing = TRUE)
  cat("Mapped root vertex table:\n")
  print(tab)
  
  g <- principal_graph(cds)[[reduction_method]]
  valid_pr_nodes <- igraph::V(g)$name
  
  cat("First few valid principal graph nodes:\n")
  print(head(valid_pr_nodes))
  
  root_pr_nodes <- intersect(names(tab), valid_pr_nodes)
  
  cat("Filtered root_pr_nodes:\n")
  print(root_pr_nodes)
  
  if (length(root_pr_nodes) == 0) {
    stop("Matched root cells exist, but projected root vertices are not valid principal graph node names after mapping.")
  }
  
  root_pr_nodes
}

main <- function() {
  opt <- parse_args(args)
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(opt$random_seed)
  
  expr <- readMM(file.path(opt$input_dir, "expr_gene_by_cell.mtx"))
  cell_md <- fread(file.path(opt$input_dir, "cell_metadata.csv"), data.table = FALSE)
  gene_md <- fread(file.path(opt$input_dir, "gene_metadata.csv"), data.table = FALSE)
  umap_df <- fread(file.path(opt$input_dir, "umap_embedding.csv"), data.table = FALSE)
  latent_df <- fread(file.path(opt$input_dir, "latent_embedding.csv"), data.table = FALSE)
  
  rownames(cell_md) <- cell_md$cell_id
  rownames(gene_md) <- gene_md$gene_short_name
  rownames(umap_df) <- umap_df[[1]]
  rownames(latent_df) <- latent_df[[1]]
  
  umap_df <- umap_df[, c("UMAP_1", "UMAP_2"), drop = FALSE]
  latent_df <- latent_df[, setdiff(colnames(latent_df), colnames(latent_df)[1]), drop = FALSE]
  
  if (nrow(gene_md) != nrow(expr)) {
    stop(sprintf("gene metadata rows (%d) != matrix rows (%d)", nrow(gene_md), nrow(expr)))
  }
  if (nrow(cell_md) != ncol(expr)) {
    stop(sprintf("cell metadata rows (%d) != matrix cols (%d)", nrow(cell_md), ncol(expr)))
  }
  
  colnames(expr) <- rownames(cell_md)
  rownames(expr) <- rownames(gene_md)
  
  cds <- new_cell_data_set(
    expression_data = expr,
    cell_metadata = cell_md,
    gene_metadata = gene_md
  )
  
  # 这里主要是为了满足 monocle3 对对象结构的要求
  cds <- preprocess_cds(cds, num_dim = min(opt$num_dim, ncol(expr) - 1))
  
  # 用你自己的 latent / UMAP 覆盖
  reducedDims(cds)$PCA <- as.matrix(latent_df[colnames(cds), , drop = FALSE])
  reducedDims(cds)$UMAP <- as.matrix(umap_df[colnames(cds), , drop = FALSE])
  
  cds <- cluster_cells(
    cds,
    reduction_method = "UMAP",
    resolution = opt$cluster_resolution
  )
  
  cds <- learn_graph(
    cds,
    use_partition = opt$use_partition,
    close_loop = FALSE
  )
  
  root_pr_nodes <- NULL
  if (!is.null(opt$stage_order_json)) {
    tmp_map <- unlist(fromJSON(opt$stage_order_json), use.names = TRUE)
    stage_order_map <- as.numeric(tmp_map)
    names(stage_order_map) <- names(tmp_map)
    
    root_pr_nodes <- get_roots_from_labels(
      cds = cds,
      label_key = opt$label_key,
      stage_order_map = stage_order_map,
      explicit_root_labels = opt$root_labels,
      reduction_method = "UMAP"
    )
  }
  
  if (!is.null(root_pr_nodes) && length(root_pr_nodes) > 0) {
    cat("Using root principal nodes for ordering:\n")
    print(root_pr_nodes)
    cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)
  } else {
    cds <- order_cells(cds)
  }
  
  pseudotime_vec <- monocle3::pseudotime(cds)
  
  # safer than unqualified partitions()/clusters()
  cd <- as.data.frame(SummarizedExperiment::colData(cds))
  
  part_col <- grep("^partition", colnames(cd), value = TRUE)
  clus_col <- grep("^cluster", colnames(cd), value = TRUE)
  
  if (length(part_col) > 0) {
    part <- cd[[part_col[1]]]
    names(part) <- rownames(cd)
  } else {
    part <- rep(NA_character_, ncol(cds))
    names(part) <- colnames(cds)
  }
  
  if (length(clus_col) > 0) {
    clus <- cd[[clus_col[1]]]
    names(clus) <- rownames(cd)
  } else {
    clus <- rep(NA_character_, ncol(cds))
    names(clus) <- colnames(cds)
  }
  
  closest_vertex_raw <- get_closest_vertex_raw(cds, reduction_method = "UMAP")
  closest_vertex <- map_closest_vertex_to_pr_nodes(cds, closest_vertex_raw, reduction_method = "UMAP")
  
  g <- principal_graph(cds)[["UMAP"]]
  edge_df <- as_data_frame(g, what = "edges")
  vertex_df <- data.frame(
    principal_node = V(g)$name,
    degree = degree(g),
    stringsAsFactors = FALSE
  )
  
  leaf_nodes <- vertex_df$principal_node[vertex_df$degree == 1]
  if (!is.null(root_pr_nodes)) {
    leaf_nodes <- setdiff(leaf_nodes, root_pr_nodes)
  }
  
  cell_out <- data.frame(
    cell_id = colnames(cds),
    pseudotime = as.numeric(pseudotime_vec[colnames(cds)]),
    partition = as.character(part[colnames(cds)]),
    cluster = as.character(clus[colnames(cds)]),
    closest_vertex = as.character(closest_vertex[colnames(cds)]),
    closest_vertex_raw = as.character(closest_vertex_raw[colnames(cds)]),
    stringsAsFactors = FALSE
  )
  
  for (nm in colnames(colData(cds))) {
    cell_out[[nm]] <- as.vector(colData(cds)[[nm]])
  }
  
  fwrite(cell_out, file.path(opt$outdir, "monocle3_cells.csv"))
  fwrite(edge_df, file.path(opt$outdir, "principal_graph_edges.csv"))
  fwrite(vertex_df, file.path(opt$outdir, "principal_graph_vertices.csv"))
  writeLines(vertex_df$principal_node, file.path(opt$outdir, "principal_nodes.txt"))
  writeLines(leaf_nodes, file.path(opt$outdir, "leaf_nodes.txt"))
  if (!is.null(root_pr_nodes)) {
    writeLines(root_pr_nodes, file.path(opt$outdir, "root_nodes.txt"))
  }
  
  saveRDS(cds, file.path(opt$outdir, "monocle3_cds.rds"))
  
  summary <- list(
    n_cells = ncol(cds),
    n_genes = nrow(cds),
    label_key = opt$label_key,
    used_root_nodes = if (is.null(root_pr_nodes)) list() else as.list(root_pr_nodes),
    n_leaf_nodes = length(leaf_nodes),
    use_partition = opt$use_partition
  )
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(opt$outdir, "monocle3_summary.json"))
  cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE))
}

main()