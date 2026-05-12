#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(igraph)
  library(pROC)
  library(reticulate)
})

use_python("/opt/anaconda3/bin/python", required = TRUE)
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  res <- list(
    monocle_dir = NULL,
    input_h5ad = NULL,
    outdir = NULL,
    label_key = NULL,
    reference_time_key = NULL,
    stage_order_json = NULL,
    future_map_json = NULL,
    future_label_key = "future_label",
    temperature = 1.0
  )
  
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    val <- if (i < length(args)) args[[i + 1]] else NULL
    
    if (key == "--monocle-dir") res$monocle_dir <- val
    if (key == "--input-h5ad") res$input_h5ad <- val
    if (key == "--outdir") res$outdir <- val
    if (key == "--label-key") res$label_key <- val
    if (key == "--reference-time-key") res$reference_time_key <- val
    if (key == "--stage-order-json") res$stage_order_json <- val
    if (key == "--future-map-json") res$future_map_json <- val
    if (key == "--future-label-key") res$future_label_key <- val
    if (key == "--temperature") res$temperature <- as.numeric(val)
    
    i <- i + 2
  }
  
  if (is.null(res$monocle_dir) || is.null(res$input_h5ad) || is.null(res$outdir)) {
    stop("Must provide --monocle-dir, --input-h5ad, and --outdir")
  }
  res
}

load_json_maybe <- function(x) {
  if (is.null(x)) return(NULL)
  x <- as.character(x)
  if (file.exists(x)) {
    return(fromJSON(txt = x, simplifyVector = TRUE))
  }
  fromJSON(txt = x, simplifyVector = TRUE)
}

find_first_existing_key <- function(container_names, candidates, kind = "key") {
  for (k in candidates) {
    if (k %in% container_names) return(k)
  }
  stop(sprintf("Could not find any valid %s. Tried: %s", kind, paste(candidates, collapse = ", ")))
}

resolve_label_key <- function(df, preferred = NULL) {
  candidates <- c()
  if (!is.null(preferred)) candidates <- c(candidates, preferred)
  candidates <- c(candidates, "celltype", "cell_type", "annotation", "annot", "labels", "leiden")
  find_first_existing_key(colnames(df), candidates, kind = "label key")
}

normalize_future_labels <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- trimws(x)
  x[x == ""] <- NA_character_
  x[tolower(x) %in% c("nan", "none", "null")] <- NA_character_
  x
}

safe_auc_binary <- function(y_true, y_score) {
  if (length(unique(y_true)) < 2) return(NA_real_)
  roc_obj <- tryCatch(
    pROC::roc(response = y_true, predictor = y_score, quiet = TRUE, direction = "<"),
    error = function(e) NULL
  )
  if (is.null(roc_obj)) return(NA_real_)
  as.numeric(pROC::auc(roc_obj))
}

interclass_pairwise_accuracy <- function(pred, ref, weighted = TRUE) {
  pred <- as.numeric(pred)
  ref <- as.numeric(ref)
  
  valid <- is.finite(pred) & is.finite(ref)
  pred <- pred[valid]
  ref <- ref[valid]
  
  correct <- 0
  total <- 0
  n <- length(pred)
  
  if (n < 2) return(NA_real_)
  
  cat(sprintf("[T3] Pairwise accuracy: %d valid cells, ~%.0f pairs\n", n, n * (n - 1) / 2))
  flush.console()
  
  for (i in seq_len(n)) {
    if (i == n) next
    for (j in seq((i + 1), n)) {
      if (ref[i] == ref[j]) next
      w <- if (weighted) abs(ref[i] - ref[j]) else 1
      is_correct <- if (ref[i] < ref[j]) pred[i] < pred[j] else pred[j] < pred[i]
      correct <- correct + w * as.numeric(is_correct)
      total <- total + w
    }
    if (i %% 500 == 0 || i == n) {
      cat(sprintf("[T3] PairAcc progress: %d / %d cells processed\n", i, n))
      flush.console()
    }
  }
  if (total > 0) correct / total else NA_real_
}

interclass_pairwise_auc <- function(pred, ref, weighted = TRUE) {
  pred <- as.numeric(pred)
  ref <- as.numeric(ref)
  
  valid <- is.finite(pred) & is.finite(ref)
  pred <- pred[valid]
  ref <- ref[valid]
  
  classes <- sort(unique(ref))
  aucs <- c()
  weights <- c()
  
  if (length(classes) < 2) return(NA_real_)
  
  total_pairs <- length(classes) * (length(classes) - 1) / 2
  done_pairs <- 0
  
  cat(sprintf("[T3] Pairwise AUC: %d classes, %d class-pairs\n", length(classes), total_pairs))
  flush.console()
  
  for (i in seq_along(classes)) {
    if (i == length(classes)) next
    a <- classes[i]
    for (b in classes[(i + 1):length(classes)]) {
      mask <- (ref == a) | (ref == b)
      y <- as.integer(ref[mask] == b)
      s <- pred[mask]
      if (length(unique(y)) >= 2) {
        auc <- safe_auc_binary(y, s)
        w <- if (weighted) abs(b - a) else 1
        aucs <- c(aucs, auc)
        weights <- c(weights, w)
      }
      done_pairs <- done_pairs + 1
      if (done_pairs %% 10 == 0 || done_pairs == total_pairs) {
        cat(sprintf("[T3] PairAUC progress: %d / %d class-pairs\n", done_pairs, total_pairs))
        flush.console()
      }
    }
  }
  
  if (length(aucs) == 0) return(NA_real_)
  weighted.mean(aucs, weights, na.rm = TRUE)
}

compute_t3_metrics <- function(df, pred_time_key = "pseudotime", reference_time_key = NULL,
                               label_key = NULL, stage_order_map = NULL, eval_mask = NULL) {
  pred <- suppressWarnings(as.numeric(df[[pred_time_key]]))
  
  if (is.null(eval_mask)) {
    eval_mask <- rep(TRUE, nrow(df))
  } else {
    eval_mask <- as.logical(eval_mask)
  }
  
  if (!is.null(reference_time_key)) {
    ref <- suppressWarnings(as.numeric(df[[reference_time_key]]))
    reference_name <- reference_time_key
  } else {
    if (is.null(label_key)) label_key <- resolve_label_key(df)
    if (is.null(stage_order_map)) {
      stop("stage_order_map must be provided when reference_time_key is NULL.")
    }
    labels <- as.character(df[[label_key]])
    ref <- as.numeric(stage_order_map[labels])
    reference_name <- sprintf("ordinal(%s)", label_key)
  }
  
  valid <- eval_mask & is.finite(pred) & is.finite(ref)
  
  cat(sprintf("[T3] Valid cells for evaluation: %d / %d\n", sum(valid), nrow(df)))
  flush.console()
  
  if (sum(valid) < 3) {
    return(data.frame(
      PairAcc = NA_real_,
      PairAUC = NA_real_,
      n_eval = sum(valid),
      pred_time_key = pred_time_key,
      reference = reference_name,
      stringsAsFactors = FALSE
    ))
  }
  
  data.frame(
    PairAcc = interclass_pairwise_accuracy(pred[valid], ref[valid], weighted = TRUE),
    PairAUC = interclass_pairwise_auc(pred[valid], ref[valid], weighted = TRUE),
    n_eval = sum(valid),
    pred_time_key = pred_time_key,
    reference = reference_name,
    stringsAsFactors = FALSE
  )
}

build_leaf_class_map <- function(cell_df, leaf_nodes, lineage_key) {
  tmp <- copy(cell_df)
  tmp[[lineage_key]] <- normalize_future_labels(tmp[[lineage_key]])
  tmp <- tmp[tmp$closest_vertex %in% leaf_nodes & !is.na(tmp[[lineage_key]]), ]
  
  if (nrow(tmp) == 0) return(character(0))
  
  leaf_class_map <- tapply(
    tmp[[lineage_key]],
    INDEX = tmp$closest_vertex,
    FUN = function(x) {
      tb <- sort(table(as.character(x)), decreasing = TRUE)
      names(tb)[1]
    }
  )
  unlist(leaf_class_map)
}

graph_distance_dict <- function(g, source) {
  d <- distances(g, v = source, to = V(g), weights = E(g)$weight)
  out <- as.numeric(d[1, ])
  names(out) <- V(g)$name
  out
}

compute_monocle_future_probabilities <- function(cell_df, edge_df, root_nodes,
                                                 lineage_key = "future_label",
                                                 temperature = 1.0) {
  cat("[T4] Building graph...\n")
  flush.console()
  
  g <- graph_from_data_frame(edge_df[, c("from", "to")], directed = FALSE)
  E(g)$weight <- 1.0
  
  deg <- degree(g)
  leaf_nodes <- names(deg)[deg == 1 & !(names(deg) %in% root_nodes)]
  
  if (length(root_nodes) == 0) stop("No root nodes were provided; cannot orient the graph for T4.")
  if (length(leaf_nodes) == 0) stop("No leaf nodes found in principal graph.")
  
  cat(sprintf("[T4] Graph nodes: %d, edges: %d, root nodes: %d, leaf nodes: %d\n",
              vcount(g), ecount(g), length(root_nodes), length(leaf_nodes)))
  flush.console()
  
  cat("[T4] Computing root distances...\n")
  flush.console()
  root_dist <- rep(Inf, length(V(g)))
  names(root_dist) <- V(g)$name
  pb_root <- txtProgressBar(min = 0, max = length(V(g)$name), style = 3)
  idx_root <- 0
  for (n in V(g)$name) {
    vals <- c()
    for (r in root_nodes) {
      d <- tryCatch(distances(g, v = r, to = n, weights = E(g)$weight)[1, 1], error = function(e) Inf)
      if (is.finite(d)) vals <- c(vals, d)
    }
    root_dist[n] <- if (length(vals) > 0) min(vals) else Inf
    idx_root <- idx_root + 1
    setTxtProgressBar(pb_root, idx_root)
  }
  close(pb_root)
  cat("\n[T4] Finished root distances.\n")
  flush.console()
  
  leaf_class_map <- build_leaf_class_map(cell_df, leaf_nodes, lineage_key)
  valid_leaf_nodes <- leaf_nodes[leaf_nodes %in% names(leaf_class_map)]
  if (length(valid_leaf_nodes) == 0) {
    stop("None of the graph leaves could be annotated from future labels. Check lineage_key / future_map / cell labels.")
  }
  
  class_order <- sort(unique(unname(leaf_class_map[valid_leaf_nodes])))
  cat(sprintf("[T4] Valid annotated leaf nodes: %d, classes: %s\n",
              length(valid_leaf_nodes), paste(class_order, collapse = ", ")))
  flush.console()
  
  cat("[T4] Caching leaf-to-node graph distances...\n")
  flush.console()
  pb_leaf <- txtProgressBar(min = 0, max = length(valid_leaf_nodes), style = 3)
  dist_cache <- vector("list", length(valid_leaf_nodes))
  names(dist_cache) <- valid_leaf_nodes
  for (ii in seq_along(valid_leaf_nodes)) {
    leaf <- valid_leaf_nodes[ii]
    dist_cache[[leaf]] <- graph_distance_dict(g, leaf)
    setTxtProgressBar(pb_leaf, ii)
  }
  close(pb_leaf)
  cat("\n[T4] Finished leaf distance cache.\n")
  flush.console()
  
  scores <- matrix(0, nrow = nrow(cell_df), ncol = length(class_order))
  colnames(scores) <- paste0("prob_", class_order)
  
  per_leaf_scores <- list(
    cell_id = character(0),
    leaf_node = character(0),
    leaf_class = character(0),
    score = numeric(0)
  )
  
  cat(sprintf("[T4] Computing future probabilities for %d cells...\n", nrow(cell_df)))
  flush.console()
  pb_cell <- txtProgressBar(min = 0, max = nrow(cell_df), style = 3)
  
  for (i in seq_len(nrow(cell_df))) {
    v <- as.character(cell_df$closest_vertex[i])
    if (v %in% V(g)$name) {
      downstream <- valid_leaf_nodes[
        is.finite(root_dist[v]) &
          is.finite(root_dist[valid_leaf_nodes]) &
          root_dist[valid_leaf_nodes] > root_dist[v]
      ]
      candidates <- if (length(downstream) > 0) downstream else valid_leaf_nodes
      
      leaf_scores <- c()
      for (leaf in candidates) {
        d <- dist_cache[[leaf]][v]
        if (!is.finite(d)) next
        leaf_scores[leaf] <- exp(-d / max(temperature, 1e-8))
      }
      
      if (length(leaf_scores) > 0) {
        class_score <- setNames(rep(0, length(class_order)), class_order)
        for (leaf in names(leaf_scores)) {
          cls <- leaf_class_map[[leaf]]
          class_score[[cls]] <- class_score[[cls]] + leaf_scores[[leaf]]
        }
        
        vec <- as.numeric(class_score[class_order])
        if (sum(vec) > 0) vec <- vec / sum(vec)
        scores[i, ] <- vec
        
        for (leaf in names(leaf_scores)) {
          per_leaf_scores$cell_id <- c(per_leaf_scores$cell_id, as.character(cell_df$cell_id[i]))
          per_leaf_scores$leaf_node <- c(per_leaf_scores$leaf_node, leaf)
          per_leaf_scores$leaf_class <- c(per_leaf_scores$leaf_class, leaf_class_map[[leaf]])
          per_leaf_scores$score <- c(per_leaf_scores$score, as.numeric(leaf_scores[[leaf]]))
        }
      }
    }
    setTxtProgressBar(pb_cell, i)
  }
  close(pb_cell)
  cat("\n[T4] Finished future probability computation.\n")
  flush.console()
  
  proba_df <- data.frame(cell_id = cell_df$cell_id, scores, check.names = FALSE, stringsAsFactors = FALSE)
  
  leaf_meta <- data.frame(
    leaf_node = valid_leaf_nodes,
    leaf_class = unname(leaf_class_map[valid_leaf_nodes]),
    root_distance = as.numeric(root_dist[valid_leaf_nodes]),
    stringsAsFactors = FALSE
  )
  
  list(
    proba_df = proba_df,
    class_order = class_order,
    leaf_meta = leaf_meta,
    per_leaf_scores = as.data.frame(per_leaf_scores, stringsAsFactors = FALSE)
  )
}

compute_t4_metrics <- function(cell_df, proba_df, class_order, lineage_key = "future_label", eval_mask = NULL) {
  y_true_str <- normalize_future_labels(cell_df[[lineage_key]])
  
  if (is.null(eval_mask)) {
    eval_mask <- rep(TRUE, nrow(cell_df))
  } else {
    eval_mask <- as.logical(eval_mask)
  }
  
  valid_classes <- y_true_str %in% class_order
  valid <- eval_mask & valid_classes & !is.na(y_true_str)
  
  cat(sprintf("[T4] Valid cells for evaluation: %d / %d\n", sum(valid), nrow(cell_df)))
  flush.console()
  
  if (sum(valid) < 3) {
    return(data.frame(
      F1 = NA_real_,
      AUROC = NA_real_,
      n_eval = sum(valid),
      class_order = paste(class_order, collapse = ","),
      stringsAsFactors = FALSE
    ))
  }
  
  y_true <- match(y_true_str[valid], class_order) - 1L
  score_cols <- paste0("prob_", class_order)
  y_score <- as.matrix(proba_df[valid, score_cols, drop = FALSE])
  y_pred <- max.col(y_score, ties.method = "first") - 1L
  
  # macro F1
  f1_per_class <- c()
  for (k in seq_along(class_order) - 1L) {
    tp <- sum(y_true == k & y_pred == k)
    fp <- sum(y_true != k & y_pred == k)
    fn <- sum(y_true == k & y_pred != k)
    precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    recall <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    f1 <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      NA_real_
    }
    f1_per_class <- c(f1_per_class, f1)
  }
  f1_macro <- mean(f1_per_class, na.rm = TRUE)
  
  # AUROC
  if (length(class_order) == 2) {
    auc <- safe_auc_binary(y_true, y_score[, 2])
  } else {
    aucs <- c()
    for (k in seq_along(class_order) - 1L) {
      y_bin <- as.integer(y_true == k)
      aucs <- c(aucs, safe_auc_binary(y_bin, y_score[, k + 1]))
    }
    auc <- mean(aucs, na.rm = TRUE)
  }
  
  data.frame(
    F1 = as.numeric(f1_macro),
    AUROC = as.numeric(auc),
    n_eval = sum(valid),
    class_order = paste(class_order, collapse = ","),
    stringsAsFactors = FALSE
  )
}

build_parser <- function() {
  parse_args(args)
}

main <- function() {
  opt <- build_parser()
  dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
  
  cat("[INFO] Loading h5ad...\n")
  flush.console()
  ad <- reticulate::import("anndata", convert = FALSE)
  adata <- ad$read_h5ad(opt$input_h5ad)
  obs <- reticulate::py_to_r(adata$obs)
  obs <- as.data.frame(obs)
  obs_names <- reticulate::py_to_r(adata$obs_names$to_list())
  obs$cell_id <- as.character(obs_names)
  cat(sprintf("[INFO] h5ad loaded: %d cells\n", nrow(obs)))
  flush.console()
  
  cat("[INFO] Loading Monocle outputs...\n")
  flush.console()
  cell_df <- fread(file.path(opt$monocle_dir, "monocle3_cells.csv"), data.table = FALSE)
  edge_df <- fread(file.path(opt$monocle_dir, "principal_graph_edges.csv"), data.table = FALSE)
  
  if (!("cell_id" %in% colnames(cell_df))) {
    stop("monocle3_cells.csv must contain cell_id")
  }
  
  cat("[INFO] Aligning cell order to h5ad...\n")
  flush.console()
  cell_df <- merge(cell_df, obs, by = "cell_id", all.x = TRUE, suffixes = c("", "_obs"))
  cell_df <- cell_df[cell_df$cell_id %in% obs$cell_id, , drop = FALSE]
  cell_df <- cell_df[match(obs$cell_id, cell_df$cell_id), , drop = FALSE]
  cat(sprintf("[INFO] Aligned cell_df: %d cells\n", nrow(cell_df)))
  flush.console()
  
  label_key <- if (!is.null(opt$label_key)) opt$label_key else resolve_label_key(cell_df)
  cat(sprintf("[INFO] Using label_key: %s\n", label_key))
  flush.console()
  
  stage_order_map <- load_json_maybe(opt$stage_order_json)
  future_map <- load_json_maybe(opt$future_map_json)
  
  if (!(opt$future_label_key %in% colnames(cell_df))) {
    if (is.null(future_map)) {
      stop(sprintf("%s not found in cell dataframe and no --future-map-json was provided.", opt$future_label_key))
    }
    raw_labels <- as.character(cell_df[[label_key]])
    cell_df[[opt$future_label_key]] <- unname(future_map[raw_labels])
  }
  
  root_nodes <- readLines(file.path(opt$monocle_dir, "root_nodes.txt"), warn = FALSE)
  root_nodes <- trimws(root_nodes)
  root_nodes <- root_nodes[nzchar(root_nodes)]
  cat(sprintf("[INFO] Loaded %d root nodes\n", length(root_nodes)))
  flush.console()
  
  eval_mask <- is.finite(suppressWarnings(as.numeric(cell_df$pseudotime)))
  cat(sprintf("[INFO] Finite pseudotime cells: %d / %d\n", sum(eval_mask), nrow(cell_df)))
  flush.console()
  
  cat("[INFO] Starting T3 computation...\n")
  flush.console()
  t3 <- compute_t3_metrics(
    df = cell_df,
    pred_time_key = "pseudotime",
    reference_time_key = opt$reference_time_key,
    label_key = label_key,
    stage_order_map = stage_order_map,
    eval_mask = eval_mask
  )
  fwrite(t3, file.path(opt$outdir, "T3_metrics.csv"))
  cat("[INFO] Finished T3. Wrote T3_metrics.csv\n")
  flush.console()
  
  cat("[INFO] Starting T4 probability computation...\n")
  flush.console()
  t4_obj <- compute_monocle_future_probabilities(
    cell_df = cell_df,
    edge_df = edge_df,
    root_nodes = root_nodes,
    lineage_key = opt$future_label_key,
    temperature = opt$temperature
  )
  
  fwrite(t4_obj$proba_df, file.path(opt$outdir, "T4_probabilities.csv"))
  fwrite(t4_obj$leaf_meta, file.path(opt$outdir, "T4_leaf_class_map.csv"))
  fwrite(t4_obj$per_leaf_scores, file.path(opt$outdir, "T4_leaf_scores_long.csv"))
  cat("[INFO] Finished T4 probabilities. Wrote T4_probabilities.csv / T4_leaf_class_map.csv / T4_leaf_scores_long.csv\n")
  flush.console()
  
  cat("[INFO] Starting T4 metrics...\n")
  flush.console()
  t4 <- compute_t4_metrics(
    cell_df = cell_df,
    proba_df = t4_obj$proba_df,
    class_order = t4_obj$class_order,
    lineage_key = opt$future_label_key,
    eval_mask = eval_mask
  )
  fwrite(t4, file.path(opt$outdir, "T4_metrics.csv"))
  cat("[INFO] Finished T4 metrics. Wrote T4_metrics.csv\n")
  flush.console()
  
  #merged <- merge(cell_df, t4_obj$proba_df, by = "cell_id", all.x = TRUE)
  #fwrite(merged, file.path(opt$outdir, "monocle3_cells_with_eval.csv"))
  cat("[INFO] Wrote monocle3_cells_with_eval.csv\n")
  flush.console()
  
  summary <- list(
    label_key = label_key,
    reference_time_key = opt$reference_time_key,
    future_label_key = opt$future_label_key,
    n_cells = nrow(cell_df),
    t3 = as.list(t3[1, ]),
    t4 = as.list(t4[1, ])
  )
  write(toJSON(summary, auto_unbox = TRUE, pretty = TRUE), file.path(opt$outdir, "summary.json"))
  cat("[INFO] Wrote summary.json\n")
  flush.console()
  
  cat(toJSON(summary, auto_unbox = TRUE, pretty = TRUE))
}

main()