#' Save png, pdf and plot object
#'
#' @param plot `ggplot2` object to save
#' @param filename Filename
#' @param root Directories to add to /plots or /data
#' @param h Height
#' @param w Width
#' @param s Scale
#'
#' @return NULL
#' @export
SavePlot <- function(
  plot,
  filename,
  root = "publication",
  h = 6,
  w = 6,
  s = 1,
  save.data = TRUE
) {
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.png"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  ggplot2::ggsave(plot = plot, filename = glue::glue("plots/{root}/{filename}.pdf"),
                  scale = s, width = w, height = h, units = "in", dpi = 300)
  if (save.data) {
    saveRDS(plot, file = glue::glue("data/{root}/plots/{filename}.rds"))
  }
  usethis::ui_done("Saved")
}

#' RemoveAxes
#'
#' Modified from Seurat::NoAxes()
#'
#' @return
#' @export
RemoveAxes <- function (..., keep.text = FALSE, keep.ticks = FALSE)
{
  blank <- element_blank()
  no_axes_theme <- theme(axis.line.x = blank,
                         axis.line.y = blank,
                         validate = TRUE, ...)
  if (!keep.text) {
    no_axes_theme <- no_axes_theme + theme(
      axis.text.x = blank,
      axis.text.y = blank,
      validate = TRUE, ...)
  }

  if (!keep.ticks) {
    no_axes_theme <- no_axes_theme + theme(
      axis.ticks.x = blank,
      axis.ticks.y = blank,
      validate = TRUE, ...)
  }

  return(no_axes_theme)
}

#' Remove backgrounds
#'
#' @param outline Keep plot outline
#' @param ...
#'
#' @return
#' @export
RemoveBackgrounds <- function(outline = FALSE, ...)
{
  if (outline) {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "black", size = 1),
                         plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.border = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.grid = element_blank(),
                         axis.line = element_blank(),
                         validate = TRUE, ...)
  } else {
    no_bg_theme <- theme(panel.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         legend.box.background = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.border = element_rect(fill = "transparent", color = "transparent", size = 0),
                         panel.grid = element_blank(),
                         validate = TRUE, ...)
  }

  return(no_bg_theme)
}


#' Space axis titles away from plot area
#'
#' @param scale Increase spacing
#' @param ...
#'
#' @return
#' @export
SpaceAxisTitles <- function(scale = 1, ...) {
  theme (
    axis.title.x = element_text(face = "plain", margin = margin(12*scale, 0, 0, 0)),
    axis.title.y = element_text(face = "plain", margin = margin(0, 12*scale, 0, 0)),
    validate = TRUE
  )
}

#' General plotting theme
#'
#' @param base_size
#' @param ...
#'
#' @return
#' @export
GeneralTheme <- function(base_size, ...) {
  theme_classic(base_size = base_size, ...) +
    ggeasy::easy_all_text_color("black") +
    theme(
      axis.line = element_blank(),
      plot.title = element_text(size =  base_size, color = "black", face = "bold", margin = margin(0,0,4,0)),
      plot.subtitle = element_text(size = base_size - 2, color = "black", margin = margin(0,0,4,0)),
      panel.background = element_rect(fill = "transparent", color = "black", size = 1),
      plot.background = element_rect(fill = "transparent", color = "transparent", size = 0),
      panel.border = element_rect(size = 1, color = "black", fill = "transparent"),
      plot.caption = element_text(hjust = 0, color = "gray40", margin = margin(12)),
      legend.title = element_text(size = base_size, face = "plain"),
      legend.text = element_text(size = base_size - 2),
      legend.background = element_rect(fill = "transparent", color = "transparent"),
      legend.box.background = element_rect(fill = "transparent", color = "transparent"),
      legend.position = "right",
      legend.justification = "top",
      legend.key.size = unit(1, "line"),
      validate = TRUE
    )
}

#' Run Scrublet
#'
#' @param seurat.obj Seurat object
#' @param batch.var Batch variable
#'
#' @return Seurat object
#' @export
Scrub <- function(seurat.obj, batch.var = "batch") {
  scrublet <- reticulate::import("scrublet")
  seurat.obj$scrublet_score <- "NA"
  seurat.obj$scrublet_label <- "NA"
  sample_df <- as.data.frame(table(seurat.obj[[batch.var]]))
  colnames(sample_df) <- c("batch", "freq")
  sample_df$batch <- as.character(sample_df$batch)
  sample_df$freq <- as.numeric(sample_df$freq)

  for (i in seq_along(1:length(sample_df$batch))) {
    freq <- as.numeric(sample_df$freq[i])
    cells <- (seurat.obj[[batch.var, drop = TRUE]] == sample_df$batch[i])
    if (freq < 100) {
      message(glue(">> Only {freq} cells, skipping doublet prediction"))
      seurat.obj[["scrublet_score"]][cells, ] <- NA
      seurat.obj[["scrublet_label"]][cells, ] <- NA
    } else {
      matrix <- as.matrix(GetAssayData(seurat.obj, slot = "counts")[, cells])
      scrublet_object <- scrublet$Scrublet(t(matrix), expected_doublet_rate = 4.6e-06*freq)
      message(glue(">> Scrublet object created for iteration {i}/{length(sample_df$batch)}"))
      scores <- scrublet_object$scrub_doublets(min_counts = 3, min_cells = 3,
                                               min_gene_variability_pctl = 85, n_prin_comps = as.integer(30), verbose = TRUE)
      message(glue(">> Identified {sum(as.vector(scores[[2]]))}/{length(scores[[2]])} cells as doublets"))
      seurat.obj[["scrublet_score"]][cells, ] <- scores[[1]]
      seurat.obj[["scrublet_label"]][cells, ] <- scores[[2]]
    }
  }
  return(seurat.obj)
}


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#' Process Seurat seurat.object through PCA
#'
#' @param seurat.obj Seurat object
#' @param n.var.features # of variable features to use
#' @param n.pcs # of PCs to use
#' @param remove.genes Genes to remove from variable features and PCA
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom Seurat NormalizeData FindVariableFeatures ScaleData RunPCA
ReduceDims <- function(
  seurat.obj,
  n.var.features = 3000,
  n.pcs = 50,
  remove.genes = NULL,
  pca.name = "pca",
  normalize.data = TRUE,
  ...
) {
  if (normalize.data) {
    seurat.obj <- NormalizeData(object = seurat.obj)
  }
  seurat.obj <- FindVariableFeatures(object = seurat.obj, selection.method = "vst", nfeatures = n.var.features)
  if (!missing(remove.genes) & !is.null(remove.genes)) {
    usethis::ui_done("Removed {length(remove.genes)} genes from variable features")
    VariableFeatures(seurat.obj) <- setdiff(VariableFeatures(seurat.obj), remove.genes)
  }
  seurat.obj <- ScaleData(object = seurat.obj, features = VariableFeatures(seurat.obj), block.size = 1000)
  usethis::ui_todo("Running PCA...")
  seurat.obj <- RunPCA(object = seurat.obj, npcs = n.pcs, features = VariableFeatures(seurat.obj),
                       verbose = FALSE, seed.use = 1, weight.by.var = TRUE, reduction.name = pca.name)
  return(seurat.obj)
}

Embed <- function(
  seurat.obj,
  reduction = "pca",
  heuristics.name = "pca_metrics",
  dims.use,
  knn.use = 20,
  store.graph = TRUE
) {
  dims_name <- glue::glue("{reduction}_dims")
  if (!missing(dims.use)) {
    seurat.obj@misc[[dims_name]] <- dims.use
  } else if (unique(seurat.obj@misc[[heuristics.name]]$tp_pc) == 101 | unique(seurat.obj@misc[[heuristics.name]]$tp_pc) == 0 | is.nan(unique(seurat.obj@misc[[heuristics.name]]$tp_pc))) {
    seurat.obj@misc[[dims_name]] <- unique(seurat.obj@misc[[heuristics.name]]$percent_cutoff_pcs)
  } else {
    seurat.obj@misc[[dims_name]] <- unique(seurat.obj@misc[[heuristics.name]]$tp_pc)
  }
  usethis::ui_info("\nUsing {seurat.obj@misc[[dims_name]]} dimensions for UMAP...")

  seurat.obj <- Seurat::RunUMAP(object = seurat.obj, reduction = reduction, reduction.name = glue::glue("{reduction}_umap"),
                                reduction.key = glue::glue("{substring(reduction, 1, 1)}UMAP_"),
                                dims = 1:seurat.obj@misc[[dims_name]], seed.use = 1)
  if (store.graph) {
    seurat.obj <- PrepareGraph(object = seurat.obj, reduction = reduction, dims = seurat.obj@misc[[dims_name]], knn = knn.use)
  }
  return(seurat.obj)
}

#' Create igraph-compatible graph and save in Seurat object
#'
#' @param object Seurat object
#' @param reduction Reduced dimension slot to pull from
#' @param dims Dimensions to use
#' @param knn Number of k-nearest neighbors to use
#'
#' @return Seurat object with igraph graph stored in `object@misc$westerlund_graph`
#' @export
PrepareGraph <- function(object, reduction, dims, knn) {
  graph.name <- glue::glue("{reduction}_snn_{knn}")
  stopifnot(ncol(object@reductions[[reduction]]) >= dims)
  object <- Seurat::FindNeighbors(object = object, k.param = knn, prune.SNN = 1/15, dims = 1:dims,
                                  reduction = reduction, graph.name = c("nn", graph.name), compute.SNN = TRUE, verbose = FALSE)
  g <- object@graphs[[graph.name]]
  attributes(g)[[1]] <- NULL
  attributes(g)$class <- "dgCMatrix"
  #adj_matrix <- Matrix::Matrix(as.matrix(object@graphs[[graph.name]]), sparse = TRUE)
  g <- igraph::graph_from_adjacency_matrix(adjmatrix = g, mode = "undirected", weighted = TRUE, add.colnames = TRUE)
  object@misc[[glue::glue("{graph.name}")]] <- g
  return(object)
}

#' Leiden graph-based clustering
#'
#' leidenbase wrapper function to parse and select results
#'
#' @param res.param Resolution
#' @param graph Graph
#' @param seed.param Seed
#' @param return.params Return cluster-level results
#' @param return.clusters Return cell-level results
#' @param num.iter Number of leidenbase iterations
#'
#' @return List containing clustering and/or parameter results
#' @export
#'
#' @examples
Cluster <- function(res.param, graph, seed.param, return.params = TRUE, return.clusters = FALSE, num.iter = 1) {

  cluster_res <- leidenbase::leiden_find_partition(graph,
                                                   partition_type = "CPMVertexPartition",
                                                   initial_membership = NULL,
                                                   edge_weights = NULL,
                                                   node_sizes = NULL,
                                                   seed = seed.param,
                                                   resolution_parameter = res.param,
                                                   num_iter = num.iter,
                                                   verbose = FALSE)

  param_res <- data.frame(
    res_param = res.param,
    quality = cluster_res[["quality"]],
    modularity = cluster_res[["modularity"]],
    significance = cluster_res[["significance"]],
    cluster_count = max(cluster_res[["membership"]]))

  if(return.params & return.clusters) {
    return(list(cluster_res, param_res))
  } else if (return.params == TRUE & return.clusters == FALSE) {
    return(param_res)
  } else if (return.params == FALSE & return.clusters == TRUE) {
    return(cluster_res)
  }
}

#' Automated clustering at max or near-max resolution
#'
#' @param seurat.obj Seurat object
#' @param graph.name Graph name
#' @param col.name Clustering column name
#' @param mod.percentile Percentile of modularity to use
#'
#' @return Seurat object
#' @export
AutoCluster <- function(seurat.obj, graph.name = "pca_snn_20", col.name = "pca_leiden", mod.percentile = 1, min.cluster = 9) {

  StartFuture()

  # scan resolutions broadly
  first_results <- ScanResolutions(g = seurat.obj@misc[[graph.name]],
                                   res.lower = 1E-10,
                                   res.upper = 1,
                                   num.res = 50,
                                   clusters.lb = 10,
                                   clusters.ub = 50,
                                   use.mod = FALSE,
                                   seed.param = 1)

  # scan again for optimal modularity
  second_results <- ScanResolutions(g = seurat.obj@misc[[graph.name]],
                                    res.lower = first_results$range[1],
                                    res.upper = first_results$range[2],
                                    num.res = 30,
                                    use.mod = TRUE,
                                    mod.percentile = 0.75,
                                    seed.param = 1)


  StopFuture()

  # identify and cluster at optimal
  param_results <- second_results$results
  max_modularity <- max(param_results$modularity)
  modularity_cutoff <- max_modularity*mod.percentile
  max_resolution <- param_results$res_param[param_results$modularity == max_modularity]
  usethis::ui_info("{sum(param_results$modularity >= modularity_cutoff)} resolutions")

  if (any(param_results$modularity >= modularity_cutoff)) {
    final_resolution_parameter <- max(param_results$res_param[param_results$modularity >= modularity_cutoff])
  } else if (!any(param_results$modularity >= modularity_cutoff)) {
    final_resolution_parameter <- param_results$res_param[param_results$modularity == max_modularity]
  }
  usethis::ui_done("Clustering @ {final_resolution_parameter} vs. maximum resolution @ {max_resolution}")

  final_results <- Cluster(res.param = final_resolution_parameter,
                           graph = seurat.obj@misc[[graph.name]],
                           seed.param = 1,
                           return.params = TRUE,
                           return.clusters = TRUE,
                           num.iter = 30)

  # group and return membership vector
  ids <- final_results[[1]]$membership
  names(ids) <- Cells(seurat.obj)
  ids <- GroupSingletons(ids, seurat.obj@graphs[[graph.name]],
                         min.size = min.cluster, group.singletons = TRUE, verbose = TRUE)
  seurat.obj <- AddMetaData(seurat.obj, ids, col.name = col.name)
  usethis::ui_done("Identified {length(unique(seurat.obj[[col.name, drop = TRUE]]))} Leiden clusters")
  Clean()
  return(seurat.obj)
}

#' Scan resolutions for appropriate ranges
#'
#' @param g igraph graph
#' @param res.lower starting lower bound
#' @param res.upper starting upper bound
#' @param num.res number of steps to take
#' @param clusters.lb cluster lower bound
#' @param clusters.ub cluster upper bound
#' @param use.mod Use modularity instead of cluster numbers
#' @param mod.percentile Range
#' @param seed.param Seed
#'
#' @return Vector of lower and upper resolution bounds based on modularity percentiles or clusters
#' @export
ScanResolutions <- function(g,
                            res.lower = 1E-10,
                            res.upper = 1,
                            num.res = 100,
                            clusters.lb = 10,
                            clusters.ub = 50,
                            use.mod = FALSE,
                            mod.percentile = 50,
                            seed.param = 1
) {

  res_params <-  signif(exp(seq(log(res.lower), log(res.upper), length.out = num.res)), 3)
  usethis::ui_todo("Scanning resolutions ({res_params[1]}, {res_params[length(res_params)]})")
  params_res <- future.apply::future_lapply(X = res_params, FUN = Cluster, graph = g, seed.param = seed.param, return.params = TRUE)
  params_res <- do.call(rbind, params_res)

  if (use.mod) {
    max_mod <- max(params_res$modularity)
    mod_threshold <- quantile(params_res$modularity[params_res$modularity > 0], mod.percentile)
    lower_mod_res <- min(params_res$res_param[params_res$modularity >= mod_threshold])
    upper_mod_res <- max(params_res$res_param[params_res$modularity >= mod_threshold])
    params <- c(lower_mod_res, upper_mod_res)
  } else {
    lower_cluster_res <- min(params_res$res_param[params_res$cluster_count >= clusters.lb])
    upper_cluster_res <- max(params_res$res_param[params_res$cluster_count <= clusters.ub])
    params <- c(lower_cluster_res, upper_cluster_res)
  }
  usethis::ui_done("Found bounds ({params[1]}, {params[2]})")
  return(list(results = params_res, range = params))
}

# Modified from Seurat
GroupSingletons <- function(ids, snn, min.size = 9, clusters.to.merge, group.singletons = TRUE, verbose = TRUE) {

  usethis::ui_info("Merging small or provided clusters")

  # identify singletons
  singletons <- c()
  singletons <- names(x = which(x = table(ids) <= min.size))
  singletons <- intersect(x = unique(x = ids), singletons)

  if (!missing(clusters.to.merge)) {
    singletons <- c(singletons, as.character(clusters.to.merge))
  }

  if (!group.singletons) {
    ids[which(ids %in% singletons)] <- "singleton"
    return(ids)
  }

  # calculate connectivity of singletons to other clusters, add singleton
  # to cluster it is most connected to

  if (!is_empty(singletons)) {
    cluster_names <- as.character(x = unique(x = ids))
    cluster_names <- setdiff(x = cluster_names, y = singletons)
    connectivity <- vector(mode = "numeric", length = length(x = cluster_names))
    names(x = connectivity) <- cluster_names
    new.ids <- ids
    for (i in singletons) {
      i.cells <- names(which(ids == i))
      for (j in cluster_names) {
        j.cells <- names(which(ids == j))
        subSNN <- snn[i.cells, j.cells]
        set.seed(1) # to match previous behavior, random seed being set in WhichCells
        if (is.object(x = subSNN)) {
          connectivity[j] <- sum(subSNN) / (nrow(x = subSNN) * ncol(x = subSNN))
        } else {
          connectivity[j] <- mean(x = subSNN)
        }
      }
      m <- max(connectivity, na.rm = T)
      mi <- which(x = connectivity == m, arr.ind = TRUE)
      closest_cluster <- sample(x = names(x = connectivity[mi]), 1)
      ids[i.cells] <- closest_cluster
    }
  }

  if (length(x = singletons) > 0 && verbose) {
    message(paste(
      length(x = singletons),
      "singletons identified.",
      length(x = unique(x = ids)),
      "final clusters."
    ))
  }

  return(ids)
}

Walktrap <- function(
  seurat.obj,
  col.name = "pca_walktrap",
  graph.name = "pca_snn_20",
  min.size = 9,
  ...
) {
  walktrap <- igraph::cluster_walktrap(seurat.obj@misc[[graph.name]], ...)
  ids <- walktrap$membership
  names(ids) <- Cells(seurat.obj)
  ids <- GroupSingletons(ids, seurat.obj@graphs[[graph.name]],
                         min.size = min.size, group.singletons = TRUE, verbose = TRUE)
  seurat.obj <- AddMetaData(seurat.obj, ids, col.name = col.name)
  usethis::ui_done("Identified {length(unique(seurat.obj[[col.name, drop = TRUE]]))} walktrap clusters")
  return(seurat.obj)
}