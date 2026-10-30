#!/usr/bin/env Rscript
# Single-cell RNA-seq
# Revised 2026-10-22. Default: reproduce and evaluate the manuscript cell populations.
# Annotation combines canonical markers, cluster DE profiles and eGFP context (Fig. S9C).
# eGFP is supporting annotation evidence, NOT a per-cell inclusion threshold.
# Upstream reporter quantification: custom FASTA/GTF + Cell Ranger (7.2.0 in Methods).
# Reference-build sequences/commands are not reconstructed by this downstream script.
# Rscript TMEM41B_scRNAseq_analysis_revised_20261022.R /absolute/project/root [--quick]
# Outputs: project/results_TMEM41B_review_20261022 (never overwrites submission inputs).
# --quick: annotation, QC, raw reporter summaries and plot reproduction only.
# Default additionally runs DE markers, SMC reclustering and parameter sensitivity.
# --rebuild: corrected-QC, reporter-excluded discovery pipeline (new clusters need review).
# All relative paths resolve from the working directory, not the script directory.
# Dependencies: Seurat (supporting leiden_method='leidenbase'), Matrix, igraph,
#       ggplot2, future, leidenbase. Packages are not installed or upgraded automatically.
# Retains reference reconstruction, Louvain comparison, multi-seed checks, and marker export.
# Generic reconstruction functions precede project-specific reproducibility functions.

cfg <- list(
  samples = data.frame(
    sample = c('sample_1', 'sample_2', 'sample_3', 'sample_4'),
    group = c('control', 'control', 'treatment', 'treatment'),
    path = file.path('data', c('sample_1', 'sample_2', 'sample_3', 'sample_4'))
  ),
  output_dir = 'results',
  species = 'mouse',                         # mouse or human; gene symbols must match
  cell_naming = 'sample_prefix',             # or index_suffix (legacy barcode_1 format)
  archive_rds = NULL,                        # Optional: archived Seurat object
  reference_clusters = NULL,                 # Optional: CSV with cell and cluster as the first two columns
  reference_embeddings = NULL,               # Optional: CSV with cell followed by two tSNE coordinates
  reference_markers = NULL,                  # Optional: TSV with cluster/gene/avg_log2FC
  use_archived_features = FALSE,             # Requires archive_rds when TRUE
  require_reference_qc_match = FALSE,        # Require exact agreement between QC-passing and reference cell sets
  qc = list(min_cells = 3, min_features = 250, feature_min = 250,
            feature_max = 7500, count_min = 100, mt_max = 25,
            mt_pattern = NULL, ribo_pattern = NULL),
  nfeatures = 2000, npcs = 50, anchor_dims = 1:20, dims = 1:40,
  seed = 42, tsne_seed = 1, tsne_perplexity = 30,
  resolutions = seq(0.3, 1.1, 0.1), seeds = 42:46,
  preferred_resolution = 0.7, baseline_resolution = 0.7,
  selection = 'strict',                      # revised marker-independent rule
  # Revised selection uses connectivity and seed stability ONLY.
  # Marker coherence and reference similarity are diagnostic outputs, never selection gates.
  min_seed_ari = 0.85, coherence_tolerance = 0.05,
  marker_panels = NULL,                      # Custom named list; NULL uses the example panels below
  annotation_min_cells = 20, annotation_margin = 0.15,
  annotation_agreement = 0.5,
  make_plots = TRUE, plot_genes = character(), raw_tsne_diagnostic = FALSE,
  marker_assays = c('RNA'),
  marker_min_pct = 0.25, marker_logfc = 0.25,
  marker_p_column = 'p_val_adj', marker_p_cutoff = 0.05, # Seurat Bonferroni, NOT BH FDR
  marker_top_n = c(10, 20, 50, 100)
)

# Mouse example panels come from the original scripts; adapt panels to the tissue and study.
# Scores provide provisional annotations; review markers, differential expression, and biology.
default_panels <- function(species) {
  if (species == 'human') return(list(
    T_cell = c('CD3D', 'CD3E', 'TRAC'), B_cell = c('CD19', 'CD79A', 'MS4A1'),
    NK_cell = c('NCR1', 'KLRD1', 'NKG7'),
    Myeloid = c('LYZ', 'LST1', 'CSF1R', 'CD68'),
    Endothelial = c('PECAM1', 'CDH5', 'VWF'),
    Fibroblast = c('DCN', 'LUM', 'COL1A1'),
    Plasma_cell = c('JCHAIN', 'MZB1', 'SDC1'),
    Erythroid = c('HBB', 'HBA1', 'ALAS2')
  ))
  list(
    Smooth_muscle = c('Myh11', 'Cnn1', 'Tagln', 'Myl9', 'Acta2'),
    Endothelial = c('Pecam1', 'Cdh5', 'Vwf', 'Tek', 'Flt1'),
    Fibroblast = c('Dcn', 'Lum', 'Pdgfra', 'Col1a1', 'Dpt'),
    Macrophage = c('Adgre1', 'Cd68', 'C1qb', 'Csf1r', 'Lyz2'),
    T_cell = c('Cd3d', 'Cd3e', 'Cd3g', 'Trac'),
    NK_cell = c('Ncr1', 'Klrd1', 'Klrb1c', 'Nkg7'),
    B_cell = c('Cd19', 'Cd79a', 'Cd79b', 'Bank1'),
    Schwann = c('Plp1', 'Cnp', 'Sox10', 'Mpz'),
    Dendritic = c('Flt3', 'Itgax', 'Cd209a', 'H2-Ab1'),
    Adipocyte = c('Adipoq', 'Cidec', 'Dgat2', 'Fabp4'),
    Erythroid = c('Hbb-bs', 'Hbb-bt', 'Hba-a1', 'Hba-a2', 'Alas2'),
    Plasma_cell = c('Jchain', 'Mzb1', 'Sdc1', 'Igkc'),
    Neutrophil = c('S100a8', 'S100a9', 'Retnlg', 'Csf3r'),
    Mast_cell = c('Kit', 'Cpa3', 'Ms4a2', 'Tpsb2'),
    ILC2_candidate = c('Il5', 'Il1rl1', 'Gata3', 'Rora')
  )
}

write_tsv <- function(data, name, cfg) {
  write.table(data, file.path(cfg$output_dir, 'tables', paste0(name, '.tsv')),
              sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA')
}
checkpoint <- function(name, cfg) file.path(cfg$output_dir, 'checkpoints', paste0(name, '.rds'))
save_plot <- function(p, name, cfg, width = 9, height = 7) {
  for (ext in c('pdf', 'png')) ggplot2::ggsave(
    file.path(cfg$output_dir, 'figures', paste0(name, '.', ext)), p,
    width = width, height = height, dpi = 150, limitsize = FALSE)
}

# ARI/NMI are label-invariant; single-cluster cases are handled explicitly.
ari <- function(a, b) {
  z <- table(a, b); n <- sum(z)
  if (n < 2) return(1)
  choose2 <- function(v) v * (v - 1) / 2
  A <- sum(choose2(rowSums(z))); B <- sum(choose2(colSums(z)))
  expected <- A * B / choose2(n); denominator <- (A + B) / 2 - expected
  if (denominator == 0) return(1)
  (sum(choose2(z)) - expected) / denominator
}
nmi <- function(a, b) {
  p <- table(a, b); p <- p / sum(p)
  pa <- rowSums(p); pb <- colSums(p); q <- outer(pa, pb); ok <- p > 0
  h <- -sum(pa[pa > 0] * log(pa[pa > 0])) - sum(pb[pb > 0] * log(pb[pb > 0]))
  if (h == 0) return(1)
  2 * sum(p[ok] * log(p[ok] / q[ok])) / h
}
separate_singletons <- function(labels) {
  labels <- as.character(labels); i <- which(labels == 'singleton')
  labels[i] <- paste0('singleton_', seq_along(i)); labels
}
connectivity <- function(graph, labels) {
  do.call(rbind, lapply(unique(labels), function(k) {
    cc <- igraph::components(igraph::induced_subgraph(graph, which(labels == k)))
    data.frame(cluster = k, n = sum(labels == k), components = cc$no,
               largest = max(cc$csize))
  }))
}
run_leiden <- function(snn, resolution, seed) {
  z <- Seurat::FindClusters(snn, algorithm = 4, leiden_method = 'leidenbase',
    resolution = resolution, n.iter = 10, random.seed = seed,
    group.singletons = FALSE, verbose = FALSE)
  separate_singletons(z[, 1])
}

read_reference <- function(cfg) {
  if (is.null(cfg$reference_clusters)) return(NULL)
  d <- read.csv(cfg$reference_clusters, colClasses = 'character', check.names = FALSE)
  if (ncol(d) < 2 || !nrow(d) || anyNA(d[, 1:2]) || anyDuplicated(d[[1]]))
    stop('Reference labels must contain unique cell IDs and nonmissing cluster labels.')
  setNames(d[[2]], d[[1]])
}

rebuild <- function(cfg, reference) {
  mt <- cfg$qc$mt_pattern
  ribo <- cfg$qc$ribo_pattern
  if (is.null(mt)) mt <- if (cfg$species == 'mouse') '^mt-' else '^MT-'
  if (is.null(ribo)) ribo <- if (cfg$species == 'mouse') '^Rp[sl]' else '^RP[SL]'
  objects <- audits <- vector('list', nrow(cfg$samples))
  for (i in seq_len(nrow(cfg$samples))) {
    s <- cfg$samples[i, ]; message('Reading sample: ', s$sample)
    counts <- Read10X(s$path)
    if (is.list(counts)) {
      if (!'Gene Expression' %in% names(counts)) stop('The 10x input does not contain Gene Expression.')
      counts <- counts[['Gene Expression']]
    }
    x <- CreateSeuratObject(counts, project = s$sample,
      min.cells = cfg$qc$min_cells, min.features = cfg$qc$min_features)
    x$sample <- s$sample; x$group <- s$group
    if (!any(grepl(mt, rownames(x)))) stop('No mitochondrial genes matched. Check species, gene IDs, or mt_pattern.')
    x$percent.mt <- PercentageFeatureSet(x, pattern = mt)
    x$percent.ribo <- PercentageFeatureSet(x, pattern = ribo)
    cell_ids <- if (cfg$cell_naming == 'index_suffix') paste0(colnames(x), '_', i) else
      paste(s$sample, colnames(x), sep = '_')
    x <- RenameCells(x, new.names = cell_ids)
    passed <- colnames(x)[x$nFeature_RNA > cfg$qc$feature_min &
      x$nFeature_RNA < cfg$qc$feature_max & x$nCount_RNA > cfg$qc$count_min &
      x$percent.mt < cfg$qc$mt_max]
    # When reference labels are supplied, rebuild that cell set and audit differences from current QC.
    keep <- if (is.null(reference)) passed else intersect(colnames(x), names(reference))
    identical_qc <- setequal(passed, keep)
    if (cfg$require_reference_qc_match && !identical_qc) stop('Current QC and reference cell sets differ for sample: ', s$sample)
    audits[[i]] <- data.frame(sample = s$sample, created = ncol(x), qc_pass = length(passed),
      retained = length(keep), qc_matches_retained = identical_qc)
    if (!length(keep)) stop('No cells retained for sample: ', s$sample)
    x <- subset(x, cells = keep)
    x <- NormalizeData(x, scale.factor = 10000, verbose = FALSE)
    # Reporters remain in RNA counts for auditing but are absent from feature modeling.
    endogenous <- setdiff(rownames(x), cfg$reporter_features)
    temp <- FindVariableFeatures(subset(x, features = endogenous), selection.method = 'vst',
      nfeatures = cfg$nfeatures, verbose = FALSE)
    VariableFeatures(x) <- VariableFeatures(temp); rm(temp)
    objects[[i]] <- x
  }
  write_tsv(do.call(rbind, audits), 'qc_audit', cfg)
  archive <- if (!is.null(cfg$archive_rds)) readRDS(cfg$archive_rds) else NULL
  features <- if (cfg$use_archived_features) {
    archive@commands$FindIntegrationAnchors@params$anchor.features
  } else SelectIntegrationFeatures(objects, nfeatures = cfg$nfeatures)
  features <- setdiff(features, cfg$reporter_features)
  if (!length(features) || !all(features %in% Reduce(intersect, lapply(objects, rownames))))
    stop('Integration features are empty or are not present in every sample.')
  set.seed(cfg$seed)
  anchors <- FindIntegrationAnchors(objects, anchor.features = features,
    dims = cfg$anchor_dims, reduction = 'cca', k.anchor = 5, k.filter = 200,
    k.score = 30, max.features = 200, nn.method = 'annoy', n.trees = 50)
  x <- IntegrateData(anchors, dims = cfg$anchor_dims, k.weight = 100, sd.weight = 1)
  DefaultAssay(x) <- 'integrated'
  x <- ScaleData(x, features = rownames(x), verbose = FALSE)
  pca_features <- if (cfg$use_archived_features)
    archive@commands$RunPCA.integrated@params$features else VariableFeatures(x)
  pca_features <- setdiff(pca_features, cfg$reporter_features)
  if (!length(pca_features) || !all(pca_features %in% rownames(x))) stop('Invalid PCA features.')
  if (cfg$npcs >= min(length(pca_features), ncol(x))) stop('Insufficient cells or features. Reduce npcs.')
  VariableFeatures(x) <- pca_features
  x <- RunPCA(x, features = pca_features, npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
  x <- RunUMAP(x, dims = cfg$dims, seed.use = cfg$seed, n.neighbors = 30,
               metric = 'cosine', min.dist = 0.3, verbose = FALSE)
  x <- RunTSNE(x, dims = cfg$dims, seed.use = cfg$tsne_seed,
               perplexity = cfg$tsne_perplexity, max_iter = 1000, verbose = FALSE)
  x <- FindNeighbors(x, dims = cfg$dims, k.param = 20, prune.SNN = 1/15,
    nn.method = 'annoy', n.trees = 50, annoy.metric = 'euclidean',
    l2.norm = FALSE, graph.name = c('integrated_nn', 'integrated_snn'))
  if (!is.null(reference)) {
    if (!setequal(colnames(x), names(reference))) stop('Some reference cells were not reconstructed. Check cell naming.')
    x$reference_cluster <- unname(reference[colnames(x)])
  }
  if (!is.null(archive) && 'pca' %in% Reductions(archive)) {
    shared <- intersect(colnames(archive), colnames(x))
    a <- Embeddings(archive, 'pca'); b <- Embeddings(x, 'pca')
    pcs <- seq_len(min(ncol(a), ncol(b), max(cfg$dims)))
    if (length(shared) >= 3) write_tsv(data.frame(PC = pcs,
      correlation = vapply(pcs, function(j) cor(a[shared, j], b[shared, j]), numeric(1))),
      'reference_pca_comparison', cfg)
  }
  DefaultAssay(x) <- 'RNA'; x <- NormalizeData(x, verbose = FALSE)
  DefaultAssay(x) <- 'integrated'
  saveRDS(x, checkpoint('precluster', cfg), compress = FALSE)
  x
}

marker_evidence <- function(x, cfg) {
  panels <- cfg$marker_panels
  if (is.null(panels)) panels <- default_panels(cfg$species)
  dat <- GetAssayData(x, assay = 'RNA', layer = 'data')
  panels <- lapply(panels, intersect, rownames(dat))
  empty <- lengths(panels) == 0
  if (any(empty)) warning('Removing panels with no matching genes: ', paste(names(panels)[empty], collapse = ', '))
  panels <- panels[!empty]
  if (length(panels) < 2) stop('At least two marker panels with matching genes are required.')
  write_tsv(do.call(rbind, lapply(names(panels), function(k)
    data.frame(celltype = k, gene = panels[[k]]))), 'annotation_marker_panels', cfg)
  # Densify only the small marker-gene subset, not the full expression matrix.
  z <- t(scale(t(as.matrix(dat[unique(unlist(panels)), , drop = FALSE]))))
  z[!is.finite(z)] <- 0
  scores <- vapply(panels, function(p) colMeans(z[p, , drop = FALSE]), numeric(ncol(x)))
  rownames(scores) <- colnames(x)
  evidence <- list(panels = panels, scores = scores,
    cell_type = colnames(scores)[max.col(scores, ties.method = 'first')])
  saveRDS(evidence, checkpoint('marker_evidence', cfg)); evidence
}

select_clusters <- function(x, evidence, cfg) {
  snn <- x[['integrated_snn']]
  graph <- igraph::graph_from_adjacency_matrix(as(snn, 'dgCMatrix'),
    mode = 'undirected', weighted = TRUE, diag = FALSE)
  baseline <- as.character(FindClusters(snn, algorithm = 1,
    resolution = cfg$baseline_resolution, random.seed = 0, n.start = 10,
    n.iter = 10, verbose = FALSE)[, 1])
  write_tsv(connectivity(graph, baseline), 'louvain_connectivity', cfg)
  reference <- x[[]]$reference_cluster
  if (!is.null(reference)) write_tsv(data.frame(ARI = ari(baseline, reference),
    NMI = nmi(baseline, reference)), 'reconstruction_agreement', cfg)
  assignments <- data.frame(cell = colnames(x), baseline_louvain = baseline)
  evaluations <- connections <- pairs <- candidates <- list()
  seed_index <- match(cfg$seed, cfg$seeds)
  for (j in seq_along(cfg$resolutions)) {
    r <- cfg$resolutions[j]; message('Evaluating Leiden resolution = ', r)
    labs <- lapply(cfg$seeds, function(s) run_leiden(snn, r, s))
    candidates[[j]] <- labs[[seed_index]]
    for (i in seq_along(labs)) assignments[[paste0('candidate_', j, '_seed_', cfg$seeds[i])]] <- labs[[i]]
    pp <- combn(seq_along(labs), 2)
    agreement <- apply(pp, 2, function(i) ari(labs[[i[1]]], labs[[i[2]]]))
    pairs[[j]] <- data.frame(resolution = r, seed1 = cfg$seeds[pp[1, ]],
                            seed2 = cfg$seeds[pp[2, ]], ARI = agreement)
    cc <- do.call(rbind, lapply(seq_along(labs), function(i) {
      d <- connectivity(graph, labs[[i]]); d$resolution <- r; d$seed <- cfg$seeds[i]; d
    }))
    connections[[j]] <- cc
    lab <- candidates[[j]]; composition <- table(lab, evidence$cell_type)
    evaluations[[j]] <- data.frame(candidate = j, resolution = r, clusters = length(unique(lab)),
      reference_ARI = if (is.null(reference)) NA_real_ else ari(lab, reference),
      reference_NMI = if (is.null(reference)) NA_real_ else nmi(lab, reference),
      baseline_ARI = ari(lab, baseline), mean_seed_ARI = mean(agreement),
      min_seed_ARI = min(agreement), disconnected_across_seeds = sum(cc$components > 1),
      marker_coherence = sum(apply(composition, 1, max)) / sum(composition))
  }
  ev <- do.call(rbind, evaluations)
  ev$eligible <- ev$disconnected_across_seeds == 0 & ev$mean_seed_ARI >= cfg$min_seed_ari
  ev$selection_strategy <- cfg$selection; ev$selected <- FALSE
  write_tsv(assignments, 'candidate_assignments', cfg)
  write_tsv(do.call(rbind, pairs), 'seed_pairwise_ari', cfg)
  write_tsv(do.call(rbind, connections), 'candidate_connectivity', cfg)
  write_tsv(ev, 'resolution_evaluation', cfg)
  if (!any(ev$eligible)) stop('No candidate passed selection. Inspect the evaluation table; thresholds are not relaxed automatically.')
  eligible <- ev[ev$eligible, ]
  ranking <- order(-eligible$mean_seed_ARI, abs(eligible$resolution - cfg$preferred_resolution), eligible$resolution)
  chosen <- eligible[ranking[1], ]; lab <- candidates[[chosen$candidate]]
  ev$selected <- ev$candidate == chosen$candidate
  write_tsv(ev, 'resolution_evaluation', cfg)
  x <- FindClusters(x, graph.name = 'integrated_snn', algorithm = 4,
    leiden_method = 'leidenbase', resolution = chosen$resolution, random.seed = cfg$seed,
    n.iter = 10, group.singletons = FALSE, cluster.name = 'leiden_raw', verbose = FALSE)
  # Rerun and verify exact same-seed assignments; relabel clusters bijectively by size.
  stopifnot(identical(lab, separate_singletons(x$leiden_raw)),
            identical(lab, run_leiden(snn, chosen$resolution, cfg$seed)))
  sizes <- sort(table(lab), decreasing = TRUE)
  mapping <- setNames(as.character(seq_along(sizes) - 1), names(sizes))
  x$leiden_raw <- lab
  x$seurat_clusters <- factor(unname(mapping[lab]), levels = unname(mapping))
  Idents(x) <- 'seurat_clusters'
  cc <- connectivity(graph, as.character(Idents(x)))
  stopifnot(all(cc$components == 1))
  write_tsv(cc, 'final_connectivity', cfg)
  write_tsv(data.frame(raw_id = names(mapping), cluster = unname(mapping)), 'cluster_id_mapping', cfg)
  write_tsv(data.frame(test = c('same_seed_exact_assignments', 'all_clusters_connected'),
    passed = TRUE), 'reproducibility_checks', cfg)
  if (!is.null(reference)) {
    mat <- table(old = reference, new = x$seurat_clusters)
    write_tsv(as.data.frame(mat), 'reference_overlap', cfg)
    cross <- do.call(rbind, lapply(colnames(mat), function(k) {
      v <- mat[, k]; j <- which.max(v)
      data.frame(cluster = k, n = sum(v), reference_cluster = names(v)[j],
        intersection = v[j], purity = v[j]/sum(v), recall = v[j]/sum(mat[j, ]),
        Jaccard = v[j]/(sum(v) + sum(mat[j, ]) - v[j]))
    }))
    write_tsv(cross, 'cluster_crosswalk', cfg)
  }
  x@misc$leiden_reanalysis <- list(resolution = chosen$resolution, config = cfg,
    algorithm = 4, leiden_method = 'leidenbase', n.iter = 10, group.singletons = FALSE)
  saveRDS(x, checkpoint('selected_leiden', cfg), compress = FALSE)
  list(object = x, evaluation = ev, assignments = assignments)
}

export_plots <- function(x, result, cfg) {
  for (red in c('umap', 'tsne')) {
    for (field in intersect(c('seurat_clusters', 'group', 'sample',
                              'reference_cluster', 'celltype_provisional'), colnames(x[[]]))) {
      p <- DimPlot(x, reduction = red, group.by = field,
        label = field == 'seurat_clusters', repel = TRUE, raster = FALSE)
      save_plot(p, paste(red, field, sep = '_'), cfg, width = 11)
    }
    # Use configured groups; number output files to avoid path characters in group names.
    for (i in seq_along(unique(x$group))) {
      gr <- unique(x$group)[i]
      p <- DimPlot(subset(x, cells = colnames(x)[x$group == gr]), reduction = red,
        label = TRUE, raster = FALSE) + ggplot2::ggtitle(gr)
      save_plot(p, paste0(red, '_group_', i), cfg)
    }
    genes <- intersect(cfg$plot_genes, rownames(x))
    for (i in seq_along(genes)) save_plot(FeaturePlot(x, features = genes[i],
      reduction = red, split.by = 'group', raster = FALSE),
      paste0(red, '_feature_', i), cfg, width = 12)
  }
  save_plot(VlnPlot(x, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt',
    'percent.ribo'), group.by = 'sample', pt.size = 0, ncol = 2), 'qc_metrics', cfg, 12, 8)
  save_plot(ElbowPlot(x, ndims = cfg$npcs), 'pca_elbow', cfg)
  ev <- result$evaluation
  metric <- if (all(is.na(ev$reference_ARI))) 'mean_seed_ARI' else 'reference_ARI'
  d <- data.frame(resolution = ev$resolution, agreement = ev[[metric]], selected = ev$selected)
  save_plot(ggplot2::ggplot(d, ggplot2::aes(resolution, agreement)) +
    ggplot2::geom_line() + ggplot2::geom_point(ggplot2::aes(color = selected)) +
    ggplot2::theme_classic() + ggplot2::labs(y = metric), 'resolution_comparison', cfg)
  # Track cells between adjacent resolutions; candidate IDs avoid floating-point column-name collisions.
  nodes <- flows <- list()
  for (j in seq_len(nrow(ev))) {
    v <- result$assignments[[paste0('candidate_', j, '_seed_', cfg$seed)]]
    lev <- sort(unique(v))
    nodes[[j]] <- data.frame(resolution = ev$resolution[j], cluster = lev,
      y = seq_along(lev), n = as.numeric(table(factor(v, levels = lev))))
    if (j > 1) {
      a <- result$assignments[[paste0('candidate_', j - 1, '_seed_', cfg$seed)]]
      tt <- as.data.frame(table(a, v)); tt <- tt[tt$Freq > 0, ]
      flows[[j]] <- data.frame(x = ev$resolution[j - 1], xend = ev$resolution[j],
        y = match(tt$a, nodes[[j - 1]]$cluster), yend = match(tt$v, lev), n = tt$Freq)
    }
  }
  p <- ggplot2::ggplot()
  if (length(flows)) p <- p + ggplot2::geom_segment(data = do.call(rbind, flows),
    ggplot2::aes(x = x, xend = xend, y = y, yend = yend, linewidth = n), alpha = 0.25)
  p <- p + ggplot2::geom_point(data = do.call(rbind, nodes),
    ggplot2::aes(resolution, y, size = n), color = 'steelblue') + ggplot2::theme_classic()
  save_plot(p, 'resolution_transitions', cfg, 12, 9)
  if ('reference_cluster' %in% colnames(x[[]])) {
    ov <- as.data.frame(table(old = x$reference_cluster, new = x$seurat_clusters))
    save_plot(ggplot2::ggplot(ov, ggplot2::aes(new, old, fill = Freq)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(trans = 'sqrt') +
      ggplot2::theme_classic(), 'reference_overlap', cfg)
  }
  if (!is.null(cfg$reference_embeddings)) {
    emb <- read.csv(cfg$reference_embeddings, row.names = 1, check.names = FALSE)
    if (ncol(emb) < 2 || !all(colnames(x) %in% rownames(emb))) stop('Reference embeddings are missing cells or coordinate columns.')
    d <- data.frame(x = emb[colnames(x), 1], y = emb[colnames(x), 2], cluster = Idents(x))
    save_plot(ggplot2::ggplot(d, ggplot2::aes(x, y, color = cluster)) +
      ggplot2::geom_point(size = 0.15) + ggplot2::theme_classic() +
      ggplot2::ggtitle('Leiden labels on historical coordinates'), 'historical_tsne', cfg)
  }
  if (cfg$raw_tsne_diagnostic) {
    raw <- CreateSeuratObject(GetAssayData(x, assay = 'RNA', layer = 'counts'), meta.data = x[[]])
    raw <- NormalizeData(raw, verbose = FALSE)
    raw <- subset(raw, features = setdiff(rownames(raw), cfg$reporter_features))
    raw <- FindVariableFeatures(raw, nfeatures = cfg$nfeatures, verbose = FALSE)
    raw <- ScaleData(raw, verbose = FALSE)
    raw <- RunPCA(raw, npcs = cfg$npcs, seed.use = cfg$seed, verbose = FALSE)
    raw <- RunTSNE(raw, dims = cfg$anchor_dims, seed.use = cfg$tsne_seed,
      perplexity = cfg$tsne_perplexity, max_iter = 1000, verbose = FALSE)
    save_plot(DimPlot(raw, reduction = 'tsne', group.by = 'sample'), 'unintegrated_tsne', cfg)
  }
}

export_markers <- function(x, cfg) {
  clusters <- levels(Idents(x)); counts <- table(Idents(x)); exported <- list()
  # Include the input object and marker parameters in the cache key to prevent stale reuse.
  key_file <- tempfile(fileext = '.rds')
  saveRDS(list(input = unname(tools::md5sum(checkpoint('selected_leiden', cfg))),
    min_pct = cfg$marker_min_pct, logfc = cfg$marker_logfc,
    p_column = cfg$marker_p_column, p_cutoff = cfg$marker_p_cutoff,
    seed = cfg$seed, test = 'wilcox', seurat = as.character(packageVersion('Seurat'))), key_file)
  key <- unname(tools::md5sum(key_file)); unlink(key_file)
  for (assay in cfg$marker_assays) {
    DefaultAssay(x) <- assay
    cache <- file.path(cfg$output_dir, 'checkpoints', paste0('markers_', assay, '_', key))
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    # Cache each cluster; run sequentially for Windows compatibility without cell downsampling.
    results <- lapply(clusters, function(k) {
      f <- file.path(cache, paste0('cluster_', k, '.rds'))
      if (file.exists(f)) return(readRDS(f))
      if (counts[k] < 3 || sum(counts) - counts[k] < 3) return(NULL)
      message('Computing markers: ', assay, ' / cluster ', k)
      m <- FindMarkers(x, ident.1 = k, only.pos = TRUE, min.pct = cfg$marker_min_pct,
        logfc.threshold = cfg$marker_logfc, test.use = 'wilcox', random.seed = cfg$seed,
        verbose = FALSE)
      if (!nrow(m)) { saveRDS(NULL, f); return(NULL) }
      if (!cfg$marker_p_column %in% names(m)) stop('Marker results lack the configured P-value column.')
      m$gene <- rownames(m); m$cluster <- rep(k, nrow(m))
      m <- m[!is.na(m[[cfg$marker_p_column]]) & m[[cfg$marker_p_column]] < cfg$marker_p_cutoff, , drop = FALSE]
      saveRDS(m, f); m
    })
    m <- do.call(rbind, results)
    if (is.null(m)) m <- data.frame(p_val = numeric(), avg_log2FC = numeric(),
      pct.1 = numeric(), pct.2 = numeric(), p_val_adj = numeric(), gene = character(), cluster = character())
    rownames(m) <- NULL
    write_tsv(m, paste0('markers_', assay, '_all'), cfg)
    for (n in cfg$marker_top_n) {
      tops <- if (!nrow(m)) m else do.call(rbind, lapply(split(m, m$cluster), function(z) {
        z <- z[order(-z$avg_log2FC, z$gene), , drop = FALSE]
        # Preserve ties at rank N, as in the original scripts; output may exceed N rows.
        if (nrow(z) > n) z <- z[z$avg_log2FC >= z$avg_log2FC[n], , drop = FALSE]
        z
      }))
      write_tsv(tops, paste0('markers_', assay, '_top', n), cfg)
    }
    exported[[assay]] <- m
  }
  # Optional integrated-marker comparison only; biological interpretation uses RNA markers.
  if (!is.null(cfg$reference_markers) && 'reference_cluster' %in% colnames(x[[]]) &&
      'integrated' %in% names(exported)) {
    old <- read.delim(cfg$reference_markers)
    if (!all(c('cluster', 'gene', 'avg_log2FC') %in% names(old))) stop('The reference marker table lacks required columns.')
    cross <- read.delim(file.path(cfg$output_dir, 'tables', 'cluster_crosswalk.tsv'), colClasses = 'character')
    overlap <- do.call(rbind, lapply(seq_len(nrow(cross)), function(i) {
      o <- old[old$cluster == cross$reference_cluster[i], ]
      m <- exported$integrated; m <- m[m$cluster == cross$cluster[i], ]
      og <- head(o$gene[order(-o$avg_log2FC)], 20); ng <- head(m$gene[order(-m$avg_log2FC)], 20)
      shared <- intersect(og, ng); total <- union(og, ng)
      data.frame(cluster = cross$cluster[i], reference_cluster = cross$reference_cluster[i],
        overlap = length(shared), Jaccard = if (length(total)) length(shared)/length(total) else NA_real_,
        shared_genes = paste(shared, collapse = ','))
    }))
    write_tsv(overlap, 'reference_marker_overlap', cfg)
  }
}

export_results <- function(result, evidence, cfg) {
  x <- result$object; DefaultAssay(x) <- 'RNA'; Idents(x) <- 'seurat_clusters'
  stopifnot(identical(rownames(evidence$scores), colnames(x)))
  annotation <- do.call(rbind, lapply(levels(Idents(x)), function(k) {
    i <- which(as.character(Idents(x)) == k)
    scores <- colMeans(evidence$scores[i, , drop = FALSE]); o <- order(scores, decreasing = TRUE)
    data.frame(cluster = k, n = length(i), marker_candidate = names(scores)[o[1]],
      second_candidate = names(scores)[o[2]], score = unname(scores[o[1]]),
      margin = unname(scores[o[1]] - scores[o[2]]),
      cell_marker_agreement = mean(evidence$cell_type[i] == names(scores)[o[1]]))
  }))
  annotation$annotation <- ifelse(annotation$n < cfg$annotation_min_cells, 'Unresolved_small_cluster',
    ifelse(annotation$margin < cfg$annotation_margin | annotation$cell_marker_agreement < cfg$annotation_agreement,
      paste0('Ambiguous_', annotation$marker_candidate), annotation$marker_candidate))
  annotation$status <- 'marker_based_provisional'
  x$celltype_provisional <- annotation$annotation[match(as.character(Idents(x)), annotation$cluster)]
  write_tsv(annotation, 'annotation_evidence', cfg)
  write_tsv(data.frame(cell = colnames(x), x[[]]), 'cell_metadata', cfg)
  prop <- as.data.frame(table(cluster = Idents(x), sample = x$sample)); names(prop)[3] <- 'cells'
  write_tsv(prop, 'cluster_by_sample_counts', cfg)
  prop$fraction <- prop$cells / as.numeric(table(x$sample)[as.character(prop$sample)])
  write_tsv(prop, 'cluster_by_sample_fractions', cfg)
  write_tsv(as.data.frame(table(celltype = x$celltype_provisional, sample = x$sample)), 'celltype_by_sample_counts', cfg)
  if (cfg$make_plots) export_plots(x, result, cfg)
  export_markers(x, cfg)
  saveRDS(x, file.path(cfg$output_dir, 'seurat_final.rds'))
  invisible(x)
}

main <- function(cfg) {
  required <- c('Seurat', 'Matrix', 'igraph', 'ggplot2', 'future', 'leidenbase')
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install the required packages first: ', paste(missing, collapse = ', '))
  suppressPackageStartupMessages(library(Seurat))
  if (!'leiden_method' %in% names(formals(getS3method('FindClusters', 'default'))))
    stop('This Seurat version does not support leiden_method. Use a version with the leidenbase backend.')
  stopifnot(cfg$species %in% c('mouse', 'human'), identical(cfg$selection, 'strict'),
    cfg$cell_naming %in% c('sample_prefix', 'index_suffix'),
    all(c('sample', 'group', 'path') %in% names(cfg$samples)), nrow(cfg$samples) >= 2,
    !anyNA(cfg$samples), !anyDuplicated(cfg$samples$sample),
    all(nzchar(cfg$samples$sample)), all(nzchar(cfg$samples$group)),
    all(dir.exists(cfg$samples$path)), length(unique(cfg$seeds)) >= 2,
    !anyDuplicated(cfg$seeds), cfg$seed %in% cfg$seeds,
    length(cfg$resolutions) >= 1, all(is.finite(cfg$resolutions)), all(cfg$resolutions > 0),
    !anyDuplicated(cfg$resolutions), length(cfg$dims) > 0, length(cfg$anchor_dims) > 0,
    all(cfg$dims >= 1), all(cfg$anchor_dims >= 1), max(cfg$dims) <= cfg$npcs,
    all(cfg$marker_assays %in% c('RNA', 'integrated')),
    cfg$marker_p_column %in% c('p_val', 'p_val_adj'), all(cfg$marker_top_n >= 1))
  if (cfg$use_archived_features && is.null(cfg$archive_rds)) stop('Reusing archived features requires archive_rds.')
  if (cfg$require_reference_qc_match && is.null(cfg$reference_clusters)) stop('Reference QC verification requires reference_clusters.')
  for (path in list(cfg$archive_rds, cfg$reference_clusters, cfg$reference_embeddings, cfg$reference_markers))
    if (!is.null(path) && !file.exists(path)) stop('Input file not found: ', path)
  # Use v3 assays to match the original reconstruction and avoid unjoined v5 layers.
  options(Seurat.object.assay.version = 'v3', future.globals.maxSize = 16 * 1024^3)
  future::plan('sequential'); set.seed(cfg$seed)
  cfg$resolutions <- sort(cfg$resolutions)
  for (d in c('tables', 'figures', 'checkpoints', 'logs'))
    dir.create(file.path(cfg$output_dir, d), recursive = TRUE, showWarnings = FALSE)
  dput(cfg, file.path(cfg$output_dir, 'logs', 'config.R'))
  on.exit(writeLines(capture.output(sessionInfo()), file.path(cfg$output_dir, 'logs', 'session_info.txt')), add = TRUE)
  reference <- read_reference(cfg)
  x <- rebuild(cfg, reference)
  evidence <- marker_evidence(x, cfg)
  result <- select_clusters(x, evidence, cfg)
  export_results(result, evidence, cfg)
  message('Analysis complete. Results: ', normalizePath(cfg$output_dir))
}



# ---------------- Revision modules ----------------
# Historical manuscript IDs are used ONLY with their archived cell identities.
# They must never be transferred by numeric ID to a newly clustered object.
cfg$reporter_features <- c('gene-GFP','gene-RFP','GFP','eGFP','EGFP','RFP','tdTomato')
revision_write <- function(d, name, out) {
  write.table(d, file.path(out,paste0(name,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE,na='NA')
}
revision_plot <- function(p,name,out,w=12,h=7) {
  ggplot2::ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,limitsize=FALSE)
}
# Descriptive cluster markers; cells are not independent mouse replicates.
# Retain all returned tests, including nonsignificant rows, plus an adjusted-P top table.
revision_markers <- function(x, field, prefix, out) {
  DefaultAssay(x)<-'RNA';Idents(x)<-field
  all<-lapply(levels(Idents(x)),function(k){
    message(prefix,' markers ',k)
    if(sum(Idents(x)==k)<3 || sum(Idents(x)!=k)<3) return(NULL)
    m<-FindMarkers(x,ident.1=k,only.pos=TRUE,min.pct=.25,logfc.threshold=.25,
      test.use='wilcox',random.seed=42,verbose=FALSE)
    m$gene<-rownames(m);m$cluster<-rep(k,nrow(m));m
  })
  m<-do.call(rbind,all);revision_write(m,paste0(prefix,'_markers_all_tests'),out)
  ok<-m[!is.na(m$p_val_adj)&m$p_val_adj<.05,,drop=FALSE]
  top<-do.call(rbind,lapply(split(ok,ok$cluster),function(z)head(z[order(-z$avg_log2FC,z$gene),],10)))
  revision_write(top,paste0(prefix,'_markers_top10_Bonferroni05'),out)
  genes<-unique(top$gene)
  if(length(genes)) {
    # Split the large panel into readable pages, retaining original cluster labels.
    for(i in seq_along(split(genes,ceiling(seq_along(genes)/30)))) {
      g<-split(genes,ceiling(seq_along(genes)/30))[[i]]
      revision_plot(DotPlot(x,features=g,assay='RNA',group.by=field)+RotatedAxis(),
        paste0(prefix,'_DE_marker_panel_',i),out,14,8)
    }
  }
  invisible(m)
}
revision_qc <- function(x,field,prefix,out) {
  cts<-GetAssayData(x,assay='RNA',layer='counts')
  mt<-grep('^mt-',rownames(cts),value=TRUE)
  if(!length(mt))stop('No mouse mitochondrial genes found.')
  d<-data.frame(cell=colnames(x),cluster=as.character(x[[field]][,1]),
    sample=x$Sample,group=x$Group,nCount=Matrix::colSums(cts),nFeature=Matrix::colSums(cts>0),
    percent_mt=100*Matrix::colSums(cts[mt,,drop=FALSE])/Matrix::colSums(cts))
  d$corrected_qc_pass<-d$nFeature>250 & d$nFeature<7500 & d$nCount>100 & d$percent_mt<25
  revision_write(d,paste0(prefix,'_cell_QC'),out)
  sets<-split(seq_len(nrow(d)),d$cluster)
  summary<-do.call(rbind,lapply(names(sets),function(k){i<-sets[[k]];z<-d[i,];n<-table(factor(z$sample,levels=unique(d$sample)))
    data.frame(cluster=k,n=nrow(z),median_UMI=median(z$nCount),median_genes=median(z$nFeature),
      median_mt=median(z$percent_mt),mt_ge25=sum(z$percent_mt>=25),qc_fail=sum(!z$corrected_qc_pass),
      samples_present=sum(n>0),samples_ge10=sum(n>=10),max_sample_fraction=max(n)/sum(n))}))
  revision_write(summary,paste0(prefix,'_cluster_QC'),out)
  tab<-as.data.frame(table(cluster=d$cluster,sample=d$sample));names(tab)[3]<-'cells'
  tab$fraction_of_cluster<-tab$cells/as.numeric(table(d$cluster)[as.character(tab$cluster)])
  tab$fraction_of_sample<-tab$cells/as.numeric(table(d$sample)[as.character(tab$sample)])
  revision_write(tab,paste0(prefix,'_cluster_sample_contribution'),out)
  for(v in c('nCount','nFeature','percent_mt')) revision_plot(
    ggplot2::ggplot(d,ggplot2::aes(x=factor(cluster,levels=sort(unique(as.integer(cluster)))),y=.data[[v]]))+
      ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::facet_wrap(~sample,ncol=2)+
      ggplot2::theme_bw()+ggplot2::labs(x='Original cluster',y=v),paste0(prefix,'_QC_',v),out,13,8)
  invisible(d)
}
# Rebuild SMC PCA and graph on the historically selected cells; do not gate on reporter.
# Historical state IDs may be reapplied only after exact cell-wise cluster agreement.
revision_smc_recluster <- function(x, archived, out) {
  s<-subset(x,cells=colnames(x)[x$celltype_customer=='Smooth muscle cell'])
  stopifnot(setequal(colnames(s),colnames(archived)))
  DefaultAssay(s)<-'integrated';set.seed(42)
  s<-RunPCA(s,npcs=10,seed.use=42,verbose=FALSE)
  s<-FindNeighbors(s,dims=1:10,k.param=20,prune.SNN=1/15,nn.method='annoy',n.trees=50,
    annoy.metric='euclidean',l2.norm=FALSE,verbose=FALSE)
  lab<-run_leiden(s[['integrated_snn']],1,42)
  lev<-names(sort(table(lab),decreasing=TRUE));mapped<-unname(setNames(as.character(seq_along(lev)-1),lev)[lab])
  truth<-as.character(archived$smc_leiden)[match(colnames(s),colnames(archived))]
  revision_write(data.frame(cell=colnames(s),recomputed=mapped,archived=truth),'SMC_recluster_crosscheck',out)
  revision_write(data.frame(ARI=ari(mapped,truth),exact=identical(mapped,truth)), 'SMC_recluster_validation',out)
  if(!identical(mapped,truth))stop('SMC rerun differs: do not transplant historical state labels.')
  invisible(s)
}
# PC/k sensitivity uses no marker score or labels in graph construction or clustering.
# Comparison to historical clusters is a diagnostic, not an optimization objective.
revision_sensitivity <- function(x,field,prefix,pcs,ks,resolutions,out) {
  DefaultAssay(x)<-'integrated';truth<-as.character(x[[field]][,1]);ans<-maps<-list();i<-0L
  for(pc in pcs)for(k in ks){
    set.seed(42);z<-FindNeighbors(x,dims=seq_len(pc),k.param=k,prune.SNN=1/15,
      nn.method='annoy',n.trees=50,annoy.metric='euclidean',l2.norm=FALSE,verbose=FALSE)
    for(r in resolutions){i<-i+1L;message(prefix,' sensitivity ',pc,' ',k,' ',r)
      labs<-lapply(42:46,function(seed)run_leiden(z[['integrated_snn']],r,seed))
      pairs<-combn(1:5,2);av<-apply(pairs,2,function(j)ari(labs[[j[1]]],labs[[j[2]]]))
      lab<-labs[[1]];tb<-table(truth,lab)
      maps[[i]]<-do.call(rbind,lapply(rownames(tb),function(old){v<-tb[old,];best<-which.max(v)
        data.frame(PC=pc,k=k,resolution=r,original_cluster=old,matched_cluster=colnames(tb)[best],
          Jaccard=v[best]/(sum(v)+sum(tb[,best])-v[best]),recall=v[best]/sum(v))}))
      ans[[i]]<-data.frame(PC=pc,k=k,resolution=r,clusters=length(unique(lab)),
        ARI_to_manuscript=ari(lab,truth),mean_seed_ARI=mean(av),min_seed_ARI=min(av))
      revision_write(do.call(rbind,ans),paste0(prefix,'_parameter_sensitivity'),out)
      revision_write(do.call(rbind,maps),paste0(prefix,'_cluster_sensitivity'),out)
    }
  }
}
revision_reporter <- function(x,s,out) {
  cts<-GetAssayData(x,assay='RNA',layer='counts');stopifnot('gene-GFP'%in%rownames(cts))
  d<-data.frame(cell=colnames(x),sample=x$Sample,group=x$Group,global_cluster=x$leiden_global,
    population=as.character(x$celltype_customer),GFP=as.numeric(cts['gene-GFP',]),
    UMI=Matrix::colSums(cts),genes=Matrix::colSums(cts>0))
  ix<-match(colnames(s),d$cell);stopifnot(!anyNA(ix));d$population[ix]<-as.character(s$celltype_smc)
  d$GFP_per10000<-10000*d$GFP/d$UMI
  revision_write(d,'reporter_per_cell',out)
  pop<-c('Contractile SMC','Fibroblast like SMC1','Fibroblast like SMC2','Fibroblast like SMC3','fibroblast')
  z<-d[d$population%in%pop,];z$population<-factor(z$population,levels=pop)
  summary<-do.call(rbind,lapply(split(z,interaction(z$sample,z$population,drop=TRUE)),function(a){
    data.frame(sample=a$sample[1],group=a$group[1],population=as.character(a$population[1]),n=nrow(a),
      median_GFP=median(a$GFP),mean_GFP=mean(a$GFP),q90_GFP=unname(quantile(a$GFP,.9)),
      pct_ge1=100*mean(a$GFP>=1),pct_ge2=100*mean(a$GFP>=2),pct_ge5=100*mean(a$GFP>=5),
      median_UMI=median(a$UMI),median_genes=median(a$genes),
      rho_GFP_UMI=suppressWarnings(cor(a$GFP,a$UMI,method='spearman')),
      rho_GFP_genes=suppressWarnings(cor(a$GFP,a$genes,method='spearman')))
  }))
  revision_write(summary,'reporter_by_sample_population',out)
  # Thresholds are descriptive sensitivity checks, not validated lineage gates.
  revision_plot(ggplot2::ggplot(z,ggplot2::aes(population,log1p(GFP),fill=population))+
    ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::facet_wrap(~sample,ncol=2)+
    ggplot2::theme_bw()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=30,hjust=1),legend.position='none')+
    ggplot2::labs(y='log(1 + raw GFP UMI)',x=NULL),'GFP_raw_distribution',out,14,9)
  for(v in c('UMI','genes')) revision_plot(ggplot2::ggplot(z,ggplot2::aes(x=.data[[v]],y=GFP,color=population))+
    ggplot2::geom_point(size=.3,alpha=.2)+ggplot2::scale_x_log10()+ggplot2::scale_y_continuous(trans='log1p')+
    ggplot2::facet_wrap(~sample,ncol=2)+ggplot2::theme_bw(),paste0('GFP_vs_',v),out,13,8)
}
# Conservative analysis labels are separate from historical display labels.
# Exclude the erythroid-enriched/mixed global cluster 14, do not silently rename it SMC.
# This fixed-label sensitivity is NOT a substitute for re-integration/reclustering after QC.
revision_review <- function(x,s,out) {
  q<-read.delim(file.path(out,'SMC_cell_QC.tsv'))
  pass<-q$corrected_qc_pass[match(colnames(s),q$cell)]
  d<-data.frame(cell=colnames(s),sample=s$Sample,group=s$Group,global_cluster=s$leiden_global,
    smc_cluster=s$smc_leiden,historical_state=s$celltype_smc,corrected_qc_pass=pass)
  d$exclude_erythroid_mixed<-as.character(d$global_cluster)=='14'
  d$analysis_label<-as.character(d$historical_state)
  d$analysis_label[d$exclude_erythroid_mixed]<-'Excluded_erythroid_enriched_or_mixed'
  d$analysis_label[!d$exclude_erythroid_mixed & as.character(d$smc_cluster)=='4']<-'SMC_high_mt_state_unresolved'
  d$analysis_label[as.character(d$smc_cluster)=='15']<-'Matrix_like_low_complexity_candidate'
  d$analysis_label[as.character(d$smc_cluster)=='17']<-'Cycling_matrix_SMC_like_candidate'
  revision_write(d,'SMC_annotation_review_and_exclusions',out)
  modes<-list(historical=rep(TRUE,nrow(d)),exclude_global14=!d$exclude_erythroid_mixed,
    corrected_QC_and_exclude14=pass & !d$exclude_erythroid_mixed,
    conservative_exclude_low_complexity=pass & !d$exclude_erythroid_mixed & !as.character(d$smc_cluster)%in%c('3','15'))
  a<-do.call(rbind,lapply(names(modes),function(mode){
    z<-d[modes[[mode]],];tb<-as.data.frame(table(sample=factor(z$sample,levels=unique(d$sample)),
      historical_state=factor(z$historical_state,levels=unique(d$historical_state))))
    names(tb)[3]<-'cells';tb$mode<-mode
    tb$fraction_within_retained_SMC<-tb$cells/as.numeric(table(z$sample)[as.character(tb$sample)])
    tb
  }))
  revision_write(a,'SMC_fixed_label_exclusion_sensitivity',out)
  cross<-as.data.frame(table(global_cluster=s$leiden_global,SMC_cluster=s$smc_leiden))
  revision_write(cross,'global_to_SMC_crosswalk',out)
}


# Final manual maps are distinct from automated provisional marker scores.
# These preserve the manuscript labels for reproducibility, not evidence that every
# historical state is equally secure. Independent DE/QC diagnostics remain visible.
revision_annotation <- function(x,s,out) {
  map<-unique(data.frame(cluster=as.character(s$smc_leiden),state=as.character(s$celltype_smc)))
  map<-map[order(as.integer(map$cluster)),]
  panels<-list('Contractile SMC'=c('Cnn1','Myh11','Tagln','Tpm1'),
    'Proliferating SMC'=c('Ccnl2','Ccnt2','Malat1'),
    'Synthetic SMC'=c('Mfap5','Vcan','Gsn'),
    'Fibroblast like SMC1'=c('Dcn','Serpinh1','Lum'),
    'Fibroblast like SMC2'=c('Dcn','Serpinh1','Lum'),
    'Fibroblast like SMC3'=c('Dcn','Serpinh1','Lum'))
  map$S9F_panel<-vapply(map$state,function(k)paste(panels[[k]],collapse=';'),character(1))
  map$annotation_method<-'manual state map; marker expression with eGFP context'
  map$additional_validation<-'RNA cluster DE, QC and sample contribution; Mki67/Top2a/Ube2c for cycling'
  revision_write(map,'SMC_state_annotation_criteria',out)
  selected<-as.character(x$celltype_customer)=='Smooth muscle cell'
  stopifnot(setequal(colnames(x)[selected],colnames(s)))
  revision_write(data.frame(cell=colnames(x),global_cluster=x$leiden_global,
    final_annotation=x$celltype_customer,selected_for_SMC=selected,
    selection_rule='final SMC annotation; no individual GFP UMI or normalized-expression gate'),
    'SMC_selection_rule',out)
  # Record the source figure separately from the quantitative reporter analysis.
  revision_write(data.frame(item=c('S9C source','reporter feature','quantification','positivity gate',
      'annotation','historical resolution selection','revised selection'),
    detail=c('240624_new_result/1/tsne_GFP.pdf','gene-GFP in RNA assay',
      'custom FASTA/GTF + Cell Ranger; version 7.2.0 reported in Supplemental Methods',
      'none for SMC extraction; raw-count thresholds are descriptive sensitivity checks',
      'marker scores are provisional; final cluster maps reflect manual review with eGFP and markers',
      'connectivity and marker-coherence eligibility; reference ARI ranking in similarity mode',
      'graph connectivity and seed stability only; markers and reference labels are diagnostics')),
    'method_provenance',out)
}
revision_S9C <- function(x,out) {
  DefaultAssay(x)<-'RNA'
  stopifnot('gene-GFP'%in%rownames(x[['RNA']]))
  # Continuous log-normalized expression, not binary lineage classification.
  revision_plot(FeaturePlot(x,features='gene-GFP',reduction='tsne',split.by='Group',
    raster=FALSE),'S9C_eGFP_continuous_expression',out,12,6)
}

revision_audit <- function(root,out,full=TRUE) {
  for(p in c('Seurat','Matrix','igraph','ggplot2','leidenbase'))if(!requireNamespace(p,quietly=TRUE))stop('Missing package ',p)
  suppressPackageStartupMessages(library(Seurat));future::plan('sequential');set.seed(42)
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  on.exit(writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt')),add=TRUE)
  x<-readRDS(file.path(root,'240630_new/3/scRNA_customer_Leiden.rds'))
  s<-readRDS(file.path(root,'240724_new/2/scRNA_smc_240724_new.rds'))
  stopifnot(ncol(x)==27592,ncol(s)==19946,length(unique(x$leiden_global))==23,length(unique(s$smc_leiden))==18)
  # Provenance, not a new assertion of biological validation.
  gm<-read.delim(file.path(root,'240630_new/source/customer_mapping.tsv'))
  sm<-read.delim(file.path(root,'240724_new/tables/customer_SMC_mapping.tsv'))
  stopifnot(all(as.character(x$celltype_customer)==gm$customer_celltype[match(x$leiden_global,gm$cluster)]),
    all(as.character(s$celltype_smc)==sm$celltype_smc[match(s$smc_leiden,sm$smc_leiden)]))
  revision_write(gm,'global_annotation_mapping',out);revision_write(sm,'SMC_annotation_mapping',out)
  revision_write(data.frame(cell=colnames(x),original_cluster=x$leiden_global,
    provisional=x$celltype_provisional,reviewed=x$celltype_reviewed,manuscript=x$celltype_customer),'annotation_provenance',out)
  revision_qc(x,'leiden_global','global',out);revision_qc(s,'smc_leiden','SMC',out)
  revision_reporter(x,s,out)
  revision_annotation(x,s,out)
  revision_S9C(x,out)
  revision_review(x,s,out)
  features<-list(integration=x@commands$FindIntegrationAnchors@params$anchor.features,
    global_PCA=rownames(Loadings(x[['pca']])),SMC_PCA=rownames(Loadings(s[['pca']])))
  audit<-do.call(rbind,lapply(names(features),function(n)data.frame(stage=n,reporter=cfg$reporter_features,
    included=cfg$reporter_features%in%features[[n]])))
  revision_write(audit,'reporter_feature_inclusion',out)
  ev<-read.delim(file.path(root,'240624_new/tables/resolution_evaluation.tsv'))
  ev$marker_free_eligible<-ev$disconnected_across_seeds==0 & ev$mean_seed_ARI>=.85
  ev$marker_free_selected<-FALSE
  good<-which(ev$marker_free_eligible)
  if(length(good)){j<-good[order(-ev$mean_seed_ARI[good],abs(ev$resolution[good]-.7))][1];ev$marker_free_selected[j]<-TRUE}
  revision_write(ev,'resolution_selection_audit',out)
  DefaultAssay(x)<-'RNA';DefaultAssay(s)<-'RNA'
  gp<-c('Myh11','Cnn1','Tagln','Acta2','Tpm1','Dcn','Timp1','Lum','Pdgfra','Col1a1','Dpt',
    'Pecam1','Cdh5','Adgre1','Cd68','Cd3d','Cd3e','Cd79a','Cd19','Plp1','Sox10','Flt3','Itgax',
    'Hbb-bs','Hba-a1','Alas2','Jchain','Mzb1','Il1rl1','Gata3','gene-GFP','gene-RFP')
  sp<-c('Myh11','Cnn1','Tagln','Tpm1','Dcn','Timp1','Serpinh1','Lum','Mfap5','Vcan','Gsn',
    'Ccnl2','Ccnt2','Malat1','Mki67','Top2a','Ube2c','Pdgfra','gene-GFP','gene-RFP')
  for(level in c('global','SMC')) {
    o<-if(level=='global')x else s;field<-if(level=='global')'leiden_global' else 'smc_leiden'
    genes<-intersect(if(level=='global')gp else sp,rownames(o[['RNA']]))
    p<-DotPlot(o,features=genes,group.by=field,assay='RNA')+RotatedAxis()
    revision_plot(p,paste0(level,'_original_cluster_markers'),out,16,9)
    revision_write(p$data,paste0(level,'_original_cluster_markers'),out)
  }
  # Reuse already-computed global RNA tests; provenance is recorded below.
  m<-read.delim(file.path(root,'240624_new/tables/markers_RNA_all.tsv'))
  revision_write(m,'global_archived_RNA_markers_p001',out)
  ok<-m[!is.na(m$p_val_adj)&m$p_val_adj<.05,]
  top<-do.call(rbind,lapply(split(ok,ok$cluster),function(z)head(z[order(-z$avg_log2FC,z$gene),],10)))
  revision_write(top,'global_markers_top10_Bonferroni05',out)
  if(full) {
  revision_markers(s,'smc_leiden','SMC',out)
  revision_smc_recluster(x,s,out)
  revision_sensitivity(x,'leiden_global','global',c(20,30,40),c(15,20,30),c(.6,.8,1),out)
  revision_sensitivity(s,'smc_leiden','SMC',c(5,10),c(15,20,30),c(.8,1,1.2),out)
  }
  writeLines('Completed archived-object evaluation; no ambient correction or new lineage validation claimed.',file.path(out,'COMPLETE.txt'))
}
# Corrected end-to-end discovery mode, separate from the manuscript-object audit.
# New annotation/SMC inclusion requires cell-level review; historical cluster IDs are not reused.
if (sys.nframe()==0L) {
  args<-commandArgs(trailingOnly=TRUE)
  positional<-args[!grepl('^--',args)]
  flags<-args[grepl('^--',args)]
  if(any(!flags%in%c('--quick','--rebuild')))stop('Valid flags: --quick, --rebuild')
  f<-sub('^--file=','',grep('^--file=',commandArgs(FALSE),value=TRUE))
  root<-normalizePath(if(length(positional))positional[1] else
    file.path(dirname(normalizePath(f)),'../..'),mustWork=TRUE)
  local_lib<-file.path(root,'240724_new/R_library')
  if(dir.exists(local_lib)).libPaths(c(local_lib,.libPaths()))
  out<-file.path(root,'results_TMEM41B_review_20261022')
  if('--rebuild'%in%flags){
    cfg$samples<-data.frame(sample=c('HMY_KO1','HMY_KO2','HMY_WT1','HMY_WT2'),group=c('KO','KO','WT','WT'))
    cfg$samples$path<-file.path(root,'240624',cfg$samples$sample)
    cfg$output_dir<-file.path(out,'reporter_excluded_discovery');cfg$reference_clusters<-NULL
    cfg$archive_rds<-NULL;cfg$use_archived_features<-FALSE;cfg$selection<-'strict'
    # New cluster IDs receive provisional labels only. Manual annotation is required
    # before applying a final SMC subset to a newly reconstructed clustering solution.
    main(cfg)
  }else revision_audit(root,out,full=!'--quick'%in%flags)
}
