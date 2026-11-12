## dispersionDiagnostics.R
## Diagnostics for whether one DESeq2 dispersion per gene, estimated from all libraries (design ~ Condition),
## is appropriate for a data set, or whether conditions differ in replicate variability enough that a
## comparison should be fitted on its own libraries (or with group-specific weights).
##
## Background: with all conditions in one model, DESeq2 estimates a single dispersion per gene from every
## library and uses it in every contrast. If the replicates of some conditions vary much more than others,
## contrasts involving the noisy conditions are tested against too small a variance and contrasts among the
## quiet conditions against too large a one. Cook's distance does not detect this: it looks for one library
## that disagrees with its group, not for a group whose libraries are all spread out.
##
## Three diagnostics, each a plot plus a table:
##   1. conditionVariance()  replicate variance of each condition relative to the data-set average
##                           (limma::voomWithQualityWeights with var.group = condition; Figure 1A of the
##                           Rizvi DE method report). 1 = average; 2 = twice the average variance.
##   2. libraryWeights()     the same model with one weight per library: separates one aberrant library
##                           (a low weight in an otherwise normal condition) from a whole noisy condition.
##   3. dispersionByCondition()  each condition's own gene-wise DESeq2 dispersion (fitted on that condition's
##                           libraries only) relative to the pooled dispersion that ~ Condition would use.
##
## How to read them: if conditions sit within about 0.5-2 x the average (plot 1) and the per-condition
## dispersion ratios are centred near the same value (plot 3), the pooled dispersion is a reasonable summary.
## Conditions well outside that range mean contrasts involving them should be examined: refit the comparison
## on its own conditions (dds[, Condition %in% c(A, B)]) and compare the calls, or use limma-voom with
## var.group. The 0.5-2 band is a rule of thumb for "look in more detail", not a test.
##
## Caveats: conditions with fewer than min_rep libraries are left out. Gene-wise dispersions from three
## libraries are noisy and biased low, so in plot 3 compare conditions with the same number of libraries
## with each other rather than reading the absolute ratio.
##
## Requires: DESeq2, limma, edgeR, ggplot2. Written 2026-10-29 for UPSCb-common (src/R).

suppressPackageStartupMessages({library(DESeq2); library(limma); library(edgeR); library(ggplot2)})

.ddCond <- function(dds, condition, min_rep) {
  grp <- droplevels(factor(colData(dds)[[condition]]))
  n <- table(grp)
  keep <- grp %in% names(n)[n >= min_rep]
  if (any(!keep)) message("Left out (fewer than ", min_rep, " libraries): ",
                          paste(names(n)[n < min_rep], collapse = ", "))
  list(dds = dds[, keep], grp = droplevels(grp[keep]))
}

.ddVoom <- function(dds, grp, var.group) {
  cnt <- counts(dds)
  y <- calcNormFactors(DGEList(cnt, group = grp, lib.size = colSums(cnt)))
  y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = TRUE]
  g <- grp; levels(g) <- make.names(levels(g))
  design <- model.matrix(~ 0 + g)
  voomWithQualityWeights(y, design, var.group = if (var.group) g else NULL)
}

## Long condition names: draw horizontally so the names are readable.
.ddLayout <- function(p, conds, ylim = NULL, long = 15) {
  p <- p + theme(plot.title.position = "plot")
  if (max(nchar(as.character(conds))) > long) {
    p + coord_flip(ylim = ylim) + scale_x_discrete(limits = rev(levels(factor(conds)))) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  } else if (!is.null(ylim)) p + coord_cartesian(ylim = ylim) else p
}

## 1. Replicate variance per condition ------------------------------------------------------------------
conditionVariance <- function(dds, condition = "Condition", min_rep = 2, colour_by = NULL,
                              band = c(0.5, 2), title = NULL) {
  x <- .ddCond(dds, condition, min_rep)
  v <- .ddVoom(x$dds, x$grp, var.group = TRUE)
  w <- tapply(v$targets$sample.weights, x$grp, mean)
  tab <- data.frame(condition = factor(names(w), levels = levels(x$grp)),
                    libraries = as.integer(table(x$grp)[names(w)]),
                    weight = as.numeric(w), relative_variance = 1 / as.numeric(w))
  tab$outside_band <- tab$relative_variance < band[1] | tab$relative_variance > band[2]
  if (!is.null(colour_by)) tab$colour <- factor(colData(x$dds)[[colour_by]][match(tab$condition, x$grp)])
  rng <- max(tab$relative_variance) / min(tab$relative_variance)
  if (is.null(title)) title <- sprintf("Replicate variance differs %.1f-fold between conditions", rng)
  p <- ggplot(tab, aes(condition, relative_variance)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = band[1], ymax = band[2], fill = "grey92") +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_segment(aes(xend = condition, y = 1, yend = relative_variance), linewidth = 0.4) +
    geom_point(size = 1.8) +
    scale_y_log10(breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), labels = c("1/8", "1/4", "1/2", "1", "2", "4", "8")) +
    labs(x = NULL, y = "replicate variance relative to average", title = title,
         subtitle = "voomWithQualityWeights, one weight per condition; grey band = 0.5-2") +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(colour_by)) p <- p + aes(colour = colour) + labs(colour = colour_by)
  p <- .ddLayout(p, tab$condition)
  list(table = tab, plot = p, fold_range = rng)
}

## 2. Per-library weights ------------------------------------------------------------------------------
libraryWeights <- function(dds, condition = "Condition", min_rep = 2, label_below = 0.5) {
  x <- .ddCond(dds, condition, min_rep)
  v <- .ddVoom(x$dds, x$grp, var.group = FALSE)
  tab <- data.frame(library = colnames(x$dds), condition = x$grp, weight = v$targets$sample.weights)
  tab$condition_median <- ave(tab$weight, tab$condition, FUN = median)
  p <- ggplot(tab, aes(condition, weight)) +
    geom_hline(yintercept = 1, colour = "grey55", linewidth = 0.3) +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = 1), size = 1.4, colour = "grey30") +
    (if (requireNamespace("ggrepel", quietly = TRUE))
       ggrepel::geom_text_repel(data = subset(tab, weight < label_below), aes(label = library), size = 2.3,
                                colour = "#D55E00", min.segment.length = 0, seed = 1)
     else geom_text(data = subset(tab, weight < label_below), aes(label = library), size = 2.3,
                    hjust = -0.2, colour = "#D55E00")) +
    scale_y_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2, 4), labels = c("0.05", "0.1", "0.2", "0.5", "1", "2", "4")) +
    labs(x = NULL, y = "library weight (log scale; 1 = average)",
         title = "Library weights: one low library, or a whole noisy condition?",
         subtitle = sprintf("voomWithQualityWeights, one weight per library; labelled if below %.1f", label_below)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  p <- .ddLayout(p, tab$condition)
  list(table = tab, plot = p)
}

## 3. Each condition's own dispersion against the pooled dispersion ------------------------------------
dispersionByCondition <- function(dds, condition = "Condition", min_rep = 3, min_mean = 10) {
  x <- .ddCond(dds, condition, min_rep)
  d <- x$dds
  colData(d)$.grp <- x$grp
  design(d) <- ~ .grp
  if (is.null(sizeFactors(d)) && is.null(normalizationFactors(d))) d <- estimateSizeFactors(d)
  d <- estimateDispersions(d, quiet = TRUE)               # what ~ Condition would use in every test
  pooled <- dispersions(d)
  nc <- counts(d, normalized = TRUE)
  res <- lapply(levels(x$grp), function(g) {
    s <- d[, x$grp == g]
    design(s) <- ~ 1
    s <- estimateDispersionsGeneEst(s, quiet = TRUE)        # this condition's libraries only
    ok <- rowMeans(nc[, x$grp == g, drop = FALSE]) >= min_mean & !is.na(pooled) &
          mcols(s)$dispGeneEst > 1e-7                       # drop genes at the lower bound
    data.frame(condition = g, libraries = sum(x$grp == g), gene = rownames(d)[ok],
               log2_ratio = log2(mcols(s)$dispGeneEst[ok] / pooled[ok]))
  })
  long <- do.call(rbind, res)
  long$condition <- factor(long$condition, levels = levels(x$grp))
  tab <- aggregate(log2_ratio ~ condition + libraries, long, median)
  names(tab)[3] <- "median_log2_ratio"
  tab$genes <- as.integer(table(long$condition)[as.character(tab$condition)])
  p <- ggplot(long, aes(condition, log2_ratio)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3, fill = "grey95") +
    labs(x = NULL, y = "log2(condition's own dispersion / pooled dispersion)",
         title = "Each condition's own dispersion relative to the pooled ~ Condition dispersion",
         subtitle = sprintf("genes with mean normalised count >= %d in the condition; compare conditions with equal numbers of libraries", min_mean)) +
    theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
  p <- .ddLayout(p, long$condition, ylim = quantile(long$log2_ratio, c(0.02, 0.98)))
  list(table = tab, long = long, plot = p)
}
