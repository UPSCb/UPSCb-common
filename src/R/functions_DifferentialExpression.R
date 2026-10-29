## HAVE some toy example to test - I mean unit test :-)

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  library(RColorBrewer)
  library(tidyverse)
})


# Functions

## Plot specific gene expression


"line_plot" <- function(dds=dds,vst=vst,gene_id=gene_id){
  message(paste("Plotting",gene_id))
  sel <- grepl(gene_id,rownames(vst))
  stopifnot(sum(sel)==1)
  
  p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                        data.frame(value=vst[sel,])),
              aes(x=Treatment,y=value,col=Treatment,group=Treatment)) +
    geom_point() + geom_smooth() +
    scale_y_continuous(name="VST expression") + 
    ggtitle(label=paste("Expression for: ",gene_id))
  
  suppressMessages(suppressWarnings(plot(p)))
  return(NULL)
}


## Get suspicious genes


"get_suspicious_genes" <- function(vst,relevant_samples, threshold=10){
  
  
  ## column_to_rownames
  vstdf <- as_tibble(vst[,unlist(relevant_samples)], rownames = "Gene")
  
  
  # Convert to long format for easier manipulation
  vstdf_long <- vstdf %>%
    pivot_longer(cols = -Gene, names_to = "sample", values_to = "expression")
  
  
  # Find the maximum length of the vectors in relevant_samples
  ## UNLIKELY TO NEED AN sapply on LENGTH - R is VECTORISED
  max(elementNROWS(relevant_samples))
  max_length <- max(sapply(relevant_samples, length))
  
  # Pad the vectors with NA to make them the same length
  ## THERE IS no need in R to initialise objects
  padded_samples <- lapply(relevant_samples, function(x) {
    length(x) <- max_length
    x
  })
  
  # Convert to tibble
  group_info <- as_tibble(padded_samples) %>%
    pivot_longer(names_to = "group", values_to = "sample", cols = everything()) %>%
    arrange(group)
  
  # Add group information to vstdf_long
  vstdf_long <- vstdf_long %>%
    left_join(group_info, by = "sample")
  
  
  # First, I replace zeroes with values equal to half the minimum detected values.
  # This high values of z-replacement is so genes with very low expression values and 1 or 2 zeroes do not come up as suspicious
  # Calculate the minimum non-zero value in the expression column
  min_non_zero <- min(vstdf_long$expression[vstdf_long$expression > 0])
  
  # Define the low value as a small fraction of the minimum non-zero value
  low_value <- min_non_zero * 0.5  # Adjust the fraction as needed
  
  
  
  # Function to calculate the required values
  ## TRY to reuse the same object (to avoid duplicating memory usage)
  results <- vstdf_long %>%
    mutate(expression = ifelse(expression == 0, low_value, expression)) %>%
    group_by(Gene, group) %>%
    summarise(
      max_expr = max(expression),
      min_expr = min(expression),
      avg_expr = mean(expression),
      max_div_avg = max(expression) / mean(expression[expression != max(expression)]),
      avg_non_min = mean(expression[expression != min(expression)]),
      min_div_avg = avg_non_min / min(expression),
      .groups = "drop"
    ) %>%
    #As genes with only 0 values were returning NA in the max_div_avg and min_div_avg columns, I replace this values with 1
    mutate_all(~replace_na(., 1))
  
  
  # In the end I only use max_div_avg to determine suspicious genes
  suspicious_genes  <- results %>%
    filter(max_div_avg > threshold) %>%
    dplyr::select(Gene) %>%
    distinct() %>%
    pull()
  
  
  return(suspicious_genes)
}


## function to find out the name of columns not identical in two dataframes


compare_columns <- function(tibble1, tibble2) {
  # Initialize an empty vector to store names of non-identical columns
  non_identical_columns <- c()
  
  # Loop through each column in tibble1
  for (col in colnames(tibble1)) {
    # Check if the column exists in tibble2
    if (col %in% colnames(tibble2)) {
      # Compare the columns
      if (!all(tibble1[[col]] == tibble2[[col]])) {
        # If not identical, add the column name to the vector
        non_identical_columns <- c(non_identical_columns, col)
      }
    } else {
      # If the column does not exist in tibble2, add the column name to the vector
      non_identical_columns <- c(non_identical_columns, col)
    }
  }
  
  return(non_identical_columns)
}




# Process a comparison of interest

process_comparison <- function(comparison, Interest_col, variables_interest, output_dir= here("data/analysis/DE"), ...) {
  result_name <- gsub(" ", "", paste(comparison, collapse="vs"))
  
  result <- apply_extract_results(
    comparison = comparison,
    variables_interest = variables_interest,
    NameConditionColumn = Interest_col,
    default_prefix = paste0("DE-", result_name),
    default_dir=output_dir,
    ...
  )
  
  up <- tibble(Gene = result$up)
  down <- tibble(Gene = result$dn)
  
  write_tsv(up, here(paste0(output_dir, "/DE-", result_name, "_UP.tsv")), col_names = FALSE)
  write_tsv(down, here(paste0(output_dir, "/DE-", result_name, "_DOWN.tsv")), col_names = FALSE)
  
  return(result)
}




## Plot information about changes in differential expression depending on the cooks cutoff:


"plot_cooks_cutoffs" <- function(results_cooks_list, lfc, thres) {
  
  ## Plot 1
  #Remove complete DESeq2 results from the object to plot
  plotting_results <- lapply(results_cooks_list, function(x) {
    x[setdiff(names(x), "res")]
  })
  
  # Extract percentile values
  percentile_list <- sapply(plotting_results, function(x) x$percentile)
  
  # Set names of results_cooks_list
  names(plotting_results) <- percentile_list
  
  #Make a tibble for plotting
  df <- t(as_tibble(plotting_results))
  
  colnames(df) <- c("DEGs", "Suspicious_DEGs", "Percentile")
  
  df <- df %>%
    as_tibble() %>% 
    unnest(cols = c(DEGs, Suspicious_DEGs, Percentile))
  
  original_DEGs <- plotting_results[["0.99"]]$DEGs
  
  
  df <- df %>% mutate(Tens_lost_DEGs = (original_DEGs - DEGs)/10)
  
  # Now we plot
  
  p1 <- ggplot(df, aes(x = Percentile)) +
    geom_line(aes(y = Suspicious_DEGs, color = "Suspicious DEGs")) +
    geom_line(aes(y = Tens_lost_DEGs, color = "Tens of lost DEGs")) +
    scale_x_reverse() +
    labs(title = "Changes in suspicious DEGs and tens of lost DEGs by cooks percentile",
         x = "Percentile",
         y = "Gene number",
         color = "Legend") +
    theme_minimal()
  
  
  print(p1)
  
  
  ## Plot 2
  
  percentile_list <- sapply(results_cooks_list, function(x) x$percentile)
  
  expression_results <- lapply(results_cooks_list, function(x) {
    x[setdiff(names(x), c("suspicious", "DEGs", "percentile"))]
  })
  
  names(expression_results) <- percentile_list
  
  expression_results <- unlist(expression_results)
  
  expression_results <- map(expression_results, ~ as_tibble(.x, rownames = "Gene"))
  
  original_DEGs <- expression_results[["0.99.res"]] %>%
    filter(log2FoldChange > lfc | log2FoldChange < -lfc) %>% 
    filter(!is.na(padj)) %>% 
    filter(padj < thres) %>% 
    dplyr::select(Gene) %>% 
    pull()
  
  
  # Define function to filter tibbles
  process_tibble <- function(tibble, original_DEGs) {
    tibble %>%
      filter(Gene %in% original_DEGs) %>% 
      filter(is.na(padj)) %>% 
      dplyr::select(baseMean) %>%
      pull()
  }
  
  # Remove expression results of 0.99 as no original DEGs will have been removed there
  expression_results[["0.99.res"]] <- NULL
  
  # Apply the function to each tibble in the list
  plottable_results <- map2(expression_results, names(expression_results), ~setNames(list(process_tibble(.x, original_DEGs)), .y))
  
  
  # Convert to a single list of vectors
  plottable_results <- lapply(plottable_results, function(x) setNames(unlist(x), NULL))
  
  
  # Convert the list to a data frame
  df <- reshape2::melt(plottable_results)
  colnames(df) <- c("baseMean", "Percentile")
  
  
  df <- df  %>%
    mutate(LogBaseMean=log10(baseMean),
           Percentile=str_replace(Percentile, ".res", ""))
  
  
  chosen_percentiles <- df %>% dplyr::select(Percentile) %>% distinct() %>% slice(seq(1, n(), by = 5)) %>% pull()
  
  df <- df %>% filter(Percentile %in% chosen_percentiles)
  
  # Visualize stuff
  
  p2 <- ggplot(df, aes(x = Percentile, y = LogBaseMean)) +
    geom_boxplot() +
    labs(title = "Box Plot of Vectors (Logarithmic Scale)",
         x = "Percentile",
         y = "Log(BaseMean) of removed DEGs")
  
  
  print(p2)
  
}



## Extract the DE results. Default cutoffs are from Schurch _et al._, RNA, 2016

## Additions / replacements for UPSCb-common src/R/functions_DifferentialExpression.R
## - extract_results() gains method = c("DESeq2", "voomQWB")
## - the DESeq2 filter defaults to the mean of the libraries compared, not the
##   median of all libraries; the old behaviour is still available as "median"
## - the independent filter is optimised at alpha = padj instead of 0.1
##   (filter = "median", filter_alpha = 0.1 reproduces the previous DESeq2 results)
## - fit_voomQWB() adds nuisance covariates from design(dds) (e.g. batch) to ~ 0 + group
## - the Cook's double step refits ~ covariates + comparison_group, dropping covariates
##   that are confounded with the two groups (previously: "not full rank" error)
## Everything else (dds, VST heatmaps, baseMean, export, volcano) stays DESeq2-based.
## Source this file AFTER functions_DifferentialExpression.R so these versions win.

suppressPackageStartupMessages({
  library(DESeq2)
  library(limma)
  library(edgeR)
  library(tidyverse)
})

## Statistic for DESeq2's independent filtering
## "comparison_mean":   mean normalised count over the libraries of the two stages compared
## "comparison_median": median normalised count over the same libraries; a single outlying
##                      library cannot lift a gene over the filter, but a gene expressed in
##                      fewer than half of these libraries (e.g. only in the smaller group)
##                      has a median near zero and is filtered out
## "mean":              DESeq2's own default (baseMean, all libraries in the model)
## "median":            previous UPSCb-common default (median over all libraries)
"get_filter_stat" <- function(dds, strategy = c("comparison_mean", "comparison_median", "mean", "median"),
                              samples = colnames(dds)) {
  strategy <- match.arg(strategy)
  nc <- counts(dds, normalized = TRUE)
  switch(strategy,
         comparison_mean = rowMeans(nc[, samples, drop = FALSE]),
         comparison_median = rowMedians(nc[, samples, drop = FALSE]),
         mean = rowMeans(nc),
         median = rowMedians(nc))
}

## results() by name or by contrast
"run_deseq_results" <- function(dds, contrast, ...) {
  if (length(contrast) == 1) {
    results(dds, name = contrast, ...)
  } else {
    results(dds, contrast = contrast, ...)
  }
}

## Fit voomQWB ONCE per series (all libraries, all stages), then reuse the fit for
## every stage-to-stage comparison. group = stage/cell of each library, in colnames(dds) order;
## it must use the same labels as the NameConditionColumn of apply_extract_results().
## covariates: NULL (default) adds every variable of design(dds) that varies within a group
## (e.g. batch in ~ v1 * v2 * v3 + batch); or give a one-sided formula, e.g. ~ batch, or ~ 1 for none.
## norm = "DESeq2" uses the DESeq2 size factors as effective library sizes (as for the
## miRNAs in the report); "TMM" uses edgeR's TMM normalisation.
"fit_voomQWB" <- function(dds, group, covariates = NULL, norm = c("TMM", "DESeq2"),
                          min_count = 10, plot = FALSE) {
  norm <- match.arg(norm)
  group <- droplevels(factor(group))
  stopifnot(length(group) == ncol(dds), !anyNA(group))
  cd <- as.data.frame(colData(dds))
  
  covars <- if (is.null(covariates)) {
    all.vars(design(dds)) %>% keep(~ varies_within(cd[[.x]], group))
  } else {
    all.vars(covariates)
  }
  
  y <- DGEList(counts(dds), group = group)
  keep <- filterByExpr(y, group = group, min.count = min_count)
  y <- y[keep, , keep.lib.sizes = FALSE]
  
  if (norm == "DESeq2") {
    sf <- sizeFactors(dds)
    if (is.null(sf)) {
      ## e.g. a tximport dds, which has gene-wise normalisation factors instead of size factors
      sf <- estimateSizeFactorsForMatrix(counts(dds))
      message("No size factors in dds; using estimateSizeFactorsForMatrix(counts(dds))")
    }
    nf <- sf / y$samples$lib.size
    y$samples$norm.factors <- nf / exp(mean(log(nf)))
  } else {
    y <- calcNormFactors(y)
  }
  
  ## one column per group (cell means) + additive covariates
  design <- model.matrix(reformulate(c("0", "group", covars)),
                         data = cd %>% mutate(group = !!group))
  colnames(design)[seq_len(nlevels(group))] <- levels(group)
  if (qr(design)$rank < ncol(design)) {
    stop("The voomQWB design ~ 0 + group", paste0(" + ", covars, collapse = ""),
         " is not full rank; a covariate may be nested in group. Set covariates explicitly.")
  }
  
  ## var.group = one weight per stage, estimated from all genes ("block" weights)
  v <- voomWithQualityWeights(y, design, var.group = group, plot = plot)
  fit <- lmFit(v, design)
  
  stage_weights <- tibble(group = group, weight = v$targets$sample.weights) %>%
    group_by(group) %>%
    summarise(weight = first(weight), n_libraries = n(), .groups = "drop") %>%
    mutate(relative_variance = 1 / weight)
  
  message(sprintf("voomQWB: %d of %d genes kept by filterByExpr; design ~ 0 + group%s",
                  sum(keep), length(keep), paste0(" + ", covars, collapse = "")))
  list(fit = fit, keep = keep, stage_weights = stage_weights, norm = norm, covariates = covars)
}

## TRUE if x takes more than one value within at least one level of group
"varies_within" <- function(x, group) {
  any(tapply(as.character(x), group, n_distinct) > 1)
}

## Test one contrast from a voomQWB fit and return a DESeq2-like DataFrame
## (same column names and order as DESeqResults) so the rest of extract_results works.
## groups = c(numerator, denominator); genes removed by filterByExpr get NA p-values,
## like genes removed by DESeq2's independent filtering.
"voomQWB_contrast_results" <- function(voom_fit, dds, groups,
                                       lfc = 0.5, lfc_test = FALSE) {
  fit <- voom_fit$fit
  stopifnot(length(groups) == 2)
  if (!all(groups %in% colnames(fit$design))) {
    stop("Groups ", paste(setdiff(groups, colnames(fit$design)), collapse = ", "),
         " not found in the voomQWB fit; fit_voomQWB(group = ) must use the labels of the comparison column")
  }
  
  cvec <- setNames(numeric(ncol(fit$design)), colnames(fit$design))
  cvec[groups] <- c(1, -1)
  fit2 <- contrasts.fit(fit, cvec)
  
  if (lfc_test) {
    fit2 <- treat(fit2, lfc = lfc)
    tt <- topTreat(fit2, coef = 1, number = Inf, sort.by = "none")
  } else {
    fit2 <- eBayes(fit2)
    tt <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
  }
  
  tested <- tt %>%
    as_tibble(rownames = "gene") %>%
    transmute(gene,
              log2FoldChange = logFC,
              lfcSE = sqrt(fit2$s2.post) * fit2$stdev.unscaled[, 1],
              stat = t,
              pvalue = P.Value,
              padj = adj.P.Val)
  
  res <- tibble(gene = rownames(dds),
                baseMean = rowMeans(counts(dds, normalized = TRUE))) %>%
    left_join(tested, by = "gene") %>%
    column_to_rownames("gene") %>%
    DataFrame()
  
  mcols(res) <- DataFrame(
    type = "results",
    description = c("mean of normalized counts for all samples",
                    sprintf("log2 fold change (limma voomQWB): %s vs %s", groups[1], groups[2]),
                    "moderated standard error",
                    if (lfc_test) sprintf("treat statistic, |log2FC| > %s", lfc) else "moderated t-statistic",
                    "p-value",
                    "BH adjusted p-values"))
  metadata(res)$method <- "voomQWB"
  metadata(res)$n_tested <- nrow(tested)
  res
}

## Extract the DE results. Default cutoffs are from Schurch _et al._, RNA, 2016
"extract_results" <- function(dds, vst, contrast, relevant_samples,
                              padj = 0.01, lfc = 0.5,
                              method = c("DESeq2", "voomQWB"),
                              voom_fit = NULL,
                              lfc_test = FALSE,
                              plot = TRUE, verbose = TRUE,
                              export = TRUE, default_dir = here("data/analysis/DE"),
                              default_prefix = "DE-",
                              labels = colnames(dds),
                              sample_sel = 1:ncol(dds),
                              cooks_percentile = 0.99,
                              expression_cutoff = 0,
                              debug = FALSE,
                              filter = c("comparison_mean", "comparison_median", "mean", "median"),
                              filter_alpha = padj,
                              double_step_filter = TRUE, ...) {
  
  method <- match.arg(method)
  filtering_strategy <- match.arg(filter)
  relevant <- unlist(relevant_samples)
  stopifnot(all(relevant %in% colnames(dds)))
  stopifnot(length(sample_sel) == ncol(vst))
  
  if (method == "voomQWB") {
    ## the stage-to-stage test comes from the series-wide voomQWB fit;
    ## the contrast is taken from the names of relevant_samples (numerator first)
    if (is.null(voom_fit)) {
      stop("method = 'voomQWB' needs voom_fit = fit_voomQWB(dds, group), fitted once per series")
    }
    res <- voomQWB_contrast_results(voom_fit, dds, groups = names(relevant_samples),
                                    lfc = lfc, lfc_test = lfc_test)
    if (verbose) {
      message(sprintf("voomQWB: %s genes tested (filterByExpr), %s vs %s",
                      metadata(res)$n_tested, names(relevant_samples)[1],
                      names(relevant_samples)[2]))
    }
  } else {
    filter_stat <- get_filter_stat(dds, filtering_strategy, relevant)
    if (verbose) {
      message(sprintf("Independent filtering on: %s, optimised at alpha = %s",
                      filtering_strategy, filter_alpha))
    }
    
    ## Setting a cooks_cutoff according to the cooks_percentile parameter
    m <- ncol(counts(dds))
    p <- ncol(model.matrix(design(dds), colData(dds)))
    cooks_cutoff <- qf(cooks_percentile, p, m - p)
    
    res <- run_deseq_results(dds, contrast,
                             filter = filter_stat,
                             alpha = filter_alpha,
                             lfcThreshold = if (lfc_test) lfc else 0,
                             cooksCutoff = cooks_cutoff)
    
    if (double_step_filter) {
      ## Refit DESeq2 on the libraries of the two groups only, as ~ comparison_group, and use
      ## it only to flag outlier libraries with Cook's distance. Other design variables (e.g. a
      ## batch) are left out on purpose: DESeq2 computes Cook's distances only for libraries
      ## whose design-matrix row is shared by at least 3 libraries, so with one library per
      ## group and batch, ~ batch + comparison_group would never flag anything.
      dds_relevant <- dds[, colnames(dds) %in% relevant]
      dds_relevant$comparison_group <- factor(
        if_else(colnames(dds_relevant) %in% relevant_samples[[1]], "num", "den"),
        levels = c("den", "num"))
      design(dds_relevant) <- ~ comparison_group
      dds_relevant <- DESeq(dds_relevant, quiet = !verbose)
      
      ## all libraries of dds_relevant are in the comparison
      filter_relevant <- get_filter_stat(dds_relevant, filtering_strategy)
      
      m_relevant <- ncol(counts(dds_relevant))
      p_relevant <- ncol(model.matrix(design(dds_relevant), colData(dds_relevant)))
      cooks_cutoff_relevant <- qf(cooks_percentile, p_relevant, m_relevant - p_relevant)
      
      res_relevant <- results(dds_relevant, contrast = c("comparison_group", "num", "den"),
                              filter = filter_relevant,
                              alpha = filter_alpha,
                              cooksCutoff = cooks_cutoff_relevant)
      
      ## Genes removed by the Cook's distance test have NA pvalue and padj
      ## but a log2FoldChange
      genes_to_remove <- as_tibble(res_relevant, rownames = "Gene") %>%
        filter(is.na(padj), is.na(pvalue), !is.na(log2FoldChange)) %>%
        pull(Gene)
      
      res$pvalue[rownames(res) %in% genes_to_remove] <- NA
      res$padj[rownames(res) %in% genes_to_remove] <- NA
    }
  }
  
  if (plot) {
    par(mar = c(5, 5, 5, 5))
    volcanoPlot(res, alpha = padj, lfc = lfc)
    par(mar = mar)
  }
  
  ## independent filtering diagnostics exist only for DESeq2
  if (method == "DESeq2") {
    if (plot) {
      plot(metadata(res)$filterNumRej,
           type = "b", ylab = "number of rejections",
           xlab = "quantiles of filter")
      lines(metadata(res)$lo.fit, col = "red")
      abline(v = metadata(res)$filterTheta)
    }
    
    if (verbose) {
      message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                      round(metadata(res)$filterThreshold, digits = 5),
                      names(metadata(res)$filterThreshold)))
      max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej), "theta"]
      message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                      round(max.theta * 100, digits = 5),
                      round(quantile(counts(dds, normalized = TRUE), probs = max.theta), digits = 5)))
    }
    
    if (plot) {
      qtl.exp <- quantile(counts(dds, normalized = TRUE), probs = metadata(res)$filterNumRej$theta)
      dat <- data.frame(thetas = metadata(res)$filterNumRej$theta,
                        qtl.exp = qtl.exp,
                        number.degs = sapply(lapply(qtl.exp, function(qe) {
                          res$padj <= padj & abs(res$log2FoldChange) >= lfc &
                            !is.na(res$padj) & res$baseMean >= qe
                        }), sum))
      if (debug) {
        plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) +
               geom_line() + geom_point() +
               scale_x_continuous("quantiles of expression") +
               scale_y_continuous("base mean expression") +
               geom_hline(yintercept=expression_cutoff,
                          linetype="dotted",col="red"))
        
        p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) +
          geom_line() + geom_point() +
          scale_x_continuous("quantiles of expression") +
          scale_y_log10("base mean expression") +
          geom_hline(yintercept=expression_cutoff,
                     linetype="dotted",col="red")
        suppressMessages(suppressWarnings(plot(p)))
        
        plot(ggplot(dat,aes(x=thetas,y=number.degs)) +
               geom_line() + geom_point() +
               geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
               scale_x_continuous("quantiles of expression") +
               scale_y_continuous("Number of DE genes"))
        
        plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) +
               geom_line() + geom_point() +
               scale_x_continuous("quantiles of expression") +
               scale_y_continuous("Cumulative number of DE genes"))
        
        plot(ggplot(data.frame(x=dat$thetas[-1],
                               y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) +
               geom_line() + geom_point() +
               scale_x_continuous("quantiles of expression") +
               scale_y_continuous("Number of DE genes per interval"))
        
        plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                               y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) +
               geom_line() + geom_point() +
               scale_x_continuous("base mean of expression") +
               scale_y_continuous("Number of DE genes per interval"))
        
        p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                               y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) +
          geom_line() + geom_point() +
          scale_x_log10("base mean of expression") +
          scale_y_continuous("Number of DE genes per interval") +
          geom_vline(xintercept=expression_cutoff,
                     linetype="dotted",col="red")
        suppressMessages(suppressWarnings(plot(p)))
      }
    }
  }
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & !is.na(res$padj) &
    res$baseMean >= expression_cutoff
  
  if (verbose) {
    message(sprintf(paste(
      ifelse(sum(sel) == 1,
             "There is %s gene that is DE",
             "There are %s genes that are DE"),
      "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s (%s)"),
      sum(sel), padj, lfc, expression_cutoff, method))
  }
  
  ## proceed only if there are DE genes
  if (sum(sel) > 0) {
    val <- rowSums(vst[sel, sample_sel, drop = FALSE]) == 0
    if (sum(val) > 0) {
      warning(sprintf(paste(
        ifelse(sum(val) == 1,
               "There is %s DE gene that has",
               "There are %s DE genes that have"),
        "no vst expression in the selected samples"), sum(val)))
      sel[sel][val] <- FALSE
    }
    
    if (export) {
      if (!dir.exists(default_dir)) {
        dir.create(default_dir, showWarnings = FALSE, recursive = TRUE, mode = "0771")
      }
      write.csv(res, file = file.path(default_dir, paste0(default_prefix, "_results.csv")), quote = FALSE)
      write.csv(res[sel, ], file.path(default_dir, paste0(default_prefix, "_genes.csv")), quote = FALSE)
    }
    
    par(mar = c(0, 5, 5, 5) + 0.1)
    if (plot & sum(sel) > 1) {
      heatmap.2(t(scale(t(vst[sel, sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X) { hclust(X, method = "ward.D2") }, dendrogram = "row", Rowv = TRUE,
                trace = "none", col = hpal, labRow = FALSE, srtCol = 45, margins = c(6, 6),
                labCol = labels[sample_sel], ...)
    }
  }
  
  return(list(all = rownames(res[sel, ]),
              up = rownames(res[sel & res$log2FoldChange > 0, ]),
              dn = rownames(res[sel & res$log2FoldChange < 0, ])))
}




## Apply extract results in a contrast-dependent manner

# variables_interest must be a vector containing all the names of column in samples containing variables across which the comparisons will run.

# For example, if you plan to run comparisons across Treatment and Time, then variables_interest shoulc be c("Treatment", "Time")

# NameConditionColumn must be the column from which you will pull the values to specify your comparison.
# For example, if the comparison is c("T0_treated", "T0_control"), then NameConditionColumn will need to specify which samples belong to each of these two groups.



"apply_extract_results" <- function(dds,
                                    samples,
                                    NameConditionColumn,
                                    comparison,
                                    variables_interest,
                                    ...)
{
  
  
  
  
  print(gsub(" ", "", paste(comparison, collapse="vs")))
  
  #Make a new column whose value is equal to NameConditionColumn
  
  samples$NewCol <- samples[[NameConditionColumn]]
  
  # Determine the correct baseline
  
  values_comparison <- samples %>% filter(NewCol==comparison[1])
  
  values_comparison <- values_comparison[variables_interest] %>% distinct()
  
  
  values_baseline <- samples %>% filter(NewCol==comparison[2])
  
  values_baseline <- values_baseline[variables_interest] %>% distinct()
  
  
  difference <- compare_columns(values_comparison, values_baseline)
  
  comparison_name <- paste0(difference, "_", pull(values_comparison[difference]), "_vs_", pull(values_baseline[difference]))
  
  print(paste0("The name of the comparison is ", comparison_name))
  print("Baseline values are:")
  
  ## my religion is NO_FOR_LOOP_EVER_in_R_and_spaghetti_monsters
  for (var in variables_interest) {
    print(values_baseline[[var]])
  }
  
  # Rerun deseq2 with the correct baseline
  
  for (variable in variables_interest) {
    if (is.factor(values_baseline[[variable]])) {
      dds[[variable]] <- relevel(dds[[variable]], ref = as.character(values_baseline[[variable]]))
    } else {
      dds[[variable]] <- relevel(dds[[variable]], ref = values_baseline[[variable]])
    }
  }
  
  dds <- DESeq(dds)
  
  vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
  vst <- assay(vsd)
  vst <- vst - min(vst)
  
  # Make list of samples that will be considered to calculate the cook distance
  relevant <- list() 
  
  relevant[[comparison[1]]] <- samples %>% filter(NewCol == comparison[1]) %>% dplyr::select(SampleID) %>% pull()
  
  relevant[[comparison[2]]] <- samples %>% filter(NewCol == comparison[2]) %>% dplyr::select(SampleID) %>% pull()
  
  
  # Apply the function and store the result in the list
  return(extract_results(dds=dds,
                         vst=vst,
                         contrast=str_replace_all(comparison_name, "-", "."),
                         relevant_samples = relevant,
                         ...
  )
  )
  
}
