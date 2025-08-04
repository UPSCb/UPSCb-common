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
  
  vstdf <- as_tibble(vst[,unlist(relevant_samples)], rownames = "Gene")
  
  
  # Convert to long format for easier manipulation
  vstdf_long <- vstdf %>%
    pivot_longer(cols = -Gene, names_to = "sample", values_to = "expression")
  
  
  # Find the maximum length of the vectors in relevant_samples
  max_length <- max(sapply(relevant_samples, length))
  
  # Pad the vectors with NA to make them the same length
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




## Process a comparison of interest

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

# Suspicious degs are also returned now, users can autonomolusly decide if they want to keep them (they must be removed from the up and dn lists if that is the case)


"extract_results" <- function(dds,vst,contrast, relevant_samples,
                              padj=0.01,lfc=0.5,
                              plot=TRUE,verbose=TRUE,
                              export=TRUE,default_dir=here("data/analysis/DE"),
                              default_prefix="DE-",
                              labels=colnames(dds),
                              sample_sel=1:ncol(dds),
                              cooks_percentile = 0.99,
                              cook_optimization = FALSE,
                              plot_cooks_threshold_effect = FALSE,
                              expression_cutoff=0,
                              debug=FALSE,filter=c("median",NULL),
                              double_step_filter=TRUE, ...){
  
  # get the filter
  # get the filter
  if(!is.null(match.arg(filter))){
    filter <- rowMedians(counts(dds,normalized=TRUE))
    filtering_strategy <- "median"
    message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
  }
  
  # Setting a cooks_cutoff according to the cooks_percentile parameter
  m <- ncol(counts(dds))
  
  p <- ncol(model.matrix(design(dds), colData(dds)))
  
  cooks_cutoff <- qf(cooks_percentile, p, m - p)
  
  
  # validation
  if(length(contrast)==1){
    res <- results(dds,name=contrast,filter = filter, cooksCutoff = cooks_cutoff)
  } else {
    res <- results(dds,contrast=contrast,filter = filter, cooksCutoff = cooks_cutoff)
  }
  
  stopifnot(length(sample_sel)==ncol(vst))
  
  if(double_step_filter==TRUE){
    # Remove samples not among the relevant_samples
    dds_relevant <- dds[, colnames(dds) %in% unlist(relevant_samples)]
    
    #Find factor columns
    factor_columns <- sapply(colData(dds_relevant), is.factor)
    factor_column_names <- names(colData(dds_relevant))[factor_columns]
    
    # Read design to remove irrelevant variable in our small subset dds
    vars <- all.vars(design(dds_relevant))
    
    # Apply droplevels to each factor column
    for (element in factor_column_names) {
      dds_relevant[[element]] <- droplevels(dds_relevant[[element]])
      if (element %in% vars & length(levels(dds_relevant[[element]])) == 1) {
        vars <- setdiff(vars, element)}
    }
    
    # Put design to the only relevant variable
    design(dds_relevant) <- as.formula(paste("~", paste(vars, collapse = " + "), sep=" "))
    
    #Repeat dds analysis
    dds_relevant <- DESeq(dds_relevant)
    
    # Set filter
    if(filtering_strategy == "median"){
      filter_relevant <- rowMedians(counts(dds_relevant,normalized=TRUE))
      message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
    }
    
    
    # Setting a cooks_cutoff for relevant samples according to the cooks_percentile parameter
    m_relevant <- ncol(counts(dds_relevant))
    
    p_relevant <- ncol(model.matrix(design(dds_relevant), colData(dds_relevant)))
    
    cooks_cutoff_relevant <- qf(cooks_percentile, p_relevant, m_relevant - p_relevant)
    
    
    # Extract the results considering only relevant samples
    if(length(contrast)==1){
      res_relevant <- results(dds_relevant,name=contrast,
                              filter = filter_relevant,
                              cooksCutoff = cooks_cutoff_relevant)
    } else {
      res_relevant <- results(dds_relevant,contrast=contrast,
                              filter = filter_relevant,
                              cooksCutoff = cooks_cutoff_relevant)
    }
    
    # Where a gene was removed by the cooks distance test, its pvalue and adjusted pvalue
    # will have been set to NA, while columns like log2FoldChange will have a value
    
    genes_to_remove <- as_tibble(res_relevant, rownames="Gene") %>%
      filter(is.na(padj) & is.na(pvalue)) %>%
      filter(!is.na(log2FoldChange)) %>%
      dplyr::select(Gene) %>%
      pull()
    
    
    res$pvalue[rownames(res) %in% genes_to_remove] <- NA
    res$padj[rownames(res) %in% genes_to_remove] <- NA
  }
  
  #Find suspicious genes
  
  
  suspicious_genes <- get_suspicious_genes(vst, relevant_samples, threshold=10)
  
  thres <- padj
  
  suspicious_degs <- as_tibble(res, rownames = "Gene") %>%
    filter(Gene %in% suspicious_genes) %>%
    filter(log2FoldChange > lfc | log2FoldChange < -lfc) %>%
    filter(!is.na(padj)) %>%
    filter(padj < thres) %>%
    dplyr::select(Gene) %>%
    pull()
  
  
  #Set maximum number of suspicious degs that I wnat to accept after the optimization phase, only if the suspicious degs are more than 10
  message(sprintf("The number of suspicious DEGs is %s",
                  length(suspicious_degs)))
  
  if (length(suspicious_degs) > 10 & plot_cooks_threshold_effect == TRUE) {
    maximum_suspicious <- length(suspicious_degs)/5
    
    
    # Checking for cooks optimization in any case, just to do the plots
    
    # Define a function to check the condition and update results
    optimize_cooks <- function(cooks_percentile) {
      cooks_cutoff <- qf(cooks_percentile, p, m - p)
      res <- if (length(contrast) == 1) {
        results(dds, name = contrast, filter = filter, cooksCutoff = cooks_cutoff)
      } else {
        results(dds, contrast = contrast, filter = filter, cooksCutoff = cooks_cutoff)
      }
      
      suspicious_DEGs <- as_tibble(res, rownames = "Gene") %>%
        filter(Gene %in% suspicious_genes) %>%
        filter(log2FoldChange > lfc | log2FoldChange < -lfc) %>%
        filter(!is.na(padj)) %>%
        filter(padj < thres) %>%
        dplyr::select(Gene) %>%
        pull() %>%
        length()
      
      DEGs <- as_tibble(res, rownames="Gene") %>%
        filter(log2FoldChange > lfc | log2FoldChange < -lfc) %>%
        filter(!is.na(padj)) %>%
        filter(padj < thres) %>%
        dplyr::select(Gene) %>%
        pull() %>%
        length()
      
      list(res = res, DEGs=DEGs, suspicious = suspicious_DEGs, percentile = cooks_percentile)
    }
    
    # Use map to iterate over cooks_percentile values
    results_cooks_list <- map(seq(0.99, 0.50, by = -0.01), optimize_cooks)
    
    # Extract the count of suspicious genes from results
    counts_suspicious <- map_dbl(results_cooks_list, ~ .x$suspicious)
    
    
    
    if (cook_optimization == TRUE) {
      # Find the index where the condition is met
      index <- which(counts_suspicious < maximum_suspicious)
      
      #Get the less stringent threshold that meets the condition, or the most stringent one if no threshold meets it
      
      if (length(index) > 0) {
        res <- results_cooks_list[[index[1]]]$res
        cooks_percentile <- results_cooks_list[[index[1]]]$percentile
      } else {
        res <- tail(results_cooks_list, 1)[[1]]$res
        cooks_percentile <- tail(results_cooks_list, 1)[[1]]$percentile
      }
      
      suspicious_degs <- as_tibble(res, rownames = "Gene") %>%
        filter(Gene %in% suspicious_genes) %>%
        filter(log2FoldChange > lfc | log2FoldChange < -lfc) %>%
        filter(!is.na(padj)) %>%
        filter(padj < thres) %>%
        dplyr::select(Gene) %>%
        pull()
      
    }
    
    
    print(paste0("The cooks percentile was ", cooks_percentile, ", corresponding to a cooks cutoff of ", qf(cooks_percentile, p, m - p)))
    
    
    # Also do plots about the cooks cutoff optimization
    
    if(plot){
      plot_cooks_cutoffs(results_cooks_list=results_cooks_list, lfc=lfc, thres=padj)
    }
    
  }
  
  if(plot){
    par(mar=c(5,5,5,5))
    volcanoPlot(res, alpha=padj, lfc=lfc)
    par(mar=mar)
  }
  
  # a look at independent filtering
  if(plot){
    plot(metadata(res)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
  }
  
  if(verbose){
    message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                    round(metadata(res)$filterThreshold,digits=5),
                    names(metadata(res)$filterThreshold)))
    
    max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
    message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                    round(max.theta*100,digits=5),
                    round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
  }
  
  if(plot){
    qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
    dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                      qtl.exp=qtl.exp,
                      number.degs=sapply(lapply(qtl.exp,function(qe){
                        res$padj <= padj & abs(res$log2FoldChange) >= lfc &
                          ! is.na(res$padj) & res$baseMean >= qe
                      }),sum))
    if(debug){
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
  
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) &
    res$baseMean >= expression_cutoff
  
  if(verbose){
    message(sprintf(paste(
      ifelse(sum(sel)==1,
             "There is %s gene that is DE",
             "There are %s genes that are DE"),
      "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
      sum(sel),padj,
      lfc,expression_cutoff))
  }
  
  # proceed only if there are DE genes
  if(sum(sel) > 0){
    val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
    if (sum(val) >0){
      warning(sprintf(paste(
        ifelse(sum(val)==1,
               "There is %s DE gene that has",
               "There are %s DE genes that have"),
        "no vst expression in the selected samples"),sum(val)))
      sel[sel][val] <- FALSE
    }
    
    if(export){
      if(!dir.exists(default_dir)){
        dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
      }
      write.csv(res,file=file.path(default_dir,paste0(default_prefix,"_results.csv")), quote=FALSE)
      write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"_genes.csv")),quote=FALSE)
    }
    par(mar=c(0,5,5,5)+0.1)
    if(plot & sum(sel)>1){
      heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X,method="ward.D2")}, dendrogram='row',  Rowv=TRUE,
                trace="none",col=hpal,labRow = FALSE, srtCol = 45, margins=c(6,6),
                labCol=labels[sample_sel],...
      )
    }
  }
  return(list(all=rownames(res[sel,]),
              up=rownames(res[sel & res$log2FoldChange > 0,]),
              dn=rownames(res[sel & res$log2FoldChange < 0,]),
              suspicious=suspicious_degs))
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
