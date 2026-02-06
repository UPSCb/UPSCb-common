
find_difference <- function(interest, control, significant_cols, condition_col, metadata) {
  # Step 1: Divide interest and control into vectors
  interest_vector <- str_split(interest, "_")[[1]]
  control_vector <- str_split(control, "_")[[1]]
  
  # Step 2: Check if any elements or combinations exist in significant_cols
  check_values <- function(vector, metadata, significant_cols) {
    combinations <- unlist(lapply(1:length(vector), function(i) combn(vector, i, paste, collapse = "_")))
    existing_values <- map(significant_cols, ~ metadata[[.x]] %in% combinations)
    names(existing_values) <- significant_cols
    existing_values
  }
  
  interest_check <- check_values(interest_vector, metadata, significant_cols)
  
  control_check <- check_values(control_vector, metadata, significant_cols)
  
  # Step 3: Find out which of the significant_cols has different values in interest and control
  different_col <- significant_cols[map_lgl(significant_cols, ~ any(interest_check[[.x]] != control_check[[.x]]))]
  
  if (length(different_col) > 1) {
    stop("Error: More than one column in formula has different values between interest and control.")
  }
  
  if (length(different_col) == 0) {
    stop("Error: The two element of the comparison do not differ in any of the columns present in the design formula")
  }
  
  
  # Step 4: Find out the value of interest and control for the significant_cols where they have different values
  interest_values = metadata[metadata[[condition_col]] == interest][, significant_cols] %>% 
    as_tibble() %>% distinct() %>% unlist() %>% as.character()
  names(interest_values) <- significant_cols
  control_values = metadata[metadata[[condition_col]] == control][, significant_cols] %>% 
                            as_tibble() %>% distinct() %>% unlist() %>% as.character()
  names(control_values) <- significant_cols
  
  result <- list(different_col, interest_values[[different_col]],
                 control_values[[different_col]],
                 interest_values, control_values)
  
  names(result) <- c("diff_col", "diff_col_interest", "diff_col_control",
                     "interest_values", "control_values")
  
  return(result)
  
}

"run_ancombc2" <- function(phylodata, comparison, terms_in_formula, condition_col, lfc=0.5, ...) {
  interest <- comparison[1]
  control <- comparison[2]
  metadata <- as.data.frame(sample_data(phylodata))
  significant_cols <- terms_in_formula
  
  ancombc_contrast <- find_difference(interest=interest, control=control,
                                      significant_cols=significant_cols,
                                      condition_col=condition_col,
                                      metadata=metadata)
  
  for (column in names(ancombc_contrast$control_values)) {
    metadata[[column]] = relevel(factor(metadata[[column]]), ref = ancombc_contrast$control_values[[column]])
  }
  
  sample_data(phylodata) <- metadata
  
  set.seed(123)
  output = ancombc2(data = phylodata,
                    ...)
  
  term_interest <- paste0(ancombc_contrast$diff_col, ancombc_contrast$diff_col_interest)
  res = output$res %>% 
    dplyr::select(taxon,
                  paste0("lfc_", term_interest),
                  paste0("se_", term_interest),
                  paste0("W_", term_interest),
                  paste0("p_", term_interest),
                  paste0("q_", term_interest),
                  paste0("diff_", term_interest),
                  paste0("passed_ss_", term_interest))
  
  # FIlter results accoridng to adjusted pvalue and ss test of ancombc2
  
  res_filtered <- res %>%
    filter(
      !!sym(paste0("diff_", term_interest)) == TRUE,
      !!sym(paste0("passed_ss_", term_interest)) == TRUE
    )
  
  # Filter results according to logarthmic fold change
  res_filtered <- res_filtered[!(res_filtered[[2]] > (-1 * lfc) & res_filtered[[2]] < lfc)]
  
  
  final = list(output, res, res_filtered)
  
  names(final) = c("all", "res", "filtered_res")
  
  return(final)
  
}







