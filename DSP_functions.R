

# Required libraries for functions
library(pheatmap)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)

subset_counts_for_lmm <- function(counts, 
                                   annotation, 
                                   subset.list){ 
  
  subset.counts <- counts
  subset.annotation <- annotation
  
  # Subset the object based on the given annotations
  for(column in names(subset.list)){ 
    
    subset.annotation <- subset.annotation %>% 
      filter(.[[column]] %in% subset.list[[column]])
    
    subset.IDs <- subset.annotation$Sample_ID
    
    subset.columns <- c("gene", subset.IDs)
    
    subset.counts <- subset.counts %>% 
      select(all_of(subset.columns))
    
    # Factor the columns with relevant annotations
    subset.annotation[[column]] <- factor(subset.annotation[[column]])
    
  }
  
  # Factor the slide column
  subset.annotation[["slide_name"]] <- factor(subset.annotation[["slide_name"]])
  
  # Create log2 counts
  subset.counts.log2 <-  subset.counts %>%
    mutate(across(where(is.numeric), log2))
  
  return(list("subset.counts" = subset.counts, 
              "subset.log.counts" = subset.counts.log2, 
              "subset.annotation" = subset.annotation))
  
}

subset_object_for_lmm <- function(object, 
                                  subset.list){ 
  
  # Set up the object to subset
  subset.object <- object
  
  # Subset the object based on the given annotations
  for(column in names(subset.list)){ 
    
    subset.indices <- pData(subset.object)[[column]] %in% subset.list[[column]]
    subset.object <- subset.object[, subset.indices]
    
    # Factor the columns with relevant annotations
    pData(subset.object)[[column]] <- factor(pData(subset.object)[[column]])
    
  }
  
  # Factor the slide column
  pData(subset.object)[["slide_name"]] <- 
    factor(pData(subset.object)[["slide_name"]])
  
  # Create log2 counts
  assayDataElement(object = subset.object, elt = "log_q") <-
    assayDataApply(subset.object, 2, FUN = log, base = 2, elt = "q_norm")
  
  assayDataElement(object = subset.object, elt = "log_raw") <-
    assayDataApply(subset.object, 2, FUN = log, base = 2, elt = "exprs")
  
  
  # Gather the log counts and annotation to return
  log.counts <- subset.object@assayData$log_q
  raw.log.counts <- subset.object@assayData$log_raw
  annotation.df <- pData(subset.object)
  
  # Replace all bad characters in column names
  annotation.df <- annotation.df %>%
    rename_all(~str_replace_all(., " ", "_"))
  
  return(list("subset.object" = subset.object, 
              "log.counts" = log.counts, 
              "raw.log.counts" = raw.log.counts, 
              "annotation" = annotation.df))
  
}


run_limma <- function(counts, 
                      annotation, 
                      include.slide, 
                      within.slide, 
                      contrast, 
                      contrast.levels){
  
  # Create the DGE object
  DGE.list <- DGEList(counts = counts, 
                      samples = annotation)
  
  if(include.slide == FALSE){ 
    # Create the LM model design
    design <- model.matrix(formula(paste0("~ 0 + ", contrast)), 
                           data = DGE.list$samples)
    
  } else {
    
    if(within.slide == TRUE){ 
      # For within slide we use a random slope in the mixed effect
      
      # Create the LM model design with slide as a mixed effect
      design <- model.matrix(formula(paste0("~ 1 + ", 
                                            contrast, 
                                            " + (1 + " , 
                                            contrast, 
                                            " | slide_name)")), 
                             data = DGE.list$samples)
      
    } else{
      # For between slide we use slide in the mixed effect, no random slope
      
      # Create the LM model design with slide as a mixed effect
      design <- model.matrix(formula(paste0("~ 1 + ", 
                                            contrast, 
                                            " + (1 | slide_name)")), 
                             data = DGE.list$samples)
    }
    
  }
  
  
  # Create the fit for the model
  fit <- lmFit(DGE.list$counts, design)
  
  # Set up the contrast
  contrast.level.ref <- paste0(contrast, contrast.levels[[1]])
  contrast.level.condition <- paste0(contrast, contrast.levels[[2]])
  
  
  contrast <- makeContrasts(paste0(contrast.level.condition, 
                                   " - ", 
                                   contrast.level.ref),
                            levels = colnames(coef(fit)))
  
  # Generate the estimate of the contrast
  contrast.estimate <- contrasts.fit(fit, contrast)
  
  # Run Empirical Bayes smoothing of standard errors 
  fit.eb <- eBayes(contrast.estimate, robust = TRUE)
  
  # Generate the results table
  results <- topTable(fit.eb, sort.by = "P", n=Inf)
  
  
  return(list("results" = results, 
              "fit" = fit.eb, 
              "design" = design))
}



run_lmm <- function(object, contrast, within.slide){
  
  if(within.slide == TRUE){
    
    # Run the linear model with random slope
    lmm.results <- mixedModelDE(object, 
                                elt = "log_q", 
                                modelFormula = formula(paste0("~ 1 + ", 
                                                              contrast, 
                                                              " + (1 + " , 
                                                              contrast, 
                                                              " | slide_name)")), 
                                groupVar = contrast, 
                                nCores = parallel::detectCores(), 
                                multiCore = TRUE)
  } else {
    
    lmm.results <- mixedModelDE(object, 
                                elt = "log_q", 
                                modelFormula = formula(paste0("~ 1 + ", 
                                                              contrast, 
                                                              " + (1 | slide_name)")), 
                                groupVar = contrast, 
                                nCores = parallel::detectCores(), 
                                multiCore = TRUE)
    
  }
  
  # Gather the results into an output table
  lmm.results.summary <- do.call(rbind, lmm.results["lsmeans", ])
  lmm.results.summary <- as.data.frame(lmm.results.summary)
  
  # use lapply in case you have multiple levels of your test factor to
  # correctly associate gene name with it's row in the results table
  lmm.results.summary$gene <- 
    unlist(lapply(colnames(lmm.results),
                  rep, nrow(lmm.results["lsmeans", ][[1]])))
  
  # Run multiple test correction
  lmm.results.summary$FDR <- p.adjust(lmm.results.summary$`Pr(>|t|)`,
                                      method = "fdr")
  
  
  # Rename columns
  lmm.results.summary$pval <- lmm.results.summary[["Pr(>|t|)"]]
  lmm.results.summary$adj_pval <- lmm.results.summary$FDR
  lmm.results.summary$logfc <- lmm.results.summary$Estimate
  
  # Format final summary data frame
  lmm.results.summary <- lmm.results.summary[, c("gene", "logfc", 
                                                 "pval", "adj_pval")]
  
  return(list("results" = lmm.results.summary, "lm.output" = lmm.results))
  
}

# Set up the MA plot table
make_MA <- function(contrast.field, 
                    condition.label, 
                    reference.label, 
                    results.df, 
                    log.counts, 
                    raw.log.counts, 
                    annotation){
  
  # Gather the sample IDs for condition and reference groups
  condition.samples <- rownames(annotation[annotation[[contrast.field]] == condition.label, ])
  reference.samples <- rownames(annotation[annotation[[contrast.field]] == reference.label, ])
  
  # Gather normalized and raw counts for both groups
  condition.counts <- as.data.frame(log.counts[, condition.samples])
  reference.counts <- as.data.frame(log.counts[, reference.samples])
  
  condition.raw.counts <- as.data.frame(raw.log.counts[, condition.samples])
  reference.raw.counts <- as.data.frame(raw.log.counts[, reference.samples])  
  
  # Get the mean log score for each gene for both 
  # normalized counts
  condition.row.order <- rownames(condition.counts)
  condition.counts <- as.data.frame(sapply(condition.counts, as.numeric))
  condition.counts$cond_mean <- rowMeans(condition.counts)
  condition.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.counts)
  reference.counts <- as.data.frame(sapply(reference.counts, as.numeric))
  reference.counts$ref_mean <- rowMeans(reference.counts)
  reference.counts$gene <- reference.row.order
  
  # raw counts
  condition.row.order <- rownames(condition.raw.counts)
  condition.raw.counts <- as.data.frame(sapply(condition.raw.counts, as.numeric))
  condition.raw.counts$cond_raw_mean <- rowMeans(condition.raw.counts)
  condition.raw.counts$gene <- condition.row.order
  
  reference.row.order <- rownames(reference.raw.counts)
  reference.raw.counts <- as.data.frame(sapply(reference.raw.counts, as.numeric))
  reference.raw.counts$ref_raw_mean <- rowMeans(reference.raw.counts)
  reference.raw.counts$gene <- reference.row.order
  
  
  # Create a new data frame of the gene and group means with M and A values
  normalized.counts <- merge(condition.counts, reference.counts, by = "gene") %>% 
    select(gene, cond_mean, ref_mean) %>% 
    mutate(M.value = cond_mean - ref_mean) %>% 
    mutate(A.value = (cond_mean + ref_mean)/2)
  
  raw.counts <- merge(condition.raw.counts, reference.raw.counts, by = "gene") %>% 
    select(gene, cond_raw_mean, ref_raw_mean) %>% 
    mutate(M.raw.value = cond_raw_mean - ref_raw_mean) %>% 
    mutate(A.raw.value = (cond_raw_mean + ref_raw_mean)/2)
  
  # Add the DE results and log counts together
  ma.plot.counts <- merge(normalized.counts, raw.counts, by = "gene")
  
  # Set the bounds for the y axix so that they are aligned
  min.y <- min(c(min(ma.plot.counts$M.value),min(ma.plot.counts$M.raw.value)))
  max.y <- max(c(max(ma.plot.counts$M.value),max(ma.plot.counts$M.raw.value)))
  
  ma.plot.norm <- ggplot(ma.plot.table, aes(x = A.value, y = M.value)) +
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=lm, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Post-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  ma.plot.raw <- ggplot(ma.plot.table, aes(x = A.raw.value, y = M.raw.value)) + 
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=lm, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Pre-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  combined.MA.plots <- arrangeGrob(ggplotGrob(ma.plot.raw), 
                                   ggplotGrob(ma.plot.norm), 
                                   nrow = 1, ncol = 2)
  
  return(combined.MA.plots)
  
  
  
  
}

run_GSEA <- function(){
  
  
  
}

make_heatmap <- function(normalized.log.counts.df = q3.norm.log.counts, 
                         de.results = NULL, 
                         top.degs = FALSE, 
                         top.variable = FALSE, 
                         logfc.column = NULL, 
                         logfc.cutoff = NULL, 
                         annotation.column, 
                         annotation.row = NULL, 
                         anno.colors, 
                         cluster.rows = TRUE, 
                         cluster.columns = TRUE, 
                         main.title, 
                         row.gaps = NULL, 
                         column.gaps = NULL, 
                         show.rownames = FALSE, 
                         show.colnames = FALSE, 
                         min.genes.to.display = 2, 
                         max.genes.to.display = 500, 
                         font.size.row = 4){
  
  
  
  if(top.degs == TRUE & top.variable == TRUE){ 
  
    stop("Set only one of top.degs or top.variable to TRUE, not both")  
    
  }
  
  #if (top.variable == TRUE){
    
  #}
  
  # Filter genes by top DEGs, if applicable
  if(top.degs == TRUE){ 
    
    # Arrange by adjusted p-value
    degs.df <- de.results %>% 
      filter(padj < 0.05) %>% 
      arrange(desc(padj))
    
    # Arrange by log FC
    degs.df <- degs.df %>% arrange(desc(logfc))
    
    if(!is.null(logfc.cutoff)){
      
      degs.df <- degs.df %>% 
        filter(.data[[logfc.column]] > logfc.cutoff | .data[[logfc.column]] < -(logfc.cutoff))
    }
    
    # Use DEGs meeting adj p-val and logfc cutoffs
    if(length(rownames(degs.df)) >= min.genes.to.display){
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [logfc (+-", logfc.cutoff, ")")
      
      main.title <- paste0(main.title, " adj p-val (<0.05)]")
      
    }    
    
    # Revert to adj p-val cutoff and no logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(padj < 0.05) %>% 
        arrange(desc(padj))
      
      print("Not enough DEGs with listed logFC cutoff, reverting to all DEGs with adj p-value < 0.05")
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [no logfc cutoff")
      
      main.title <- paste0(main.title, " adj p-val (<0.05)]")
      
    }
    
    # Revert to nonadj p-val cutoff and logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(pval < 0.05) %>% 
        arrange(desc(padj)) %>% 
        filter(.data[[logfc.column]] > logfc.cutoff | .data[[logfc.column]] < -(logfc.cutoff))
      
      print("Not enough DEGs with adj p-val < 0.05, reverting to all DEGs with p-value < 0.05")
      
      main.title <- paste0(main.title, " [logfc (+-", logfc.cutoff, ")")
      
      main.title <- paste0(main.title, " p-val (<0.05, no adj)]")
      
    }
    
    # Revert to only p-val if no DEGs with p-val cutoff and logfc cutoff
    if(length(rownames(degs.df)) < min.genes.to.display){
      
      degs.df <- de.results %>% 
        filter(pval < 0.05) %>% 
        arrange(desc(pval))
      
      print("Not enough DEGs with adj p-val < 0.05, reverting to NON-adjusted p-value < 0.05")
      
      # Set up main title with logfc cutoff
      main.title <- paste0(main.title, " [no logfc cutoff ")
      
      main.title <- paste0(main.title, " p-val (<0.05, no adj)]")
      
    }
    
    # If there are more then 500 DEGs, trim down to top 500
    if(length(rownames(degs.df)) > max.genes.to.display){
      degs.df <- degs.df %>% slice(1:max.genes.to.display)
    }
    
    # Grab the list of DEGs
    degs.list <- degs.df$gene
    
    # Subset the counts df for the DEGs and order based on the DEGs list
    counts <- normalized.log.counts.df[rownames(normalized.log.counts.df) %in% degs.list, ]
    counts <- counts[match(degs.list, rownames(counts)), ]
    
  } else {
    
    counts <- normalized.log.counts.df
    
  }
  
  # Arrange by annotations if no column clustering
  if(cluster.columns == FALSE){
    
    # First arrange the annotation by the annotation groups
    anno.col.names <- colnames(annotation.column)
    
    for(col in anno.col.names){
      
      #annotation.column[[col]] <- as.factor(annotation.column[[col]])
      
      annotation.column <- annotation.column %>% 
        arrange(.data[[col]])
      
    }
    
    # Next match the counts file to the row order of the annotation
    counts <- counts[, rownames(annotation.column), drop = FALSE]
    
  }
  
  heatmap.plot <- pheatmap(counts, 
                           main = main.title, 
                           show_rownames = show.rownames, 
                           scale = "row",   
                           show_colnames = show.colnames,
                           border_color = NA, 
                           cluster_rows = cluster.rows, 
                           cluster_cols = cluster.columns, 
                           clustering_method = "average", 
                           clustering_distance_rows = "correlation", 
                           clustering_distance_cols = "correlation", 
                           color = colorRampPalette(c("blue", "white", "red"))(120), 
                           annotation_row = annotation.row, 
                           annotation_col = annotation.column,  
                           annotation_colors = anno.colors, 
                           gaps_row = row.gaps, 
                           gaps_col = column.gaps, 
                           fontsize_row = font.size.row)
  
  
  
  return(heatmap.plot)
  
}


calculate_signal2noise <- function(){
  
  
}


normalize_counts <- function() {}

gsea_preranked_list <- function(contrast.field, 
                                contrast.levels, 
                                annotation, 
                                log.counts){
  
  # Gather the signal to noise ratio for GSEA ranking
  # Default method for ranking genes from GSEA manual:
  # https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideTEXT.htm#_Metrics_for_Ranking
  
  # Contrast level A is the "condition" (positive when calculating fold change)
  contrast.A.annotation <- annotation %>% 
    filter(!!sym(contrast.field) == contrast.levels[1])
  
  contrast.A.sampleIDs <- rownames(contrast.A.annotation)
  
  contrast.A.counts <- as.data.frame(log.counts) %>% 
    select(all_of(contrast.A.sampleIDs))
  
  contrast.A.counts$gene <- rownames(contrast.A.counts)
  
  # Contrast level B is the "reference" (negative when calculating fold change)
  
  contrast.B.annotation <- annotation %>% 
    filter(!!sym(contrast.field) == contrast.levels[2])
  
  contrast.B.sampleIDs <- rownames(contrast.B.annotation)
  
  contrast.B.counts <- as.data.frame(log.counts) %>% 
    select(all_of(contrast.B.sampleIDs))
  
  contrast.B.counts$gene <- rownames(contrast.B.counts)
  
  # Add a column to each contrast level for the mean and standard deviation
  contrast.A.counts <- contrast.A.counts %>% 
    mutate(mean.A = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.A = apply(select_if(., is.numeric), 1, sd))
  
  contrast.B.counts <- contrast.B.counts %>% 
    mutate(mean.B = rowMeans(select_if(., is.numeric))) %>%  
    mutate(stdev.B = apply(select_if(., is.numeric), 1, sd))
  
  GSEA.preanked.df <- merge(contrast.A.counts, contrast.B.counts, by = "gene")
  
  GSEA.preanked.df <- GSEA.preanked.df %>% 
    mutate(signal2noise = (mean.A - mean.B)/(stdev.A + stdev.B)) %>% 
    arrange(desc(signal2noise)) %>% 
    select(c(gene, mean.A, mean.B, stdev.A, stdev.B, signal2noise))
  
  return(GSEA.preanked.df)
  
}

make_volcano <- function(lmm.results, 
                         title, 
                         legend.title, 
                         x.axis.title, 
                         fc.limit = 1, 
                         pos.label.limit = 1, 
                         neg.label.limit = -1, 
                         custom.gene.labels = NULL){ 
  
  ## Make a volcano plot for the comparison
  
  # Create a column for direction of DEGs
  lmm.results$de_direction <- "NONE"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc > fc.limit] <- "UP"
  lmm.results$de_direction[lmm.results$padj < 0.05 & 
                             lmm.results$logfc < -fc.limit] <- "DOWN"
  
  # Create a label for DEGs based on label limits
  lmm.results$deglabel <- ifelse((lmm.results$logfc > pos.label.limit | 
                                    lmm.results$logfc < neg.label.limit) & 
                                   lmm.results$padj < 0.05, 
                                 lmm.results$gene,
                                 NA
  )
  
  # Create a label for DEGs
  if(is.null(custom.gene.labels)){
    
    lmm.results$deglabel <- ifelse(lmm.results$de_direction == "NONE", 
                                   NA, 
                                   lmm.results$gene)
    
  } else {
    
    lmm.results$deglabel <- ifelse(lmm.results$gene %in% custom.gene.labels, 
                                   lmm.results$gene, 
                                   NA)
    
  }
  
  # Compute the scale for the volcano x-axis
  log2.scale <- max(abs(lmm.results$logfc))
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("steelblue4", "grey", "violetred4")
  names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
  
  # Make the volcano plot based on custom gene labels
  if(is.null(custom.gene.labels)){
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
    
    volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                                   y = -log10(padj), 
                                                   col = de_direction, 
                                                   label = deglabel)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors) + 
      geom_text_repel(max.overlaps = Inf) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5))
    
  } else {
    
    # Label the custom genes depending on significance
    lmm.results <- lmm.results %>% 
      mutate(custom.label = ifelse(!is.na(deglabel) & de_direction == "NONE", 
                                   "BLACK", 
                                   ifelse(!is.na(deglabel) & de_direction != "NONE", 
                                          de_direction, 
                                          "NONE")))
    
    contrast.level.colors <- c("steelblue4", "grey", "violetred4", "black")
    names(contrast.level.colors) <- c("DOWN", "NONE", "UP", "BLACK")
    
    lmm.results.labeled <- lmm.results %>%
      filter(custom.label != "NONE")
    
    lmm.results.unlabeled <- lmm.results %>% 
      filter(custom.label == "NONE")
    
    
    volcano.plot <- ggplot() + 
      geom_point(data = lmm.results.unlabeled, aes(x = logfc, 
                                                   y = -log10(padj), 
                                                   col = custom.label, 
                                                   alpha = 0.5)) + 
      geom_point(data = lmm.results.labeled, aes(x = logfc, 
                                                 y = -log10(padj), 
                                                 col = custom.label, 
                                                 alpha = 1)) +
      geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
      geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
      labs(x = x.axis.title,
           y = "-log10 adjusted p-value", 
           title = title) + 
      geom_point(size = 2) +
      scale_color_manual(legend.title, 
                         values = contrast.level.colors, 
                         breaks = c("DOWN", "UP")) + 
      geom_text_repel(data = lmm.results.labeled,
                      aes(x = logfc, 
                          y = -log10(padj), 
                          label = deglabel, 
                          col = custom.label), 
                      max.overlaps = Inf, 
                      size = 6, 
                      show.legend = FALSE) + 
      xlim(-log2.scale-1, log2.scale+1) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      scale_alpha_identity(guide = "none")
    
  }
  
  return(list("volcano.plot" = volcano.plot))
  
}

make_dual_volcano <- function(lmm.results.1, 
                              lmm.results.2,
                              lmm.results.1.label,
                              lmm.results.2.label, 
                              title, 
                              legend.title, 
                              x.axis.title, 
                              fc.limit = 1, 
                              pos.label.limit = 1, 
                              neg.label.limit = -1, 
                              custom.gene.labels = NULL){
  
  # A combined list of the two results
  lmm.results.list <- list()
  lmm.results.list[[lmm.results.1.label]] <- lmm.results.1
  lmm.results.list[[lmm.results.2.label]] <- lmm.results.2
  
  # Create a column for direction of DEGs
  for(lmm.results.label in names(lmm.results.list)){
    
    lmm.results <- lmm.results.list[[lmm.results.label]]
    
    lmm.results$de_direction <- "NONE"
    lmm.results$de_direction[lmm.results$padj < 0.05 & 
                               lmm.results$logfc > fc.limit] <- "UP"
    lmm.results$de_direction[lmm.results$padj < 0.05 & 
                               lmm.results$logfc < -fc.limit] <- "DOWN"
    
    #lmm.results$de_direction <- factor(lmm.results$de_direction, 
    #                                   levels = "UP", "NONE", "DOWN")
    
    # Create a label for DEGs based on label limits
    lmm.results$deglabel <- ifelse((lmm.results$logfc > pos.label.limit | 
                                      lmm.results$logfc < neg.label.limit) & 
                                     lmm.results$padj < 0.05, 
                                   lmm.results$gene,
                                   NA)
    
    # Convert to -log10 p-value
    lmm.results$neg.log10.pval <- -log10(lmm.results$padj)
    
    # Add an analysis label 
    lmm.results$analysis.label <- lmm.results.label
    
    lmm.results.list[[lmm.results.label]] <- lmm.results
    
  }
  
  # Convert -log10 p-value to neg for mirrored second analysis
  lmm.results.list[[lmm.results.2.label]]$neg.log10.pval <- -(lmm.results.list[[lmm.results.2.label]]$neg.log10.pval)
  
  # Combine the two analyses into a master df
  lmm.results.combine <- bind_rows(lmm.results.list[[lmm.results.1.label]], 
                                   lmm.results.list[[lmm.results.2.label]])
  
  # Establish the limits and breaks for x and y axes
  log2.scale <- ceiling(max(abs(lmm.results$logfc)))
  pval.scale <- ceiling(max(abs(lmm.results.combine$neg.log10.pval)))

  y.axis.breaks <- seq(-pval.scale, pval.scale, by = 1)
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("violetred4", "grey", "steelblue4")
  names(contrast.level.colors) <- c("UP", "NONE", "DOWN")
  
  
  
  # Make the plot
  dual.volcano.plot <- ggplot(data = lmm.results.combine, 
                              aes(x = logfc, 
                                  y = neg.log10.pval, 
                                  col = de_direction, 
                                  label = deglabel)) +
    geom_vline(xintercept = c(-fc.limit, fc.limit), 
               col = "darkgray", 
               linetype = 'dashed') + 
    geom_vline(xintercept = 0, 
               col = "black") + 
    geom_hline(yintercept = c(-log10(0.05), -(-log10(0.05))), 
               col = "darkgray", 
               linetype = 'dashed') + 
    geom_hline(yintercept = 0, 
               col = "black") + 
    xlim(-log2.scale - 1, log2.scale + 1) + 
    scale_y_continuous(
      limits = c(-pval.scale, pval.scale),
      breaks = y.axis.breaks,
      labels = abs(y.axis.breaks)) + 
    labs(x = x.axis.title,
         y = "-log10 adjusted p-value", 
         title = title) + 
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(legend.title, 
                       values = contrast.level.colors) + 
    geom_text_repel(max.overlaps = Inf, show.legend = FALSE) + 
    theme(plot.title = element_text(hjust = 0.5), 
          ) + 
    theme_linedraw()
  
  #dual.volcano.labeled <- ggdraw() +
  #  draw_plot(dual.volcano.plot, x = 0.02, width = 0.95) +  # Move the plot slightly right
  #  draw_text(lmm.results.1.label, x = 0.01, y = 0.6, angle = 90, size = 14, hjust = 0) +
  #  draw_text(lmm.results.2.label, x = 0.01, y = 0.2, angle = 90, size = 14, hjust = 0)
  
  return(dual.volcano.plot)

}