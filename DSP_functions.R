

# Required libraries for functions
library(pheatmap)

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


make_volcano <- function(lmm.results, 
                         title, 
                         legend.title, 
                         x.axis.title, 
                         fc.limit = 1, 
                         pos.label.limit = 1, 
                         neg.label.limit = -1){ 
  
  ## Make a volcano plot for the comparison
  
  # Define the columns for the volcano plot data
  #logfc.column.name <- paste0("logFC_", comparison)
  #padj.column.name <- paste0("adj.pval", comparison)
  
  #results$logfc <- results[[logfc.column.name]]
  #results$padj <- results[[padj.column.name]]
  
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
  
  # Compute the scale for the volcano x-axis
  log2.scale <- max(abs(lmm.results$logfc))
  
  # Establish the color scheme for the volcano plot
  contrast.level.colors <- c("steelblue4", "grey", "violetred4")
  names(contrast.level.colors) <- c("DOWN", "NONE", "UP")
  
  # Make the volcano plot
  volcano.plot <- ggplot(data = lmm.results, aes(x = logfc, 
                                                 y = -log10(padj), 
                                                 col = de_direction, 
                                                 label = deglabel)) +
    geom_vline(xintercept = c(-fc.limit, fc.limit), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    xlim(-7.5, 7.5) + 
    labs(x = x.axis.title,
         y = "-log10 adjusted p-value", 
         title = title) + 
    geom_point(size = 2) +
    scale_color_manual(legend.title, 
                       values = contrast.level.colors) + 
    geom_text_repel(max.overlaps = Inf) + 
    xlim(-log2.scale-1, log2.scale+1) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  return(list("volcano.plot" = volcano.plot))
  
}

region.types <- c("tumor", "vessel")

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

make_heatmap <- function(normalized.log.counts.df, 
                         de.results, 
                         top.degs, 
                         logfc.column = NULL, 
                         logfc.cutoff = NULL, 
                         annotation.column, 
                         annotation.row = NULL, 
                         anno.colors, 
                         cluster.rows = FALSE, 
                         cluster.columns = FALSE, 
                         main.title, 
                         row.gaps = NULL, 
                         column.gaps = NULL, 
                         show.rownames = FALSE, 
                         show.colnames = FALSE){
  
  
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
    
    # Revert to only p-value correction if no DEGs with logFC cutoff
    if(length(rownames(degs.df)) < 2){
      
      degs.df <- de.results %>% 
        filter(padj < 0.05) %>% 
        arrange(desc(padj))
      
      print("Not enough DEGs with listed logFC cutoff, reverting to all DEGs with adj p-value < 0.05")
      
    }
    
    # If there are more then 500 DEGs, trim down to top 500
    if(length(rownames(degs.df)) > 500){
      degs.df <- degs.df %>% slice(1:500)
    }
    
    # Grab the list of DEGs
    degs.list <- degs.df$gene
    
    # Subset the counts df for the DEGs and order based on the DEGs list
    counts <- normalized.log.counts.df[rownames(normalized.log.counts.df) %in% degs.list, ]
    counts <- counts[match(degs.list, rownames(counts)), ]
    
  } else {
    
    counts <- normalized.log.counts.df
    
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
                           fontsize_row = 4)
  
  
  return(heatmap.plot)
  
}


calculate_signal2noise <- function(){
  
  
}


normalize_counts <- function(object, norm.type, facet.annotation) {
  
  if(class(object)[1] != "NanoStringGeoMxSet"){
    stop(paste0("Error: You have the wrong data class, must be NanoStringGeoMxSet" ))
  }
  
  # run reductions
  color.variable <- Value <- Statistic <- NegProbe <- Q3 <- Annotation <- NULL
  
  # Start Function
  neg.probes<- "NegProbe-WTX"
  ann.of.interest <- facet.annotation
  
  stat.data <- base::data.frame(row.names = colnames(exprs(object)),
                                AOI = colnames(exprs(object)),
                                Annotation = Biobase::pData(object)[, ann.of.interest],
                                Q3 = unlist(apply(exprs(object), 2,
                                                  quantile, 0.75, na.rm = TRUE)),
                                NegProbe = exprs(object)[neg.probes, ])
  
  stat.data.m <- melt(stat.data, measures.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  stat.data.mean <- stat.data.m %>% 
    mutate(group = paste0(Annotation, Statistic)) %>% 
    group_by(group) %>% 
    mutate(group_mean = mean(Value)) %>% 
    ungroup() %>% 
    select(Annotation, Statistic, group_mean) %>% 
    distinct()
  
  distribution.plot <- ggplot(stat.data.m, aes(x=Value, 
                                               color=Statistic, 
                                               fill=Statistic)) + 
    geom_density(alpha=0.6) +
    geom_vline(data=stat.data.mean, aes(xintercept=group_mean, color=Statistic),
               linetype="dashed") +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00")) + 
    scale_x_continuous(limits = c(0, max(stat.data.m$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) +  
    facet_wrap(~Annotation, nrow = 1) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Statistic", 
         fill = "Statistic") +
    theme_bw()
  
  #scale_x_continuous(trans = "log2") + 
  #scale_y_continuous(trans = "log2") +
  
  q3.neg.plot <- ggplot(stat.data,
                 aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(alpha = 0.5, intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point(alpha = 0.5) + 
    geom_smooth(method = "loess", 
                se = FALSE, 
                linetype = "dashed", 
                alpha = 0.5) + 
    theme_bw() + 
    theme(aspect.ratio = 1) +
    labs(title = "Q3 versus Negative Mean", 
         x = "Negative Probe GeoMean per AOI", 
         y = "Q3 of all Probes per AOI ") +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2")
  
  plt3 <- ggplot(stat.data,
                 aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point() + theme_bw() +
    scale_x_continuous(trans = "log2") + 
    scale_y_continuous(trans = "log2") +
    theme(aspect.ratio = 1) +
    labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")
  
  btm.row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                       rel_widths = c(0.43,0.57))
  multi.plot <- plot_grid(plt1, btm.row, ncol = 1, labels = c("A", ""))
  
  if(norm == "q3"){
    # Q3 norm (75th percentile) for WTA/CTA  with or without custom spike-ins
    object <- normalize(object,
                        norm_method = "quant", 
                        desiredQuantile = .75,
                        toElt = "q_norm")
    
    # The raw counts boxplot
    transform1.raw<- exprs(object[,1:10])
    transform2.raw<- as.data.frame(transform1.raw)
    transform3.raw<- melt(transform2.raw)
    ggboxplot.raw <- ggplot(transform3.raw, aes(variable, value)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot(fill="#2CA02C") +
      scale_y_log10() +
      xlab("Segment") + 
      ylab("Counts, Raw") +
      ggtitle("Q3 Norm Counts") +
      scale_x_discrete(labels=c(1:10))
    
    # The normalized counts boxplot
    transform1.norm<- assayDataElement(object[,1:10], elt = "q_norm")
    transform2.norm<- as.data.frame(transform1.norm)
    transform3.norm<- melt(transform2.norm)
    ggboxplot.norm <- ggplot(transform3.norm, aes(variable, value)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot(fill="#2CA02C") +
      scale_y_log10() +
      xlab("Segment") + 
      ylab("Counts, Q3 Normalized") +
      ggtitle("Quant Norm Counts") +
      scale_x_discrete(labels=c(1:10))
  }
  if(norm == "Q3"){
    stop(paste0("Error: Q3 needs to be q3" ))
  }
  if(norm == "quantile"){
    stop(paste0("Error: quantile needs to be q3" ))
  }
  if(norm == "Quantile"){
    stop(paste0("Error: Quantile needs to be q3" ))
  }
  if(norm == "quant"){
    stop(paste0("Error: quant needs to be q3" ))
  }
  
  if(norm == "neg"){
    # Background normalization for WTA/CTA without custom spike-in
    object <- normalize(object,
                        norm_method = "neg", 
                        fromElt = "exprs",
                        toElt = "neg_norm")
    
    # The raw counts boxplot
    transform1.raw<- exprs(object[,1:10])
    transform2.raw<- as.data.frame(transform1.raw)
    transform3.raw<- melt(transform2.raw)
    ggboxplot.raw <- ggplot(transform3.raw, aes(variable, value)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot(fill="#FF7F0E") +
      scale_y_log10() +
      xlab("Segment") + 
      ylab("Counts, Raw") +
      ggtitle("Neg Norm Counts") +
      scale_x_discrete(labels=c(1:10))
    
    # The normalized counts boxplot
    transform1.norm<- assayDataElement(object[,1:10], elt = "neg_norm")
    transform2.norm<- as.data.frame(transform1.norm)
    transform3.norm<- melt(transform2.norm)
    ggboxplot.norm <- ggplot(transform3.norm, aes(variable, value)) +
      stat_boxplot(geom = "errorbar") +
      geom_boxplot(fill="#FF7F0E") +
      scale_y_log10() +
      xlab("Segment") + 
      ylab("Counts, Neg. Normalized") +
      ggtitle("Neg Norm Counts") +
      scale_x_discrete(labels=c(1:10))
  }
  if(norm == "Neg"){
    stop(paste0("Error: Neg needs to be neg" ))
  }
  if(norm == "negative"){
    stop(paste0("Error: negative needs to be neg" ))
  }
  if(norm == "Negative"){
    stop(paste0("Error: Negative needs to be neg" ))
  }
  
  return(list("multi.plot" = multi.plot, "boxplot.raw" = ggboxplot.raw, "boxplot.norm" = ggboxplot.norm, "object" = object))
}

run_lme4_lmm <- function(){

  
  
}