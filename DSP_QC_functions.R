
initialize_object <- function(dcc.files,
                              pkc.files,
                              annotation.file,
                              annotation.sheet.name = "Template",
                              sample.id.field.name = "Sample_ID",
                              roi.field.name = "roi",
                              panel.field.name = "panel",
                              slide.field.name = "slide name", 
                              class.field.name = "class", 
                              region.field.name = "region", 
                              segment.field.name = "segment",
                              area.field.name = "area",
                              nuclei.field.name = "nuclei", 
                              segment.id.length = 4){
    
    # load all input data into a GeoMX object
    object <-
      readNanoStringGeoMxSet(
        dccFiles = dcc.files,
        pkcFiles = pkc.files,
        phenoDataFile = annotation.file,
        phenoDataSheet = annotation.sheet.name,
        phenoDataDccColName = sample.id.field.name, 
        experimentDataColNames = panel.field.name
      )
    
    # Check the column names for required fields exist in the annotation
    
    required.field.names = c(slide.field.name, 
                             class.field.name, 
                             region.field.name, 
                             segment.field.name, 
                             roi.field.name)
    given.field.names = colnames(sData(object))
    
    # Check each of the required fields for correct naming
    for (field in required.field.names) {
      if (!(field %in% given.field.names)) {
        stop(
          paste0(
            field,
            " is not found in the annotation sheet field names.\n"
          )
        )
      }
    }
    
    # Check for the optional fields
    optional.field.names = c("area", "nuclei")
    for (field in optional.field.names) {
      if (!(field %in% given.field.names)) {
        warning(
          paste0(
            field,
            " is not found in the annotation and will not be considered \n"
          )
        )
      }
    }
    
    # Rename all of the required columns based on user parameters in data
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == slide.field.name] = "slide_name"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == class.field.name] = "class"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == region.field.name] = "region"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == segment.field.name] = "segment"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == roi.field.name] = "roi"
    
    # Rename all of the required columns based on user parameters in metadata
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == slide.field.name] = "slide_name"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == class.field.name] = "class"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == region.field.name] = "region"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == segment.field.name] = "segment"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == roi.field.name] = "roi"
    
    # Rename optional columns if they are present
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == area.field.name] = "area"
    colnames(object@phenoData@data)[colnames(object@phenoData@data) == nuclei.field.name] = "nuclei"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == area.field.name] = "area"
    rownames(object@phenoData@varMetadata)[rownames(object@phenoData@varMetadata) == nuclei.field.name] = "nuclei" 
    
    # Reformat to remove spaces and dashes in the main annotation columns
    annotation.columns <- c("class", "region", "segment", "slide_name")
    
    for(column in annotation.columns){
      pData(object)[[column]] <- gsub("\\s+", "", pData(object)[[column]])
      pData(object)[[column]] <- gsub("-", "", pData(object)[[column]])
    }
    
    # Establish the segment specific IDs
    pData(object)$segmentID <- paste0(substr(pData(object)$class, 1, segment.id.length),
                                      "|",
                                      substr(pData(object)$region, 1, segment.id.length),
                                      "|",
                                      substr(pData(object)$segment, 1, segment.id.length),
                                      "|", 
                                      substr(pData(object)$slide_name, 1, segment.id.length), 
                                      "|", 
                                      sData(object)$roi)

    
    return(object)
  
}

# Set up the MA plot table
make_MA <- function(contrast.field, 
                    condition.label, 
                    reference.label, 
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
  
  ma.plot.norm <- ggplot(ma.plot.counts, aes(x = A.value, y = M.value)) +
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=loess, col="steelblue1") + 
    geom_hline(yintercept = 0, lty = "dashed") + 
    labs(x = "Average log expression",
         y = paste0("log(", condition.label, ") - log(", reference.label, ")"), 
         title = "Post-normalization") + 
    ylim(min.y, max.y) + 
    theme_classic()
  
  ma.plot.raw <- ggplot(ma.plot.counts, aes(x = A.raw.value, y = M.raw.value)) + 
    geom_point(alpha = 0.5, col = "black") + 
    geom_smooth(method=loess, col="steelblue1") + 
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


plot_distribution <- function(object, annotation.fields){
  
  # run reductions
  color.variable <- Value <- Statistic <- NegProbe <- Q3 <- Annotation <- NULL
  
  # Start Function
  neg.probes<- "NegProbe-WTX"
  
  # Set up a list of annotation fields and values
  annotation.list <- list()
  for(field in annotation.fields){
    annotation.list[[field]] <- unique(pData(object)[[field]])
  }
  
  count.data <- t(exprs(object))
  
  annotation.data <- pData(object)
  
  stat.data <- base::data.frame(row.names = colnames(exprs(object)),
                                AOI = colnames(exprs(object)),
                                Annotation = Biobase::pData(object)[, annotation.fields],
                                Q3 = unlist(apply(exprs(object), 2,
                                                  quantile, 0.75, na.rm = TRUE)),
                                NegProbe = exprs(object)[neg.probes, ])
  
  
  
  stat.data <- stat.data %>% 
    mutate(sig2noise = Q3 / NegProbe)
  
  
  stat.data.melt <- melt(stat.data, measures.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  stat.data.melt <- melt(stat.data, 
                         measure.vars = annotation.fields, 
                         variable.name = "field", 
                         value.name = "annotation")
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Annotation, 
                                               fill=Annotation)) + 
    geom_density(alpha=0.6) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Annotation", 
         fill = "Annotation") +
    theme_bw()
  
  #stat.data.mean <- stat.data.m %>% 
  #  mutate(group = paste0(Annotation, Statistic)) %>% 
  #  group_by(group) %>% 
  #  mutate(group_mean = mean(Value)) %>% 
  #  ungroup() %>% 
  #  select(Annotation, Statistic, group_mean) %>% 
  #  distinct()
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Statistic, 
                                               fill=Statistic)) + 
    geom_density(alpha=0.6) +
    geom_vline(data=stat.data.mean, aes(xintercept=group_mean, color=Statistic),
               linetype="dashed") +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00")) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) +  
    facet_wrap(~Annotation, nrow = 1) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Statistic", 
         fill = "Statistic") +
    theme_bw()
  
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
  
  stat.data.melt <- melt(stat.data, measures.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  stat.data.mean <- stat.data.melt %>% 
    mutate(group = paste0(Annotation, Statistic)) %>% 
    group_by(group) %>% 
    mutate(group_mean = mean(Value)) %>% 
    ungroup() %>% 
    select(Annotation, Statistic, group_mean) %>% 
    distinct()
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Statistic, 
                                               fill=Statistic)) + 
    geom_density(alpha=0.6) +
    geom_vline(data=stat.data.mean, aes(xintercept=group_mean, color=Statistic),
               linetype="dashed") +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00")) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) +  
    facet_wrap(~Annotation, nrow = 1) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Statistic", 
         fill = "Statistic") +
    theme_bw()
  
  distribution.plot <- ggplot(stat.data.melt, aes(x=Value, 
                                               color=Annotation, 
                                               fill=Annotation)) + 
    geom_density(alpha=0.6) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Annotation", 
         fill = "Annotation") +
    theme_bw()
  
  #scale_x_continuous(trans = "log2") + 
  #scale_y_continuous(trans = "log2") +
  
  q3.neg.plot <- ggplot(stat.data,
                        aes(x = NegProbe, y = Q3, color = Annotation)) +
    geom_abline(alpha = 0.5, intercept = 0, slope = 1, lty = "solid", color = "darkgray") +
    geom_point(alpha = 0.3) + 
    geom_smooth(method = "loess", 
                se = FALSE, 
                linetype = "longdash", 
                alpha = 0.2) + 
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
    object.norm <- normalize(object,
                        norm_method = "quant", 
                        desiredQuantile = .75,
                        toElt = "q_norm")
    
    # The raw counts boxplot
    #transform1.raw<- exprs(object[,1:10])
    #transform2.raw<- as.data.frame(transform1.raw)
    #transform3.raw<- melt(transform2.raw)
    #ggboxplot.raw <- ggplot(transform3.raw, aes(variable, value)) +
    #  stat_boxplot(geom = "errorbar") +
    #  geom_boxplot(fill="#2CA02C") +
    #  scale_y_log10() +
    #  xlab("Segment") + 
    #  ylab("Counts, Raw") +
    #  ggtitle("Q3 Norm Counts") +
    #  scale_x_discrete(labels=c(1:10))
    
    # The normalized counts boxplot
    #transform1.norm<- assayDataElement(object[,1:10], elt = "q_norm")
    #transform2.norm<- as.data.frame(transform1.norm)
    #transform3.norm<- melt(transform2.norm)
    #ggboxplot.norm <- ggplot(transform3.norm, aes(variable, value)) +
    #  stat_boxplot(geom = "errorbar") +
    #  geom_boxplot(fill="#2CA02C") +
    #  scale_y_log10() +
    #  xlab("Segment") + 
    #  ylab("Counts, Q3 Normalized") +
    #  ggtitle("Quant Norm Counts") +
    #  scale_x_discrete(labels=c(1:10))
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
  
  stat.data.norm <- base::data.frame(row.names = colnames(object.norm@assayData$q_norm),
                                AOI = colnames(object.norm@assayData$q_norm),
                                Annotation = Biobase::pData(object.norm)[, ann.of.interest],
                                Q3 = unlist(apply(object.norm@assayData$q_norm, 2,
                                                  quantile, 0.75, na.rm = TRUE)),
                                NegProbe = object.norm@assayData$q_norm[neg.probes, ])
  
  stat.data.norm.m <- melt(stat.data.norm, measures.vars = c("Q3", "NegProbe"),
                      variable.name = "Statistic", value.name = "Value")
  
  stat.data.norm.mean <- stat.data.norm.m %>% 
    mutate(group = paste0(Annotation, Statistic)) %>% 
    group_by(group) %>% 
    mutate(group_mean = mean(Value)) %>% 
    ungroup() %>% 
    select(Annotation, Statistic, group_mean) %>% 
    distinct()
  
  distribution.plot.norm <- ggplot(stat.data.norm.m, aes(x=Value, 
                                               color=Statistic, 
                                               fill=Statistic)) + 
    geom_density(alpha=0.6) +
    geom_vline(data=stat.data.mean, aes(xintercept=group_mean, color=Statistic),
               linetype="dashed") +
    scale_color_manual(values = c("#56B4E9", "#E69F00")) +
    scale_fill_manual(values=c("#56B4E9", "#E69F00")) + 
    scale_x_continuous(limits = c(0, max(stat.data.melt$Value) + 10), 
                       expand = expansion(mult = c(0, 0))) +  
    facet_wrap(~Annotation, nrow = 1) + 
    labs(title=" Distribution per AOI of All Probes vs Negative", 
         x="Probe Counts per AOI", 
         y = "Density from AOI Count", 
         color = "Statistic", 
         fill = "Statistic") +
    theme_bw()
  
  
  
  multi.plot <- plot_grid(distribution.plot, 
                          distribution.plot.norm, 
                          ncol = 1)
  
  return(list("multi.plot" = multi.plot, "boxplot.raw" = ggboxplot.raw, "boxplot.norm" = ggboxplot.norm, "object" = object))
}

top_variable_heatmap <- function(log2.counts, 
                                 top.x.genes = 500, 
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
  
  # create Coefficient of Variation (CV) function and apply to the log counts
  calc_CV <- function(x) {sd(x) / mean(x)}
  cv.df <- data.frame(CV = apply(log2.counts, 1, calc_CV))
  
  # Take the top X most variable genes by CV score
  cv.df.top <- cv.df %>% arrange(desc(CV)) %>% slice(1:top.x.genes)
  
  # Get the list of top CV genes
  top.cv.gene.list <- rownames(cv.df.top)
  
  # Subset the counts for the top CV genes
  top.cv.heatmap.counts <- log2.counts[rownames(log2.counts) %in% top.cv.gene.list, ]
  
  # Order the counts by top CV
  top.cv.heatmap.counts <- top.cv.heatmap.counts[match(top.cv.gene.list, rownames(top.cv.heatmap.counts)), ]
  
  # Subset the annotation and arrange the order
  annotation.column.fields <- names(anno.colors)
  
  annotation.row.order <- gsub("\\.dcc", "", rownames(annotation.column))
  
  # Order the samples in counts the same as the annotation
  top.cv.heatmap.counts <- top.cv.heatmap.counts[, annotation.row.order]
  
  heatmap.plot <- pheatmap(top.cv.heatmap.counts, 
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

plot_umap <- function(log.counts, 
                      annotation, 
                      group.field, 
                      roi.field, 
                      slide.field){
  
  # Set up the counts and order by sample ID
  log.counts.transpose <- as.data.frame(t(log.counts))
  log.counts.transpose <- log.counts.transpose[order(rownames(log.counts.transpose)), ]
  
  # Order the annotation by sample ID
  annotation <- annotation[order(rownames(annotation)), ]
  
  # Run 2D UMAP and select PCs
  umap <- umap(log.counts.transpose, 
               n_components = 2, 
               random_state = 15) 
  layout <- umap[["layout"]] 
  layout <- data.frame(layout) 
  
  # Merge the annotation and UMAP
  layout$sampleID <- rownames(layout)
  annotation$sampleID <- rownames(annotation)
  umap.df <- merge(layout, annotation, by = "sampleID") 
  
  # Use the correct column names in mutate and select
  umap.df <- umap.df %>% 
    mutate(segmentID = paste({{ roi.field }}, {{ slide.field }}, sep = "|")) %>% 
    select(segmentID, X1, X2, {{ group.field }})
  
  # Create the UMAP plot
  umap.plot <- ggplot(umap.df, 
                         aes(x = X1, 
                             y = X2, 
                             color = !!sym(group.field), 
                             fill = !!sym(group.field))) +
    geom_point() + 
    geom_encircle(inherit.aes = TRUE, 
                  alpha = 0.2)
  
  return(umap.plot)
}

gene_detect_plot <- function(object, 
                             facet.column = NULL, 
                             loq.mat = NULL){
  
  # Create the plot for the all genes
  gene.stacked.bar.plot.total <- ggplot(fData(object),
                                        aes(x = DetectionThreshold)) +
    geom_bar(aes(fill = Module)) +
    geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
    theme_bw() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    labs(x = "Gene Detection Rate (Detected AOIs/Total AOIs)",
         y = "Genes, #",
         fill = "Probe Set")
  
  
  # If a facet has been selected also make a faceted bar plot
  if(!is.null(facet.column)) {
    
    # Gather the facet annotation information
    annotation.data <- pData(object)
    facet.values <- unique(annotation.data[[facet.column]])
    
    # A master df to hold all feature (gene) detection for facet values
    feature.detect.facet.df <- data.frame(feature = rownames(fData(object)))
    
    
    # Gather the IDs for each facet value
    for(value in facet.values){
      
      # Gather the sample IDs for only the current facet value
      value.df <- annotation.data %>% 
        filter(!!sym(facet.column) == value)
      
      value.IDs <- rownames(value.df)
      
      total.AOIs <- length(value.IDs)
      
      # Gather the detection per gene for value Sample IDs
      loq.mat.value <- loq.mat[, value.IDs]
      
      # Compute the detection for each feature
      value.feature.df <- data.frame(feature = rownames(fData(object)))
      
      value.feature.df[[value]] <- 100*(rowSums(loq.mat.value, na.rm = TRUE)/total.AOIs)
      
      # Add the detection per feature for this value to the master df
      feature.detect.facet.df <- merge(feature.detect.facet.df, 
                                       value.feature.df, 
                                       by = "feature")
    }
    
    # Melt the feature detect facet df for easier ggplot faceting
    
    facet.df.melt <- feature.detect.facet.df %>% 
      pivot_longer(cols = -feature, 
                   names_to = "class", 
                   values_to = "detection")
    
    # Create bins for the boxplot
    detection.bins <- c("0", 
                        "<1", 
                        "1-5", 
                        "5-10", 
                        "10-20", 
                        "20-30", 
                        "30-40", 
                        "40-50", 
                        ">50")
    
    # Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
    facet.df.melt$detection_bin <- 
      cut(facet.df.melt$detection,
          breaks = c(-1, 0, 1, 5, 10, 20, 30, 40, 50, 100),
          labels = detection.bins)
    
    facet.table <- table(facet.df.melt$detection_bin,
                         facet.df.melt$class)
    
    max.count.facet <- max(facet.table)
    
    gene.stacked.bar.plot.facet <- ggplot(facet.df.melt,
                                          aes(x = detection_bin, 
                                              fill = class)) +
      geom_bar(position = "dodge") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                         breaks = seq(0, max(max.count.facet), by = 500)) +
      labs(x = "Gene Detection Rate (Detected AOIs/Total AOIs)",
           y = "Number of Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  }
  
  
  return(list("total.plot" = gene.stacked.bar.plot.total, 
                 "facet.plot" = gene.stacked.bar.plot.facet, 
                 "facet.table" = facet.table))
}
