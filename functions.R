
  
"%ni%" <- Negate("%in%")
source("packages.R")



add_genes_annotations <- function(vector_of_genes, attribute_known) {
  
  # Download homo sapiens genes ensembl database
  
  # 
  coordinates <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), filters = attribute_known, values = vector_of_genes$ID, mart = ensembl)
  genes_annoted <- merge(x = vector_of_genes, y = coordinates, by.x = "ID", by.y = attribute_known, all.x = TRUE)
  return(genes_annoted)
}




Differential_analysis <- function(Focused_variable, DATA){
  design.pairs <- function(levels) {
    n <- length(levels)
    design <- matrix(0,n,choose(n,2))
    rownames(design) <- levels
    colnames(design) <- 1:choose(n,2)
    k <- 0
    for (i in 1:(n - 1))
      for (j in (i + 1):n) {
        k <- k + 1
        design[i,k] <- 1
        design[j,k] <- -1
        colnames(design)[k] <- paste(levels[i], "-", levels[j],sep = "")
      }
    design
  }
  design <- model.matrix(~0 + Focused_variable)
  contr.matrix <- design.pairs(levels(factor(Focused_variable)))
  colnames(design) <- rownames(contr.matrix)   
  Fit <- lmFit(DATA, design) %>%
    contrasts.fit(., contr.matrix) %>%
    eBayes(., trend = TRUE)
  
  FitList <- list()
  for (i in 1:ncol(contr.matrix)) {
    FitList[[i]] <- topTable(Fit, coef = i, adjust.method = "BH", number = nrow(DATA)) %>%
      mutate(ID = rownames(.))
    
    message(paste0(i, " done"))
    
  }
  names(FitList) <- colnames(contr.matrix)
  return(FitList)
}




add_genes_coordinates <- function(vector_of_genes) {
  gene_name_annotation <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), filters = 'ensembl_gene_id', values = vector_of_genes$ID, mart = ensembl)
  genes_annoted <- merge(x = vector_of_genes, y = gene_name_annotation, by.x = "ID", by.y = "ensembl_gene_id", all.x = TRUE)
  return(genes_annoted)
}




prepare_pchic <- function(cell_lines = "all", minimum_interaction = 5){
  load("~/PCHIC/pchic.RData")
  if (cell_lines == "all") {
    cell_lines = c("Mon", "Mac0", "Mac1", "Mac2", "Neu", "MK", "EP", "Ery", "FoeT", "nCD4", "tCD4", "aCD4", "naCD4", "nCD8", "tCD8", "nB", "tB")
  }
  pchic <- data.frame(pchic[rowSums(pchic[,cell_lines] >= minimum_interaction) >= 1, 1:10]) %>% na.omit(.)
  colnames(pchic)[c(1:5, 6:10)] <- rep(c("chr", "start", "end", "ID", "Name"), 2)
  return(pchic)
}




Create_pchic_Grange <- function(pchic){
  PCHiC_bed <- unique(rbind(pchic[, c(1:3, 5)], pchic[, c(6:8, 10)]))
  PCHiC_GRange <- GRanges(
    seqnames = PCHiC_bed$chr,
    IRanges(start = PCHiC_bed$start, end = PCHiC_bed$end),
    Gene_Pchic = PCHiC_bed$Name,
    start_fragment = PCHiC_bed$start,
    end_fragment = PCHiC_bed$end
  )
  PCHiC_GRange$ID <- paste(PCHiC_bed$chr, PCHiC_bed$start, sep = "_")
  return(PCHiC_GRange)
}



Gene_name_finding <- function(gene){
  Blueprint_gene_name <- Blueprint$gene_names %>% unique(.)
  gene_referenced <- c(Blueprint_gene_name) %>%
    unlist(.) %>%
    unique(.) %>%
    .[str_detect(., gene)]
  return(gene_referenced)
}


Look_at_differential_methylation <- function(Phenotype, DATA_met){
  TS <- Phenotype$Phenotype
  TS <- factor(TS, levels = levels(as.factor(Phenotype$Phenotype)))
  design <- model.matrix(~ 0 + TS)
  colnames(design) <- levels(TS)
  fit <- lmFit(DATA_met, design)
  
  cont.matrix <- makeContrasts(
    IDHm_vs_CD34 = CD34 - IDHm,
    IDHm_vs_WT = IDHm - WT,
    WT_vs_CD34 = WT - CD34,
    FLT3m_vs_WT = FLT3m - WT,
    levels = design
  )
  fit.contrast <- contrasts.fit(fit, cont.matrix)
  efit.contrast <- eBayes(fit.contrast)
  return(efit.contrast)
}



Filter_comparison_methylation <- function(comparison, DATA_met, efit.contrast=FALSE, Phenotype=0){
  if (type(efit.contrast)!="list"){
    efit.contrast <- Look_at_differential_methylation(Phenotype, DATA_met)
  }
  filtered_DMP <- data.frame(logFC = efit.contrast[["coefficients"]][,comparison], pvalue = efit.contrast[["p.value"]][,comparison]) %>%
    dplyr::filter(., logFC > 0.3 | logFC < -0.3) %>%
    dplyr::filter(., pvalue < 0.01) %>%
    merge(., DATA_met, by.x = 0, by.y = 0)
  rownames(filtered_DMP) <- filtered_DMP$Row.names
  filtered_DMP <- filtered_DMP[,-1]
  return(filtered_DMP)
}






Focus_global_differencial_methylation <- function(match_hit_CpGs_Blueprint_promoter,
                                                  Illumina_annotation_promoter,
                                                  Phenotype,
                                                  DATA_methylation,
                                                  comparison,
                                                  fill_colour,
                                                  CpGs_DM,
                                                  save = FALSE){
  
  CpGs_of_Interest <- Find_CpGs_of_Interest(gene_cpgs_association = match_hit_CpGs_Blueprint_promoter, cpg_annotation = Illumina_annotation_promoter)  %>%
    .[. %in% CpGs_DM]
  
  Methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, DATA_methylation)
  
  nb_cpgs <- count_cpgs(Methylation)
  
  if(nb_cpgs==0){
    print(paste0("No Differentially methylated cpgs found in promoter of ", gene))
    return(FALSE)
  }
  
  Aggregated_CpGs <- Create_DATA_jitter(Methylation, Phenotype)
  
  Beta_value <- Create_DATA_boxplot(Methylation, Phenotype)
  
  text_cpg <- paste0("Number of cpgs\nin promoter:\n", nb_cpgs)
  
  title <- make_title(diff = TRUE,
                      comparison = comparison)
  
  filename <- Create_filename(diff = TRUE,
                              comparison = comparison)
  
  nb_sample <- Nb_sample_function(Phenotype)
  
  nb_phenot <- length(levels(as.factor(Phenotype$Phenotype)))
  
  text_cpg <- create_text_cpg(diff = TRUE, rownames(Methylation))
  
  Create_graph(DATA_boxplot = Beta_value,
               DATA_jitter = Aggregated_CpGs,
               Phenotype = Phenotype,
               nb_phenot = nb_phenot,
               title = title,
               text_cpg = text_cpg,
               nb_sample = nb_sample,
               fill_colour = fill_colour,
               filename = filename,
               save = save)
}




Focus_global_methylation <- function(match_hit_CpGs_Blueprint_promoter,
                                     Illumina_annotation_promoter,
                                     Phenotype,
                                     DATA_methylation,
                                     comparison,
                                     fill_colour,
                                     Change_directory="",
                                     save = FALSE,
                                     CD34_DATA = FALSE){
  
  CpGs_of_Interest <- Find_CpGs_of_Interest(gene_cpgs_association = match_hit_CpGs_Blueprint_promoter, cpg_annotation = Illumina_annotation_promoter)
  
  Methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, DATA_methylation)
  
  nb_cpgs <- count_cpgs(Methylation)
  
  Aggregated_CpGs <- Create_DATA_jitter(Methylation, Phenotype)
  
  Beta_value <- Create_DATA_boxplot(Methylation, Phenotype)
  
  if (type(CD34_DATA) != "logical") {
    Phenotype_CD_34 <- data.frame("Phenotype" = rep("CD34", length(colnames(CD34_DATA))))
    CD34_methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, CD34_DATA)
    print(colnames(CD34_methylation))
    Beta_value_CD34 <- Create_DATA_boxplot(CD34_methylation, data.frame("Phenotype" = rep("CD34", length(colnames(CD34_methylation)))))
    Beta_value <- rbind(Beta_value, Beta_value_CD34)
    Aggregated_CpGs_CD34 <- Create_DATA_jitter(CD34_methylation, data.frame("Phenotype" = rep("CD34", length(colnames(CD34_methylation)))))
    Aggregated_CpGs <- rbind(Aggregated_CpGs, Aggregated_CpGs_CD34)
    Phenotype <- data.frame("Phenotype" = Phenotype$Phenotype) %>%
      rbind(., Phenotype_CD_34)
  }
  
  text_cpg <- paste0("Number of cpgs:\n", nb_cpgs)
  
  title <- "Global methylation"
  
  filename <- Create_filename(diff = TRUE,
                              comparison = comparison, 
                              Change_directory = Change_directory)
  
  nb_sample <- Nb_sample_function(Phenotype)
  
  nb_phenot <- length(levels(as.factor(Phenotype$Phenotype)))
  
  Create_graph(DATA_boxplot = Beta_value,
               DATA_jitter = Aggregated_CpGs,
               Phenotype = Phenotype,
               nb_phenot = nb_phenot,
               title = title,
               text_cpg = text_cpg,
               nb_sample = nb_sample,
               fill_colour = fill_colour,
               filename = filename,
               save = save,
               comparison = comparison)
}




Focus_Gene_differential_methylation <- function(gene,
                                                match_hit_CpGs_Blueprint_promoter,
                                                Illumina_annotation_promoter,
                                                Phenotype,
                                                DATA_methylation,
                                                CpGs_DM,
                                                comparison,
                                                fill_colour,
                                                save = FALSE){
  
  CpGs_of_Interest <- Find_CpGs_of_Interest(gene, match_hit_CpGs_Blueprint_promoter, Illumina_annotation_promoter) %>%
    .[. %in% CpGs_DM]
  
  
  Methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, DATA_methylation)
  
  nb_cpgs <- count_cpgs(Methylation)
  
  if(nb_cpgs==0){
    print(paste0("No Differentially methylated cpgs found in promoter of ", gene))
    return(FALSE)
  }
  
  Aggregated_CpGs <- Create_DATA_jitter(Methylation, Phenotype)
  
  Beta_value <- Create_DATA_boxplot(Methylation, Phenotype)
  
  text_cpg <- paste0("Number of cpgs\nin the promoter of\n", gene, ":\n", nb_cpgs)
  
  title <- make_title(diff = TRUE,
                      gene = gene,
                      comparison = comparison)
  
  filename <- Create_filename(diff = TRUE,
                              gene = gene,
                              comparison = comparison)
  
  nb_sample <- Nb_sample_function(Phenotype)
  
  nb_phenot <- length(levels(as.factor(Phenotype$Phenotype)))
  
  text_cpg <- create_text_cpg(diff = TRUE, rownames(Methylation), gene)
  
  Create_graph(DATA_boxplot = Beta_value,
               DATA_jitter = Aggregated_CpGs,
               Phenotype = Phenotype,
               nb_phenot = nb_phenot,
               title = title,
               text_cpg = text_cpg,
               nb_sample = nb_sample,
               fill_colour = fill_colour,
               filename = filename, 
               save = save)
}




Focus_Gene_global_methylation <- function(gene,
                                          match_hit_CpGs_Blueprint_promoter,
                                          Illumina_annotation_promoter,
                                          Phenotype,
                                          DATA_methylation,
                                          comparison,
                                          fill_colour,
                                          Change_directory="", 
                                          save = FALSE) {
  
  CpGs_of_Interest <- Find_CpGs_of_Interest(gene, match_hit_CpGs_Blueprint_promoter, Illumina_annotation_promoter)
  
  Methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, DATA_methylation)
  
  nb_cpgs <- count_cpgs(Methylation)
  
  Aggregated_CpGs <- Create_DATA_jitter(Methylation, Phenotype)
  
  Beta_value <- Create_DATA_boxplot(Methylation, Phenotype)
  
  text_cpg <- paste0("Number of cpgs\nin the promoter of\n", gene, ":\n", nb_cpgs)
  
  title <- make_title(diff = FALSE,
                      gene = gene,
                      comparison = comparison)
  
  filename <- Create_filename(diff = FALSE,
                              gene = gene,
                              comparison = comparison,
                              Change_directory = Change_directory)
  
  nb_sample <- Nb_sample_function(Phenotype)
  
  nb_phenot <- length(levels(as.factor(Phenotype$Phenotype)))
  
  text_cpg <- create_text_cpg(diff = FALSE, rownames(Methylation), gene)
  
  Create_graph(DATA_boxplot = Beta_value,
               DATA_jitter = Aggregated_CpGs,
               Phenotype = Phenotype,
               nb_phenot = nb_phenot,
               title = title,
               text_cpg = text_cpg,
               nb_sample = nb_sample,
               fill_colour = fill_colour,
               comparison = comparison,
               filename = filename,
               save = save)
}



Focus_gene_promoter_neighborhood_methylation <- function(gene,
                                                         match_hit_CpGs_Blueprint_promoter,
                                                         Illumina_annotation_promoter,
                                                         Phenotype,
                                                         DATA_methylation,
                                                         comparison,
                                                         fill_colour,
                                                         chromatin_annotation,
                                                         chromatin_network, 
                                                         gene_annotation,
                                                         Change_directory = "",
                                                         save = FALSE) {
  
  CpGs_of_Interest <- Find_CpGs_of_Interest_in_neighborhood_of_gene(gene, match_hit_CpGs_Blueprint_promoter, Illumina_annotation_promoter, chromatin_annotation, chromatin_network, gene_annotation, DATA_methylation)
  
  Methylation <- Add_Betavalue_to_cpgs_of_interest(CpGs_of_Interest, DATA_methylation)
  
  nb_cpgs <- count_cpgs(Methylation)
  
  if(nb_cpgs==0){
    print(paste0("No Cpgs found in neighbor of promoter of ", gene))
    return(FALSE)
  }
  
  Aggregated_CpGs <- Create_DATA_jitter(Methylation, Phenotype)
  
  Beta_value <- Create_DATA_boxplot(Methylation, Phenotype)
  
  text_cpg <- paste0("Number of cpgs\nin the promoter of\n", gene, ":\n", nb_cpgs)
  
  title <- make_title(diff = FALSE,
                      gene = gene,
                      comparison = comparison)
  
  filename <- Create_filename(diff = FALSE,
                              gene = gene,
                              comparison = comparison,
                              neighbor = TRUE,
                              Change_directory = Change_directory)
  
  nb_sample <- Nb_sample_function(Phenotype)
  
  nb_phenot <- length(levels(as.factor(Phenotype$Phenotype)))
  
  text_cpg <- create_text_cpg(diff = FALSE, rownames(Methylation), gene)
  
  Create_graph(DATA_boxplot = Beta_value,
               DATA_jitter = Aggregated_CpGs,
               Phenotype = Phenotype,
               nb_phenot = nb_phenot,
               title = title,
               text_cpg = text_cpg,
               nb_sample = nb_sample,
               fill_colour = fill_colour,
               comparison = comparison,
               filename = filename,
               save = save)
}





Find_CpGs_of_Interest_in_neighborhood_of_gene <- function(gene, gene_cpgs_association, cpg_annotation, chromatin_annotation, chromatin_network, gene_annotation, DATA_methylation){
  Gene_promoter_cpg <- Find_CpGs_of_Interest(gene, gene_cpgs_association, cpg_annotation)
  
  
  focused_gene_fragment <- chromatin_annotation %>%
    dplyr::filter(., CpG %in% Gene_promoter_cpg)
  
  neighborhood <- chromatin_network %>%
    dplyr::filter(., IDbait %in% focused_gene_fragment$ID | IDoe %in% focused_gene_fragment$ID) %>%
    .[,c(1:3,11, 5:8,12, 10)]
  colnames(colnames(neighborhood) <- rep(c("chr", "start", "end", "ID", "Gene_name"), 2))
  
  Network_nodes <- rbind(neighborhood[,c(1:5)], neighborhood[,c(6:10)]) %>%
    merge(., chromatin_annotation, by.x = "ID", by.y = "ID") %>%
    dplyr::select(., c("ID", "CpG")) %>%
    merge(., gene_annotation, by.x = "CpG", by.y = "CpG") %>%
    dplyr::filter(., Blueprint_gene_names != gene & Illumina_Gene_name != gene) %>%
    dplyr::select(., CpG)
  return(Network_nodes$CpG)
}






Find_CpGs_of_Interest <- function(gene = "", gene_cpgs_association, cpg_annotation){
  if(gene == ""){
    CpGs <- dplyr::select(gene_cpgs_association, CpG)
    return(CpGs$CpG)
  }
  focus_gene_cpgs_association <- dplyr::filter(gene_cpgs_association, Blueprint_gene_names == gene) %>%
    dplyr::select(., CpG)
  CpGs <- focus_gene_cpgs_association$CpG
  if(length(focus_gene_cpgs_association$CpG)==0){
    focus_gene_cpgs_association <- dplyr::filter(cpg_annotation, UCSC_RefGene_Name == gene) %>%
      dplyr::select(., Name)
    CpGs <- focus_gene_cpgs_association$Name
  }
  return(CpGs)
}



Add_Betavalue_to_cpgs_of_interest <- function(CpGs, DATA_methylation){
  DATA_methylation <- DATA_methylation %>%
    .[rownames(.) %in% CpGs,]
  return(DATA_methylation)
}




count_cpgs <- function(DATA_met){
  return(length(rownames(DATA_met)))
}




Create_DATA_jitter <- function(Methylation_values, Phenotype){
  Aggregated <- data.frame("Beta_values" = Methylation_values %>%
                             as.matrix(.) %>%
                             colMedians(.),
                           "Phenotype" = Phenotype$Phenotype) %>%
    unique(.)
  return(Aggregated)
}




Create_DATA_boxplot <- function(Methylation_values, Phenotype){
  Beta_value <- data.frame("Beta_values" = c(Methylation_values) %>%
                             unlist(.) %>%
                             na.omit(.),
                           "Phenotype" = rep(Phenotype$Phenotype, each = length(rownames(Methylation_values)))) %>%
    unique(.)
  return(Beta_value)
}




make_title <- function(diff = FALSE, gene = "", comparison){
  if (diff){
    title <- paste0("Differential")
    if(gene != ""){
      title <- paste0(title, " ", gene, " promoter methylation between ", comparison[1], " and ", comparison[2])
    }else{
      title <- paste0("Global differential promoter methylation between ", comparison[1], " and ", comparison[2])
    }
  }else{
    title <- paste0("Global")
    if(gene != ""){
      title <- paste0(title, " ", gene, " promoter methylation")
    }else{
      title <- paste0(title, " promoter methylation")
    }
  }
  return(title)
}




Create_filename <- function(diff, gene = "", comparison, neighbor = FALSE, Change_directory = ""){
  if (diff){
    filename <- paste0("Differential_methylation/")
    if(gene != ""){
      filename <- paste0(filename, gene, " promoter methylation between ", comparison[1], " and ", comparison[2], ".png")
    }else{
      filename <- paste0("Global differential promoter methylation between ", comparison[1], " and ", comparison[2], ".png")
    }
  }else{
    filename <- paste0("Global_methylation/")
    if(gene != ""){
      filename <- paste0(filename, gene, " promoter methylation.png")
    }else{
      filename <- paste0(filename, "Promoter global methylation.png")
    }
  }
  if(neighbor){
    filename <- paste0("../Results/",Change_directory,"Neighbor_methylation/", filename)
  }else{
    filename <- paste0("../Results/",Change_directory,"Specific_gene_focus/", filename)
  }
  return(filename)
}




create_text_cpg <- function(diff = FALSE, cpgs, gene = ""){
  if(diff == TRUE){
    text <- "Number of CpGs\ndifferentially\nmethylated"
    if(gene != ""){
      text <- paste0(text, " in\nthe promoter\nof ", gene, ": ")
    }else{
      text <- paste0(text, " in\npromoter: ")
    }
  }else{
    text <- paste0("Number of CpGs in\n")
    if (gene != ""){
      text <- paste0(text, " the promoter of\n", gene, ": ")
    }else{
      text <- paste0(text, "promoter :")
    }
  }
  text <- paste0(text, length(cpgs))
  return(text)
}






Nb_sample_function <- function(Phenotype){
  levels_pheno <- levels(as.factor(Phenotype$Phenotype))
  nb_sample <- list()
  for (pheno in levels_pheno){
    nb_sample[[pheno]] <- length(Phenotype[Phenotype$Phenotype == pheno, 1])
  }
  return(nb_sample)
}




Create_graph <- function(DATA_boxplot,
                         DATA_jitter,
                         Phenotype,
                         nb_phenot,
                         title,
                         text_cpg,
                         nb_sample,
                         fill_colour,
                         comparison=c(),
                         filename,
                         save){
  xlim_graph_up <- nb_phenot
  if(nb_phenot==4){
    xpos_annotations <- 5.25
  }else if(nb_phenot == 3){
    xpos_annotations <- 4.125
  }else if(nb_phenot == 2){
    xpos_annotations <- 3
  }else{
    xpos_annotations <- 2
  }
  plot <- ggplot(DATA_boxplot, aes(y=Beta_values, x = Phenotype, fill = Phenotype))+
    geom_boxplot(alpha=0.25, outlier.size = -1)+
    geom_jitter(DATA_jitter, inherit.aes = FALSE, mapping = aes(y = Beta_values, x = Phenotype), width = 0.25, alpha = 0.5, colour = "darkred")+
    coord_cartesian(xlim = c(1, xlim_graph_up), ylim = c(0,1), clip = "off")+
    scale_fill_manual(values = fill_colour)+
    annotate("text", x = 1, y = -0.11, label = paste0("n = ", nb_sample[[1]]))+
    annotate("text", x = 2, y = -0.11, label = paste0("n = ", nb_sample[[2]]))+
    #annotate("text", x = 3, y = -0.11, label = paste0("n = ", nb_sample[[3]]))+
    #annotate("text", x = 4, y = -0.11, label = paste0("n = ", nb_sample[[4]]))+
    # annotate("text", x = 5, y = -0.11, label = paste0("n = ", nb_sample[[5]]))+
    # annotate("text", x = 6, y = -0.11, label = paste0("n = ", nb_sample[[6]]))+
    annotate("text", x = xpos_annotations, y = 0.25, label = text_cpg) +
    ggtitle(title)+
    xlab("")+
    theme(plot.title = element_text(vjust = 2))+
    ylab("Beta-value of methylation")+
    scale_shape_manual(values=c(0,1,2,5,6,7))+
    geom_signif(comparisons = list(comparison), map_signif_level=TRUE)+  
    theme(plot.margin=unit(c(0.2,1.5,0.5,0.2),"cm"))
  if(save){
    ggsave(
      filename = filename,
      plot = last_plot(),
      device = NULL,
      path = NULL,
      scale = 1,
      width = NA,
      height = NA,
      units = c("in", "cm", "mm"),
      dpi = 300,
      limitsize = TRUE
    )
  }
  return(plot)
}


Focus_high_variance_CpGs <- function(DATA_met, percentage = 1) {
  CpGs_variance <- rowVars(DATA_met) %>% as.data.frame(.)
  
  colnames(CpGs_variance) <- c("Variance")
  CpGs_variance$CpGs <- rownames(CpGs_variance)
  
  CpGs_variance <- CpGs_variance[order(-CpGs_variance$Variance),]
  
  Top_CpGs <- head(CpGs_variance, round(length(CpGs_variance[,1])*(percentage/100)))
  
  DATA_met_high_variance <- DATA_met %>%
    .[Top_CpGs$CpGs,]
  return(DATA_met_high_variance)
}




Focus_low_variance_CpGs <- function(DATA_met, percentage = 1) {
  CpGs_variance <- rowVars(DATA_met) %>% as.data.frame(.)
  
  colnames(CpGs_variance) <- c("Variance")
  CpGs_variance$CpGs <- rownames(CpGs_variance)
  
  CpGs_variance <- CpGs_variance[order(CpGs_variance$Variance),]
  
  Top_CpGs <- head(CpGs_variance, round(length(CpGs_variance[,1])*(percentage/100)))
  
  DATA_met_low_variance <- DATA_met %>%
    .[Top_CpGs$CpGs,]
  return(DATA_met_low_variance)
}




Make_heatmap <- function(DATA, Phenotype, method = "pearson", title, annotation_color, kmeans_k = NA, cuttree = NA) {
  annotation_for_heatmap <- data.frame(Phenotype = Phenotype$Phenotype)
  rownames(annotation_for_heatmap) <- colnames(DATA)
  corr <- rcorr(as.matrix(DATA), type = method)$r
  colnames(corr) <- colnames(DATA)
  title <- paste0(title, " ", method)
  heatmap <- pheatmap(corr, 
                      color = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(100),
                      annotation_col = annotation_for_heatmap,
                      annotation_colors = annotation_color,
                      legend = TRUE,
                      treeheight_row = 20,
                      main = title, 
                      fontsize = 10,
                      cutree_cols = cuttree
  )
  return(heatmap)
}




Look_at_gene_with_CpGs_in_promoter <- function(list_of_cpgs, Gene_cpgs_annotation_Promoter){
  Specific_Cpgs_gene_annotation <- Gene_cpgs_annotation_Promoter %>%
    dplyr::filter(., CpG %in% list_of_cpgs)
  Gene_found <- unique(Specific_Cpgs_gene_annotation$Blueprint_gene_names)
  
  return(Gene_found)
}




Look_at_genes_connected_to_CpGs <- function(list_of_cpgs,
                                            overlap_cpgs_pchic,
                                            pchic,
                                            overlap_pchic_genes_promoter, 
                                            Gene_cpgs_annotation_Promoter) {
  
  Specific_CpGs_pchic_fragments <- overlap_cpgs_pchic %>%
    dplyr::filter(., CpG %in% list_of_cpgs)
  Pchic_of_interest <- pchic %>%
    dplyr::filter(., IDbait %in% Specific_CpGs_pchic_fragments$ID | IDoe %in% Specific_CpGs_pchic_fragments$ID)
  Pchic_of_interest$ID_bait <- Pchic_of_interest$IDbait
  Pchic_of_interest$ID_oe <- Pchic_of_interest$IDoe
  Pchic_of_interest <- dplyr::select(Pchic_of_interest, c(1:10))
  colnames(Pchic_of_interest) <- rep(c("chr", "start", "end", "ID", "Gene_name"), 2)
  Fragment_of_interest <- unique(rbind(Pchic_of_interest[, c(1:4)], Pchic_of_interest[, c(6:9)]))
  overlap_pchic_genes_promoter_of_interest <- overlap_pchic_genes_promoter %>%
    dplyr::filter(., ID %in% Fragment_of_interest$ID)
  Gene_with_cpgs_in_prom <- Look_at_gene_with_CpGs_in_promoter(list_of_cpgs, Gene_cpgs_annotation_Promoter)
  Genes_connected_to_CpGs_of_Interest <- overlap_pchic_genes_promoter_of_interest %>%
    dplyr::filter(., Blueprint_gene_names %ni% Gene_with_cpgs_in_prom) %>%
    dplyr::select(., Blueprint_gene_names) %>%
    unique(.)
  
  return(Genes_connected_to_CpGs_of_Interest$Blueprint_gene_names)
}



volcanoplot_methylation <- function(DATA, gene_association = match_hit_CpGs_Blueprint_EPIC, title){
  enhancedvolcano_data1 <- data.frame(logFC = DATA[, "logFC"],
                                      pvalue = DATA[, "P.Value"],
                                      cpgs = DATA[, "ID"])
  enhancedvolcano_data1 <- merge(x = enhancedvolcano_data1, y = gene_association, by.x = "cpgs", by.y = "CpG")
  EnhancedVolcano(toptable = enhancedvolcano_data1, 
                  lab = enhancedvolcano_data1$Blueprint_gene_names, 
                  x = "logFC", 
                  y = "pvalue",
                  FCcutoff = 0.1,
                  pCutoff = 0.0001,
                  title = title,
                  subtitle = NA,
                  legendPosition = "right",
                  subtitleLabSize = 0,
                  legendLabSize = 10,
                  ylim = c(0,10)
  )
}

volcanoplot_gene_expression <- function(DATA, title, xlim = 10, ylim = 5){
  gene <- DATA[,"ID"] %>% 
    str_split(., "[|]") %>% 
    lapply(., function(.){.[1]}) %>% 
    unlist(.)
  
  enhancedvolcano_data1 <- data.frame(logFC = DATA[, "logFC"],
                                      pvalue = DATA[, "P.Value"],
                                      gene = gene)
  
  EnhancedVolcano(toptable = enhancedvolcano_data1, 
                  lab = enhancedvolcano_data1$gene, 
                  x = "logFC", 
                  y = "pvalue",
                  FCcutoff = 2.5,
                  pCutoff = 0.01,
                  title = title,
                  subtitle = NA,
                  legendPosition = "right",
                  subtitleLabSize = 0,
                  legendLabSize = 10,
                  ylim = c(0, ylim),
                  xlim = c(-xlim, xlim)
  )
}


Filter_gene_DM_in_promoter <- function(DMPs, Gene_cpgs_annotation_Promoter){
  DMPs_and_genes <- merge(DMPs, Gene_cpgs_annotation_Promoter, by.x = "ID", by.y = "CpG") %>%
    dplyr::select(., "Blueprint_gene_names", "logFC")
  
  return(DMPs_and_genes)
}

Genes_DMR_analysis <- function(DMR_analysis, gene_universe = gene_universe_450K){
  res <- list()
  DMR_Grange <- GRanges(
    seqnames = DMR_analysis$BumphunterDMR$seqnames,
    ranges = IRanges(start = DMR_analysis$BumphunterDMR$start, end = DMR_analysis$BumphunterDMR$end),
    Value = DMR_analysis$BumphunterDMR$value,
    Length = DMR_analysis$BumphunterDMR$L,
    P.value = DMR_analysis$BumphunterDMR$p.value
  )
  message("================ Overlapping between DMR results and Blueprint annotation ==================")
  overlaps <- findOverlaps(Blueprint_Granges, DMR_Grange)
  match_hit <- data.frame(mcols(Blueprint_Granges[queryHits(overlaps),]),
                          data.frame(mcols(DMR_Grange[subjectHits(overlaps),]))) %>%
    dplyr::filter(., type == "P" & P.value < 0.01)
  message(paste0(length(rownames(match_hit)), " fragments found"))
  Genes_hyper <- match_hit %>%
    dplyr::filter(., Value < 0) %>%
    dplyr::select(., Blueprint_gene_names) %>%
    .$Blueprint_gene_names %>%
    unique(.)
  message(paste0(length(Genes_hyper), " genes hypermethylated"))
  Genes_hypo <- match_hit %>%
    dplyr::filter(., Value > 0) %>%
    dplyr::select(., Blueprint_gene_names) %>%
    .$Blueprint_gene_names %>%
    unique(.)
  message(paste0(length(Genes_hypo), " genes hypomethylated"))
  
  message("=================== Gene ontology analysis ======================")
  Genes_hypermetylated_ego <- enrichGO(
    gene = Genes_hyper,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = gene_universe
  )
  
  Genes_hypometylated_ego <- enrichGO(
    gene = Genes_hypo,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = gene_universe
  )
  
  res[["DMR overlapping"]] <- match_hit
  res[["Genes_hypermethylated"]] <- Genes_hyper
  res[["Genes_hypomethylated"]] <- Genes_hypo
  res[["GO_hypermeth"]] <- Genes_hypermetylated_ego
  res[["GO_hypometh"]] <- Genes_hypometylated_ego
  message("=========== Enhancer analysis =============")
  overlaps <- findOverlaps(PCHiC_GRange, DMR_Grange)
  match_hit <- data.frame(mcols(PCHiC_GRange[queryHits(overlaps),]),
                          data.frame(mcols(DMR_Grange[subjectHits(overlaps),])))
  
  message(paste0(length(match_hit$ID), " fragments found"))
  match_hit_hyper <- match_hit %>%
    dplyr::filter(., Value < 0)
  
  match_hit_hypo <- match_hit %>%
    dplyr::filter(., Value > 0)
  
  neighbor_hyper <- look_at_neighborhood(match_hit_hyper, pchic)
  
  neighbor_hypo <- look_at_neighborhood(match_hit_hypo, pchic)
  
  neighbor_genes_hyper <- matchit_Blueprint_Pchic %>%
    dplyr::filter(., ID %in% neighbor_hyper$ID)
  
  neighbor_genes_hypo <- matchit_Blueprint_Pchic %>%
    dplyr::filter(., ID %in% neighbor_hypo$ID)
  
  genes_enh_hyper <- neighbor_genes_hyper$Blueprint_gene_names %>% unique(.)
  genes_enh_hypo <- neighbor_genes_hypo$Blueprint_gene_names %>% unique(.)
  message(paste0(length(genes_enh_hyper), " genes hypermethylated"))
  message(paste0(length(genes_enh_hypo), " genes hypomethylated"))
  
  message("=================== Gene ontology analysis ======================")

  
  Genes_hypermetylated_ego <- enrichGO(
    gene = genes_enh_hyper,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = gene_universe
  )
  
  Genes_hypometylated_ego <- enrichGO(
    gene = genes_enh_hypo,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = gene_universe
  )
  Genes_DMR_enh_ego <- c(genes_enh_hypo, genes_enh_hyper) %>% unique(.) %>% enrichGO(.,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = gene_universe
  )
  
  res[["neighbor_genes_hyper"]] <- neighbor_genes_hyper
  res[["neighbor_genes_hypo"]] <- neighbor_genes_hypo
  res[["Genes_enh_hyper"]] <- genes_enh_hyper
  res[["Genes_enh_hypo"]] <- genes_enh_hypo
  res[["GO_enh_hypermeth"]] <- Genes_hypermetylated_ego
  res[["GO_enh_hypometh"]] <- Genes_hypometylated_ego
  res[["GO_DMR_enh"]] <- Genes_DMR_enh_ego
  
  return(res)
}


look_at_neighborhood <- function(match_hit_DMR, pchic){
  neighbor_bait <- pchic %>%
    dplyr::filter(., IDbait %in% match_hit_DMR$ID)
  neighbor_bait$ID_oe <- neighbor_bait$IDoe
  neighbor_bait <- neighbor_bait %>%
    dplyr::select(., c(6:10)) %>%
    unique(.)
  colnames(neighbor_bait) <- c("chr", "start", "end", "ID", "Gene_name")
  
  neighbor_oe <- pchic %>%
    dplyr::filter(., IDoe %in% match_hit_DMR$ID)
  neighbor_oe$ID_bait <- neighbor_bait$IDbait
  neighbor_oe <- neighbor_oe %>%
    dplyr::select(., c(1:5)) %>%
    unique(.)
  colnames(neighbor_oe) <- c("chr", "start", "end", "ID", "Gene_name")
  
  neighbor <- rbind(neighbor_bait, neighbor_oe) %>%
    unique(.)
  
  return(neighbor)
}


T_test_on_methylation_promoter <- function(BMIQ, match_hit_CpGs_Blueprint_promoter = match_hit_CpGs_Blueprint_promoter_450, Phenotype, comparison){
  message("Starting ....")
  match_hit_CpGs_Blueprint_promoter <- match_hit_CpGs_Blueprint_promoter[,c(1,4)]
  message(" == Merging Data with Blueprint annotation == ")
  BMIQ_DMP_analysis <- merge(BMIQ, match_hit_CpGs_Blueprint_promoter, by.x = 0, by.y = "CpG") %>% unique(.) %>% dplyr::filter(., Blueprint_gene_names != "")
  Gene_promoter_methylation_state <- split(BMIQ_DMP_analysis, BMIQ_DMP_analysis$Blueprint_gene_names)
  message(" == Calculating Average CpGs Bvalue per Genes == ")
  Gene_promoter_global_methylation <- lapply(Gene_promoter_methylation_state, function(x){
    nbcolumn <- length(colnames(x)) - 1
    if(length(x$Row.names) >1){
      tmp <- x[,c(2:nbcolumn)] %>% as.matrix(.) %>% colMeans2(.) %>% data.frame(.) %>% t(.)
      colnames(tmp) <- colnames(x)[c(2:nbcolumn)]
      tmp <- as.data.frame(tmp)
    }else{
      tmp <- x[,c(2:nbcolumn)]
    }
    row.names(tmp) <- "Means_Bvalues"
    tmp
  })
  
  A <- Phenotype$Phenotype == comparison[1]
  B <- Phenotype$Phenotype == comparison[2]
  
  message(" == t.test on each gene through all samples == ")
  
  test_t.student <- lapply(Gene_promoter_global_methylation, function(x){
    t.test(x[,A], x[,B])
  })
  
  Genes_significally_hyper <- list.filter(test_t.student, p.value < 0.05 & (estimate[1] - estimate[2]) > 0.1)
  Genes_significally_hypo <- list.filter(test_t.student, p.value < 0.05 & (estimate[2] - estimate[1]) > 0.1)
  
  Genes_hyper <- names(Genes_significally_hyper) %>% unique(.)
  Genes_hypo <- names(Genes_significally_hypo) %>% unique(.)
  Gene_universe <- names(test_t.student) %>% unique(.)
  
  message(paste0(length(Genes_hyper), " genes significantly hypermethylated found !"))
  message(paste0(length(Genes_hypo), " genes significantly hypomethylated found !"))
  
  res <- list()
  res[["test_t.student"]] <- test_t.student
  
  message(" == Gene ontology analysis == ")
  
  Genes_hyper_ego <- enrichGO(
    gene = Genes_hyper,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = Gene_universe
  )
  
  Genes_hypo_ego <- enrichGO(
    gene = Genes_hypo,
    keyType = "SYMBOL",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",
    pAdjustMethod = "none", 
    universe = Gene_universe
  ) 
  
  Genes_hypermet_rank <- list.select(Genes_significally_hyper, c(p.value, estimate[1], estimate[2])) %>%
    lapply(., function(gene){
      data.frame("p.value" = unlist(gene)[[1]], "mean of x" = unlist(gene)["mean of x"], "mean of y" = unlist(gene)["mean of y"], "delta" = unlist(gene)["mean of x"] - unlist(gene)["mean of y"])
    }) %>%
    bind_rows(., .id = "Gene_name") 
  
  Genes_hypomet_rank <- list.select(Genes_significally_hypo, c(p.value, estimate[1], estimate[2])) %>%
    lapply(., function(gene){
      data.frame("p.value" = unlist(gene)[[1]], "mean of x" = unlist(gene)["mean of x"], "mean of y" = unlist(gene)["mean of y"], "delta" = unlist(gene)["mean of x"] - unlist(gene)["mean of y"])
    }) %>%
    bind_rows(., .id = "Gene_name")
  
  Genes_methylation_analysis <- list.select(test_t.student, c(p.value, estimate[1], estimate[2])) %>%
    lapply(., function(gene){
      data.frame("p.value" = unlist(gene)[[1]], "mean of x" = unlist(gene)["mean of x"], "mean of y" = unlist(gene)["mean of y"], "delta" = unlist(gene)["mean of x"] - unlist(gene)["mean of y"])
    }) %>%
    bind_rows(., .id = "Gene_name")
  
  res[["Genes_methylation_analysis"]] <- Genes_methylation_analysis
  res[["Genes_significally_Hyper"]] <- Genes_hypermet_rank
  res[["Genes_significally_Hypo"]] <- Genes_hypomet_rank
  res[["Genes_universe"]] <- Gene_universe
  res[["GO_hyper"]] <- Genes_hyper_ego
  res[["GO_hypo"]] <- Genes_hypo_ego
  
  return(res)
}

Visualize_Ttest_results <- function(Ttest_Results, dataname){
  write.csv(Ttest_Results$Genes_hypermet, file = paste0(Gene_lists_folder, dataname, "_T_Hypermet.test.csv"), row.names = FALSE)
  write.csv(Ttest_Results$Genes_hypomet, file = paste0(Gene_lists_folder,  dataname, "_T_Hypomet.test.csv"), row.names = FALSE)
  write.csv(Ttest_Results$Genes_universe, file = paste0(Gene_lists_folder,  dataname, "_gene_universe.csv"), row.names = FALSE)
  if("GO_hyper" %in% names(Ttest_Results)){
    write.csv(Ttest_Results$GO_hyper@result, file = paste0(GO_folder,  dataname, "_Ttest_hyper.csv"))
    auto_stress_pathways <- Ttest_Results$GO_hyper@result %>%
      dplyr::filter(., str_detect(Description, "metabolism") | str_detect(Description, "stress"))
    write.csv(auto_stress_pathways, file = paste0(GO_folder, dataname, "_stress_metabolism_pathway_hyper.csv"))
    message(" GO hypermeth written ")
    }
  if("GO_hypo" %in% names(Ttest_Results)){
    write.csv(Ttest_Results$GO_hypo@result, file = paste0(GO_folder,  dataname, "_Ttest_hypo.csv"))
    auto_stress_pathways <- Ttest_Results$GO_hypo@result %>%
      dplyr::filter(., str_detect(Description, "metabolism") | str_detect(Description, "stress"))
    write.csv(auto_stress_pathways, file = paste0(GO_folder, dataname, "_stress_metabolism_pathway_hypo.csv"))
    message(" GO hypometh written ")
    }
}

T_test_on_methylation_enhancer <- function(Mean_CpGs_per_fragment, comparisons, Fragment_connected_per_gene){
  res <- list()
  A <- comparisons[[1]]
  B <- comparisons[[2]]
  message("T test analysis")
  T_test <- lapply(Mean_CpGs_per_fragment, function(x){
    t.test(x[,A], x[,B])
  })
  
  res[["T_Test"]] <- T_test
  
  T_test_analysis <- lapply(T_test, function(x){
    if(x[["p.value"]] < 0.05){
      if(x[["estimate"]]["mean of x"] - x[["estimate"]]["mean of y"] > 0.05){
        delta <- "Hyper"
      }else if(x[["estimate"]]["mean of x"] - x[["estimate"]]["mean of y"] < -0.05){
        delta <- "Hypo"
      }else{
        delta <- "NULL"
      }
    }else{
      delta <- "NULL"
    }
    tmp <- data.frame("p.value" = x[["p.value"]], "delta" = delta, "estimate" = x[["estimate"]]["mean of x"] - x[["estimate"]]["mean of y"])
    tmp
  })
  
  tmp <- names(T_test_analysis)
  
  T_test_analysis <- rbindlist(T_test_analysis) %>% as.data.frame(.)
  rownames(T_test_analysis) <- tmp
  res[["T_test_analysis"]] <- T_test_analysis
  
  
  res[["Enhancer"]] <- lapply(Fragment_connected_per_gene, function(x){
    x[x %in% names(Mean_CpGs_per_fragment)]
  })

  res[["Enhancer"]] <- res[["Enhancer"]][lapply(res[["Enhancer"]],length)>0]
  res[["Enhancer"]] <- lapply(res[["Enhancer"]], function(x){
    T_test_analysis[x,]
  })
  Gene_universe <- names(res[["Enhancer"]]) %>% unique(.)

  res[["Enhancer_summary"]] <- lapply(res[["Enhancer"]], function(x){
    nb_enhancer_hyper <- length(x[,"delta"][x[,"delta"] == "Hyper"])
    nb_enhancer_hypo <- length(x[,"delta"][x[,"delta"] == "Hypo"])
    nb_enhancer_unchanged <- length(x[,"delta"][x[,"delta"] == "NULL"])
    data.frame("nb_enhancer_hyper" = nb_enhancer_hyper, "nb_enhancer_hypo" = nb_enhancer_hypo, "nb_enhancer_unchanged" = nb_enhancer_unchanged, "total_enhancers" = nb_enhancer_hyper + nb_enhancer_hypo + nb_enhancer_unchanged)
  })

  res[["Genes_with_hyper_enhancers"]] <- list.filter(res[["Enhancer_summary"]], nb_enhancer_hyper > 0)

  res[["Genes_list_of_Genes_with_hyper_enhancers"]] <- res[["Genes_with_hyper_enhancers"]][list.order(res[["Genes_with_hyper_enhancers"]], (nb_enhancer_hyper))]

  res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]] <- rbindlist(res[["Genes_list_of_Genes_with_hyper_enhancers"]]) %>%
    as.data.frame(.)
  rownames(res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]) <- names(res[["Genes_list_of_Genes_with_hyper_enhancers"]])

  res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$percent_hyper <- (res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$nb_enhancer_hyper/res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$total_enhancers)*100
  res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$DM <- res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$total_enhancers - res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]$nb_enhancer_unchanged

  res[["GO"]] <- enrichGO(rownames(res[["Genes_list_of_Genes_with_hyper_enhancers_annoted"]]),
                                                        keyType = "SYMBOL",
                                                        OrgDb = "org.Hs.eg.db",
                                                        ont = "BP",
                                                        pAdjustMethod = "none",
                                                        universe = Gene_universe)
  return(res)
}

####Miguel Functions for Viper

# Function to convert a dorothea table to viper regulons

# Function extracted from dorothea code
# https://github.com/saezlab/dorothea/blob/master/R/helpers.R#L17
dorothea2viper_regulons <- function(df) {
  regulon_list <- split(df, df$tf)
  viper_regulons <- lapply(regulon_list, function(regulon) {
    tfmode <- stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  
  return(viper_regulons)
}


# Function to convert viper regulons to a dorothea table

# Function extracted from
# https://github.com/saezlab/ConservedFootprints/blob/master/src/dorothea_analysis.R#L126
viper_regulons2dorothea <- function(r) {
  res <- r %>%
    map_df(
      .f = function(i) {
        tf_target <- i$tfmode %>%
          enframe(name = "target", value = "mor") %>%
          mutate(likelihood = i$likelihood)
      },
      .id = "tf"
    )
  return(res)
}

# Function to convert dorothea database and gene expression to aracne regulons

# We need to create a network file where the first TF is the regulon with every target with a true confidence score
# in our case from Dorothea it always be 1
# Finally use it in aracne2regulon function from viper package
dorothea2aracne2viper_regulons <- function(dorothea, exprs_m) {
  dorothea_aggregation_tf <- dorothea %>%
    select(tf, target) %>%
    group_by(tf) %>%
    summarise(targets = str_c(target, collapse = ";"))
  tmp_file <- tempfile()
  for (i in 1:nrow(dorothea_aggregation_tf)) {
    tf_targets <- str_split(dorothea_aggregation_tf$targets[i], ";")[[1]]
    row <- c(dorothea_aggregation_tf$tf[i], unlist(mapply(c, tf_targets, rep(1, length(tf_targets)), SIMPLIFY = F)))
    cat(str_c(row, collapse = "\t"), "\n", file = tmp_file, append = T)
  }
  aracne_regulons <- aracne2regulon(tmp_file, exprs_m, format = "adj", verbose = F)
  file.remove(tmp_file)
  return(aracne_regulons)
}


# Main function to run msviper

run_msviper <- function(exprs_m, dorothea, use_aracne, ref, treat, ref_name, treat_name, minsize, ges.filter) {
  # First we need to generate the phenotype table (AnnotatedDataFrame)
  conditions <- rep("NA", ncol(exprs_m))
  conditions[ref] <- ref_name
  conditions[treat] <- treat_name
  names(conditions) <- colnames(exprs_m)
  conditions <- conditions[which(conditions != "NA")]
  
  phenotype <- data.frame(condition = factor(conditions))
  rownames(phenotype) <- names(conditions)
  
  phenoData <- new("AnnotatedDataFrame", data = phenotype)
  
  exprs_m <- exprs_m[,which(colnames(exprs_m) %in% rownames(phenotype))] %>% as.matrix()
  
  # Create Expression set from phenotyble table and expression matrix
  dset_viper <- ExpressionSet(assayData = exprs_m, phenoData = phenoData)
  dset_viper$sampleID <- factor(colnames(exprs_m))
  
  # Aracne can be used to estimate the mor instead using the -1, 1 from dorothea
  regulons <- NULL
  if (use_aracne) {
    regulons <- dorothea2aracne2viper_regulons(dorothea, dset_viper)
  } else {
    regulons <- dorothea2viper_regulons(dorothea)
  }
  
  # We need to create the statistics signature from the conditions
  signature <- rowTtest(dset_viper, "condition", treat_name, ref_name)
  statistics_signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) * sign(signature$statistic))[, 1]
  # Generate the null model with bootstrapping (1000 iterations)
  nullmodel <- ttestNull(dset_viper, "condition", treat_name, ref_name, per = 1000, repos = T, verbose = F)
  # Run msviper using the statistics signature, the regulons converted from dorothea table, the null model the minSize of regulon and the ges.filter
  mrs <- msviper(ges = statistics_signature, regulon = regulons, nullmodel = nullmodel, minsize = minsize, ges.filter = ges.filter, verbose = F)
  # Convert the msviper regulons to dorothea
  dorothea_mrs_regulons <- viper_regulons2dorothea(mrs$regulon) %>%
    mutate(state = ifelse(mor > 0, "activation", "inhibition"))
  # Generate a table with the TFs, the regulon size, the NES score, the pval and the pval.fdr
  mrs_table <- tibble(TF = names(mrs$es$p.value), size = mrs$es$size, nes = mrs$es$nes, pval = mrs$es$p.value, pval.fdr = p.adjust(mrs$es$p.value, method = "fdr")) %>% arrange(pval)
  
  list(mrs_table = mrs_table, mrs = mrs, regulons = dorothea_mrs_regulons)
}


# Extra function to generate the cytoscape networok from msviper result


mrs2cytoscape <- function(mrs,full.path) {
  
  all_nodes <- unique(c(mrs$regulons$tf, mrs$regulons$target))
  tnodes <- tibble(TF = all_nodes)
  all_nodes_metadata <- right_join(mrs$mrs_table, tnodes, by = "TF")
  regulons_network <- graph.data.frame(mrs$regulons, directed = T, vertices = all_nodes_metadata)
  deleteAllNetworks()
  createNetworkFromIgraph(regulons_network, "regulons_network")
  # setVisualStyle(cytoscape_id_network, 'default')
  # setVisualStyle("default")
  
  my_style <- "my_style"
  
  
  
  
  
  createVisualStyle(my_style, list())
  setNodeColorDefault("#D3D3D3", style.name = my_style)
  blue_white_red <- c("#0000FF", "#FFFFFF", "#FF0000")
  setNodeColorMapping("nes", c(min(V(regulons_network)$nes, na.rm = T), mean(V(regulons_network)$nes, na.rm = T), max(V(regulons_network)$nes, na.rm = T)), blue_white_red, style.name = my_style)
  
  setEdgeTargetArrowShapeMapping("state", c("activation", "inhibition"), c("DELTA", "T"), style.name = my_style)
  
  setEdgeColorMapping("mor", c(min(E(regulons_network)$mor, na.rm = T), mean(E(regulons_network)$mor, na.rm = T), max(E(regulons_network)$mor, na.rm = T)), blue_white_red, style.name = my_style)
  setNodeLabelMapping('id'
  )
  setVisualStyle("my_style")
  
  createColumnFilter(filter.name='null', column='pval', 0.05, 'GREATER_THAN', network = regulons_network)
  applyFilter('null', hide=T, network = regulons_network)
  
  exportImage(full.path, 'SVG', zoom=200) 
  
}


##Function for TF activities between two conditions

tf_activities_function<- function(treatment,control,path){
  
  
  data(dorothea_mm, package = "dorothea")
  
  
  regulons = dorothea_mm %>%
    filter(confidence %in% c("A", "B"))
  
  
  
  # In my case ref is LOW CDA expression
  ref_name = control
  # Versus HIGH CDA
  treat_name= treatment
  
  # Expression matrix with all samples together
  exprs_m= counts(RNAseqDatBatchCorrection$dds,normalized=T)[,RNAseqDatBatchCorrection$colData$condition%in%c(treat_name,ref_name)]
  # ref to indicate the control samples
  ref= which(RNAseqDatBatchCorrection$colData$condition[RNAseqDatBatchCorrection$colData$condition%in%c(treat_name,ref_name)]%in%c(ref_name))
  
  
  
  msviper_results <- run_msviper(exprs_m, regulons, use_aracne = T, ref, ref_name, treat_name, minsize = 4, ges.filter = T)
  
  #assay(RNAseqDatBatchCorrection$vsd_mat)
  tf_activities <- run_viper(exprs_m, regulons, 
                             options =  list(method = "scale", minsize = 4,eset.filter = FALSE, cores = 1, verbose = FALSE))
  
  
  tf_activities=t(tf_activities) %>% as.data.frame()
  
  
  tf_activities=aggregate(x = tf_activities,                # Specify data column
                          by = list(RNAseqDatBatchCorrection$colData$condition[match(rownames(tf_activities),RNAseqDatBatchCorrection$colData$Sample)]),              # Specify group indicator
                          FUN = mean) %>% column_to_rownames(var = "Group.1")
  
  tf_activities=as.data.frame(t(tf_activities))
  
  tf_activities$pval=  msviper_results$mrs_table$pval[match(rownames(tf_activities),msviper_results$mrs_table$TF)]
  
  tf_activities[is.na(tf_activities)]=1
  
  tf_activities$pval=ifelse(tf_activities$pval<0.05,"S","NS")
  
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  
  pvalue_col_fun = circlize::colorRamp2(unique(tf_activities$pval), c("red", "blue")) 
  
  
  ha = rowAnnotation(pval = tf_activities$pval,annotation_legend_param = (pval = list(
    title = "pval",
    at = c("S", "NS"),
    labels = c("S", "NS")
  )))
  
  
  ht= Heatmap(as.matrix(tf_activities %>% select(-pval)), col = viridis::inferno(50), column_title = paste0("Transcription Factor activities ",treat_name," vs ",ref_name), row_gap = unit(5, "mm"),column_gap =unit(5, "mm"),
              left_annotation = rowAnnotation(pval = tf_activities$pval,annotation_legend_param = (pval = list(
                title = "pval",
                at = c("S", "NS"),
                labels = c("S", "NS")
              ))))
  
  
  mrs2cytoscape(mrs=msviper_results,full.path = paste0("GitHub/Epigenomic_integration/Results/",treat_name,"_",ref_name))
  
  return(ht) 
}

Gene_set_analysis <- function(Analysis, organism = org.Hs.eg.db){
  gene_list <- Analysis$logFC
  names(gene_list) <- Analysis$ID
  gene_list<-na.omit(gene_list)
  
  gene_list <- sort(gene_list, decreasing = T)

  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "SYMBOL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  return(gse)
}




KEGG_pathway_analysis <- function(Analysis, organism){
  
  gene_list <- Analysis$logFC
  names(gene_list) <- Analysis$ID
  
  ids<-bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
  # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
  dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = Analysis[Analysis$ID %in% dedup_ids$SYMBOL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2 = merge(df2, dedup_ids, by.x = "ID", by.y = "SYMBOL")
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$logFC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$ENTREZID
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
  
  kegg_organism <- "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 nPerm        = 10000,
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pAdjustMethod = "none",
                 keyType       = "ncbi-geneid")
  return(kk2)
}


Make_factor <- function(samples, 
                        samplesheet, 
                        index_of_samplesheet, 
                        index_name, 
                        First_level, 
                        Second_level, 
                        condition_A_name, 
                        condition_B_name, 
                        mutation_to_ignore = c("ABCDEF"), 
                        column_of_mutation = 0,
                        comparing_mutation = F){
  if(comparing_mutation){
    First_level_samples <- samplesheet[which(str_detect(samplesheet[,index_of_samplesheet], pattern = First_level)), index_name]
    Second_level_samples <- samplesheet[which(!str_detect(samplesheet[,index_of_samplesheet], pattern = First_level)), index_name]
  }else{
    if(column_of_mutation != 0){
      First_level_samples <- samplesheet[which(samplesheet[,index_of_samplesheet] %in% First_level & 
                                                 !duplicated(str_split(samplesheet[,column_of_mutation], pattern=","), mutation_to_ignore)), index_name]
      Second_level_samples <- samplesheet[which(samplesheet[,index_of_samplesheet] %in% Second_level & 
                                                  !duplicated(str_split(samplesheet[,column_of_mutation], pattern=","), mutation_to_ignore)), index_name]
    }else{
    First_level_samples <- samplesheet[which(samplesheet[,index_of_samplesheet] %in% First_level), index_name]
    Second_level_samples <- samplesheet[which(samplesheet[,index_of_samplesheet] %in% Second_level), index_name]
    }
  }
  res <- ifelse(samples %in% First_level_samples, condition_A_name, ifelse(samples %in% Second_level_samples, condition_B_name, "NA"))
  return(res)
}
