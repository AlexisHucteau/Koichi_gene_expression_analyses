"%ni%" <- Negate("%in%")
source("packages.R")













































Make_factor <- function(Samplesheet = Clinical_patient_data, 
                        Samples_names, 
                        Mutations_to_ignore = 0, 
                        Clinical_outcome_A, 
                        Clinical_name_A,
                        Clinical_outcome_B,  
                        Clinical_name_B,
                        Clinical_outcome_C, 
                        Clinical_name_C){
  # Function made for Clinical_patient_data 
  # Create a factor that can be used for Differential_analysis function
  # Samplesheet = Clinical_patient_data
  # Mutations_to_ignore: A vector of mutations that have to be taken into account (type 0 no mutations to ignore)
  # Clinical_outcome_A: A vector of best response corresponding to the phenotype A
  # Clinical_outcome_B: A vector of best response corresponding to the phenotype B
  # Clinical_outcome_C: A vector of best response corresponding to the phenotype C
  # Baseline_sample: A logical variable indicating whether Baseline samples are taken or not
  # Relapse_sample: A logical variable indicating whether Relapse samples are taken or not
  
  # Phenotype_A: The name of the first phenotype that have to be compared to
  # Phenotype _B: The name of the second phenotype that have to be compared to
  # Clinical_outcome_comparison: A logical variable indicating whether clinical outcome are taken into account
  # Baseline: 
  # Relapse: A logical variable indicating whether Relapse samples are taken or not
  if(typeof(Mutations_to_ignore) != "double"){
    Mutations_samples <- Samplesheet[which(duplicated(str_split(Samplesheet$mutations, pattern=","), Mutations_to_ignore)),] %>% 
      c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
      na.omit()
    Mutations_factor <- factor(ifelse(Samples_names %in% Mutations_samples, "Mut", "WT"))
  }else{
    Mutations_factor <- factor(rep("", length(Samples_names)))
  }
  
  Clinical_outcome_A <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_A),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  Clinical_outcome_B <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_B),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  Clinical_outcome_C <- Samplesheet[which(Samplesheet$Best_response %in% Clinical_outcome_C),] %>% 
    c(.$Baseline_RNAseq_data, .$Relapse_RNAseq_data) %>% 
    na.omit()
  
  Clinical_outcome <- factor(ifelse(Samples_names %in% Clinical_outcome_A, Clinical_name_A,
                                    ifelse(Samples_names %in% Clinical_outcome_B, Clinical_name_B, 
                                           ifelse(Samples_names %in% Clinical_outcome_C, Clinical_name_C, ""))))
  Sample_timing <- factor(ifelse(Samples_names %in% Samplesheet$Baseline_RNAseq_data, "B", "REL"))
  if(typeof(Mutations_to_ignore) != "double"){
    Final_factor <- paste(Mutations_factor, Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }else{
    Final_factor <- paste(Clinical_outcome, Sample_timing, sep = ".") %>% as.factor()
  }
  
  return(Final_factor)
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









run_msviper <- function(exprs_m, dorothea, use_aracne, ref_samples_all_data, ref_name, treat_samples, treat_name, minsize, ges.filter) {
  # First we need to generate the phenotype table (AnnotatedDataFrame)
  
  factor_merged <- ref_samples_all_data | treat_samples
  
  
  ref_samples <- ref_samples_all_data[factor_merged]
  exprs_m <- exprs_m[,factor_merged]
  conditions <- rep(ref_name, ncol(exprs_m))
  conditions[!ref_samples] <- treat_name
  phenotype <- data.frame(condition = factor(conditions))
  phenotype$condition[!ref_samples] <- treat_name
  rownames(phenotype) <- colnames(exprs_m)
  phenoData <- new("AnnotatedDataFrame", data = phenotype)
  # Create Expression set from phenotyble table and expression matrix
  dset_viper <- ExpressionSet(assayData = as.matrix(exprs_m), phenoData = phenoData)
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




Do_MS_viper_analysis <- function(exprs_m,
                                 dorothea, 
                                 ref_samples, 
                                 ref_name, 
                                 treat_samples, 
                                 treat_name){
  MS_vip <- run_msviper(exprs_m, dorothea, use_aracne=T, ref_samples, ref_name, treat_samples, treat_name, minsize=4, ges.filter=T)
  combi <- msviperCombinatorial(MS_vip$mrs, regulators = 25, verbose = FALSE)
  Synergy <- msviperSynergy(combi, verbose = FALSE)
  message(paste0(treat_name, " vs ", ref_name, " DONE!"))
  return(list("MSVIPER" = MS_vip, 
              "Synergy" = Synergy))
}
