"%ni%" <- Negate("%in%")
source("GitHub/Koichi_gene_expression_git/Koichi_gene_expression_analyses/packages_methylation.R")

DMR_analysis <- function(DNAmeth, Factor_A, Factor_A_name, Factor_B, Factor_B_name, Phenotype){
  Samples <- colnames(DNAmeth)
  Factors <- sapply(Samples, function(sample){
    if(sample %ni% Phenotype$Baseline_Sample){
      "NA"
    }else{
      patient_response <- Phenotype[which(Phenotype$Baseline_Sample == sample), "Best_response"]
      ifelse(patient_response %in% Factor_A, Factor_A_name, ifelse(patient_response %in% Factor_B, Factor_B_name, "NA"))
    }
  }) %>% unlist()
  DNAmeth <- DNAmeth[,c(Factors != "NA" & Factors != "")]
  Factors <- Factors[Factors != "NA" & Factors != ""]
  res <- champ.DMR(beta = DNAmeth, pheno = Factors, arraytype = "EPIC", cores = 8)
  return(res)
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


Make_factor_DNA_meth <- function(samples, Phenotype){
  res <- sapply(samples, function(s){
    if(s %in% Phenotype$Baseline_Sample){
      p <- "B"
      resp <- Phenotype[which(Phenotype$Baseline_Sample == s), "Best_response"][1]
      if(resp %in% c("CR", "CRi")){
        p <- paste(p, "R", sep = "_")
      }else if(resp %in% c("SD", "PD")){
        p <- paste(p, "NR", sep = "_")
      }else{
        p <- paste(p, "OR", sep = "_")
      }
    }else{
      p <- "PostT"
      resp <- Phenotype[which(Phenotype$Post_treatment_sample == s), "Best_response"][1]
      if(resp %in% c("CR", "CRi")){
        p <- paste(p, "R", sep = "_")
      }else if(resp %in% c("SD", "PD")){
        p <- paste(p, "NR", sep = "_")
      }else{
        p <- paste(p, "OR", sep = "_")
      }
    }
    p
  })
  res <- unname(res)
  return(res)
}


TF_target_methylation_analysis <- function(TF, reg, pch, cpgs_pchic_overlap, DMP){
  Targets <- reg[which(reg$tf == TF), "target"] %>% c() %>% unlist() %>% as.list()
  names(Targets) <- reg[which(reg$tf == TF), "target"] %>% c() %>% unlist()
  Targets <- lapply(names(Targets), function(tar){
    Target_methylation_analysis(tar, pch, cpgs_pchic_overlap, DMP)
  })
  names(Targets) <- reg[which(reg$tf == TF), "target"] %>% c() %>% unlist()
  return(Targets)
}

Target_methylation_analysis <- function(t, pch, cpgs_pchic_overlap, DMP){
  fragments <- pch[which(str_detect(pch$Name_bait, t)),]
  fragments_prom <- fragments[,"IDbait"]
  fragments_oe <- fragments[,"IDoe"]
  
  Promoter <- sapply(fragments_prom, function(prom){
    Cpgs <- cpgs_pchic_overlap[which(cpgs_pchic_overlap$ID == prom), "CpG"]
    ncpgs <- length(Cpgs)
    ndiffcpgs <- DMP[which(rownames(DMP) %in% Cpgs),]
    nrow(ndiffcpgs)/ncpgs
  }) %>% mean()
  
  Enhancer <- sapply(fragments_oe, function(oe){
    Cpgs <- cpgs_pchic_overlap[which(cpgs_pchic_overlap$ID == oe), "CpG"]
    ncpgs <- length(Cpgs)
    if(ncpgs == 0){
      0
    }else{
      ndiffcpgs <- sum(rownames(DMP) %in% Cpgs)
      ndiffcpgs/ncpgs
    }
  }) %>% mean()
  res1 <- data.frame("Promoter" = Promoter, "Enhancer" = Enhancer)
  return(res1)
}

test <- TF_target_methylation_analysis("RXRA", regulons, Pchic_data[["pchic"]], Overlaps_list[["CPGS_PCHIC"]], DMP_analysis[["NR_R_filtered"]])

Target_methylation_analysis("RUNX1", Pchic_data[["pchic"]], Overlaps_list[["CPGS_PCHIC"]], DMP_analysis[["NR_R_filtered"]])




