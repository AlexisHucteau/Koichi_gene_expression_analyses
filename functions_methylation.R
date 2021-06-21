"%ni%" <- Negate("%in%")
source("packages_methylation.R")

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
