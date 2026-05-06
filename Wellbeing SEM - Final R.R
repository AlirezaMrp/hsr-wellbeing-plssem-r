############################################################
# PLS-SEM Analysis: HSR, Well-Being and Business Outcomes
############################################################
# Required packages
packages <- c(
  "tidyverse",
  "seminr",
  "psych",
  "ggplot2",
  "DiagrammeR"
)

installed <- packages %in% rownames(installed.packages())

if (any(!installed)) {
  install.packages(packages[!installed])
}

lapply(packages, library, character.only = TRUE)

# Load data
hsr_data <- read.csv("data/hsr_wellbeing_survey_anonymised.csv",
                     stringsAsFactors = FALSE)


# Measurement model
hsr_mm <- constructs(
  composite("Service_Quality",          multi_items("SQ",  1:5)),
  composite("Travel_Experience_Affect", multi_items("TEA", 1:4)),
  composite("Cognitive_Evaluation",     multi_items("CET", 1:4)),
  composite("Subjective_Well_Being",    multi_items("SWB", 1:4)),
  composite("Customer_Loyalty",         multi_items("CL",  1:4)),
  composite("Willingness_to_Pay",       multi_items("WTP", 1:3))
)

# Structural model with two added paths: CET -> CL and TEA -> CL
hsr_sm <- relationships(
  paths(
    from = c("Service_Quality",
             "Travel_Experience_Affect",
             "Cognitive_Evaluation"),
    to = "Subjective_Well_Being"
  ),
  paths(
    from = "Subjective_Well_Being",
    to = c("Customer_Loyalty",
           "Willingness_to_Pay")
  ),
  paths(
    from = "Service_Quality",
    to = c("Travel_Experience_Affect",
           "Cognitive_Evaluation")
  ),
  paths(
    from = "Travel_Experience_Affect",
    to = "Cognitive_Evaluation"
  ),
  paths(
    from = c("Cognitive_Evaluation",
             "Travel_Experience_Affect"),
    to = "Customer_Loyalty"
  )
)

# Estimate PLS-SEM model
hsr_pls <- estimate_pls(
  data = hsr_data,
  measurement_model = hsr_mm,
  structural_model  = hsr_sm
)

# Bootstrap
set.seed(123)
boot_hsr <- bootstrap_model(hsr_pls, nboot = 5000)


############################################################
# Main Outputs
############################################################

# Full model summary
summary(hsr_pls)

# Bootstrapped results: paths, loadings, HTMT and total paths
summary(boot_hsr)

# Reliability and validity: rhoC and AVE
rhoC_AVE(hsr_pls)

# rhoA, alpha, rhoC and AVE are also shown in:
summary(hsr_pls)


############################################################
# HTMT Table
############################################################

construct_items <- list(
  Service_Quality          = paste0("SQ",  1:5),
  Travel_Experience_Affect = paste0("TEA", 1:4),
  Cognitive_Evaluation     = paste0("CET", 1:4),
  Subjective_Well_Being    = paste0("SWB", 1:4),
  Customer_Loyalty         = paste0("CL",  1:4),
  Willingness_to_Pay       = paste0("WTP", 1:3)
)

compute_htmt <- function(dataframe, constructs_list) {
  all_items <- unlist(constructs_list)
  cor_mat <- abs(cor(dataframe[, all_items], use = "pairwise.complete.obs"))
  
  construct_names <- names(constructs_list)
  n <- length(construct_names)
  
  htmt_mat <- matrix(
    NA,
    nrow = n,
    ncol = n,
    dimnames = list(construct_names, construct_names)
  )
  
  for (i in 1:n) {
    htmt_mat[i, i] <- 1
    
    for (j in i:n) {
      if (i != j) {
        items_i <- constructs_list[[i]]
        items_j <- constructs_list[[j]]
        
        heterotrait_mean <- mean(cor_mat[items_i, items_j])
        
        monotrait_i <- cor_mat[items_i, items_i]
        monotrait_i <- monotrait_i[lower.tri(monotrait_i)]
        
        monotrait_j <- cor_mat[items_j, items_j]
        monotrait_j <- monotrait_j[lower.tri(monotrait_j)]
        
        htmt_value <- heterotrait_mean /
          sqrt(mean(monotrait_i) * mean(monotrait_j))
        
        htmt_mat[i, j] <- htmt_value
        htmt_mat[j, i] <- htmt_value
      }
    }
  }
  
  round(htmt_mat, 3)
}

htmt_table <- compute_htmt(hsr_data, construct_items)
htmt_table


############################################################
# R2 Values
############################################################

r2_values <- summary(hsr_pls)$rSquared
r2_values


############################################################
# f2 Effect Sizes
############################################################


fSquared(hsr_pls, iv = "Service_Quality", dv = "Subjective_Well_Being")
fSquared(hsr_pls, iv = "Travel_Experience_Affect", dv = "Subjective_Well_Being")
fSquared(hsr_pls, iv = "Cognitive_Evaluation", dv = "Subjective_Well_Being")

fSquared(hsr_pls, iv = "Subjective_Well_Being", dv = "Customer_Loyalty")
fSquared(hsr_pls, iv = "Subjective_Well_Being", dv = "Willingness_to_Pay")


############################################################
# VIF from Construct Scores
############################################################

scores <- as.data.frame(hsr_pls$construct_scores)

calculate_vif <- function(data, dependent, predictors) {
  vif_results <- data.frame(
    Dependent = character(),
    Predictor = character(),
    VIF = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (pred in predictors) {
    other_preds <- predictors[predictors != pred]
    
    if (length(other_preds) > 0) {
      formula_vif <- as.formula(
        paste(pred, "~", paste(other_preds, collapse = " + "))
      )
      
      vif_model <- lm(formula_vif, data = data)
      r2 <- summary(vif_model)$r.squared
      vif <- 1 / (1 - r2)
    } else {
      vif <- NA
    }
    
    vif_results <- rbind(
      vif_results,
      data.frame(
        Dependent = dependent,
        Predictor = pred,
        VIF = round(vif, 3)
      )
    )
  }
  
  return(vif_results)
}

vif_swb <- calculate_vif(
  scores,
  dependent = "Subjective_Well_Being",
  predictors = c("Service_Quality",
                 "Travel_Experience_Affect",
                 "Cognitive_Evaluation")
)

vif_cet <- calculate_vif(
  scores,
  dependent = "Cognitive_Evaluation",
  predictors = c("Service_Quality",
                 "Travel_Experience_Affect")
)

vif_cl <- calculate_vif(
  scores,
  dependent = "Customer_Loyalty",
  predictors = c("Subjective_Well_Being",
                 "Cognitive_Evaluation",
                 "Travel_Experience_Affect")
)

vif_wtp <- calculate_vif(
  scores,
  dependent = "Willingness_to_Pay",
  predictors = c("Subjective_Well_Being")
)

vif_table <- rbind(vif_swb, vif_cet, vif_cl, vif_wtp)
vif_table


############################################################
# Direct, Total and Total Indirect Effects
############################################################

# Direct effects
direct_effects <- hsr_pls$path_coef
direct_effects

# Total effects and bootstrapped total paths
# Note: in this SEMinR version, total effects are shown inside summary(boot_hsr)
summary(boot_hsr)

# f2 effect sizes
fSquared(hsr_pls, iv = "Service_Quality", dv = "Subjective_Well_Being")
fSquared(hsr_pls, iv = "Travel_Experience_Affect", dv = "Subjective_Well_Being")
fSquared(hsr_pls, iv = "Cognitive_Evaluation", dv = "Subjective_Well_Being")
fSquared(hsr_pls, iv = "Travel_Experience_Affect", dv = "Customer_Loyalty")
fSquared(hsr_pls, iv = "Cognitive_Evaluation", dv = "Customer_Loyalty")
fSquared(hsr_pls, iv = "Subjective_Well_Being", dv = "Customer_Loyalty")
fSquared(hsr_pls, iv = "Subjective_Well_Being", dv = "Willingness_to_Pay")

############################################################
# Specific Indirect Effects
############################################################

# Function to extract a path coefficient safely
get_path <- function(path_matrix, from, to) {
  if (from %in% rownames(path_matrix) && to %in% colnames(path_matrix)) {
    value <- path_matrix[from, to]
    if (is.na(value)) return(0)
    return(value)
  } else {
    return(0)
  }
}

# Define specific indirect effects manually according to the model
specific_indirect_effects <- data.frame(
  Effect = c(
    "SQ -> TEA -> SWB",
    "SQ -> CET -> SWB",
    "SQ -> TEA -> CET -> SWB",
    "SQ -> TEA -> CL",
    "SQ -> CET -> CL",
    "SQ -> TEA -> CET -> CL",
    "TEA -> CET -> SWB",
    "TEA -> SWB -> CL",
    "TEA -> CET -> SWB -> CL",
    "CET -> SWB -> CL",
    "SQ -> SWB -> CL",
    "SQ -> TEA -> SWB -> CL",
    "SQ -> CET -> SWB -> CL",
    "SQ -> TEA -> CET -> SWB -> CL",
    "SQ -> SWB -> WTP",
    "TEA -> SWB -> WTP",
    "CET -> SWB -> WTP",
    "SQ -> TEA -> SWB -> WTP",
    "SQ -> CET -> SWB -> WTP",
    "SQ -> TEA -> CET -> SWB -> WTP"
  ),
  
  Indirect_Effect = c(
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Subjective_Well_Being"),
    
    get_path(direct_effects, "Service_Quality", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Customer_Loyalty"),
    
    get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being"),
    
    get_path(direct_effects, "Travel_Experience_Affect", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Customer_Loyalty"),
    
    get_path(direct_effects, "Service_Quality", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay"),
    
    get_path(direct_effects, "Travel_Experience_Affect", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay"),
    
    get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay"),
    
    get_path(direct_effects, "Service_Quality", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay"),
    
    get_path(direct_effects, "Service_Quality", "Travel_Experience_Affect") *
      get_path(direct_effects, "Travel_Experience_Affect", "Cognitive_Evaluation") *
      get_path(direct_effects, "Cognitive_Evaluation", "Subjective_Well_Being") *
      get_path(direct_effects, "Subjective_Well_Being", "Willingness_to_Pay")
  )
)

specific_indirect_effects$Indirect_Effect <- round(
  specific_indirect_effects$Indirect_Effect,
  3
)

specific_indirect_effects


############################################################
# Outer Loadings
############################################################

outer_loadings <- hsr_pls$outer_loadings
outer_loadings


############################################################
# SEM Diagram
############################################################

graph_hsr <- dot_graph(hsr_pls)
grViz(graph_hsr)


############################################################
# Export Results
############################################################

write.csv(rhoC_AVE(hsr_pls), "reliability_ave_results.csv")
write.csv(htmt_table, "htmt_results.csv")

r2_values <- data.frame(summary(hsr_pls)$rSquared)
write.csv(r2_values, "r2_values.csv")

f2_values <- data.frame(
  SQ_SWB = fSquared(hsr_pls, iv="Service_Quality", dv="Subjective_Well_Being"),
  TEA_SWB = fSquared(hsr_pls, iv="Travel_Experience_Affect", dv="Subjective_Well_Being"),
  CET_SWB = fSquared(hsr_pls, iv="Cognitive_Evaluation", dv="Subjective_Well_Being"),
  TEA_CL = fSquared(hsr_pls, iv="Travel_Experience_Affect", dv="Customer_Loyalty"),
  CET_CL = fSquared(hsr_pls, iv="Cognitive_Evaluation", dv="Customer_Loyalty"),
  SWB_CL = fSquared(hsr_pls, iv="Subjective_Well_Being", dv="Customer_Loyalty"),
  SWB_WTP = fSquared(hsr_pls, iv="Subjective_Well_Being", dv="Willingness_to_Pay")
)


write.csv(hsr_pls$path_coef, "direct_effects.csv")
write.csv(hsr_pls$outer_loadings, "outer_loadings.csv")
