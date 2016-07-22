#' Main OrthoMCL Analysis
#' 
#' Main function for analyzing the statistical association of PDG (phylogenetic distribution group) presence with phenotype data
#' @param mcl_data output of format_afterOrtho --list of 2 things-- 1: binary matrix indicating the presence / absence of genes in each COG and 2: vector of names of COGs
#' @param pheno_data R object with column names of the following variables
#' @param model Linear Model with gene presence as fixed effect (lm),Linear Mixed Effect models with gene presence as fixed effect and additional variables specified as: one random effect (lmeR1); two independent random effects (lmeR2ind); two random effects with rndm2 nested in rndm1 (lmeR2nest); or two independent random effects with one additional fixed effect (lmeF2), Wilcox Test with gene presence as fixed effect (wx), Survival Test with support for multi core design
#' @param species_name Column name in pheno_data containing 4-letter species designations
#' @param resp Column name in pheno_data containing response variable
#' @param fix2 Column name in pheno_data containing second fixed effect
#' @param rndm1 Column name in pheno_data containing first random variable
#' @param rndm2 Column name in pheno_data containing second random variable
#' @param multi (can only be used with survival tests) Number of cores
#' @param time (can only be used with survival tests) Column name in pheno_data containing first time
#' @param event (can only be used with survival tests) Column name in pheno_data containing event
#' @param time2 (can only be used with survival tests) Column name in pheno_data containing econd time
#' @param startnum number of test to start on
#' @param stopnum number of test to stop on
#' @return A matrix with the following columns: COG, p-values, Bonferroni corrected p-values, mean phenotype of COG-containing taxa, mean pheotype of COG-lacking taxa, taxa included in COG, taxa not included in COG
#' @references Some sort of reference
#' @examples 
#' #Linear Model
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'lm',
#'  'Treatment', resp='RespVar')
#'
#'
#' #Linear Mixed Effect with one random effect
#' \dontrun{
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'lmeR1',
#' 'Treatment', resp='RespVar', rndm1='Experiment')
#' }
#'
#' #Linear Mixed Effect with two independent random effects
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'lmeR2ind',
#'  'Treatment', resp='RespVar', rndm1='Experiment', rndm2='Vial')
#'
#'
#' #Linear Mixed Effect with rndm2 nested in rndm1
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'lmeR2nest',
#'  'Treatment',  resp='RespVar', rndm1='Experiment', rndm2='Vial')
#'
#'
#' #Linear Mixed Effect with two independent random effects and one additional fixed effect
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'lmeF2',
#'  'Treatment', resp='RespVar', fix2='Treatment', rndm1='Experiment', rndm2='Vial')
#'
#'
#' #Wilcox Test
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, pheno_data, 'wx',
#'  'Treatment', resp='RespVar')
#'
#'# 5 minutes!
#' #Survival with two independent random effects and one additional fixed effect, run on multiple cores
#' \dontrun{
#' mcl_mtrx <- analyze_OrthoMCL(after_ortho_format, starv_pheno_data, 'TRT', model='survmulti',
#'  time='t2', event='event', rndm1='EXP', rndm2='VIAL', multi=1)
#' }
#'
#' #5 minutes!
#' #Survival with two independent random effects and one additional fixed effect,
#' #including drops on multi cores
#' \dontrun{
#' mtrx <- analyze_OrthoMCL(after_ortho_format, starv_pheno_data, 'TRT', model='survmulticensor',
#'  time='t1', time2='t2', event='event', rndm1='EXP', rndm2='VIAL', fix2='BACLO', multi=1)
#'  }
#' #to be appended with surv_append_matrix
#' @importFrom foreach %dopar%
#' @importFrom multcomp glht mcp
#' @export




analyze_OrthoMCL <- function(mcl_data, pheno_data, model, species_name, resp = NULL, fix2 = NULL, rndm1 = NULL, rndm2 = NULL, 
    multi = 1, time = NULL, event = NULL, time2 = NULL, startnum = 1, stopnum = "end") {    
    
    cat("Importing Data\n")
  
    pa_mtrx <- t(mcl_data$pa_matrix)
    colnames(pa_mtrx) <- NULL
    haplo_names <- row.names(mcl_data$pa_matrix)
    
    if (model == "lm") 
      mtrx <- analyze.f(pa_mtrx, haplo_names, pheno_data, species_name, resp)
    else if (model == "lmeR1") 
      mtrx <- analyze.fr(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1)
    else if (model == "lmeR2ind") 
        mtrx <- analyze.frr.plus(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1, rndm2)
    else if (model == "lmeR2nest") 
      mtrx <- analyze.frr.div(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1, rndm2)
    else if (model == "lmeF2") 
        mtrx <- analyze.ffrr(pa_mtrx, haplo_names, pheno_data, species_name, resp, fix2, rndm1, rndm2)  
    else if (model == "wx") 
        mtrx <- analyze.wilcox(pa_mtrx, haplo_names, pheno_data, species_name, resp)
    else if (model == "survmulti") 
        mtrx <- analyze.surv.multi(pa_mtrx, haplo_names, pheno_data, species_name, time, event, rndm1, 
            rndm2, multi, startnum, stopnum) 
    else if (model == "survmulticensor") 
        mtrx <- analyze.surv.censor.multi(pa_mtrx, haplo_names, pheno_data, species_name, time, time2, 
            event, rndm1, rndm2, fix2, multi, startnum, stopnum) 
    else cat("Error: Could not find a correct match for your model declaration\n")
    
    return(mtrx)
}

analyze.f <- function(pa_mtrx, haplo_names, tx, species_name, resp_var) {
  
  #library(multcomp)
      
  cat("Merging Files\n")
  haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
  
  count <- 0
  for (i in 1:length(haplo_names)) {
    count <- count + length(unlist(strsplit(haplo_names[i],',')))      
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow = count, ncol = 7)
  
  sub <- list()
  sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]

  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for (i in 1:(num_pdg)) {
    
    ### getting column names for experimental fixed effect
    name <- paste("V", i, sep = "")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
    
    ### fitting the linear model
    lm <- try(stats::lm(resp_var ~ fix, data = sub), T)
    if (class(lm) != "try-error") {
      l1 <- try(glht(lm, mcp(fix = "Tukey")), T)
      l2 <- try(summary(l1), T)
      pval <- try(l2$test$pvalues[1], T)
      
      ### calculating meta-data
      mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
      mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
      mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
      mean_contain <- mean_calc2[2, 2]
      mean_missing <- mean_calc2[1, 2]
      
      taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
      colnames(taxa) <- c("taxa_name", "presence")
      
      taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
      taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
      
      taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
      taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
      
      try(pval_corrected <- pval * num_pdg, T)
      
      try(if (pval_corrected > 1) {
        pval_corrected <- 1
      }, T)
      
      
      
      
      ### adding data to output matrix
      for (j in unlist(strsplit(haplo_names[i],','))) {
        
        if (class(pval) != "try-error") {
          output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, 
                                 taxa_contain, taxa_missing)
          out_cnt <- out_cnt + 1
        }
        
      }
    }
    ### printing progress
    if (((i - 1)%%100) == 0) 
      cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output)) != 7, ]
  colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
                              "taxa_miss")
  
  return(output_clean)
}

analyze.fr <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1) {
  
  #library(lme4)
  #library(multcomp)
      
  cat("Merging Files\n")
  haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
  
  count <- 0
  for (i in 1:length(haplo_names)) {
    count <- count + length(unlist(strsplit(haplo_names[i],',')))      
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow = count, ncol = 7)
  
  sub <- list()
  sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
  sub$rndm <- haplo_tx[, which(colnames(haplo_tx) == rndm1)]
  
  sub$rndm <- factor(sub$rndm)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for (i in 1:(num_pdg)) {
    
    ### getting column names for experimental fixed effect
    name <- paste("V", i, sep = "")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
    
    ### fitting the linear mixed model
    lmm <- try(lme4::lmer(resp_var ~ fix + (1 | rndm), data = sub), T)
    if (class(lmm) != "try-error") {
      l1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
      l2 <- try(summary(l1), T)
      pval <- try(l2$test$pvalues[1], T)
      
      ### calculating meta-data
      mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
      mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
      mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
      mean_contain <- mean_calc2[2, 2]
      mean_missing <- mean_calc2[1, 2]
      
      taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
      colnames(taxa) <- c("taxa_name", "presence")
      
      taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
      taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
      
      taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
      taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
      
      try(pval_corrected <- pval * num_pdg, T)
      
      try(if (pval_corrected > 1) {
        pval_corrected <- 1
      }, T)
      
      
      
      
      ### adding data to output matrix
      for (j in unlist(strsplit(haplo_names[i],','))) {
        
        if (class(pval) != "try-error") {
          output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, 
                                 taxa_contain, taxa_missing)
          out_cnt <- out_cnt + 1
        }
      }
    }
    ### printing progress
    if (((i - 1)%%100) == 0) 
      cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output)) != 7, ]
  colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
                              "taxa_miss")
  
  return(output_clean)
}

analyze.frr.plus <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1, rndm2) {
  
#   library(lme4)
#   library(multcomp)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
  
  count <- 0
  for (i in 1:length(haplo_names)) {
    count <- count + length(unlist(strsplit(haplo_names[i],',')))    
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow = count, ncol = 7)
  
  ### create a good matrix changed it to a data frame (2/10/16)
  sub <- cbind(haplo_tx[, which(colnames(haplo_tx) == resp_var)], haplo_tx[, which(colnames(haplo_tx) == rndm1)], 
               haplo_tx[, which(colnames(haplo_tx) == rndm2)])
  sub <- as.data.frame(sub)
  colnames(sub) <- c("resp_var", "rndm1", "rndm2")
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  
  for (i in 1:(num_pdg)) {
    
    ### getting column names for experimental fixed effect
    name <- paste("V", i, sep = "")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
    
    ### fitting the linear mixed model
    lmm <- try(lme4::lmer(resp_var ~ fix + (1 | rndm1) + (1 | rndm2), data = sub), T)
    if (class(lmm) != "try-error") {
      l1 <- try(glht(lmm, mcp(fix = "Tukey")),T)
      l2 <- try(summary(l1), T)
      pval <- try(l2$test$pvalues[1], T) 
      
      ### calculating meta-data
      mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
      mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
      mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
      mean_contain <- mean_calc2[2, 2]
      mean_missing <- mean_calc2[1, 2]
      
      taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
      colnames(taxa) <- c("taxa_name", "presence")
      
      
      taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
      taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
      
      taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
      taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
      
      try(pval_corrected <- pval * num_pdg, T)
      
      try(if (pval_corrected > 1) {
        pval_corrected <- 1
      },T)
      
      
      
      
      ### adding data to output matrix
      for (j in unlist(strsplit(haplo_names[i],','))) {
        if (class(pval) != "try-error") {
          output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, 
                                 taxa_contain, taxa_missing)
          out_cnt <- out_cnt + 1
        }
      }
    }
    ### printing progress
    if (((i - 1)%%100) == 0) 
      cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output)) != 7, ]
  colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
                              "taxa_miss")
  
  return(output_clean)
}

analyze.frr.div <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1, rndm2) {
  
#   library(lme4)
#   library(multcomp)
#   
  cat("Merging Files\n")
  haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
  
  count <- 0
  for (i in 1:length(haplo_names)) {
    count <- count + length(unlist(strsplit(haplo_names[i],',')))      
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow = count, ncol = 7)
  
  ### create a matrix
  sub <- cbind(haplo_tx[, which(colnames(haplo_tx) == resp_var)], haplo_tx[, which(colnames(haplo_tx) == rndm1)], 
               haplo_tx[, which(colnames(haplo_tx) == rndm2)])
  sub <- as.data.frame(sub)
  colnames(sub) <- c("resp_var", "rndm1", "rndm2")  # connected to lines 337-338
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for (i in 1:(num_pdg)) {
    
    ### getting column names for experimental fixed effect
    name <- paste("V", i, sep = "")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
    
    ### fitting the linear mixed model
    lmm <- try(lme4::lmer(resp_var ~ fix + (1 | rndm1/rndm2), data = sub), T)
    if (class(lmm) != "try-error") {
      l1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
      l2 <- try(summary(l1), T)
      pval <- try(l2$test$pvalues[1], T)
      
      try(pval_corrected <- pval * num_pdg, silent = T)
      
      try(if (pval_corrected > 1) {
        pval_corrected <- 1
      }, T)
      
      ### calculating meta-data
      mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
      mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
      mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
      mean_contain <- mean_calc2[2, 2]
      mean_missing <- mean_calc2[1, 2]
      
      taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
      colnames(taxa) <- c("taxa_name", "presence")
      
      taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
      taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
      
      taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
      taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
      
      
      
      
      ### adding data to output matrix
      for (j in unlist(strsplit(haplo_names[i],','))) {
        if (class(pval) != "try-error") {
          output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, 
                                 taxa_contain, taxa_missing)
          out_cnt <- out_cnt + 1
        }
      }
    }  
    ### printing progress
    if (((i - 1)%%100) == 0) 
      cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output)) != 7, ]
  colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
                              "taxa_miss")
  
  return(output_clean)
}

analyze.ffrr <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, fix_var2, rndm1, rndm2) {
    
#     library(lme4)
#     library(multcomp)
     
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
      count <- count + length(unlist(strsplit(haplo_names[i],',')))      
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 9)
    
    sub <- cbind(haplo_tx[, which(colnames(haplo_tx) == resp_var)], haplo_tx[, which(colnames(haplo_tx) == fix_var2)], 
        haplo_tx[, which(colnames(haplo_tx) == rndm1)], haplo_tx[, which(colnames(haplo_tx) == rndm2)])
    sub <- as.data.frame(sub)
    colnames(sub) <- c("resp_var", "fix2", "rndm1", "rndm2")
    
    sub$rndm1 <- factor(sub$rndm1)
    sub$rndm2 <- factor(sub$rndm2)
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear mixed model
        lmm <- try(lme4::lmer(resp_var ~ fix + fix2 + (1 | rndm1) + (1 | rndm2), data = sub), T)
        if (class(lmm) != "try-error") {
          l1.1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
          l1.2 <- try(summary(l1.1), T)
          pval1 <- try(l1.2$test$pvalues[1], T)
          
          l2.1 <- try(glht(lmm, mcp(fix2 = "Tukey")), T)
          l2.2 <- try(summary(l2.1), T)
          pval2 <- "variable not a factor"
          try(pval2 <- l2.2$test$pvalues[1], T)
          
          pval2_corrected <- "ibid"
          try(pval2_corrected <- pval2 * num_pdg, T)
          
          try(pval1_corrected <- pval1 * num_pdg, T)
          
          
          try(if (pval1_corrected > 1) {
              pval1_corrected <- 1
          },T)
          
          
          ### calculating meta-data
          mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
          mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
          mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
          mean_contain <- mean_calc2[2, 2]
          mean_missing <- mean_calc2[1, 2]
          
          taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
          colnames(taxa) <- c("taxa_name", "presence")
          
          taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
          taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
          
          taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
          taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
          
  
  
          
          ### adding data to output matrix
          for (j in unlist(strsplit(haplo_names[i],','))) {
            if (class(pval1) != "try-error") {
              output[out_cnt, ] <- c(j, pval1, pval1_corrected, pval2, pval2_corrected, mean_contain, mean_missing, 
                                     taxa_contain, taxa_missing)
              out_cnt <- out_cnt + 1
            }
          }
        }
        ### printing progress
        if (((i - 1)%%100) == 0) 
            cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 9, ]
    colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "pval2", "pval2_corrected", "mean_COGContain", "mean_COGLack", 
        "taxa_contain", "taxa_miss")
        
    return(output_clean)
}

analyze.wilcox <- function(pa_mtrx, haplo_names, tx, species_name, resp_var) {
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
      count <- count + length(unlist(strsplit(haplo_names[i],',')))      
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 7)
    
    sub <- list()
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear model
        wix <- try(stats::wilcox.test(resp_var ~ fix, data = sub), T)
        pval <- try(wix$p.value, T)
        
        ### calculating meta-data
        mean_calc <- stats::aggregate(resp_var ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
        mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
        mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
        mean_contain <- mean_calc2[2, 2]
        mean_missing <- mean_calc2[1, 2]
        
        taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
        colnames(taxa) <- c("taxa_name", "presence")
        
        taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
        taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
        
        taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
        taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
        
        try(pval_corrected <- pval * num_pdg, T)
        
        try(if (pval_corrected > 1) {
            pval_corrected <- 1
        },T)
        


        
        ### adding data to output matrix
        for (j in unlist(strsplit(haplo_names[i],','))) {
          
          if (class(pval) != "try-error") {
            output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, 
                                   taxa_contain, taxa_missing)
            out_cnt <- out_cnt + 1
          }
        
        }
        
        ### printing progress
        if (((i - 1)%%100) == 0) 
            cat(paste(round(((i - 1)/num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
        "taxa_miss")
        
    return(output_clean)
}

### from analyze_surve_0.0.2.R
analyze.surv.multi <- function(pa_mtrx, haplo_names, tx, species_name, time, event, rndm1, rndm2, multi, 
    startnum, stopnum) {
#     
#     library(survival)
#     library(parallel)
#     library(foreach)
#     library(doParallel)
#     
    ### merge the binary matrix with phenotype data; tx = phenotypes; pa_mtrx=binary matrix produced by
    ### 'parse_orthologGroups'
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    
    ### set the number of tests to peform
    num_pdg <- length(haplo_names) - 1
    
    cat("Creating output directory in", getwd(),'\n\n')
    suppressWarnings(dir.create("outputs"))
    
    ### determine the number of COGs per PDG - needed for the subsequent loop
    count <- 0
    for (i in 1:length(haplo_names)) {
      count <- count + length(unlist(strsplit(haplo_names[i],',')))      
    }
    
    cat("Running Analysis. There is no 'percent complete' update. If you want to estimate how much time it will take, run for a short time with the 'surv' model, estimate time it will take in that model on 1 core. Time with the multiple cores ~ (time on 1 core / # cores)*2 \n")
    
    ### make a new phenotype matrix with only the specified data and with assigned names columns (makes the model run better)
    sub <- cbind(haplo_tx[, which(colnames(haplo_tx) == time)], haplo_tx[, which(colnames(haplo_tx) == event)], haplo_tx[, 
        which(colnames(haplo_tx) == species_name)], haplo_tx[, which(colnames(haplo_tx) == rndm1)], haplo_tx[, which(colnames(haplo_tx) == 
        rndm2)])  # values in matrix \t###
    sub <- data.frame(sub)  #make into a dataframe, compatible with subsequent analyses
    colnames(sub) <- c("time", "event", "trt", "rndm1", "rndm2")  # define column names
    sub$S <- survival::Surv(sub$time, sub$event)  # create the response variable, a survival model
    sub$rndm1 <- as.numeric(as.character(sub$rndm1))  # set the random effects as factors
    sub$rndm2 <- as.numeric(as.character(sub$rndm2))  # set the random effects as factors
    sub$rndm1b <- factor(sub$rndm1, labels = c(1:length(table(list(sub$rndm1)))))
    sub$rndm2b <- factor(sub$rndm2, labels = c(1:length(table(list(sub$rndm2)))))
    sub$rndm1 <- as.numeric(as.character(sub$rndm1b))  # set the random effects as factors
    sub$rndm2 <- as.numeric(as.character(sub$rndm2b))  # set the random effects as factors
    
    ### specify the start and stop positions.  stop number
    if (stopnum == "end") {
        stopnum <- as.numeric(as.character(length(haplo_names)))
    } else {
        stopnum <- as.numeric(as.character(stopnum))
    }
    
    # start number
    startnum <- as.numeric(as.character(startnum))
    
    ### set up the parallelization
    cl <- parallel::makeCluster(multi)  # of nodes
    doParallel::registerDoParallel(cl)
    
    outmulti <- c(rep(1, 7), unlist(foreach::foreach(i = iterators::icount(stopnum - startnum + 1)) %dopar% {
        
        
#         library(coxme)
#         library(multcomp)
#         library(survival)
#         
        i = i + startnum - 1
        
        ### start over with fresh values each iteration
        try(rm(output, name, lmm, l1.1, pval1, mean_contain, mean_missing, taxa_contain, taxa_missing), silent = T)
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$fix <- factor(sub$fix)
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the survival model
        lmm <- try(coxme::coxme(S ~ fix + (1 | trt) + (1 | rndm1/rndm2), data = sub), T)
        if (class(lmm) != "try-error") {
          l1.1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
          pval1 <- try(summary(l1.1)$test$pvalues[1], T)
          
          ### calculating meta-data means
          mean_calc <- stats::aggregate(time ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
          mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
          mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
          mean_contain <- mean_calc2[2, 2]
          mean_missing <- mean_calc2[1, 2]
          
          # taxa
          taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub, mean)
          colnames(taxa) <- c("taxa_name", "presence")
          taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
          taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
          taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
          taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
          
          output <- c(rep(NA, 7))
          
          try(pval1_corrected <- pval1 * num_pdg, T)
          
          try(if (pval1_corrected > 1) {
              pval1_corrected <- 1
          },T)
          
          
          ### adding data to output matrix
          for (j in unlist(strsplit(haplo_names[i],','))) {
              rm(output)
              output <- try(c(j, pval1, pval1_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing), T)
              #write(j, file = "data.csv") MAYBE FROM JOHNNY TESTS? /
              write(output, file = paste("outputs/", j, ".csv", sep = ""))
          }
          as.vector(output)  
      }
    }))
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    
    output_clean <- matrix(data = outmulti, ncol = 7, byrow = T)
    output_clean <- output_clean[-1,]
    
    colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
        "taxa_miss")
        
    return(output_clean)
}

analyze.surv.censor.multi <- function(pa_mtrx, haplo_names, tx, species_name, time, time2, event, rndm1, rndm2, 
    fix2, multi, startnum, stopnum) {
#     
#     library(survival)
#     library(parallel)
#     library(foreach)
#     library(doParallel)
#     library(coxme)
#     library(multcomp)
#     
    ### merge the binary matrix with phenotype data: tx = phenotypes; pa_mtrx=binary matrix produced by
    ### 'parse_orthologGroups'
    
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    
    ### number of phylogenetic distribution groups (PDG)
    num_pdg <- length(haplo_names) - 1
    
    ### directory for output files
    cat("Creating output directory in", getwd(),'\n\n')
    suppressWarnings(dir.create("outputs"))
    
    ### determine the number of COGs per PDG - needed for the subsequent loop
    count <- 0
    for (i in 1:length(haplo_names)) {
      count <- count + length(unlist(strsplit(haplo_names[i],',')))      
    }
    
    cat("Running Analysis. There is no 'percent complete' update. If you want to estimate how much time it will take, run for a short time with the 'surv' model, estimate time it will take in that model on 1 core. Time with the multiple cores ~ (time on 1 core / # cores)*2 \n")
    
    ### make a new phenotype matrix with only the specified data and with assigned names columns (makes the model run better)
    sub3 <- cbind(haplo_tx[, which(colnames(haplo_tx) == time)], haplo_tx[, which(colnames(haplo_tx) == time2)],
                  haplo_tx[,which(colnames(haplo_tx) == event)], haplo_tx[, which(colnames(haplo_tx) == species_name)],
                  haplo_tx[, which(colnames(haplo_tx) == rndm1)], haplo_tx[, which(colnames(haplo_tx) == rndm2)],
                  haplo_tx[, which(colnames(haplo_tx) == fix2)])  # values in matrix \t###
    sub3 <- data.frame(sub3)  #make into a dataframe, compatible with subsequent analyses
    colnames(sub3) <- c("time1", "time2", "event", "trt", "rndm1", "rndm2", "fix2")  # define column names
    sub3$S <- survival::Surv(as.numeric(as.character(sub3$time1)), as.numeric(as.character(sub3$time2)), as.numeric(as.character(sub3$event)))  # create the response variable, a survival model
    sub3$rndm1 <- as.numeric(as.character(sub3$rndm1))  # set the random effects as numeric
    sub3$v <- factor(sub3$rndm2, labels = c(1:length(table(list(sub3$rndm2)))))  # 
    sub3$v2 <- as.numeric(sub3$v)  # set the random effects as numeric
    sub3$fix2 <- factor(sub3$fix2)  # set the fixed effect as a factor 
    
    
    ### specify the start and stop positions.  stop number
    if (stopnum == "end") {
        stopnum <- as.numeric(as.character(length(haplo_names)))
    } else {
        stopnum <- as.numeric(as.character(stopnum))
    }
    
    ### start number
    startnum <- as.numeric(as.character(startnum))
    
    ### set up the parallelization
    cl <- parallel::makeCluster(multi)  # of nodes
    doParallel::registerDoParallel(cl)
    
    outmulti <- c(rep(1, 7), unlist(foreach::foreach(i = iterators::icount(stopnum - startnum + 1)) %dopar% {
        
#         library(coxme)
#         library(multcomp)
#         library(survival)
#         
        i = i + startnum - 1
        
        ### start over with fresh values each iteration
        try(rm(output, name, lmm, l1.1, pval1, mean_contain, mean_missing, taxa_contain, taxa_missing), silent = T)
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub3$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub3$fix <- factor(sub3$fix)
        sub3$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the survival model
        lmm <- try(coxme::coxme(S ~ fix + fix2 + (1 | trt) + (1 | rndm1/rndm2), data = sub3), T)
        if (class(lmm) != "try-error") {
          l1.1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
          pval1 <- try(summary(l1.1)$test$pvalues[1], T)  # should this be 2?
          
          
          #summary(l1.1)$test$pvalues[1], file = "data.csv")  MAYBE FROM JOHNNY TESTS? /
          #write(summary(l1.1)$test$pvalues[2], file = "data.csv")
          
          
          
          ### calculating meta-data means
          mean_calc <- stats::aggregate(time2 ~ fix, sub3, mean)  #matrix with mean_contain and mean_missing
          mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
          mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing = F), ]
          mean_contain <- mean_calc2[2, 2]
          mean_missing <- mean_calc2[1, 2]
          # taxa
          taxa <- stats::aggregate(as.numeric(as.character(fix)) ~ species, sub3, mean)
          colnames(taxa) <- c("taxa_name", "presence")
          taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence == 1))
          taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse = "|")
          taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence == 0))
          taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse = "|")
          
          output <- c(rep(NA, 7))
          
          try(pval1_corrected <- pval1 * num_pdg, T)
          
          try(if (pval1_corrected > 1) {
              pval1_corrected <- 1
          },T)
          
          for (j in unlist(strsplit(haplo_names[i],','))) {
              rm(output)
              output <- try(c(j, pval1, pval1_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing), T)
              #write(i, file = "data.csv") #MAYBE FROM JOHNNY TESTS? /
              write(output, file = paste("outputs/", j,".csv", sep = ""))
          }
          
          as.vector(output)
      }
    }))
    cat("Finished!\n")
    
    cat("Cleaning Final Output\n\n")
    
    output_clean <- matrix(data = outmulti, ncol = 7, byrow = T)
    output_clean <- output_clean[-1,]
    
    colnames(output_clean) <- c("COG", "pval1", "corrected_pval1", "mean_COGContain", "mean_COGLack", "taxa_contain", 
                                "taxa_miss")    
    
    cat(paste("Wrote small output files to:", getwd(), "\n"))
    cat("Outputs should be concatenated with surv_append_matrix")
    
    return(output_clean)
    
}

