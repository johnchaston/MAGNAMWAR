#' Main OrthoMCL Analysis
#' 
#' Main function for analyzing the statistical association of OG (orthologous group) presence with phenotype data
#' @param mcl_data output of FormatAfterOrtho; a list of matrices; (1) a presence/absence matrix of taxa per OG, (2) a list of the specific protein ids within each OG
#' @param pheno_data a data frame of phenotypic data with specific column names used to specify response variable as well as other fixed and random effects
#' @param model linear model with gene presence as fixed effect (lm), linear mixed mffect models with gene presence as fixed effect and additional variables specified as: one random effect (lmeR1); two independent random effects (lmeR2ind); two random effects with rndm2 nested in rndm1 (lmeR2nest); or two independent random effects with one additional fixed effect (lmeF2), Wilcox Test with gene presence as fixed effect (wx), Survival Tests with support for multi core design: with two random effects (survmulti), and with two times as well as an additional fixed variable (survmulticensor)
#' @param species_name Column name in pheno_data containing 4-letter species designations
#' @param resp Column name in pheno_data containing response variable
#' @param fix2 Column name in pheno_data containing second fixed effect
#' @param rndm1 Column name in pheno_data containing first random variable
#' @param rndm2 Column name in pheno_data containing second random variable
#' @param multi (can only be used with survival tests) Number of cores
#' @param time (can only be used with survival tests) Column name in pheno_data containing first time
#' @param event (can only be used with survival tests) Column name in pheno_data containing event
#' @param time2 (can only be used with survival tests) Column name in pheno_data containing second time
#' @param startnum number of test to start on
#' @param stopnum number of test to stop on
#' @param output_dir (if using survival tests) directory where small output files will be placed before using SurvAppendMatrix. Must specify a directory if choosing to output small files, else only written as a matrix
#' @param sig_digits amount of digits to display for p-values and means of data; default to NULL (no rounding)
#' @param princ_coord the number of principle coordinates to be included in model as fixed effects (1, 2, or 3), if a decimal is specified, as many principal coordinates as are needed to account for that percentage of the variance will be included in the analysis
#' @return A matrix with the following columns: OG, p-values, Bonferroni corrected p-values, mean phenotype of OG-containing taxa, mean pheotype of OG-lacking taxa, taxa included in OG, taxa not included in OG
#' @examples 
#' #Linear Model
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'lm',
#'  'Treatment', resp='RespVar')
#'
#'
#' # the rest of the examples are not run for time's sake
#' #Linear Mixed Effect with one random effect
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'lmeR1',
#' 'Treatment', resp='RespVar', rndm1='Experiment')
#' }
#'
#' 
#' #Linear Mixed Effect with two independent random effects
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'lmeR2ind',
#'  'Treatment', resp='RespVar', rndm1='Experiment', rndm2='Vial')
#' }
#' 
#' 
#' #Linear Mixed Effect with rndm2 nested in rndm1
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'lmeR2nest',
#'  'Treatment',  resp='RespVar', rndm1='Experiment', rndm2='Vial')
#' 
#'
#' 
#' #Linear Mixed Effect with two independent random effects and one additional fixed effect
#' \donttest{
#' mcl_mtrx3 <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'lmeF2',
#'  'Treatment', resp='RespVar', fix2='Treatment', rndm1='Experiment', rndm2='Vial', princ_coord = 4)
#' }
#' 
#' 
#' #Wilcoxon Test
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, pheno_data, 'wx',
#'  'Treatment', resp='RespVar')
#' }
#'
#' # ~ 5 minutes
#' #Survival with two independent random effects, run on multiple cores
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, starv_pheno_data, 'TRT', model='survmulti',
#'  time='t2', event='event', rndm1='EXP', rndm2='VIAL', multi=1)
#' }
#'
#' # ~ 5 minutes
#' #Survival with two independent random effects and one additional fixed effect,
#' #including drops on multi cores
#' \donttest{
#' mcl_mtrx <- AnalyzeOrthoMCL(after_ortho_format, starv_pheno_data, 'TRT', model='survmulticensor',
#'  time='t1', time2='t2', event='event', rndm1='EXP', rndm2='VIAL', fix2='BACLO', multi=1)
#'  }
#' #to be appended with SurvAppendMatrix
#' @importFrom foreach %dopar%
#' @importFrom multcomp glht mcp
#' @importFrom stats as.formula prcomp var
#' @export

AnalyzeOrthoMCL <- function(mcl_data, pheno_data, model, species_name, resp = NULL, fix2 = NULL, rndm1 = NULL, rndm2 = NULL, 
    multi = 1, time = NULL, event = NULL, time2 = NULL, startnum = 1, stopnum = "end", output_dir = NULL, sig_digits = NULL, princ_coord = 0) {
    
    cat("Importing Data\n")
    
    pa_mtrx <- t(mcl_data$pa_matrix)
    colnames(pa_mtrx) <- NULL
    haplo_names <- row.names(mcl_data$pa_matrix)
    test <- NULL
    pr_coor_mtx <- NULL
    
    if (princ_coord != 0) {
      pa <- t(mcl_data$pa_matrix)
      
      # turn to numeric
      x <- mapply(pa, FUN=as.numeric)
    
      rows <- dimnames(pa)[[1]]
      cols <- dimnames(pa)[[2]]
      
      m <- matrix(data=x, ncol=length(cols), nrow=length(rows))
      
      rownames(m) <- rows
      colnames(m) <- cols
      
      nm <- m[ , apply(m, 2, var) != 0]
      
      prm <- prcomp(nm, scale=T)
      pr_coor_mtx = prm$x
      
      if (princ_coord < 1) {
        x <- summary(prm)
        importance <- x$importance[2,]
        count <- 1
        total_var <- 0
        test <- ""
        while (total_var < princ_coord) {
          test = paste(test, names(importance)[count], sep=" + ") 
          total_var = total_var + importance[count]
          count = count + 1
        }
      }
    }
    
    
      
    
    # error checking and model selection
    
    if (model == "lm") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp)
        }
        mtrx <- analyze.f(pa_mtrx, haplo_names, pheno_data, species_name, resp, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "lmeR1") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data)) || !(rndm1 %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp, "\n\tRandom1 Variable Column Name: ", rndm1)
        }
        mtrx <- analyze.fr(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "lmeR2ind") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data)) || !(rndm1 %in% colnames(pheno_data)) || 
            !(rndm2 %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp, "\n\tRandom1 Variable Column Name: ", rndm1, "\n\tRandom2 Variable Column Name: ", rndm2)
        }
        mtrx <- analyze.frr.plus(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "lmeR2nest") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data)) || !(rndm1 %in% colnames(pheno_data)) || !(rndm2 %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp, "\n\tRandom1 Variable Column Name: ", rndm1, "\n\tRandom2 Variable Column Name: ", rndm2)
        }
        mtrx <- analyze.frr.div(pa_mtrx, haplo_names, pheno_data, species_name, resp, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "lmeF2") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data)) || !(rndm1 %in% colnames(pheno_data)) || !(rndm2 %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp, "\n\tRandom1 Variable Column Name: ", rndm1, "\n\tRandom2 Variable Column Name: ", rndm2, "\n\tFix2 Variable Column Name: ", 
                fix2)
        }
        mtrx <- analyze.ffrr(pa_mtrx, haplo_names, pheno_data, species_name, resp, fix2, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "wx") {
        if (!(species_name %in% colnames(pheno_data)) || !(resp %in% colnames(pheno_data))) {
            stop("Invalid column names specified for phenotype data file\n\tSpecies Column Name: ", species_name, "\n\tResponse Variable Column Name: ", 
                resp)
        }
        mtrx <- analyze.wilcox(pa_mtrx, haplo_names, pheno_data, species_name, resp, sig_digits)
    } else if (model == "survmulti") {
        mtrx <- analyze.surv.multi(pa_mtrx, haplo_names, pheno_data, species_name, time, event, rndm1, rndm2, multi, 
            startnum, stopnum, output_dir, sig_digits, princ_coord, pr_coor_mtx, test)
    } else if (model == "survmulticensor") {
        mtrx <- analyze.surv.censor.multi(pa_mtrx, haplo_names, pheno_data, species_name, time, time2, event, rndm1, 
            rndm2, fix2, multi, startnum, stopnum, output_dir, sig_digits, princ_coord, pr_coor_mtx, test)
    } else stop("Error: Could not find a correct match for your model declaration\n")
    
    
    
    return(mtrx)
}

analyze.f <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, sig_digits, princ_coord, pr_coor_mtx, test) {
    
    # library(multcomp)
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 7)
    
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    
    
    
    if (princ_coord == 0) {
      my_model = resp_var ~ fix
    } else if (princ_coord == 1) {
      my_model = resp_var ~ fix + PC1
    } else if (princ_coord == 2) {
      my_model = resp_var ~ fix + PC1 + PC2
    } else if (princ_coord == 3) {
      my_model = resp_var ~ fix + PC1 + PC2 + PC3
    } else if (princ_coord < 1) {
      my_model = as.formula(paste("resp_var ~ fix", test, sep = ''))
    }
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear model
        lm <- try(stats::lm(my_model, data = sub), T)
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
            for (j in unlist(strsplit(haplo_names[i], ","))) {
                
                if (class(pval) != "try-error") {
                  output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing)
                  out_cnt <- out_cnt + 1
                }
                
            }
        }
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    
    if (!is.null(sig_digits)) {
        options(digits = sig_digits)
        ### significant digits rounding
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    return(output_clean)
}

analyze.fr <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1, sig_digits, princ_coord, pr_coor_mtx, test) {
    
    # library(lme4) library(multcomp)
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 7)
    
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    sub$rndm1 <- haplo_tx[, which(colnames(haplo_tx) == rndm1)]
    
    sub$rndm1 <- factor(sub$rndm1)
    
    if (princ_coord == 0) {
      my_model = resp_var ~ fix + (1 | rndm1)
    } else if (princ_coord == 1) {
      my_model = resp_var ~ fix + (1 | rndm1) + PC1
    } else if (princ_coord == 2) {
      my_model = resp_var ~ fix + (1 | rndm1) + PC1 + PC2
    } else if (princ_coord == 3) {
      my_model = resp_var ~ fix + (1 | rndm1) + PC1 + PC2 + PC3
    } else if (princ_coord < 1) {
      my_model = as.formula(paste("resp_var ~ fix", test, " + (1 | rndm1)", sep = ''))
    }
    
    
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear mixed model
        lmm <- try(lme4::lmer(my_model, data = sub), T)
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
            for (j in unlist(strsplit(haplo_names[i], ","))) {
                
                if (class(pval) != "try-error") {
                  output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing)
                  out_cnt <- out_cnt + 1
                }
            }
        }
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        ### significant digits rounding
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    return(output_clean)
}

analyze.frr.plus <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test) {
    
    # library(lme4) library(multcomp)
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 7)
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    sub$rndm1 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm1)])
    sub$rndm2 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm2)])
  
    if (princ_coord == 0) {
      my_model = resp_var ~ fix + (1 | rndm1) + (1 | rndm2)
    } else if (princ_coord == 1) {
      my_model = resp_var ~ fix + (1 | rndm1) + (1 | rndm2) + PC1
    } else if (princ_coord == 2) {
      my_model = resp_var ~ fix + (1 | rndm1) + (1 | rndm2) + PC1 + PC2
    } else if (princ_coord == 3) {
      my_model = resp_var ~ fix + (1 | rndm1) + (1 | rndm2) + PC1 + PC2 + PC3
    } else if (princ_coord < 1) {
      my_model = as.formula(paste("resp_var ~ fix", test, " + (1 | rndm1) + (1 | rndm2)", sep = ''))
    }
    
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        
        ### fitting the linear mixed model
        lmm <- try(lme4::lmer(my_model, data = sub), T)
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
            for (j in unlist(strsplit(haplo_names[i], ","))) {
                if (class(pval) != "try-error") {
                  output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing)
                  out_cnt <- out_cnt + 1
                }
            }
        }
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        ### significant digits rounding
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    
    return(output_clean)
}

analyze.frr.div <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test) {
    
    # library(lme4) library(multcomp)
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 7)
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    sub$rndm1 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm1)])
    sub$rndm2 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm2)])
    
    if (princ_coord == 0) {
      my_model = resp_var ~ fix + (1 | rndm1/rndm2)
    } else if (princ_coord == 1) {
      my_model = resp_var ~ fix + (1 | rndm1/rndm2) + PC1
    } else if (princ_coord == 2) {
      my_model = resp_var ~ fix + (1 | rndm1/rndm2) + PC1 + PC2
    } else if (princ_coord == 3) {
      my_model = resp_var ~ fix + (1 | rndm1/rndm2) + PC1 + PC2 + PC3
    } else if (princ_coord < 1) {
      my_model = as.formula(paste("resp_var ~ fix", test, " + (1 | rndm1/rndm2)", sep = ''))
    }
    
    
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear mixed model
        lmm <- try(lme4::lmer(my_model, data = sub), T)
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
            for (j in unlist(strsplit(haplo_names[i], ","))) {
                if (class(pval) != "try-error") {
                  output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing)
                  out_cnt <- out_cnt + 1
                }
            }
        }
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        ### significant digits rounding
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    
    return(output_clean)
}

analyze.ffrr <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, fix_var2, rndm1, rndm2, sig_digits, princ_coord, pr_coor_mtx, test) {
    
    # library(lme4) library(multcomp)
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    cat("Running Analysis\n")
    output <- matrix(, nrow = count, ncol = 9)
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$resp_var <- haplo_tx[, which(colnames(haplo_tx) == resp_var)]
    sub$fix2 <- factor(haplo_tx[, which(colnames(haplo_tx) == fix_var2)])
    sub$rndm1 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm1)])
    sub$rndm2 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm2)])
    
    if (princ_coord == 0) {
      my_model = resp_var ~ fix + fix2 + (1 | rndm1) + (1 | rndm2)
    } else if (princ_coord == 1) {
      my_model = resp_var ~ fix + fix2 + (1 | rndm1) + (1 | rndm2) + PC1
    } else if (princ_coord == 2) {
      my_model = resp_var ~ fix + fix2 + (1 | rndm1) + (1 | rndm2) + PC1 + PC2
    } else if (princ_coord == 3) {
      my_model = resp_var ~ fix + fix2 + (1 | rndm1) + (1 | rndm2) + PC1 + PC2 + PC3
    } else if (princ_coord < 1) {
      my_model = as.formula(paste("resp_var ~ fix + fix2", test, " + (1 | rndm1) + (1 | rndm2)", sep = ''))
    }
    
    
    cat("\tPercent Complete:\n\t")
    out_cnt <- 1
    for (i in 1:(num_pdg)) {
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the linear mixed model
        lmm <- try(lme4::lmer(my_model, data = sub), T)
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
            for (j in unlist(strsplit(haplo_names[i], ","))) {
                if (class(pval1) != "try-error") {
                  output[out_cnt, ] <- c(j, pval1, pval1_corrected, pval2, pval2_corrected, mean_contain, mean_missing, 
                    taxa_contain, taxa_missing)
                  out_cnt <- out_cnt + 1
                }
            }
        }
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 9, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "pval2", "pval2_corrected", "mean_OGContain", "mean_OGLack", 
        "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
        output_clean[, 6] <- as.character(format(as.numeric(output_clean[, 6]), sig_digits), scientific = T)
        output_clean[, 7] <- as.character(format(as.numeric(output_clean[, 7]), sig_digits), scientific = T)
    }
    
    return(output_clean)
}

analyze.wilcox <- function(pa_mtrx, haplo_names, tx, species_name, resp_var, sig_digits) {
    
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
    
    num_pdg <- dim(pa_mtrx)[2]  #number of phylogenetic distribution groups
    
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
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
        }, T)
        
        
        
        
        ### adding data to output matrix
        for (j in unlist(strsplit(haplo_names[i], ","))) {
            
            if (class(pval) != "try-error") {
                output[out_cnt, ] <- c(j, pval, pval_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing)
                out_cnt <- out_cnt + 1
            }
            
        }
        
        ### printing progress
        if (((i - 1) %% 100) == 0) 
            cat(paste(round(((i - 1) / num_pdg * 100), digits = 2), "%__", sep = ""))
    }
    cat("Finished!\n")
    
    cat("Cleaning Final Output")
    output_clean <- output[rowSums(is.na(output)) != 7, ]
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        options(digits = sig_digits)
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    return(output_clean)
}

#mtrx <- analyze.surv.multi(pa_mtrx, haplo_names, pheno_data, species_name, time, event, rndm1, rndm2, multi, startnum, stopnum, output_dir, sig_digits, princ_coord, pr_coor_mtx, test)

#mcl_mtrx1 <- AnalyzeOrthoMCL(after_ortho_format, starv_pheno_data, 'TRT', model='survmulti', time='t2', event='event', rndm1='EXP', rndm2='VIAL', multi=1, princ_coord = 1)

### from analyze_surve_0.0.2.R
analyze.surv.multi <- function(pa_mtrx, haplo_names, tx, species_name, time, event, rndm1, rndm2, multi, startnum, stopnum, 
    output_dir, sig_digits, princ_coord, pr_coor_mtx, test) {
    # library(survival) library(parallel) library(foreach) library(doParallel) merge the binary matrix with phenotype
    # data; tx = phenotypes; pa_mtrx=binary matrix produced by 'FormatAfterOrtho'
    cat("Merging Files\n")
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    
    ### set the number of tests to peform
    num_pdg <- length(haplo_names) - 1
    
    if (!is.null(output_dir)) {
      ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
    }
    
    ### determine the number of OGs per PDG - needed for the subsequent loop
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    # cat('Running Analysis. There is no 'percent complete' update. If you want to estimate how much time it will take,
    # run for a short time with the 'surv' model, estimate time it will take in that model on 1 core. Time with the
    # multiple cores ~ (time on 1 core / # cores)*2 \n')
    cat("Running Analysis. There is no 'percent complete' update.\n")
    
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$time <- haplo_tx[, which(colnames(haplo_tx) == time)]
    sub$event <- haplo_tx[, which(colnames(haplo_tx) == event)]
    sub$S <- survival::Surv(sub$time, sub$event)  # create the response variable, a survival model
    sub$trt <- factor(haplo_tx[, which(colnames(haplo_tx) == species_name)])
    sub$rndm2 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm2)])
    sub$rndm1 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm1)])

    sub$rndm1 <- factor(sub$rndm1, labels = c(1:length(table(list(sub$rndm1)))))
    sub$rndm2 <- factor(sub$rndm2, labels = c(1:length(table(list(sub$rndm2)))))
    
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
    
    smallfile <- F
    if (!is.null(output_dir)) {
      smallfile <- T
    }
    
    outmulti <- c(rep(1, 7), unlist(foreach::foreach(i = iterators::icount(stopnum - startnum + 1)) %dopar% {
        
        
        # library(coxme) library(multcomp) library(survival)
        i <- i + startnum - 1
        
        ### start over with fresh values each iteration
        try(rm(output, name, lmm, l1.1, pval1, mean_contain, mean_missing, taxa_contain, taxa_missing), silent = T)
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$fix <- factor(sub$fix)
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the survival model
        if (princ_coord == 0) {
          my_model = S ~ fix + (1 | trt) + (1 | rndm1 / rndm2)
        } else if (princ_coord == 1) {
          my_model = S ~ fix + (1 | trt) + (1 | rndm1 / rndm2) + PC1
        } else if (princ_coord == 2) {
          my_model = S ~ fix + (1 | trt) + (1 | rndm1 / rndm2) + PC1 + PC2
        } else if (princ_coord == 3) {
          my_model = S ~ fix + (1 | trt) + (1 | rndm1 / rndm2) + PC1 + PC2 + PC3
        } else if (princ_coord < 1) {
          my_model = as.formula(paste("S ~ fix", test, " + (1 | trt) + (1 | rndm1 / rndm2)", sep = ''))
        }
        
        lmm <- try(coxme::coxme(my_model, data = sub), T)
        cat(lmm$means)
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
            }, T)
            
            
            ### adding data to output matrix
            for (j in unlist(strsplit(haplo_names[i], ","))) {
              
                if (class(pval1) != "try-error") {
                
                  output <- try(c(j, pval1, pval1_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing), T)
                  if (smallfile) {
                    write(output, file = paste(output_dir, j, ".csv", sep = ""))
                  }
                }
            }
            as.vector(output)
        }
    }))
    cat("Finished!\n")
    if (smallfile) 
        cat("Output small files to", output_dir, "\n")
    
    cat("Cleaning Matrix Final Output")
    output_clean <- matrix(data = outmulti, ncol = 7, byrow = T)
    output_clean <- output_clean[-1, ]
    
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    ### significant digits rounding
    if (!is.null(sig_digits)) {
        options(digits = sig_digits)
        ### significant digits rounding
        output_clean[, 2] <- as.character(format(as.numeric(output_clean[, 2]), sig_digits), scientific = T)
        output_clean[, 3] <- as.character(format(as.numeric(output_clean[, 3]), sig_digits), scientific = T)
        output_clean[, 4] <- as.character(format(as.numeric(output_clean[, 4]), sig_digits), scientific = T)
        output_clean[, 5] <- as.character(format(as.numeric(output_clean[, 5]), sig_digits), scientific = T)
    }
    return(output_clean)
}



analyze.surv.censor.multi <- function(pa_mtrx, haplo_names, tx, species_name, time, time2, event, rndm1, rndm2, fix2, 
    multi, startnum, stopnum, output_dir, sig_digits, princ_coord, pr_coor_mtx, test) {
    # library(survival) library(parallel) library(foreach) library(doParallel) library(coxme) library(multcomp) merge the
    # binary matrix with phenotype data: tx = phenotypes; pa_mtrx=binary matrix produced by 'parse_orthologGroups'
    
    haplo_tx <- merge(tx, pa_mtrx, by.y = "row.names", by.x = species_name, all = F)
    
    ### number of phylogenetic distribution groups (PDG)
    num_pdg <- length(haplo_names) - 1
    
    if (!is.null(output_dir)) {
      ifelse(!dir.exists(output_dir), dir.create(output_dir), FALSE)
    }
    
    ### determine the number of OGs per PDG - needed for the subsequent loop
    count <- 0
    for (i in 1:length(haplo_names)) {
        count <- count + length(unlist(strsplit(haplo_names[i], ",")))
    }
    
    # cat('Running Analysis. There is no 'percent complete' update. If you want to estimate how much time it will take,
    # run for a short time with the 'surv' model, estimate time it will take in that model on 1 core. Time with the
    # multiple cores ~ (time on 1 core / # cores)*2 \n')
    cat("Running Analysis. There is no 'percent complete' update.\n")
    
    sub <- list()
    
    if (princ_coord != 0) {
      haplo_tx <- merge(haplo_tx, pr_coor_mtx, by.x = species_name, by.y = "row.names", all = F)
      if (princ_coord < 1) {
        vec <- strsplit(test, split = " \\+ ")
        vec <- vec[[1]][-1]
        for (pc in vec){
          sub <- c(sub, list(haplo_tx[, which(colnames(haplo_tx) == pc)]))
        }
        names(sub) <- vec
      } else {
        sub$PC1 <- haplo_tx[, which(colnames(haplo_tx) == "PC1")]
        sub$PC2 <- haplo_tx[, which(colnames(haplo_tx) == "PC2")]
        sub$PC3 <- haplo_tx[, which(colnames(haplo_tx) == "PC3")]
      }
    }
    
    sub$time1 <- haplo_tx[, which(colnames(haplo_tx) == time)]
    sub$time2 <- haplo_tx[, which(colnames(haplo_tx) == time2)]
    sub$event <- haplo_tx[, which(colnames(haplo_tx) == event)]
    sub$trt <- factor(haplo_tx[, which(colnames(haplo_tx) == species_name)])
    sub$rndm2 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm2)])
    sub$rndm1 <- factor(haplo_tx[, which(colnames(haplo_tx) == rndm1)])
    sub$fix2 <- factor(haplo_tx[, which(colnames(haplo_tx) == fix2)])
    
    sub$rndm1 <- factor(sub$rndm1, labels = c(1:length(table(list(sub$rndm1)))))
    sub$rndm2 <- factor(sub$rndm2, labels = c(1:length(table(list(sub$rndm2)))))
    sub$fix2 <- factor(sub$fix2)  # set the fixed effect as a factor 
    
    sub$S <- survival::Surv(as.numeric(as.character(sub$time1)), as.numeric(as.character(sub$time2)), as.numeric(as.character(sub$event)))  # create the response variable, a survival model

    
    #sub$rndm1 <- as.numeric(as.character(sub$rndm1))  # set the random effects as numeric
    #sub$rndm2 <- as.numeric(as.character(sub$rndm2))  # set the random effects as numeric
    
    #sub$v <- factor(sub$rndm2, labels = c(1:length(table(list(sub$rndm2)))))  # 
    #sub$v2 <- as.numeric(sub$v)  # set the random effects as numeric
    
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
    
    smallfile <- F
    if (!is.null(output_dir)) {
      smallfile <- T
    }
    
    
    
    
    
    
    
    
    
    
    
    
    #!
    
    
    ### make a new phenotype matrix with only the specified data and with assigned names columns (makes the model run
    ### better)
    #sub3 <- cbind(haplo_tx[, which(colnames(haplo_tx) == time)], haplo_tx[, which(colnames(haplo_tx) == time2)], haplo_tx[, 
    #    which(colnames(haplo_tx) == event)], haplo_tx[, which(colnames(haplo_tx) == species_name)], haplo_tx[, which(colnames(haplo_tx) == 
    #    rndm1)], haplo_tx[, which(colnames(haplo_tx) == rndm2)], haplo_tx[, which(colnames(haplo_tx) == fix2)])  # values in matrix \t###
    #sub3 <- data.frame(sub3)  #make into a dataframe, compatible with subsequent analyses
    #colnames(sub3) <- c("time1", "time2", "event", "trt", "rndm1", "rndm2", "fix2")  # define column names
    #sub3$S <- survival::Surv(as.numeric(as.character(sub3$time1)), as.numeric(as.character(sub3$time2)), as.numeric(as.character(sub3$event)))  # create the response variable, a survival model
    #sub3$rndm1 <- as.numeric(as.character(sub3$rndm1))  # set the random effects as numeric
    #sub3$v <- factor(sub3$rndm2, labels = c(1:length(table(list(sub3$rndm2)))))  # 
    #sub3$v2 <- as.numeric(sub3$v)  # set the random effects as numeric
    #sub3$fix2 <- factor(sub3$fix2)  # set the fixed effect as a factor 
    
    
    ### specify the start and stop positions.  stop number
    #if (stopnum == "end") {
    #    stopnum <- as.numeric(as.character(length(haplo_names)))
    #} else {
    #    stopnum <- as.numeric(as.character(stopnum))
    #}
    
    ### start number
   # startnum <- as.numeric(as.character(startnum))
    
    ### set up the parallelization
    #cl <- parallel::makeCluster(multi)  # of nodes
    #doParallel::registerDoParallel(cl)
    
    outmulti <- c(rep(1, 7), unlist(foreach::foreach(i = iterators::icount(stopnum - startnum + 1)) %dopar% {
        
        # library(coxme) library(multcomp) library(survival)
        i <- i + startnum - 1
        
        ### start over with fresh values each iteration
        try(rm(output, name, lmm, l1.1, pval1, mean_contain, mean_missing, taxa_contain, taxa_missing), silent = T)
        
        ### getting column names for experimental fixed effect
        name <- paste("V", i, sep = "")
        sub$fix <- haplo_tx[, which(colnames(haplo_tx) == name)]
        sub$fix <- factor(sub$fix)
        sub$species <- haplo_tx[, which(colnames(haplo_tx) == species_name)]
        
        ### fitting the survival model
        if (princ_coord == 0) {
          my_model = S ~ fix + fix2 + (1 | trt) + (1 | rndm1 / rndm2)
        } else if (princ_coord == 1) {
          my_model = S ~ fix + fix2 + (1 | trt) + (1 | rndm1 / rndm2) + PC1
        } else if (princ_coord == 2) {
          my_model = S ~ fix + fix2 + (1 | trt) + (1 | rndm1 / rndm2) + PC1 + PC2
        } else if (princ_coord == 3) {
          my_model = S ~ fix + fix2 + (1 | trt) + (1 | rndm1 / rndm2) + PC1 + PC2 + PC3
        } else if (princ_coord < 1) {
          my_model = as.formula(paste("S ~ fix", test, " + fix2 + (1 | trt) + (1 | rndm1 / rndm2)", sep = ''))
        }
        # S ~ fix + fix2 + (1 | trt) + (1 | rndm1 / rndm2)
        lmm <- try(coxme::coxme(my_model, data = sub), T)
        if (class(lmm) != "try-error") {
            l1.1 <- try(glht(lmm, mcp(fix = "Tukey")), T)
            pval1 <- try(summary(l1.1)$test$pvalues[1], T)  # should this be 2?
            
            
            ### calculating meta-data means
            mean_calc <- stats::aggregate(time2 ~ fix, sub, mean)  #matrix with mean_contain and mean_missing
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
            }, T)
            
            for (j in unlist(strsplit(haplo_names[i], ","))) {

              if (class(pval1) != "try-error") {
                output <- try(c(j, pval1, pval1_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing), T)
                # write(i, file = 'data.csv') #MAYBE FROM JOHNNY TESTS? /
                if (smallfile) {
                  write(output, file = paste(output_dir, j, ".csv", sep = ""))
                }
              }
            }
            
            as.vector(output)
        }
    }))
    cat("Finished!\n")
    if (smallfile) 
        cat("Output small files to", output_dir, "\n")
        cat("Outputs should be concatenated with SurvAppendMatrix")
    
    
    cat("Cleaning Final Output\n\n")
    
    output_clean <- matrix(data = outmulti, ncol = 7, byrow = T)
    output_clean <- output_clean[-1, ]
    
    colnames(output_clean) <- c("OG", "pval1", "corrected_pval1", "mean_OGContain", "mean_OGLack", "taxa_contain", "taxa_miss")
    
    return(output_clean)
    
}
