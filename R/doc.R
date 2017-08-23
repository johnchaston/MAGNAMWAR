# Documenting Data Sets

#' Formatted output of OrthoMCL.
#'
#' A list created by inputting the output of OrthoMCL clusters
#' into the FormatAfterOrtho function.
#'
#' @format List of 2: (1) presence absence matrix, (2) protein ids:
#' \describe{
#'   \item{pa_matrix}{matrix showing taxa presence/absence in OG}
#'   \item{proteins}{matrix listing protein_id contained in each OG}
#' }
"after_ortho_format"


#' Formatted output of OrthoMCL.
#'
#' A list created by inputting the output of OrthoMCL clusters
#' into the FormatAfterOrtho function.
#'
#' @format List of 2: (1) presence absence matrix, (2) protein ids:
#' \describe{
#'   \item{pa_matrix}{matrix showing taxa presence/absence in OG}
#'   \item{proteins}{matrix listing protein_id contained in each OG}
#' }
"after_ortho_format_grps"

#' Final output of join_repset.
#'
#' A data frame containing the final results of statistical analysis with 
#' protein ids, annotations, and sequences added.
#'
#' @format A data frame with 17 rows and 11 variables:
#' \describe{
#'   \item{OG}{taxa cluster id, as defined by OrthoMCL}
#'   \item{pval1}{p-value, based on presence absence}
#'   \item{corrected_pval1}{Bonferroni p-value, corrected by number of tests}
#'   \item{mean_OGContain}{mean of all taxa phenotypes in that OG}
#'   \item{mean_OGLack}{mean of all taxa phenotypes not in that OG}
#'   \item{taxa_contain}{taxa in that cluster}
#'   \item{taxa_miss}{taxa not in that cluster}
#'   \item{rep_taxon}{randomly selected representative taxa from the cluster}
#'   \item{rep_id}{protein id, from randomly selected representative taxa}
#'   \item{rep_annot}{fasta annotation, from randomly selected representative taxa}
#'   \item{rep_seq}{AA sequence, from randomly selected representative taxa}
#' }
"joined_mtrx"

#' Final output of join_repset.
#'
#' A data frame containing the final results of statistical analysis with 
#' protein ids, annotations, and sequences added.
#'
#' @format A data frame with 10 rows and 11 variables:
#' \describe{
#'   \item{OG}{taxa cluster id, as defined by OrthoMCL}
#'   \item{pval1}{p-value, based on presence absence}
#'   \item{corrected_pval1}{Bonferroni p-value, corrected by number of tests}
#'   \item{mean_OGContain}{mean of all taxa phenotypes in that OG}
#'   \item{mean_OGLack}{mean of all taxa phenotypes not in that OG}
#'   \item{taxa_contain}{taxa in that cluster}
#'   \item{taxa_miss}{taxa not in that cluster}
#'   \item{rep_taxon}{randomly selected representative taxa from the cluster}
#'   \item{rep_id}{protein id, from randomly selected representative taxa}
#'   \item{rep_annot}{fasta annotation, from randomly selected representative taxa}
#'   \item{rep_seq}{AA sequence, from randomly selected representative taxa}
#' }
"joined_mtrx_grps"


#' Final output of AnalyzeOrthoMCL
#'
#' A matrix containing the final results of statistical analysis.
#'
#' @format A matrix with 17 rows and 7 variables:
#' \describe{
#'   \item{OG}{taxa cluster id, as defined by OrthoMCL}
#'   \item{pval1}{p-value, based on presence absence}
#'   \item{corrected_pval1}{Bonferroni p-value, corrected by number of tests}
#'   \item{mean_OGContain}{mean of all taxa phenotypes in that OG}
#'   \item{mean_OGLack}{mean of all taxa phenotypes not in that OG}
#'   \item{taxa_contain}{taxa in that cluster}
#'   \item{taxa_miss}{taxa not in that cluster}
#' }
"mcl_mtrx"

#' Final output of AnalyzeOrthoMCL
#'
#' A matrix containing the final results of statistical analysis.
#'
#' @format A matrix with 10 rows and 7 variables:
#' \describe{
#'   \item{OG}{taxa cluster id, as defined by OrthoMCL}
#'   \item{pval1}{p-value, based on presence absence}
#'   \item{corrected_pval1}{Bonferroni p-value, corrected by number of tests}
#'   \item{mean_OGContain}{mean of all taxa phenotypes in that OG}
#'   \item{mean_OGLack}{mean of all taxa phenotypes not in that OG}
#'   \item{taxa_contain}{taxa in that cluster}
#'   \item{taxa_miss}{taxa not in that cluster}
#' }
"mcl_mtrx_grps"


#' Triglyceride (TAG) content of fruit flies dataset.
#'
#' A subset of the TAG content of fruit flies, collected in the Chaston Lab,
#' to be used as a brief example for tests in AnalyzeOrthoMCL.
#'
#' @format A data frame with 586 rows and 4 variables:
#' \describe{
#'   \item{Treatment}{4-letter taxa designation of associated bacteria}
#'   \item{RespVar}{response variable, TAG content}
#'   \item{Vial}{random effect variable, vial number of flies}
#'   \item{Experiment}{random effect variable, experiment number of flies}
#' }
"pheno_data"

#' Starvation rate of fruit flies dataset.
#'
#' A subset of the Starvation rate of fruit flies, collected in the Chaston Lab,
#' to be used as a brief example for survival tests in AnalyzeOrthoMCL.
#'
#' @format A matrix with 543 rows and 7 variables:
#' \describe{
#'   \item{EXP}{random effect variable, experiment number of flies}
#'   \item{VIAL}{random effect variable, vial number of flies}
#'   \item{BACLO}{fixed effect variable, loss of bacteria in flies}
#'   \item{TRT}{4-letter taxa designation of associated bacteria}
#'   \item{t1}{time 1}
#'   \item{t2}{time 2}
#'   \item{event}{event}
#' }
"starv_pheno_data"
