#' Print analyzed matrix
#' 
#' Writes a tab separated version of the analyzed OrthoMCL data with or without the joined representative sequences
#' @param mtrx Matrix derived from AnalyzeOrthoMCL
#' @param filename File name to save final output
#' @return The path to the written file
#' @examples
#' WriteMCL(mcl_mtrx, 'matrix.tsv')
#' #mcl_mtrx previously derived from AnalyzeOrthoMCL() or join_repset()
#' @export
WriteMCL <- function(mtrx, filename) {

    write.csv(mtrx, filename, quote = F, row.names = F)

    wd <- getwd()
    file <- paste(wd, filename, sep = "/")
    cat("wrote matrix to", file, "\n")
}
