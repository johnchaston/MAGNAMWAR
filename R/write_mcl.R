#' Print analyzed matrix
#' 
#' Writes a tab separated version of the analyzed OrthoMCL data with or without the joined representative sequences
#' @param mtrx Matrix derived from analyze_OrthoMCL
#' @param filename File name to save final output
#' @return The path to the written file
#' @examples
#' write_mcl(mcl_mtrx, 'matrix.tsv')
#' #mcl_mtrx previously derived from analyze_OrthoMCL() or join_repset()
#' @export
write_mcl <- function(mtrx, filename) {

    write.csv(mtrx, filename, quote = F, row.names = F)

    wd <- getwd()
    file <- paste(wd, filename, sep = "/")
    cat("wrote matrix to", file, "\n")
}
