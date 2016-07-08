#' Append Survival Test Outputs
#' 
#' Function used to append all .csv files that are outputted from analyze_OrthoMCL into one matrix.
#' @param work_dir the directory where the output files of analyze_OrthoMCL are located
#' @param out_name file name of outputted matrix
#' @param out_dir the directory where the outputted matrix is placed
#' @return A csv file containing a matrix with the following columns: COG, p-values, Bonferroni corrected p-values, mean phenotype of COG-containing taxa, mean pheotype of COG-lacking taxa, taxa included in COG, taxa not included in COG
#' @references Some sort of reference
#' @examples 
#' file <- system.file('sample_data', 'outputs', package='MAGNAMWAR')
#' directory <- paste(file, '/', sep = '')
#' surv_append_matrix(directory)
#'
#'
#' @export


surv_append_matrix <- function(work_dir, out_name = "surv_matrix.csv", out_dir = NULL) {
    
    if (is.null(out_dir)) {
        out_dir = work_dir
    }
    
    if (getwd() != work_dir) {
        setwd(work_dir)
    }
    
    filenames <- unlist(list.files(pattern = ".csv", full.names = F))
    
    out_starve <- matrix(c(rep(1, 7)), ncol = 7)
    for (i in filenames) {
        data <- read.csv(paste(i, sep = ""), header = F)
        if (as.character(data[1]) == "Error in p_val1 * num_pdg : non-numeric argument to binary operator") {
            break
        } else {
            for (k in 1:ncol(data)) {
                final_data <- c(as.character(data[1, k]), as.character(data[2, 1]), as.character(data[3, 1]), as.character(data[4, 
                  1]), as.character(data[5, 1]), as.character(data[6, 1]), as.character(data[7, 1]))
                out_starve <- rbind(out_starve, final_data)
            }
        }
    }
    out_starve <- out_starve[-1, ]
    
    row.names(out_starve) = NULL
    
    colnames(out_starve) = c("COG", "p-val", "corrected_p-val", "mean_COGContain", "mean_COGLack", "taxa_contain", "taxa_miss")
    
    write.csv(out_starve, out_name)
    
    wd <- getwd()
    cat(paste("Wrote:", out_name, "to", wd))
    
    
}
