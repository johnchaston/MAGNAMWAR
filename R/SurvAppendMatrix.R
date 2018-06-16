#' Append Survival Test Outputs
#' 
#' Function used to append all .csv files that are outputted from AnalyzeOrthoMCL into one matrix.
#' @param work_dir the directory where the output files of AnalyzeOrthoMCL are located
#' @param out_name file name of outputted matrix
#' @param out_dir the directory where the outputted matrix is placed
#' @return A csv file containing a matrix with the following columns: OG, p-values, Bonferroni corrected p-values, mean phenotype of OG-containing taxa, mean pheotype of OG-lacking taxa, taxa included in OG, taxa not included in OG
#' @examples 
#' 
#' \dontrun{
#' file <- system.file('extdata', 'outputs', package='MAGNAMWAR')
#' directory <- paste(file, '/', sep = '')
#' SurvAppendMatrix(directory)
#' }
#'
#' @export


SurvAppendMatrix <- function(work_dir, out_name = "surv_matrix.csv",
                               out_dir = NULL) {

    if (is.null(out_dir)) {
        out_dir <- work_dir
    }

    orig_directory <- getwd()

    if (getwd() != work_dir) {
        setwd(work_dir)
    }

    filenames <- unlist(list.files(pattern = ".csv", full.names = F))

    out_starve <- matrix(c(rep(1, 7)), ncol = 7)
    for (i in filenames) {
        data <- read.csv(paste(i, sep = ""), header = F)
        if (grepl("Error", (data[1]))) {
            break
        } else {
            for (k in 1:ncol(data)) {
                final_data <- c(as.character(data[1, k]),
                                as.character(data[2, 1]),
                                as.character(data[3, 1]),
                                as.character(data[4, 1]),
                                as.character(data[5, 1]),
                                as.character(data[6, 1]),
                                as.character(data[7, 1]))
                out_starve <- rbind(out_starve, final_data)
            }
        }
    }
    out_starve <- out_starve[-1, ]

    row.names(out_starve) <- NULL

    colnames(out_starve) <- c("OG", "p-val", "corrected_p-val",
                              "mean_OGContain", "mean_OGLack",
                              "taxa_contain", "taxa_miss")

    write.csv(out_starve, out_name)

    wd <- getwd()
    cat(paste("Wrote:", out_name, "to", wd))

    setwd(orig_directory)

}
