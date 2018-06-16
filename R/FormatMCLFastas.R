#' Format all raw GenBank fastas to single OrthoMCL compatible fasta file
#' 
#' Creates the composite fasta file for use in running OrthoMCL and/or submitting to www.orthomcl.org
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have '.fasta' file descriptor
#' @param genbnk_id (Only necessary for the deprecated version of fasta headers) The index of the sequence ID in the GenBank pipe-separated annotation line (default: 4)
#' @return Returns nothing, but prints the path to the final OrthoMCL compatible fasta file
#' @examples
#' \dontrun{
#' dir <- system.file('extdata', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' formatted_file <- FormatMCLFastas(dir)
#' }
#' @export

FormatMCLFastas <- function(fa_dir, genbnk_id = 4) {

    filename <- "MCLformatted_all.fasta"
    outfile <- paste(c(fa_dir, filename), collapse = "")

    files <- dir(fa_dir, pattern = ".fasta")
    if (filename %in% files) {
        file.remove(outfile)
    }

    info_vec <- vector()
    seq_vec <- vector()

    cat("Writing out to", outfile, "\n")

    for (i in 1:length(files)) {

        if (files[i] != filename & files[i] != "repseq.fasta") {

            abs_path <- paste(c(fa_dir, files[i]), collapse = "")
            id <- strsplit(files[i], split = "\\.")

            reps_fa <- seqinr::read.fasta(file = abs_path, as.string = T,
                                          forceDNAtolower = F, seqonly = F,
                                          strip.desc = T)

            for (j in 1:length(reps_fa)) {

                info <- seqinr::getAnnot(reps_fa[[j]])
                prot_id <- strsplit(info, split = " ")

                # files
                mcl_info <- paste(c(id[[1]][1], prot_id[[1]][1]),
                                  collapse = "|")

                # Check for duplicate protein ids

                if (!(mcl_info %in% info_vec)) {
                  info_vec <- c(info_vec, mcl_info)
                  seq_vec <- c(seq_vec, reps_fa[[j]][1])
                  seqinr::write.fasta(reps_fa[[j]][1], mcl_info,
                                      outfile, open = "a")
                } else {
                  cat("Duplicate protein id found: ", mcl_info, "\n")
                }

            }
        }
    }

    cat("finished.\n")
}
