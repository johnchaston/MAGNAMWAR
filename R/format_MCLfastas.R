#' Format all raw GenBank fastas to single OrthoMCL compatible fasta file
#' 
#' Creates the composite fasta file for use in running OrthoMCL and/or submitting to www.orthomcl.org
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have '.fasta' file descriptor
#' @param genbnk_id (Only necessary for the deprecated version of fasta headers) The index of the sequence ID in the GenBank pipe-separated annotation line (default: 4)
#' @return The path to the final OrthoMCL compatible fasta file
#' @examples
#' 
#' # Not run ~ directory structure depends on system
#' dir <- system.file('sample_data', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' formatted_file <- format_MCLfastas(dir)
#' @export

format_MCLfastas <- function(fa_dir, genbnk_id = 4) {
    
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
            
            # KIND OF ANNOYING cat(abs_path) cat('\n')
            
            id <- strsplit(files[i], split = "\\.")
            
            reps_fa <- seqinr::read.fasta(file = abs_path, as.string = T, forceDNAtolower = F, seqonly = F, strip.desc = T)
            
            for (j in 1:length(reps_fa)) {
                
                info <- seqinr::getAnnot(reps_fa[[j]])
                prot_id <- strsplit(info, split = " ")
                # other_info <- paste(prot_id[[1]][2:length(prot_id[[1]])], collapse = ' ') #how to get the annotation off new fasta
                # files
                mcl_info <- paste(c(id[[1]][1], prot_id[[1]][1]), collapse = "|")
                
                # DEPRECATED WITH NEW FASTA HEADERS
                
                # info <- strsplit(getAnnot(reps_fa[[j]]), split='\\|') mcl_info <- paste(c(id[[1]][1], info[[1]][genbnk_id]),
                # collapse='|')
                
                # Check for duplicate protein ids
                
                if (!(mcl_info %in% info_vec)) {
                  info_vec <- c(info_vec, mcl_info)
                  seq_vec <- c(seq_vec, reps_fa[[j]][1])
                  seqinr::write.fasta(reps_fa[[j]][1], mcl_info, outfile, open = "a")
                } else {
                  cat("Duplicate protein id found: ", mcl_info, "\n")
                }
                
                # DEPRECATED IF : Ability to save duplicate protein ids in this form: NP_12345_1 flawed because what if more than 2 of
                # same id?
                
                # else { new_id <- paste(info[[1]][genbnk_id],'_1', sep='') mcl_info <- paste(id[[1]][1], new_id, sep='|') info_vec <-
                # c(info_vec, mcl_info) seq_vec <- c(seq_vec,reps_fa[[j]][1]) }
                
            }
        }
    }
    
    cat("finished.\n")
}
