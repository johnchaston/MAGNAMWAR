#' Pick Representative Sequences
#' 
#' Randomly picks a representative sequence from GenBank fasta files for every OrthoMCL COG
#' @param mcl_data output of format_afterOrtho() --list of 2 things-- 1: binary matrix indicating the presence / absence of genes in each COG and 2: vector of names of COGs
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have '.fasta' file descriptor
#' @param fastaformat options: new & old; defaults to old; takes care of NCBI updated fasta headers
#' @return The path to the representative sequence file in fasta format
#' @examples
#' dir <- system.file('sample_data', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' repseqs <- pick_repseq(after_ortho_format, dir, fastaformat = "old")
#' @export


pick_repseq <- function(mcl_data, fa_dir, fastaformat = "new") {
    
    ### must feed it a formatted mcl_data file from format_afterOrtho
    mcl_data <- mcl_data$proteins

    files <- dir(fa_dir, pattern = ".fasta")
    files <- files[!files %in% "MCLformatted_all.fasta"]
    
    if (getwd() != fa_dir) {
        setwd(fa_dir)
    }
    
    for (k in 1:length(files)) {
        l = strsplit(files[k], split = ".", fixed = T)
        var2 <- seqinr::read.fasta(file = files[k], seqtype = "AA", as.string = T, forceDNAtolower = F, seqonly = F, strip.desc = T)
        var3 <- matrix(rep("NA", length(var2) * 3), ncol = 3)
        colnames(var3) = c("V1", "V2", "V3")
        for (j in 1:length(var2)) {
            
            #how to get the annotation off new fasta
            if (fastaformat == "new") {
              prot_id <- strsplit(seqinr::getAnnot(var2[[j]]), split = " ")
              other_info <- paste(prot_id[[1]][2:length(prot_id[[1]])], collapse = ' ') 
              var3[j, ] <- c(prot_id, other_info, var2[[j]][1])
              
            } else if (fastaformat == "old") {
              info <- strsplit(strsplit(seqinr::getAnnot(var2[[j]]), "[", T)[[1]][1], "|", T)
              var3[j, ] <- c(info[[1]][4], info[[1]][5], var2[[j]][1])  
              
            }
        }
        var3 <- data.frame(var3)
        assign(l[[1]][1], var3)
    }
    setwd("..")
    rm(var2, var3)
    
    cat("\npicking and writing representative sequence for PDG:\n")
    
    ### run the first time
    i = 1
    total_seq <- sum(mcl_data[, i] != "")
    rndm_seq <- sample(1:total_seq, 1)
    rep_info <- strsplit(mcl_data[rndm_seq, i], split = "\\|")
    
    var2 <- get(rep_info[[1]][1])
    var4 <- data.frame(rep_info[[1]][2])
    names(var4) <- c("V1")
    var3 <- merge(var2, var4, by = "V1", all = F)
    var5 <- as.character(var3[1, 3])
    out_reps <- paste(colnames(mcl_data)[i], "\t", rep_info[[1]][1], "\t", rep_info[[1]][2], "\tNA", sep = "")
    try(out_reps <- paste(colnames(mcl_data)[i], "\t", rep_info[[1]][1], "\t", rep_info[[1]][2], "\t", var3[1, 2], 
        sep = ""), T)
    out_reps <- rbind(out_reps, var5)
    rm(var2, var3, var4, var5, total_seq, rndm_seq, rep_info)
    
    
    for (i in 2:dim(mcl_data)[2]) {
        
        total_seq <- sum(mcl_data[, i] != "")
        rndm_seq <- sample(1:total_seq, 1)
        rep_info <- strsplit(mcl_data[rndm_seq, i], split = "\\|")
        
        var2 <- get(rep_info[[1]][1])
        var4 <- data.frame(rep_info[[1]][2])
        names(var4) <- c("V1")
        var3 <- merge(var2, var4, by = "V1", all = F)  # Merging the taxa protein_id with information
        var5 <- as.character(var3[1, 3])
        
        annotations <- paste(colnames(mcl_data)[i], "\t", rep_info[[1]][1], "\t", rep_info[[1]][2], "\tNA", sep = "")
        try(annotations <- paste(colnames(mcl_data)[i], "\t", rep_info[[1]][1], "\t", rep_info[[1]][2], "\t", var3[1, 
            2], sep = ""), T)
        
        out_reps <- rbind(out_reps, annotations)
        out_reps <- rbind(out_reps, var5)
        
        if (((i - 1)%%(round(dim(mcl_data)[2]/10, digits = 0))) == 0) {
            cat(paste(round(((i - 1)/dim(mcl_data)[2] * 100), digits = 2), "%__", sep = ""))
        }
        rm(var2, var3, var4, var5, annotations, total_seq, rndm_seq, rep_info)
        
    }
    
    cat("\nFinished!\n")

    return(out_reps)
    
}
