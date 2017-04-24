#' Join Representative Sequences
#' 
#' Joins the OrthoMCL output matrix to representative sequences
#' @param mcl_data output of format_afterOrtho() --list of 2 things-- 1: binary matrix indicating the presence / absence of genes in each OG and 2: vector of names of OGs
#' @param fa_dir Path to the directory where all raw GenBank files are stored. Note, all file names must be changed to a 4-letter code representing each species and have '.fasta' file descriptor
#' @param fastaformat options: new & old; defaults to old; takes care of NCBI updated fasta headers
#' @param mcl_mtrx OrthoMCL output matrix from analyze_OrthoMCL()
#' @return Returns the original OrthoMCL output matrix with additional columns: representative sequence taxon, representative sequence id, representative sequence annotation, representative sequence 
#' @examples
#' 
#' dir <- system.file('extdata', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' joined_mtrx_grps <- join_repseq(after_ortho_format_grps, dir, mcl_mtrx_grps, fastaformat = 'old')
#' 
#' @export

join_repseq <- function(mcl_data, fa_dir, mcl_mtrx, fastaformat = "new") {

    ### must feed it a formatted mcl_data file from format_afterOrtho
    mcl_data <- mcl_data$proteins

    files <- dir(fa_dir, pattern = ".fasta")
    files <- files[!files %in% "MCLformatted_all.fasta"]

    orig_directory <- getwd()

    if (getwd() != fa_dir) {
        setwd(fa_dir)
    }

    for (k in 1:length(files)) {
        l <- strsplit(files[k], split = ".", fixed = T)
        var2 <- seqinr::read.fasta(file = files[k], seqtype = "AA",
                                   as.string = T, forceDNAtolower = F,
                                   seqonly = F, strip.desc = T)
        var3 <- matrix(rep("NA", length(var2) * 3), ncol = 3)
        colnames(var3) <- c("V1", "V2", "V3")
        for (j in 1:length(var2)) {

            # how to get the annotation off new fasta
            if (fastaformat == "new") {
                prot_id <- strsplit(seqinr::getAnnot(var2[[j]]), split = " ")
                other_info <- paste(prot_id[[1]][2:length(prot_id[[1]])],
                                    collapse = " ")

                # ERROR CHECKING FOR WRONG FASTA
                if (grepl("\\|", prot_id[[1]][1])) {
                  stop("Problem with reading fasta files, must use old version
                       of fasta: add parameter ' fastaformat=\"old\" '")
                }

                var3[j, ] <- c(prot_id, other_info, var2[[j]][1])

            } else if (fastaformat == "old") {
                info <- strsplit(strsplit(seqinr::getAnnot(var2[[j]]),
                                          "[", T)[[1]][1], "|", T)
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
    i <- 1
    total_seq <- sum(mcl_data[, i] != "")
    rndm_seq <- sample(1:total_seq, 1)
    rep_info <- strsplit(mcl_data[rndm_seq, i], split = "\\|")

    var2 <- get(rep_info[[1]][1])
    var4 <- data.frame(rep_info[[1]][2])
    names(var4) <- c("V1")
    var3 <- merge(var2, var4, by = "V1", all = F)
    var5 <- as.character(var3[1, 3])
    out_reps <- paste(colnames(mcl_data)[i], "\t",
                      rep_info[[1]][1], "\t", rep_info[[1]][2],
                      "\tNA", sep = "")
    try(out_reps <- paste(colnames(mcl_data)[i], "\t",
                          rep_info[[1]][1], "\t", rep_info[[1]][2],
                          "\t", var3[1, 2], sep = ""), T)
    out_reps <- rbind(out_reps, var5)
    rm(var2, var3, var4, var5, total_seq, rndm_seq, rep_info)

    for (i in 2:dim(mcl_data)[2]) {

        total_seq <- sum(mcl_data[, i] != "")

        if (total_seq != 0) {

            rndm_seq <- sample(1:total_seq, 1)
            rep_info <- strsplit(mcl_data[rndm_seq, i], split = "\\|")

            var2 <- get(rep_info[[1]][1])
            var4 <- data.frame(rep_info[[1]][2])
            names(var4) <- c("V1")
            # Merging the taxa protein_id with information
            var3 <- merge(var2, var4, by = "V1", all = F)
            var5 <- as.character(var3[1, 3])

            annotations <- paste(colnames(mcl_data)[i], "\t",
                                 rep_info[[1]][1], "\t", rep_info[[1]][2],
                                 "\tNA", sep = "")
            try(annotations <- paste(colnames(mcl_data)[i], "\t",
                                     rep_info[[1]][1], "\t", rep_info[[1]][2],
                                     "\t", var3[1, 2], sep = ""), T)

            out_reps <- rbind(out_reps, annotations)
            out_reps <- rbind(out_reps, var5)

            if (((i - 1) %% (round(dim(mcl_data)[2] / 10, digits = 0))) == 0) {
                cat(paste(round(((i - 1) / dim(mcl_data)[2] * 100),
                                digits = 2), "%__", sep = ""))
            }
            rm(var2, var3, var4, var5, annotations, total_seq,
               rndm_seq, rep_info)
        }
    }


    fa_mtrx <- matrix(nrow = length(out_reps) / 2, ncol = 5)
    colnames(fa_mtrx) <- c("OG", "rep_taxon", "rep_id", "rep_annot", "rep_seq")
    fa_out_reps <- lapply(out_reps, function(x) sub(",", "", x))

    count <- 1
    for (i in seq(1, length(out_reps), by = 2)) {
        info <- strsplit(out_reps[[i]], split = "\t")
        fa_mtrx[count, ] <- c(info[[1]][1], info[[1]][2], info[[1]][3],
                              info[[1]][4], out_reps[[i + 1]])
        count <- count + 1
    }

    mcl_reps <- merge(mcl_mtrx, fa_mtrx, by = "OG", all = F)

    setwd(orig_directory)

    return(mcl_reps)
}
