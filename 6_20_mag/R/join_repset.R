#' Join Representative Sequences
#' 
#' Joins the OrthoMCL output matrix to representative sequences
#' @param reps_file output of pick_repseq()
#' @param mcl_mtrx OrthoMCL output matrix from analyze_OrthoMCL()
#' @return Returns the original OrthoMCL output matrix with additional columns: representative sequence taxon, representative sequence id, representative sequence annotation, representative sequence 
#' @examples
#' joined_mtrx_grps <- join_repset(repseqs_grps, mcl_mtrx_grps)
#' #mcl_mtrx1 previously derived from analyze_OrthoMCL()
#' #repseqs1 previously derived from pick_repseq()
#' @export
join_repset <- function(reps_file, mcl_mtrx) {
    
    # library('seqinr')
    
    fa_mtrx <- matrix(nrow = length(reps_file), ncol = 5)
    colnames(fa_mtrx) <- c("COG", "rep_taxon", "rep_id", "rep_annot", "rep_seq")
    
    for (i in 1:length(reps_file)) {
        info <- strsplit(getAnnot(reps_file[[i]]), split = "\t")
        fa_mtrx[i, ] <- c(info[[1]][1], info[[1]][2], info[[1]][3], info[[1]][4], reps_file[[i]][1])
    }
    
    mcl_reps <- merge(mcl_mtrx, fa_mtrx, by = "COG", all = F)
    
    return(mcl_reps)
}
