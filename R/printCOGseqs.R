#' Print COG Sequences
#' 
#' Print all protein sequences and annotations in a given COG
#' @param COG name of COG
#' @param after_ortho output from format_afterOrtho
#' @param fasta_dir directory to fastas
#' @param out_dir complete path to output directory
#' @param outfile name of file that will be written to
#' @return A fasta file with all protein sequnces and ids for a given COG
#' @examples 
#' 
#' # Not run ~ directory structure depends on system
#' \dontrun{
#' COG <- 'OG5_126968'
#' dir <- system.file('sample_data', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' 
#' printCOGseqs(after_ortho_format, COG, dir)
#' }
#' 
#' @export

printCOGseqs <- function(after_ortho, COG, fasta_dir, out_dir = NULL, outfile = "none") {
    
    COG_proteins <- after_ortho$proteins[, COG]
    COG_proteins <- COG_proteins[COG_proteins != ""]
    
    if (getwd() != fasta_dir) {
        setwd(fasta_dir)
    }
    
    if (is.null(out_dir)) {
        out_dir = fasta_dir
    }
    
    files <- dir(fasta_dir)
    files <- files[!files %in% "MCLformatted_all.fasta"]
    
    COG_proteins = t(as.data.frame(strsplit(as.character(COG_proteins), split = "\\|")))
    row.names(COG_proteins) <- COG_proteins[, 2]
    
    myfiles <- lapply(files, function(x) seqinr::read.fasta(x, seqtype = "AA", as.string = T))
    
    if (outfile == "none") {
        outfile = paste(COG, "seqs.fasta", sep = "")
    }
    
    for (i in 1:length(COG_proteins[, 1])) {
        
        taxa_fn <- paste(COG_proteins[i, 1], ".fasta", sep = "")
        numfile <- match(taxa_fn, files)
        num_prot <- grep(COG_proteins[i, 2], seqinr::getName(myfiles[[numfile]]))
        
        if (getwd() != out_dir) {
            setwd(out_dir)
        }
        seqinr::write.fasta(seqinr::getSequence(myfiles[[numfile]][num_prot[1]]),
                    sub(">", "", seqinr::getAnnot(myfiles[[numfile]][num_prot[1]])), 
                    outfile, open = "a")
        
    }
}

