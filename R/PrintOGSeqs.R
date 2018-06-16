#' Print OG Sequences
#' 
#' Print all protein sequences and annotations in a given OG
#' @param OG name of OG
#' @param after_ortho output from FormatAfterOrtho
#' @param fasta_dir directory to fastas
#' @param out_dir complete path to output directory
#' @param outfile name of file that will be written to
#' @return A fasta file with all protein sequences and ids for a given OG
#' @examples 
#' 
#' \dontrun{
#' OG <- 'OG5_126968'
#' dir <- system.file('extdata', 'fasta_dir', package='MAGNAMWAR')
#' dir <- paste(dir,'/',sep='')
#' 
#' PrintOGSeqs(after_ortho_format, OG, dir)
#' }
#' 
#' @export

PrintOGSeqs <- function(after_ortho, OG, fasta_dir,
                        out_dir = NULL, outfile = "none") {

    OG_proteins <- after_ortho$proteins[, OG]
    OG_proteins <- OG_proteins[OG_proteins != ""]

    orig_directory <- getwd()

    if (getwd() != fasta_dir) {
        setwd(fasta_dir)
    }

    if (is.null(out_dir)) {
        out_dir <- fasta_dir
    }

    files <- dir(fasta_dir)
    files <- files[!files %in% "MCLformatted_all.fasta"]

    OG_proteins <- t(as.data.frame(strsplit(as.character(OG_proteins),
                                            split = "\\|")))
    row.names(OG_proteins) <- OG_proteins[, 2]

    myfiles <- lapply(files, function(x) seqinr::read.fasta(x, seqtype = "AA",
                                                            as.string = T))

    if (outfile == "none") {
        outfile <- paste(OG, "seqs.fasta", sep = "")
    }

    for (i in 1:length(OG_proteins[, 1])) {

        taxa_fn <- paste(OG_proteins[i, 1], ".fasta", sep = "")
        numfile <- match(taxa_fn, files)
        num_prot <- grep(OG_proteins[i, 2], seqinr::getName(myfiles[[numfile]]))

        if (getwd() != out_dir) {
            setwd(out_dir)
        }
        seqinr::write.fasta(seqinr::getSequence(myfiles[[numfile]][num_prot[1]]),
                            sub(">", "", seqinr::getAnnot(myfiles[[numfile]][num_prot[1]])),
            outfile, open = "a")
    
    }

    setwd(orig_directory)
}
