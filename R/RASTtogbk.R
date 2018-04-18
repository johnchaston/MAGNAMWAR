#' Write RAST files to Genbank formats OrthoMCL Analysis
#' 
#' Useful for reformating RAST files to GBK format
#' @param input_fasta path to input fasta file
#' @param input_reference path to a .csv file; it should be downloaded from RAST as excel format, saved as a .csv (saved as the tab-delimited version has compatibility problems)
#' @param out_name_path name and path of the file to write to
#' @examples 
#' 
#' \donttest{
#' lfrc_fasta <- system.file('extdata', 'RASTtoGBK//lfrc.fasta', package='MAGNAMWAR')
#' lfrc_reference <- system.file('extdata', 'RASTtoGBK//lfrc_lookup.csv', package='MAGNAMWAR')
#' lfrc_path <- system.file('extdata', 'RASTtoGBK//lfrc_out.fasta', package='MAGNAMWAR')
#'
#' RASTtoGBK(lfrc_fasta,lfrc_reference,lfrc_path)
#' }
#' 
#' @import seqinr
#' @export

RASTtoGBK <- function(input_fasta, input_reference, out_name_path) {
    # library(seqinr) in_fasta is a read.fasta product of the package (seqinr)
    in_fasta <- seqinr::read.fasta(paste(input_fasta, sep = ""),
                                   seqtype = "AA", as.string = T,
                                   strip.desc = T)

    ## reference file is a .csv file; it should be downloaded
    ## from RAST as excel format,
    ## saved as a .csv, and read in as a .csv
    in_rast <- read.csv(paste(input_reference, sep = ""), header = T)

    ## take only columns 2 and 8 from the input reference file
    in_rast <- data.frame(in_rast[, c(2, 8)])  # specify column name

    suppressWarnings(try(rm(out)))  ## clear the output variable
    ## do the steps below, the same as is in the loop but on
    ## the first line only, to eliminate an empty space at the top
    ## of the file when doing rbind the first time

    ## STEP1 get the annotation
    anot <- data.frame(seqinr::getAnnot(in_fasta[[1]]))

    colnames(anot) <- "ids"
    ## STEP2 merge the annotation to the reference table to extract
    ## only the line with the annotation from the reference table
    anot_rast <- merge(anot, in_rast, by.x = "ids", by.y = "feature_id",
                       all = F)

    ## STEP3A create the variable to store the unique ID,
    ## split off from 'fig': strsplit(...)  STEP3B create the variable
    ## of the annotation as.character(droplevels(anot_rast[1,2]))
    wout <- paste(">gi|NA|ref|",
                  strsplit(as.character(droplevels(anot_rast[1, 1])),
                           split = "|", fixed = T)[[1]][2],
                  "|", as.character(droplevels(anot_rast[1, 2])), sep = "")

    ## STEP4 write out the annotation
    out <- rbind(wout, seqinr::getSequence(in_fasta[[1]], as.string = T))


    ## now run those in a loop for every sequence in the database
    for (i in 2:length(in_fasta)) {
        anot <- data.frame(seqinr::getAnnot(in_fasta[[i]]))
        colnames(anot) <- "ids"

        anot_rast <- merge(anot, in_rast, by.x = "ids", by.y = "feature_id",
                           all = F)
        out <- rbind(out, paste(">gi|NA|ref|",
                                strsplit(as.character(droplevels(anot_rast[1, 1])),
                                         split = "|", fixed = T)[[1]][2],
            "|", as.character(droplevels(anot_rast[1, 2])), sep = ""))
        out <- rbind(out, seqinr::getSequence(in_fasta[[i]], as.string = T))
    }

    ## return this product
    fin_out <- unlist(out)
    write.table(fin_out, paste(out_name_path, sep = ""), quote = F,
                row.names = F, col.names = F)
    cat("Finished writing out", out_name_path)
}
