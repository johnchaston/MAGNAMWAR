#' Format file from output of OrthoMCL algorithm before use in analyzeOrthoMCL
#' 
#' After running OrthoMCL and/or submitting to www.orthomcl.org, formats the output file to be used in analyzeOrthoMCL
#' @param file Path to the orthoMCL output file
#' @param format Defaults to 'ortho' for output from orthomcl.org. Other option is 'groups' for output from local run of OrthoMCL software.
#' @return value
#' @examples
#' file <- system.file('sample_data', 'orthologGroups.txt', package='MAGNAMWAR')
#' after_ortho_format <- format_afterOrtho(file)
#' 
#' file_grps <- system.file('sample_data', 'groups_example_r.txt', package='MAGNAMWAR')
#' after_ortho_format_grps <- format_afterOrtho(file_grps, format = 'groups')
#' 
#' @export

format_afterOrtho <- function(file, format = "ortho") {
    cat(file)
    
    cat("\n\nreading in OrthoMCL data...")
    if (format == "ortho") {
        
        input <- read.table(file, header = F, sep = "\t")
        if (is.null(input$V2)) 
            stop("Probably incorrect format designation: Specified \"ortho\" (online) and might have meant \"groups\" (local).")
        linput <- input[order(input$V2), ]
        un <- unique(input[c("V1", "V2")])
        
        ### make matrix
        taxa_protein_COG <- matrix(, ncol = 0, nrow = max(table(list(un$V2))))
        cat("   done\ncreating list of COG proteins...")
        
        for (i in names(table(list(un$V2)))) {
            
            j <- droplevels(subset(un, un$V2 == i))
            group <- unique(as.character(j$V1))
            length(group) <- dim(taxa_protein_COG)[1]
            taxa_protein_COG <- cbind(taxa_protein_COG, group)
            colnames(taxa_protein_COG)[dim(taxa_protein_COG)[2]] <- i
            
        }
        
        ### set NA values to empty string
        taxa_protein_COG[is.na(taxa_protein_COG)] <- ""
        
        ### cleaning steps
        delcols <- c()
        for (i in 1:dim(taxa_protein_COG)[2]) {
            rows <- sum(taxa_protein_COG[, i] != "")
            
            ### remove no group values
            if (colnames(taxa_protein_COG)[i] == "NO_GROUP") {
                delcols[length(delcols) + 1] <- i
            }
            
            ### remove singleton groups
            if (rows <= 1) {
                delcols[length(delcols) + 1] <- i
            }
            
            ### remove species-only groups
            code <- strsplit(taxa_protein_COG[1, i], split = "\\|")[[1]][1]
            addcol <- T
            for (j in 1:rows) {
                test <- strsplit(taxa_protein_COG[j, i], split = "\\|")[[1]][1]
                if (code != test) {
                  addcol <- F
                  break
                }
            }
            if (addcol & !(i %in% delcols)) {
                delcols[length(delcols) + 1] <- i
            }
        }
        taxa_protein_COG <- taxa_protein_COG[, -delcols]
        
        
        taxa <- unique(substring(as.character(un$V1), 1, 4))
        taxa <- sort(taxa)
        
        ### NO PROTEINS INCLUDED
        short <- apply(taxa_protein_COG, c(1, 2), FUN = function(x) substring(as.character(x), 1, 4))
        
        
    } else if (format == "groups") {
        
        taxa_protein_COG <- as.matrix(read.delim(file, header = F, sep = " ", row.names = 1))
        rownames(taxa_protein_COG) <- sub(":", "", rownames(taxa_protein_COG))
        taxa_protein_COG <- as.matrix(t(taxa_protein_COG))
        
        cat("   done\ncreating list of COG proteins...")
        
        ### find all taxa
        short <- apply(taxa_protein_COG, 1:2, function(x) substring(x, 1, 4))
        taxa <- unique(c(short))
        taxa <- sort(taxa)
        taxa <- taxa[-1]
        
    } else {
        cat("format declaration invalid")
    }
    
    cat("  done\n")
    
    cat("creating presence absence matrix...")
    
    ### Create Presence Absence Matrix
    pa_mtrx <- matrix(nrow = 0, ncol = length(taxa))
    colnames(pa_mtrx) <- taxa
    class(pa_mtrx) <- "character"
    
    pa_rows = c()
    
    for (j in 1:dim(taxa_protein_COG)[2]) {
        # for each column COG
        pa_string <- ""
        taxa_in_COG <- c(unique(short[, j]))
        for (i in colnames(pa_mtrx)) {
            # for each taxa
            if (i %in% taxa_in_COG) {
                pa_string <- paste(pa_string, "1")
            } else {
                pa_string <- paste(pa_string, "0")
            }
        }
        
        ## check for repeat PDGs
        if (!(pa_string %in% pa_rows)) {
            pa_rows = c(pa_rows, pa_string)
            ind <- which(pa_rows == pa_string)
            names(pa_rows)[ind] <- colnames(taxa_protein_COG)[j]
        } else {
            ind <- which(pa_rows == pa_string)
            names(pa_rows)[ind] <- paste(names(pa_rows)[ind], colnames(taxa_protein_COG)[j], sep = ",")
        }
    }
    
    pa_rows <- lapply(pa_rows, function(x) substring(as.character(x), 2, nchar(pa_rows[1])))
    pa_rows <- strsplit(unlist(pa_rows), split = " ")
    pa_rows <- as.vector(pa_rows)
    
    for (i in 1:length(pa_rows)) {
        pa_mtrx <- rbind(pa_mtrx, pa_rows[[i]])
    }
    rownames(pa_mtrx) <- names(pa_rows)
    cat("\tdone\n\nExported list with Presence/Absence Matrix & list of COG Proteins")
    return(list(pa_matrix = pa_mtrx, proteins = taxa_protein_COG))
    
}
