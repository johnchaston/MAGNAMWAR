#' Plot of a PDG and Data with Standard Error Bars
#' 
#' Bar plot of PDG vs phenotype data with presence of taxa in PDG indicated by color
#' @param data R object of phenotype data
#' @param mcl_data mcl matrix (analyze_OrthoMCL output)
#' @param species_colname name of column in data file with taxa designations
#' @param data_colname name of column in data file with data observations
#' @param GRP optional parameter, a string with the name of chosen group (COG) to be colored
#' @param xlab string to label barplot's x axis
#' @param ylab string to label barplot's y axis
#' @param tree optional parameter (defaults to NULL) Path to tree file, orders the taxa by phylogenetic distribution, else it defaults to alphabetical
#' @param order vector with order of taxa names for across the x axis (defaults to alpha ordering)
#' @return a barplot with taxa vs phenotypic data complete with standard error bars
#' @references Some sort of reference
#' @examples 
#' #dev.off()
#' pdgplot(pheno_data, mcl_mtrx, 'OG5_126778', 'Treatment', 'RespVar', ylimit=12)
#' @export
#' 

pdgplot <- function(data, mcl_matrix, GRP = "NONE", species_colname, data_colname, xlab = "Taxa", ylab = "Data", ylimit = NULL, 
    tree = NULL, order = NULL) {
    
    grep <- mcl_matrix[grep(GRP, mcl_matrix[, 1]), ]
    l <- unlist(strsplit(as.character(unlist(grep[6])), split = "\\|"))
    if (length(l) == 0) {
        warning("Invalid COG inputted... Printing without specific coloration\n")
    }
    
    ### Calculating st dev and means
    data <- data[, c(species_colname, data_colname)]
    colnames(data) <- c("one", "two")
    
    ccc <- aggregate(data$one, data["one"], paste, collapse = " ")
    ddd <- aggregate(data$two, data["two"], paste, collapse = " ")
    
    means.sem <- ddply(data, c("one"), summarise, mean = mean(two), sem = sd(two)/sqrt(length(two)))
    means.sem <- transform(means.sem, lower = mean - sem, upper = mean + sem)
    mean <- means.sem[, 2]
    stDevs <- means.sem[, 3]
    
    names(mean) <- means.sem$one
    names(stDevs) <- means.sem$one
    
    x <- mean
    
    ### ordering by tree, order vector, or alphabetical
    
    if (!is.null(tree)) {
        # order by tree
        tree <- read.tree(tree)
        
        x <- x[tree$tip.label]
        stDevs <- stDevs[tree$tip.label]
    } else {
        
        if (!is.null(order)) {
            # order by order vector
            
            x <- x[order]
            stDevs <- stDevs[order]
            
            x <- x[!is.na(x)]
            stDevs <- stDevs[!is.na(stDevs)]
            
        } else {
            # order by alpha
            
            sorted <- sort(names(x))
            x <- x[sorted]
            stDevs <- stDevs[sorted]
            
        }
        
    }
    
    if (is.null(ylimit)) {
        ylimit = max(x) + max(x)/4
    }
    
    ### coloring for graph
    
    col <- names(x) %in% l
    col <- sub("TRUE", "palegreen", col)
    col <- sub("FALSE", "gray", col)
    
    ### plotting
    bp <- barplot(x, axes = T, axisnames = T, space = (0.2), ylim = c(0, ylimit), las = 2, col = col, names.arg = names(x))
    
    segments(bp, x - stDevs, bp, x + stDevs, lwd = 2)
    segments(bp - 0.2, x - stDevs, bp + 0.2, x - stDevs, lwd = 2)
    segments(bp - 0.2, x + stDevs, bp + 0.2, x + stDevs, lwd = 2)
    
    
    title(xlab = xlab, ylab = ylab, main = GRP, font.main = 2)
    
}

# *************************************************************************************************************************************


#' Number of PDGs vs COGs/PDG
#' 
#' Barplot that indicates the number of PDGs vs COGs(clustered orthologous groups) in a PDG
#' @param mcl_data format_afterOrtho output
#' @param num an integer indicating where the x axis should end and be compiled
#' @return a barplot with a height determined by the second column and the first column abbreviated to accomodate visual spacing
#' @references Some sort of reference
#' @examples 
#' pdg_v_cog(after_ortho_format_grps,2)
#' #dev.off() #reset margins
#' @export

pdg_v_cog <- function(mcl_data, num = 40, ...) {
    
    cgs <- mcl_data[[1]]
    names <- row.names(cgs)
    
    
    ### Number of COGS
    
    mcl_data <- as.data.frame(cgs)
    mcl_data = cbind(mcl_data, X.NAMES = names)
    
    mcl_data$X.NAMES <- as.character(mcl_data$X.NAMES)
    chk <- strsplit(mcl_data$X.NAMES, ",")
    COGS <- lapply(chk, function(x) length(x))
    COGS <- as.numeric(COGS)
    
    
    ### Number of COGs/PDGs
    
    graph <- as.data.frame(table(COGS), stringsAsFactors = FALSE)
    graph$COGS <- as.numeric(graph$COGS)
    
    sequential <- 1:graph[length(graph$COGS), 1]
    
    
    ### add zeros to columns without any PDGs
    
    zeros <- setdiff(sequential, graph$COGS)
    
    if (length(zeros) != 0) {
        zeros2 <- as.data.frame(zeros)
        zeros2$Freq <- 0
        
        colnames(zeros2) <- c("COGS", "Freq")
        zeros2$COGS <- as.numeric(zeros2$COGS)
        
        sums <- as.data.frame(graph, stringsAsFactors = FALSE)
        sums <- rbind(sums, zeros2)
        
    } else {
        sums <- as.data.frame(graph, stringsAsFactors = FALSE)
        
    }
    
    
    sums <- sums[order(sums[, 1]), ]
    
    ### Sum of last PDGs that have over a given amount of COGS (num)
    
    if (nrow(sums) > num) {
        # if the number of COG limit provided is smaller than the largest number of COGS in PDG
        
        sum <- (colSums(sums[c(num:nrow(sums)), ])[2])
        
        ### Compiling all PDGs with COGs > num and appending to plot
        plot <- sums[1:num, ]
        rownames(plot) <- NULL
        plot[, 1] <- as.numeric(plot[, 1])
        plot[num, 1] <- paste(num, "+", sep = "")
        plot[num, 2] <- sum
        
        r_axis <- round(plot[2, 2], digits = -2)
        plot[1, 2] <- r_axis + 5
        
    } else {
        
        plot <- sums
        rownames(plot) <- NULL
        r_axis <- max(plot$Freq) + 5
        warning(paste("Number provided for max number of COGs in PDG larger than largest number of COGs in PDG:\n", num, 
            " > ", nrow(sums)))
        
    }
    
    ### Plotting dev.off() # reset margins for labels
    barplot <- barplot(plot$Freq, names.arg = plot$COGS, ylim = c(0, r_axis), xlab = expression(paste("COGs PDG"^"-1")), 
        ylab = "PDGs")
    axis(1, at = barplot, labels = plot$COGS)
    
    ### Fitting label on top of bar
    par(new = T, mar = c(2, 2, 2, 2))
    text(x = barplot[1], y = plot[1, 2], label = paste(plot[1, 2]), pos = 3)
    
}


# *************************************************************************************************************************************


#' Phylogenetic Tree with Attached Bar Plot and Standard Error Bars
#' 
#' Presents data for each taxa including standard error bars next to a phylogenetic tree.
#' @param phy Path to tree file
#' @param data R object of phenotype data
#' @param mcl_data mcl matrix (analyze_OrthoMCL output)
#' @param species_colname name of column in data file with taxa designations
#' @param data_colname name of column in data file with data observations
#' @param color optional parameter, (defaults to NULL) assign colors to individual taxa by providing file (format: Taxa | Color)
#' @param GRP optional parameter, (defaults to NULL) a string with the names of chosen group to be colored
#' @param xlabel string to label barplot's x axis
#' @param ... argument to be passed from other methods such as parameters from barplot() function
#' @return A phylogenetic tree with a barplot of the data (with standard error bars) provided matched by taxa.
#' @references Some sort of reference
#' @examples 
#' file <- system.file('sample_data', 'muscle_tree2.dnd', package='MAGNAMWAR')
#' phydataerror(file, pheno_data, mcl_mtrx, species_colname = 'Treatment', data_colname = 'RespVar', GRP='OG5_126778', xlabel='TAG Content')
#' #dev.off() #reset margins and align bars
#' @export


phydataerror <- function(phy, data, mcl_matrix, species_colname, data_colname, color = NULL, GRP = NULL, xlabel = "xlabel", 
    ...) {
    # dev.off() #must reset dev (reset margins) for every run of the function to correctly line up bars and branches
    
    ### building tree
    
    phy = read.tree(phy)
    
    ### .matchDataPhylo from ape library
    
    .matchDataPhylo <- function(x, phy) {
        msg <- "'x' has no (row)names: data are assumed to be in the same order than the tips of the tree"
        labs <- phy$tip.label
        if (is.vector(x)) {
            # also for lists
            if (is.null(names(x))) 
                warning(msg) else x <- x[labs]
        } else {
            if (is.null(rownames(x))) 
                warning(msg) else x <- x[labs, ]
        }
        x
    }
    
    ### Calculating means and sd
    data <- data[, c(species_colname, data_colname)]
    colnames(data) <- c("one", "two")
    ccc <- aggregate(data$one, data["one"], paste, collapse = " ")
    ddd <- aggregate(data$two, data["two"], paste, collapse = " ")
    
    means.sem <- ddply(data, c("one"), summarise, mean = mean(two), sem = sd(two)/sqrt(length(two)))
    means.sem <- transform(means.sem, lower = mean - sem, upper = mean + sem)
    mean <- means.sem[, 2]
    stDevs <- means.sem[, 3]
    names(mean) <- means.sem$one
    names(stDevs) <- means.sem$one
    
    #### plot bar and phylo
    plottree <- plot(phy, show.tip.label = T, x.lim = 0.5, cex = 0.9, label.offset = 0.002)
    
    par(new = T, mar = c(4.7, 4, 3.7, 2))  #line up bars
    offset = 1
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- length(phy$tip.label)
    one2n <- seq_len(n)
    x <- .matchDataPhylo(mean, phy)
    stDevs <- .matchDataPhylo(stDevs, phy)
    x0 <- max(lastPP$xx[one2n]) + offset + 19
    x1 <- x0 + x
    y1 <- lastPP$yy[one2n]
    o <- order(y1)
    x <- if (is.vector(x)) 
        x[o] else x[o, ]
    
    ### Deal with color
    
    if (is.null(color)) {
        
        if (!is.null(GRP)) {
            grep <- mcl_matrix[grep(GRP, mcl_matrix[, 1]), ]
            l <- unlist(strsplit(grep[6], split = "\\|"))
            suppressWarnings(if (is.na(l)) {
                cat("Invalid COG inputted... Printing without specific coloration\n")
            })
            
            coldata <- as.data.frame(x)
            coldata$X1 <- row.names(coldata)
            
            coldata[, 3] <- coldata$X1 %in% l
            
            ### order to match phylo
            vec <- structure(coldata$V3, names = coldata$X1)
            tt <- .matchDataPhylo(vec, phy)
            
            ### assign colors
            tt <- as.character(tt)
            col <- sub("TRUE", "palegreen", tt)
            col <- sub("FALSE", "gray", col)
            
        } else {
            col = "gray80"
        }
    } else {
        
        coldata <- as.data.frame(x)
        coldata$Taxa <- row.names(coldata)
        color <- read.table(color, header = T)
        
        col <- merge(coldata, color)
        
        vec <- structure(col$Color, names = col$Taxa)
        vec <- as.vector(vec)
        names(vec) <- col$Taxa
        tt <- .matchDataPhylo(vec, phy)
        col <- tt
    }
    
    
    if (!is.null(dim(x))) 
        x <- t(x)
    bp <- barplot(x, width = 1, add = F, horiz = TRUE, offset = x0, axes = FALSE, axisnames = FALSE, space = (0.2), xlim = c(0, 
        45), col = col, ...)
    
    ### error bars
    segments(x - stDevs + x0, bp, x + stDevs + x0, bp, lwd = 2)
    segments(x - stDevs + x0, bp - 0.2, x - stDevs + x0, bp + 0.2, lwd = 2)
    segments(x + stDevs + x0, bp - 0.2, x + stDevs + x0, bp + 0.2, lwd = 2)
    px <- pretty(c(0, x + 2))
    axis(1, px + x0, labels = px, line = 0)
    title(xlab = xlabel)
}


# *************************************************************************************************************************************



#' QQPlot
#' 
#' Makes a qqplot of the p-values obtained through analyze_OrthoMCL
#' @param mcl_mtrx matrix generated by analyze_OrthoMCL
#' @return a qqplot of the p-values obtained through analyze_OrthoMCL
#' @references Some sore of reference
#' @examples 
#' qqplotter(mcl_mtrx)
#' @export
#' 
qqplotter <- function(mcl_mtrx) {
    
    
    pvector <- as.numeric(mcl_mtrx[, 2])
    
    pvector <- pvector[!is.na(pvector) & pvector < 1 & pvector > 0]
    
    o = -log10(sort(pvector, decreasing = F))
    
    e = -log10(ppoints(length(pvector)))
    
    plot(e, o, pch = 19, cex = 0.2, xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))), 
        xlim = c(0, max(e)), ylim = c(0, max(o)))
    
    abline(0, 1, col = "red")
}


# *************************************************************************************************************************************


#' Manhattan Plot of All Taxa
#' 
#' Manhattan plot that graphs all p-values for taxa.
#' @param mcl_data format_afterOrtho output
#' @param mcl_mtrx output of analyze_OrthoMCL()
#' @param equation of line of significance, defaults to -log10((.05)/dim(pdgs)[1])
#' @return a manhattan plot
#' @references Some sort of reference
#' @examples 
#' manhat_grp(after_ortho_format, mcl_mtrx)
#' @export

manhat_grp <- function(mcl_data, mcl_mtrx, tree = NULL) {
    
    mcl_data <- mcl_data[[2]]
    
    list_name <- unlist(lapply(colnames(mcl_data), FUN = function(x) rep(x, nrow(mcl_data))))
    
    taxa_prot <- strsplit(as.character(unique(mcl_data)), split = "\\|")
    list_name_parse <- list_name[lapply(taxa_prot, length) > 0]
    taxa_prot <- taxa_prot[lapply(taxa_prot, length) > 0]
    names(taxa_prot) <- list_name_parse
    
    n <- length(taxa_prot[[1]])
    DF <- t(structure(taxa_prot, row.names = c(NA, -n), class = "data.frame"))
    DF <- cbind(DF, row.names(DF))
    DF <- as.data.frame(DF)
    colnames(DF) <- c("Taxa", "Protein ID", "Gene")
    row.names(DF) <- NULL
    
    pvalue <- data.frame(Gene = mcl_mtrx[, "COG"], PValue = (mcl_mtrx[, "corrected_p-val1"]))
    spl <- strsplit(as.character(pvalue$Gene), split = ",")
    pvalue <- data.frame(Gene = unlist(spl), PValue = rep(pvalue$PValue, sapply(spl, length)))
    
    cat("ordering by protein id within taxa...\n")
    finalkey <- merge(pvalue, DF, all = F)
    colnames(finalkey) <- c("Gene", "PValue", "TAXA", "protein_id")
    
    
    ### pulling taxa names
    if (!is.null(tree)) {
        tree_lab <- read.tree(tree)
        taxa_names <- tree_lab$tip.label
    } else {
        taxa_names <- unique(finalkey[, 3])
    }
    taxsub <- subset(finalkey, grepl(paste(taxa_names, collapse = "|"), finalkey[, 3]))
    taxsub2 <- arrange(taxsub, TAXA, protein_id)
    
    ### PROBLEM: WHAT TO DO WHEN finalkey DOESN'T HAVE THE TAXA INCLUDED IN tree?  ex. efOG v efog
    
    names(taxa_names) <- 1:length(taxa_names)
    taxa_names <- as.data.frame(taxa_names)
    taxa_names$CHR <- row.names(taxa_names)
    colnames(taxa_names) <- c("TAXA", "CHR")
    chk <- merge(taxsub2, taxa_names)
    chk$CHR <- as.numeric(chk$CHR)
    
    ### BP positioning
    names(chk$TAXA) <- NULL
    names(chk$protein_id) <- NULL
    chk <- ddply(chk, .(TAXA), transform, BP = seq_along(TAXA))
    chk2 <- data.frame(lapply(chk, as.character), stringsAsFactors = F)
    chk3 <- chk2[chk2[, 3] > 0, ]
    
    ### plotting
    plot <- data.frame(SNP = chk3$Gene, BP = as.numeric(chk3$BP), P = as.numeric(as.character(chk3$PValue)), CHR = as.numeric(chk3$CHR), 
        TAXA = chk3$TAXA)
    plot <- na.omit(plot)
    cat("creating manhattan plot...\n")
    manhattan(plot, chrlabs = as.character(taxa_names$TAXA), las = 2, xlab = "", cex = 0.25, suggestiveline = -log10((0.05)/dim(mcl_mtrx)[1]), 
        genomewideline = FALSE)
    cat("finished.\n")
}

