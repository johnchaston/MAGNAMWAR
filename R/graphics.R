#' Plot of a PDG and Data with Standard Error Bars
#' 
#' Bar plot of PDG vs phenotype data with presence of taxa in PDG indicated by color
#' @param data R object of phenotype data
#' @param mcl_matrix analyze_OrthoMCL output
#' @param species_colname name of column in phenotypic data file with taxa designations
#' @param data_colname name of column in phenotypic data file with data observations
#' @param OG optional parameter, a string with the name of chosen group (OG) to be colored
#' @param xlab string to label barplot's x axis
#' @param ylab string to label barplot's y axis
#' @param ylimit optional parameter to limit y axis
#' @param tree optional parameter (defaults to NULL) Path to tree file, orders the taxa by phylogenetic distribution, else it defaults to alphabetical
#' @param order vector with order of taxa names for across the x axis (defaults to alpha ordering)
#' @param main_title string for title of the plot (defaults to OG)
#' @return a barplot with taxa vs phenotypic data complete with standard error bars
#' @examples
#' # Not run - changes options
#' \dontrun{ 
#' pdgplot(pheno_data, mcl_mtrx, 'OG5_126778', 'Treatment', 'RespVar', ylimit=12)
#' }
#' @export
#' 

pdgplot <- function(data, mcl_matrix, OG = "NONE", species_colname,
                    data_colname, xlab = "Taxa", ylab = "Data",
                    ylimit = NULL, tree = NULL, order = NULL,
                    main_title = NULL) {


    if (species_colname %in% colnames(data) && data_colname %in% colnames(data)) {
    } else {
        stop("Invalid column names specified for phenotype data file
             \n\tSpecies Column Name: ", species_colname,
             "\n\tData Column Name: ", data_colname)
    }

    grep <- mcl_matrix[grep(paste("^", OG, "$", sep = ""), mcl_matrix[, 1]), ]
    if (OG == "NONE") {
        warning("No OG inputted...
                 Printing without specific coloration\n")
    } else if (anyNA(grep[6])) {
        warning("Invalid OG inputted...
                 Printing without specific coloration\n")
        l <- unlist(strsplit(as.character(unlist(grep[6])), split = "\\|"))
    } else {
        l <- unlist(strsplit(as.character(unlist(grep[6])), split = "\\|"))
    }

    ### Calculating st dev and means
    data <- data[, c(species_colname, data_colname)]
    colnames(data) <- c("one", "two")

    two <- "two"
    sem <- "sem"

    ccc <- aggregate(data$one, data["one"], paste, collapse = " ")
    ddd <- aggregate(data$two, data["two"], paste, collapse = " ")

    means.sem <- plyr::ddply(data, c("one"), dplyr::summarise,
                             mean = mean(two), sem = sd(two) / sqrt(length(two)))
    means.sem <- transform(means.sem, lower = mean - sem, upper = mean + sem)
    mean <- means.sem[, 2]
    stdevs <- means.sem[, 3]

    names(mean) <- means.sem$one
    names(stdevs) <- means.sem$one

    x <- mean

    ### ordering by tree, order vector, or alphabetical

    if (!is.null(tree)) {
        # order by tree
        tree <- ape::read.tree(tree)

        x <- x[tree$tip.label]
        stdevs <- stdevs[tree$tip.label]
    } else {

        if (!is.null(order)) {
            # order by order vector

            x <- x[order]
            stdevs <- stdevs[order]

            x <- x[!is.na(x)]
            stdevs <- stdevs[!is.na(stdevs)]
        } else {
            # order by alpha
            sorted <- sort(names(x))
            x <- x[sorted]
            stdevs <- stdevs[sorted]

        }
    }

    if (is.null(ylimit)) {
        ylimit <- max(x) + max(x) / 4
    }

    ### coloring for graph

    col <- names(x) %in% l
    col <- sub("TRUE", "palegreen", col)
    col <- sub("FALSE", "gray", col)

    ### plotting
    bp <- barplot(x, axes = T, axisnames = T, space = (0.2),
                  ylim = c(0, ylimit), las = 2,
                  col = col, names.arg = names(x))

    segments(bp, x - stdevs, bp, x + stdevs, lwd = 2)
    segments(bp - 0.2, x - stdevs, bp + 0.2, x - stdevs, lwd = 2)
    segments(bp - 0.2, x + stdevs, bp + 0.2, x + stdevs, lwd = 2)

    if (is.null(main_title)) {

        if (OG != "NONE" && !anyNA(grep[6])) {
            main_title <- OG
            title(xlab = xlab, ylab = ylab, main = main_title, font.main = 2)
        } else {
            title(xlab = xlab, ylab = ylab, font.main = 2)
        }

    } else {
        title(xlab = xlab, ylab = ylab, main = main_title, font.main = 2)
    }
}

# --------------------------


#' Number of PDGs vs OGs/PDG
#' 
#' Barplot that indicates the number of PDGs vs OGs(clustered orthologous groups) in a PDG
#' @param mcl_data format_afterOrtho output
#' @param num an integer indicating where the x axis should end and be compiled
#' @param ... args to be passed to barplot
#' @return a barplot with a height determined by the second column and the first column abbreviated to accomodate visual spacing
#' @examples 
#' # Not run - changes options
#' \dontrun{
#' pdg_v_OG(after_ortho_format_grps,2)
#' }
#' @export

pdg_v_OG <- function(mcl_data, num = 40, ...) {

    cgs <- mcl_data[[1]]
    names <- row.names(cgs)

    ### Number of OGS

    mcl_data <- as.data.frame(cgs)
    mcl_data <- cbind(mcl_data, X.NAMES = names)

    mcl_data$X.NAMES <- as.character(mcl_data$X.NAMES)
    chk <- strsplit(mcl_data$X.NAMES, ",")
    OGS <- lapply(chk, function(x) length(x))
    OGS <- as.numeric(OGS)

    ### Number of OGs/PDGs

    graph <- as.data.frame(table(OGS), stringsAsFactors = FALSE)
    graph$OGS <- as.numeric(graph$OGS)

    sequential <- 1:graph[length(graph$OGS), 1]


    ### add zeros to columns without any PDGs

    zeros <- setdiff(sequential, graph$OGS)

    if (length(zeros) != 0) {
        zeros2 <- as.data.frame(zeros)
        zeros2$Freq <- 0

        colnames(zeros2) <- c("OGS", "Freq")
        zeros2$OGS <- as.numeric(zeros2$OGS)

        sums <- as.data.frame(graph, stringsAsFactors = FALSE)
        sums <- rbind(sums, zeros2)

    } else {
        sums <- as.data.frame(graph, stringsAsFactors = FALSE)
    }
    
    sums <- sums[order(sums[, 1]), ]

    ### Sum of last PDGs that have over a given amount of OGS (num)

    if (nrow(sums) > num) {
        # if the number of OG limit provided is smaller than the
        # largest number of OGS in PDG

        sum <- (colSums(sums[c(num:nrow(sums)), ])[2])

        ### Compiling all PDGs with OGs > num and appending to plot
        plot <- sums[1:num, ]
        rownames(plot) <- NULL
        plot[, 1] <- as.numeric(plot[, 1])
        plot[num, 1] <- paste(num, "+", sep = "")
        plot[num, 2] <- sum

        r_axis <- plot[1, 2] + 5

    } else {

        plot <- sums
        rownames(plot) <- NULL
        r_axis <- max(plot$Freq) + 5
        warning(paste("Number provided for max number of OGs in PDG larger
                      than largest number of OGs in PDG:\n",
                      num, " > ", nrow(sums)))
    }

    barplot <- barplot(plot$Freq, names.arg = plot$OGS, ylim = c(-1, r_axis),
                       xlab = expression(paste("OGs PDG"^"-1")),
                       ylab = "PDGs", ...)
    axis(1, at = barplot, labels = plot$OGS, las = 1)

    ### Fitting label on top of bar
    orig_par <- par()$mar
    par(new = T, mar = c(2, 2, 2, 2))
    text(x = barplot[1], y = plot[1, 2], label = paste(plot[1, 2]), pos = 3)

    par(mar = orig_par, no.readonly = T)
}


# --------------------------


#' Phylogenetic Tree with Attached Bar Plot and Standard Error Bars
#' 
#' Presents data for each taxa including standard error bars next to a phylogenetic tree.
#' @param phy Path to tree file
#' @param data R object of phenotype data
#' @param mcl_matrix analyze_OrthoMCL output
#' @param species_colname name of column in data file with taxa designations
#' @param data_colname name of column in data file with data observations
#' @param color optional parameter, (defaults to NULL) assign colors to individual taxa by providing file (format: Taxa | Color)
#' @param OG optional parameter, (defaults to NULL) a string with the names of chosen group to be colored
#' @param xlabel string to label barplot's x axis
#' @param ... argument to be passed from other methods such as parameters from barplot() function
#' @return A phylogenetic tree with a barplot of the data (with standard error bars) provided matched by taxa.
#' @references Some sort of reference
#' @examples 
#' 
#' # Not run - changes options
#' \dontrun{
#' file <- system.file('extdata', 'muscle_tree2.dnd', package='MAGNAMWAR')
#' phydataerror(file, pheno_data, mcl_mtrx, species_colname = 'Treatment', data_colname = 'RespVar',
#'  OG='OG5_126778', xlabel='TAG Content')
#' }
#' @importFrom graphics abline axis barplot par plot segments text title
#' @importFrom stats aggregate na.omit ppoints sd
#' @importFrom utils read.csv read.delim read.table write.csv write.table
#' @export

phydataerror <- function(phy, data, mcl_matrix, species_colname, data_colname,
                         color = NULL, OG = NULL, xlabel = "xlabel", ...) {

    ### building tree

    orig_par <- par()$mar      # make a copy of current settings
  
    phy <- ape::read.tree(phy)

    if (xlabel == "xlabel") {
        xlabel <- data_colname
    }
    ### Calculating means and sd
    data <- data[, c(species_colname, data_colname)]
    colnames(data) <- c("one", "two")
    ccc <- aggregate(data$one, data["one"], paste, collapse = " ")
    ddd <- aggregate(data$two, data["two"], paste, collapse = " ")

    two <- "two"
    sem <- "sem"

    means.sem <- plyr::ddply(data, c("one"), dplyr::summarise,
                             mean = mean(two),
                             sem = sd(two)/sqrt(length(two)))
    means.sem <- transform(means.sem, lower = mean - sem, upper = mean + sem)
    mean <- means.sem[, 2]
    stdevs <- means.sem[, 3]
    names(mean) <- means.sem$one
    names(stdevs) <- means.sem$one

    #### plot bar and phylo
    plottree <- plot(phy, show.tip.label = T, x.lim = 0.5,
                     cex = 0.9, label.offset = 0.002)

    par(new = T, mar = c(4.7, 4, 3.7, 2))  #line up bars
    offset <- 1
    lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
    n <- length(phy$tip.label)
    leng <- seq_len(n)
    labels <- phy$tip.label
    x <- mean[labels]
    stdevs <- stdevs[labels]
    x0 <- max(lastPP$xx[leng]) + offset + 19
    x1 <- x0 + x
    y1 <- lastPP$yy[leng]
    o <- order(y1)
    x <- if (is.vector(x))
        x[o] else x[o, ]

    ### Deal with color

    if (is.null(color)) {

        if (!is.null(OG)) {
            grep <- mcl_matrix[grep(OG, mcl_matrix[, 1]), ]
            l <- unlist(strsplit(grep[6], split = "\\|"))
            suppressWarnings(if (is.na(l)) {
                cat("Invalid OG inputted... Printing without
                    specific coloration\n")
            })

            coldata <- as.data.frame(x)
            coldata$X1 <- row.names(coldata)

            coldata[, 3] <- coldata$X1 %in% l

            ### order to match phylo
            vec <- structure(coldata$V3, names = coldata$X1)
            tt <- vec[labels]

            ### assign colors
            tt <- as.character(tt)
            col <- sub("TRUE", "palegreen", tt)
            col <- sub("FALSE", "gray", col)

        } else {
            col <- "gray80"
        }
    } else {

        coldata <- as.data.frame(x)
        coldata$Taxa <- row.names(coldata)
        color <- read.table(color, header = T)

        col <- merge(coldata, color)

        vec <- structure(col$Color, names = col$Taxa)
        vec <- as.vector(vec)
        names(vec) <- col$Taxa
        tt <- vec[labels]
        col <- tt
    }

    if (!is.null(dim(x)))
        x <- t(x)
    bp <- barplot(x, width = 1, add = F, horiz = TRUE, offset = x0,
                  axes = FALSE, axisnames = FALSE, space = (0.2),
                  xlim = c(0, 45), col = col, ...)

    ### error bars
    segments(x - stdevs + x0, bp, x + stdevs + x0, bp, lwd = 2)
    segments(x - stdevs + x0, bp - 0.2, x - stdevs + x0, bp + 0.2, lwd = 2)
    segments(x + stdevs + x0, bp - 0.2, x + stdevs + x0, bp + 0.2, lwd = 2)
    px <- pretty(c(0, x + 2))
    axis(1, px + x0, labels = px, line = 0)
    title(xlab = xlabel)
    par(mar = orig_par, no.readonly = T)   # restore original settings
    
}


# --------------------------



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

    o <- -log10(sort(pvector, decreasing = F))

    e <- -log10(ppoints(length(pvector)))

    plot(e, o, pch = 19, cex = 0.2,
         xlab = expression(Expected ~ ~-log[10](italic(p))),
         ylab = expression(Observed ~ ~-log[10](italic(p))),
         xlim = c(0, max(e)),
         ylim = c(0, max(o)))

    abline(0, 1, col = "red")
}


# --------------------------


#' Manhattan Plot of All Taxa
#' 
#' Manhattan plot that graphs all p-values for taxa.
#' @param mcl_data format_afterOrtho output
#' @param mcl_mtrx output of analyze_OrthoMCL()
#' @param tree tree file optional, used for ordering taxa along x axis
#' @return a manhattan plot
#' @references Some sort of reference
#' @examples 
#' manhat_grp(after_ortho_format, mcl_mtrx)
#' 
#' #@param equation of line of significance, defaults to -log10((.05)/dim(pdgs)[1])
#' @importFrom plyr as.quoted .
#' @importFrom ape .PlotPhyloEnv
#' @export

manhat_grp <- function(mcl_data, mcl_mtrx, tree = NULL) {

    mcl_data <- mcl_data[[2]]

    list_name <- unlist(lapply(colnames(mcl_data),
                               FUN = function(x) rep(x, nrow(mcl_data))))

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

    pvalue <- data.frame(Gene = mcl_mtrx[, "OG"], PValue = (mcl_mtrx[, 3]))
    spl <- strsplit(as.character(pvalue$Gene), split = ",")
    pvalue <- data.frame(Gene = unlist(spl),
                         PValue = rep(pvalue$PValue, sapply(spl, length)))

    cat("ordering by protein id within taxa...\n")
    finalkey <- merge(pvalue, DF, all = F)
    colnames(finalkey) <- c("Gene", "PValue", "TAXA", "protein_id")

    TAXA <- "TAXA"

    ### pulling taxa names
    if (!is.null(tree)) {
        tree_lab <- ape::read.tree(tree)
        taxa_names <- tree_lab$tip.label
    } else {
        taxa_names <- unique(finalkey[, 3])
    }

    taxsub <- subset(finalkey, grepl(paste(taxa_names, collapse = "|"),
                                     finalkey[, 3]))
    taxsub$TAXA <- factor(taxsub$TAXA)

    taxsub2 <- taxsub[order(taxsub$TAXA, taxsub$protein_id), ]
    row.names(taxsub2) <- NULL

    ### WHAT TO DO WHEN finalkey NOT HAVE TAXA INCLUDED IN tree?ex. efOG v efog

    names(taxa_names) <- 1:length(taxa_names)
    taxa_names <- as.data.frame(taxa_names)
    taxa_names$CHR <- row.names(taxa_names)
    colnames(taxa_names) <- c("TAXA", "CHR")
    chk <- merge(taxsub2, taxa_names)
    chk$CHR <- as.numeric(chk$CHR)

    ### BP positioning
    names(chk$TAXA) <- NULL
    names(chk$protein_id) <- NULL
    chk <- plyr::ddply(chk, .(TAXA), transform, BP = seq_along(TAXA))
    chk2 <- data.frame(lapply(chk, as.character), stringsAsFactors = F)
    chk3 <- chk2[chk2[, 3] > 0, ]

    ### plotting
    plot <- data.frame(SNP = chk3$Gene,
                       BP = as.numeric(chk3$BP),
                       P = as.numeric(as.character(chk3$PValue)),
                       CHR = as.numeric(chk3$CHR),
                       TAXA = chk3$TAXA)
    plot <- na.omit(plot)
    cat("creating manhattan plot...\n")

    final_taxa_names <- unique(plot$TAXA)
    fin_taxa <- taxa_names$TAXA[taxa_names$TAXA %in% final_taxa_names]
    fin_taxa <- factor(fin_taxa)

    qqman::manhattan(plot, xlab = "", chrlabs = as.character(fin_taxa),
                     las = 2, genomewideline = FALSE,
                     suggestiveline = -log10((0.05) / dim(mcl_mtrx)[1]))
    cat("finished.\n")
}
