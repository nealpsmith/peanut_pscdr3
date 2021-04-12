blankPlot <- ggplot()+geom_blank(aes(1,1))+
    theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank()
    )
# modification of vis.top.proportions
vis.top.proportions <- function (.data, .head = c(10, 100, 1000, 10000, 30000, 1e+05, 3e+05, 1e+06), .col = "Read.count", .x_title = "Sample", .title = "Summary proportion of the top N clones", .colors = "Dark2", .legend = TRUE)
{
    if (has.class(.data, "data.frame")) {
        .data <- list(Sample = .data)
    }
    res <- sapply(.head, function(h) top.proportion(.data, h,
                                                    .col))
    tmp <- res
    if (is.null(dim(tmp))) {
        tmp <- t(as.matrix(tmp))
        res <- t(as.matrix(res))
    }
    for (i in 2:ncol(res)) {
        tmp[, i] <- res[, i] - res[, i - 1]
    }
    res <- tmp
    colnames(res) <- paste0("[", c(1, .head[-length(.head)] +
                                       1), ":", .head, ")")
    res <- as.data.frame(res)
    res$People <- factor(row.names(res), levels = row.names(res))
    res <- melt(res)
    g <- ggplot() +
        geom_bar(aes(x = People, y = value, fill = variable), data = res, stat = "identity", position = "stack", colour = "black", show.legend = .legend)  +
        scale_fill_brewer(palette=.colors) +
        theme_linedraw() +
        theme(axis.text.x = element_blank()) +
        ylab("Clonal proportion") +
        xlab(.x_title) +
        ggtitle(.title) +
        guides(fill = guide_legend("Top N clones"))
    g
}

# cleanup function
cleanup_alleles <- function(y){
    z <- list()
    keep.names <- names(y)
    for (step in 1:length(y)) {
        y[[step]]$V.gene[y[[step]]$V.gene == ""] <- lapply(y[[step]]$V.ties[y[[step]]$V.gene == ""], function(x) {unlist(strsplit(x, ","))[1]})
        y[[step]]$V.gene[y[[step]]$V.gene == "unresolved"] <- lapply(y[[step]]$V.ties[y[[step]]$V.gene == "unresolved"], function(x) {unlist(strsplit(x, ","))[1]})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("TCR", "TR", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("^[A-Z]*?0", "TRBV", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("-0", "-", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("30-.", "30", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("TRBV2-.", "TRBV2", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("13-.", "13", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("14-.", "14", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("15-.", "15", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("16-.", "16", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("18-.", "18", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("19-.", "19", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("27-.", "27", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("28-.", "28", x)})
        y[[step]]$V.gene <- lapply(y[[step]]$V.gene, function(x) {gsub("30-.", "30", x)})
        y[[step]]$V.gene <- as.character(y[[step]]$V.gene)
        print(names(y[step]))
        z[step] <- y[step]
    }
    names(z) <- keep.names
    z
}

# modification of parse functions
parse.adaptive <- function(.filename, .nuc.seq = 'rearrangement', .aa.seq = 'amino_acid', .reads = 'templates', .barcodes = NA, .vgenes = 'v_gene', .vties = 'v_gene_ties', .jgenes = 'j_gene', .dgenes = 'd_gene', .vend = NA, .jstart = 'j_index', .dalignments = NA, .vd.insertions = 'n1_insertions', .dj.insertions = 'n2_insertions', .total.insertions = NA, .frame = 'frame_type', .skip = 0, .sep = "\t")
{
    .make.names <- function(.char) {
        if (is.na(.char[1])) {
            NA
        }
        else {
            make.names(.char)
        }
    }
    .nuc.seq <- .make.names(.nuc.seq)
    .aa.seq <- .make.names(.aa.seq)
    .reads <- .make.names(.reads)
    .barcodes <- .make.names(.barcodes)
    .vgenes <- .make.names(.vgenes)
    .vties <- .make.names(.vties)
    .jgenes <- .make.names(.jgenes)
    .dgenes <- .make.names(.dgenes)
    .vend <- .make.names(.vend)
    .jstart <- .make.names(.jstart)
    .dalignments <- .make.names(.dalignments)
    .vd.insertions <- .make.names(.vd.insertions)
    .dj.insertions <- .make.names(.dj.insertions)
    .total.insertions <- .make.names(.total.insertions)
    .frame <- .make.names(.frame)
    f <- file(.filename, "r")
    l <- readLines(f, 1)
    if (length(grep("MiTCRFullExportV1.1", l, fixed = T))) {
        .skip <- 1
    }
    if (length(strsplit(l, "-", T)[[1]]) == 3) {
        if (strsplit(l, "-", T)[[1]][2] == "header") {
            .reads <- "count"
            .barcodes <- "count"
            .skip <- 1
        }
    }
    close(f)
    table.colnames <- make.names(read.table(gzfile(.filename), sep = .sep, skip = .skip, nrows = 1, stringsAsFactors = F, strip.white = T, comment.char = "", quote = "")[1, ])
    swlist <- list("character", "character", "integer", "integer", "character", "character", "character", "character", "integer", "integer", "integer", "integer", "integer", "integer", "character")
    names(swlist) <- c(.nuc.seq, .aa.seq, .reads, .barcodes, .vgenes, .vties, .jgenes, .dgenes, .vend, .jstart, .dalignments, .vd.insertions, .dj.insertions, .total.insertions, .frame)
    swlist <- c(swlist, "NULL")
    col.classes <- unlist(sapply(table.colnames, function(x) {
        do.call(switch, c(x, swlist))
    }, USE.NAMES = F))
    suppressWarnings(df <- read.table(file = gzfile(.filename),
                                      header = T, colClasses = col.classes, sep = .sep, skip = .skip,
                                      strip.white = T, comment.char = "", quote = ""))
    df$Read.proportion <- df[, make.names(.reads)]/sum(df[, make.names(.reads)])
    .read.prop <- "Read.proportion"
    if (is.na(.barcodes)) {
        .barcodes <- "Umi.count"
        df$Umi.count <- NA
        df$Umi.proportion <- NA
    }
    else {
        df$Umi.proportion <- df[, make.names(.barcodes)]/sum(df[,
                                                                make.names(.barcodes)])
    }
    .umi.prop <- "Umi.proportion"
    if (is.na(.aa.seq)) {
        df$CDR3.amino.acid.sequence <- bunch.translate(df$CDR3.nucleotide.sequence)
        .aa.seq <- "CDR3.amino.acid.sequence"
    }
    recomb_type = "Undeterm"
    if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRA",
                                                     "TRAV", "TRGV", "IGKV", "IGLV"))) {
        recomb_type = "VJ"
    }
    else if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRB",
                                                          "TRBV", "TRDV", "IGHV"))) {
        recomb_type = "VDJ"
    }
    if (!(.vd.insertions %in% table.colnames)) {
        .vd.insertions <- "VD.insertions"
        if (!is.na(.vend) && !is.na(.dalignments)) {
            if (recomb_type == "VJ") {
                df$VD.insertions <- -1
            }
            else if (recomb_type == "VDJ") {
                df$VD.insertions <- df[[.dalignments1]] - df[[.vend]] -
                    1
                df$VD.insertions[df[[.dalignments1]] == -1] <- -1
                df$VD.insertions[df[[.vend]] == -1] <- -1
            }
            else {
                df$VD.insertions <- -1
            }
        }
        else {
            df$VD.insertions <- -1
            df$V.end <- -1
            df$D5.end <- -1
            df$D3.end <- -1
            .vend <- "V.end"
            .dalignments <- c("D5.end", "D3.end")
        }
    }
    if (!(.dj.insertions %in% table.colnames)) {
        .dj.insertions <- "DJ.insertions"
        if (!is.na(.jstart) && !is.na(.dalignments)) {
            if (recomb_type == "VJ") {
                df$DJ.insertions <- -1
            }
            else if (recomb_type == "VDJ") {
                df$DJ.insertions <- df[[.jstart]] - df[[.dalignments2]] -
                    1
                df$DJ.insertions[df[[.dalignments2]] == -1] <- -1
                df$DJ.insertions[df[[.jstart]] == -1] <- -1
            }
            else {
                df$DJ.insertions <- -1
            }
        }
        else {
            df$DJ.insertions <- -1
            df$J.start <- -1
            df$D5.end <- -1
            df$D3.end <- -1
            .jstart <- "J.start"
            .dalignments <- c("D5.end", "D3.end")
        }
    }
    if (!(.total.insertions %in% table.colnames)) {
        .total.insertions <- "Total.insertions"
        df$Total.insertions <- -1
        if (recomb_type == "VJ") {
            df$Total.insertions <- df[[.jstart]] - df[[.vend]] -
                1
            df$Total.insertions[df[[.vend]] == -1] <- -1
            df$Total.insertions[df[[.jstart]] == -1] <- -1
        }
        else if (recomb_type == "VDJ") {
            df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
            df$Total.insertions[df$Total.insertions < 0] <- -1
        }
    }
    if (is.na(.dgenes)) {
        df$D.gene <- ""
        .dgenes <- "D.gene"
    }
    df <- df[, make.names(c(.barcodes, .umi.prop, .reads, .read.prop, .nuc.seq, .aa.seq, .vgenes, .vties, .jgenes, .dgenes,.jstart, .vd.insertions, .dj.insertions, .frame))]
    colnames(df) <- c("Umi.count", "Umi.proportion", "Read.count", "Read.proportion", "CDR3.nucleotide.sequence", "CDR3.amino.acid.sequence", "V.gene", "V.ties", "J.gene", "D.gene", "J.start", "VD.insertions", "DJ.insertions", "Frame")
    df <- subset(df, Frame == "In") #remove unproductive reads
    df$Read.proportion <- df$Read.count / sum(df$Read.count) #adjust proportion accordingly
    df # return data frame
}

parse.adaptive.file.list <- function(.filenames, .namelist = NA)
{
    .remove.ext <- function(.str) {
        gsub(pattern = ".*/|[.].*$", replacement = "", x = .str)
    }
    .filenames <- as.list(.filenames)
    datalist <- list()
    for (i in 1:length(.filenames)) {
        cat(i, "/", length(.filenames), "  Parsing \"", .filenames[[i]],
            "\"\n\t", sep = "")
        datalist[[i]] <- parse.adaptive(.filenames[[i]])
        cat("Done. Cloneset with", nrow(datalist[[i]]), "clonotypes.\n")
        flush.console()
    }
    if (is.na(.namelist)) {
        namelist <- lapply(X = .filenames, FUN = .remove.ext)
        names(datalist) <- unlist(namelist)
    }
    datalist
}

parse.adaptive.folder <- function(.folderpath)
{
    parse.adaptive.file.list(list.files(.folderpath, pattern = "*.tsv", full.names = T))
}

vis.count.len <- function(.data, .ncol = 3, .name = "", .col = "Read.count") {
    if (has.class(.data, "list")) {
        return(do.call(grid.arrange, c(lapply(1:length(.data), function(i) vis.count.len(.data[[i]], .col = .col, .name = names(.data)[i])), ncol = .ncol)))
        }
    tmp <- aggregate(as.formula(paste0(.col, " ~ nchar(CDR3.amino.acid.sequence)")), .data, sum)

    names(tmp) <- c("Lengths", "Count")
    ggplot(tmp) +
        geom_bar(aes(x = Lengths, y = Count, fill = Count), stat = "identity", color = "black") +
        ggtitle(.name) +
        theme_linedraw()
}

# modification of parse functions for reading the control data set
parse.adaptive.hip <- function(.filename, .nuc.seq = 'nucleotide', .aa.seq = 'aminoAcid', .reads = 'count (reads)', .barcodes = NA, .vgenes = 'vGeneName', .vties = 'vFamilyTies', .jgenes = 'jGeneName', .dgenes = 'dGeneName', .vend = NA, .jstart = 'jIndex', .dalignments = NA, .vd.insertions = 'n1Insertion', .dj.insertions = 'n2Insertion', .total.insertions = NA, .skip = 0, .sep = "\t")
{
    .make.names <- function(.char) {
        if (is.na(.char[1])) {
            NA
        }
        else {
            make.names(.char)
        }
    }
    .nuc.seq <- .make.names(.nuc.seq)
    .aa.seq <- .make.names(.aa.seq)
    .reads <- .make.names(.reads)
    .barcodes <- .make.names(.barcodes)
    .vgenes <- .make.names(.vgenes)
    .vties <- .make.names(.vties)
    .jgenes <- .make.names(.jgenes)
    .dgenes <- .make.names(.dgenes)
    .vend <- .make.names(.vend)
    .jstart <- .make.names(.jstart)
    .dalignments <- .make.names(.dalignments)
    .vd.insertions <- .make.names(.vd.insertions)
    .dj.insertions <- .make.names(.dj.insertions)
    .total.insertions <- .make.names(.total.insertions)
    f <- file(.filename, "r")
    l <- readLines(f, 1)
    if (length(grep("MiTCRFullExportV1.1", l, fixed = T))) {
        .skip <- 1
    }
    if (length(strsplit(l, "-", T)[[1]]) == 3) {
        if (strsplit(l, "-", T)[[1]][2] == "header") {
            .reads <- "count"
            .barcodes <- "count"
            .skip <- 1
        }
    }
    close(f)
    table.colnames <- make.names(read.table(gzfile(.filename), sep = .sep, skip = .skip, nrows = 1, stringsAsFactors = F, strip.white = T, comment.char = "", quote = "")[1, ])
    swlist <- list("character", "character", "integer", "integer", "character", "character", "character", "character", "integer", "integer", "integer", "integer", "integer", "integer")
    names(swlist) <- c(.nuc.seq, .aa.seq, .reads, .barcodes, .vgenes, .vties, .jgenes, .dgenes, .vend, .jstart, .dalignments, .vd.insertions, .dj.insertions, .total.insertions)
    swlist <- c(swlist, "NULL")
    col.classes <- unlist(sapply(table.colnames, function(x) {
        do.call(switch, c(x, swlist))
    }, USE.NAMES = F))
    suppressWarnings(df <- read.table(file = gzfile(.filename),
                                      header = T, colClasses = col.classes, sep = .sep, skip = .skip,
                                      strip.white = T, comment.char = "", quote = ""))
    df$Read.proportion <- df[, make.names(.reads)]/sum(df[, make.names(.reads)])
    .read.prop <- "Read.proportion"
    if (is.na(.barcodes)) {
        .barcodes <- "Umi.count"
        df$Umi.count <- NA
        df$Umi.proportion <- NA
    }
    else {
        df$Umi.proportion <- df[, make.names(.barcodes)]/sum(df[,
                                                                make.names(.barcodes)])
    }
    .umi.prop <- "Umi.proportion"
    if (is.na(.aa.seq)) {
        df$CDR3.amino.acid.sequence <- bunch.translate(df$CDR3.nucleotide.sequence)
        .aa.seq <- "CDR3.amino.acid.sequence"
    }
    recomb_type = "Undeterm"
    if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRA",
                                                     "TRAV", "TRGV", "IGKV", "IGLV"))) {
        recomb_type = "VJ"
    }
    else if (sum(substr(head(df)[[.vgenes]], 1, 4) %in% c("TCRB",
                                                          "TRBV", "TRDV", "IGHV"))) {
        recomb_type = "VDJ"
    }
    if (!(.vd.insertions %in% table.colnames)) {
        .vd.insertions <- "VD.insertions"
        if (!is.na(.vend) && !is.na(.dalignments)) {
            if (recomb_type == "VJ") {
                df$VD.insertions <- -1
            }
            else if (recomb_type == "VDJ") {
                df$VD.insertions <- df[[.dalignments1]] - df[[.vend]] -
                    1
                df$VD.insertions[df[[.dalignments1]] == -1] <- -1
                df$VD.insertions[df[[.vend]] == -1] <- -1
            }
            else {
                df$VD.insertions <- -1
            }
        }
        else {
            df$VD.insertions <- -1
            df$V.end <- -1
            df$D5.end <- -1
            df$D3.end <- -1
            .vend <- "V.end"
            .dalignments <- c("D5.end", "D3.end")
        }
    }
    if (!(.dj.insertions %in% table.colnames)) {
        .dj.insertions <- "DJ.insertions"
        if (!is.na(.jstart) && !is.na(.dalignments)) {
            if (recomb_type == "VJ") {
                df$DJ.insertions <- -1
            }
            else if (recomb_type == "VDJ") {
                df$DJ.insertions <- df[[.jstart]] - df[[.dalignments2]] -
                    1
                df$DJ.insertions[df[[.dalignments2]] == -1] <- -1
                df$DJ.insertions[df[[.jstart]] == -1] <- -1
            }
            else {
                df$DJ.insertions <- -1
            }
        }
        else {
            df$DJ.insertions <- -1
            df$J.start <- -1
            df$D5.end <- -1
            df$D3.end <- -1
            .jstart <- "J.start"
            .dalignments <- c("D5.end", "D3.end")
        }
    }
    if (!(.total.insertions %in% table.colnames)) {
        .total.insertions <- "Total.insertions"
        df$Total.insertions <- -1
        if (recomb_type == "VJ") {
            df$Total.insertions <- df[[.jstart]] - df[[.vend]] -
                1
            df$Total.insertions[df[[.vend]] == -1] <- -1
            df$Total.insertions[df[[.jstart]] == -1] <- -1
        }
        else if (recomb_type == "VDJ") {
            df$Total.insertions <- df[[.vd.insertions]] + df[[.dj.insertions]]
            df$Total.insertions[df$Total.insertions < 0] <- -1
        }
    }
    if (is.na(.dgenes)) {
        df$D.gene <- ""
        .dgenes <- "D.gene"
    }
    df <- df[, make.names(c(.barcodes, .umi.prop, .reads, .read.prop, .nuc.seq, .aa.seq, .vgenes, .vties, .jgenes, .dgenes,.jstart, .vd.insertions, .dj.insertions))]
    colnames(df) <- c("Umi.count", "Umi.proportion", "Read.count", "Read.proportion", "CDR3.nucleotide.sequence", "CDR3.amino.acid.sequence", "V.gene", "V.ties", "J.gene", "D.gene", "J.start", "VD.insertions", "DJ.insertions")
    df # return data frame
}

parse.adaptive.hip.file.list <- function(.filenames, .namelist = NA)
{
    .remove.ext <- function(.str) {
        gsub(pattern = ".*/|[.].*$", replacement = "", x = .str)
    }
    .filenames <- as.list(.filenames)
    datalist <- list()
    for (i in 1:length(.filenames)) {
        cat(i, "/", length(.filenames), "  Parsing \"", .filenames[[i]],
            "\"\n\t", sep = "")
        datalist[[i]] <- parse.adaptive.hip(.filenames[[i]])
        cat("Done. Cloneset with", nrow(datalist[[i]]), "clonotypes.\n")
        flush.console()
    }
    if (is.na(.namelist)) {
        namelist <- lapply(X = .filenames, FUN = .remove.ext)
        names(datalist) <- unlist(namelist)
    }
    datalist
}

parse.adaptive.hip.folder <- function(.folderpath)
{
    parse.adaptive.hip.file.list(list.files(.folderpath, full.names = T))
}

vis.count.len <- function(.data, .ncol = 3, .name = "", .col = "Read.count") {
    if (has.class(.data, "list")) {
        return(do.call(grid.arrange, c(lapply(1:length(.data), function(i) vis.count.len(.data[[i]], .col = .col, .name = names(.data)[i])), ncol = .ncol)))
    }
    tmp <- aggregate(as.formula(paste0(.col, " ~ nchar(CDR3.amino.acid.sequence)")), .data, sum)

    names(tmp) <- c("Lengths", "Count")
    ggplot(tmp) +
        geom_bar(aes(x = Lengths, y = Count, fill = Count), stat = "identity", color = "black") +
        ggtitle(.name) +
        theme_linedraw()
}
