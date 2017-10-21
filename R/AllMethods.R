# initialize ####

setMethod("initialize", "Transcriptogram",
    function(.Object, association, ordering,
        radius, status) {
        .Object@association = association
        .Object@ordering = ordering
        .Object@status = status
        .Object@radius = radius
        .Object
    })

# setRadius ####

#' @rdname setRadius-method

setMethod("setRadius", "Transcriptogram",
    function(object, radius) {
        object@radius = transcriptogramerCheck("radius",
            radius)
        message("done!")
        object
    })

# orderingProperties ####

#' @rdname orderingProperties-method

setMethod("orderingProperties", "Transcriptogram",
    function(object, nCores = 1L) {
        if (is.na(object@status)) {
            stop("argument of class Transcriptogram - needs preprocessing!")
        }
        nCores <- transcriptogramerCheck("nCores", nCores)
        message("calculating node properties... step 1 of 2")
        message("** this may take some time...")
        nodeDegrees <- table(object@association$p1)
        nodeDegrees <- as.data.frame(nodeDegrees)
        colnames(nodeDegrees)[2] <- "nodeDegree"
        ord <- merge(object@ordering, nodeDegrees,
            by.x = "Protein", by.y = "Var1")
        ord <- ord[order(ord$Position), ]
        rownames(ord) <- NULL
        rm(nodeDegrees)
        g <- igraph::graph.data.frame(d = object@association,
            directed = FALSE)
        ord$nodeTriangles <- vapply(ord$Protein,
            function(x) {
                as.integer(igraph::count_triangles(g,
                  vids = x))
            }, integer(1))
        rm(g)
        ord$nodeClustering <- mapply(function(x,
            y) {
            if (y < 2) {
                return(0)
            } else {
                return(x/(y * (y - 1)/2))
            }
        }, ord$nodeTriangles, ord$nodeDegree)
        min <- ord$Position[1]
        max <- ord$Position[nrow(ord)]
        ntasks <- nrow(ord)
        pb <- progress::progress_bar$new(format = "running [:bar] :percent elapsed time :elapsed",
            total = ntasks, clear = FALSE,
            width = 60)
        adjlist <- tapply(object@association$p1,
            object@association$p2, unique)
        message("applying sliding window and mounting resulting ",
            "data... step 2 of 2")
        message("** this may take some time...")
        cl <- snow::makeSOCKcluster(nCores)
        on.exit(snow::stopCluster(cl))
        doSNOW::registerDoSNOW(cl)
        progress <- function() {
            pb$tick()
        }
        opts <- list(progress = progress)
        i <- NULL
        data <- foreach::foreach(i = seq.int(1, ntasks),
            .combine = "rbind", .options.snow = opts) %dopar%
            {
                pos <- ord[i, "Position"]
                l1 <- pos - object@radius
                l2 <- pos + object@radius
                temp <- data.frame()
                if (l1 >= min && l2 <= max) {
                  temp <- ord[which(ord$Position >=
                    l1 & ord$Position <=
                    l2), -c(2, 4)]
                } else if (l1 < min) {
                  temp <- ord[which(ord$Position <=
                    l2), -c(2, 4)]
                  temp <- rbind(temp, ord[which(ord$Position >=
                    (max + 1 + (l1 - min))),
                    -c(2, 4)])
                } else if (l2 > max) {
                  temp <- ord[which(ord$Position >=
                    l1), -c(2, 4)]
                  temp <- rbind(temp, ord[which(ord$Position <=
                    (l2%%max + min - 1)),
                    -c(2, 4)])
                }
                resultConnectivity <- mean(temp$nodeDegree)
                resultClustering <- mean(temp$nodeClustering)
                linksWindow <- sum(vapply(temp$Protein,
                  function(x) {
                    sum(temp$Protein %in%
                      adjlist[[x]])
                  }, integer(1)))
                linksTotal <- sum(temp$nodeDegree)
                resultModularity <- linksWindow/linksTotal
                return(data.frame(windowConnectivity = resultConnectivity,
                  windowClustering = resultClustering,
                  windowModularity = resultModularity))
            }
        message("done!")
        return(cbind(ord, data))
    })

# connectivityProperties ####

#' @rdname connectivityProperties-method

setMethod("connectivityProperties", "Transcriptogram",
    function(object) {
        if (is.na(object@status)) {
            stop("argument of class Transcriptogram - needs preprocessing!")
        }
        message("calculating graph properties... step 1 of 2")
        nodes <- table(object@association$p1)
        nodes <- as.data.frame(nodes)
        colnames(nodes) <- c("protein", "nodeDegree")
        g <- igraph::graph.data.frame(d = object@association,
            directed = FALSE)
        nodes$nodeTriangles <- vapply(nodes$protein,
            function(x) {
                as.integer(igraph::count_triangles(g,
                  vids = x))
            }, integer(1))
        rm(g)
        nodes$nodeClustering <- mapply(function(x,
            y) {
            if (y < 2) {
                return(0)
            } else {
                return(x/(y * (y - 1)/2))
            }
        }, nodes$nodeTriangles, nodes$nodeDegree)
        adjlist <- tapply(object@association$p1,
            object@association$p2, unique)
        nodes$nodeAssortativity <- vapply(nodes$protein,
            function(x) {
                return(mean(nodes[which(nodes$protein %in%
                  adjlist[[x]]), "nodeDegree"]))
            }, numeric(1))
        rm(adjlist)
        message("mounting resulting data... step 2 of 2")
        result <- unique(nodes$nodeDegree)
        result <- sort(result)
        n <- nrow(nodes)
        aux <- vapply(result, function(x) {
            temp <- nodes[which(nodes$nodeDegree ==
                x), c("nodeClustering", "nodeAssortativity")]
            return(c(nrow(temp)/n, mean(temp[,
                "nodeAssortativity"]), mean(temp[,
                "nodeClustering"])))
        }, vector("numeric", 3))
        message("done!")
        return(result <- data.frame(k = result,
            pk = aux[1, ], ak = aux[2, ],
            ck = aux[3, ]))
    })

# transcriptogramS1 ####

#' @rdname transcriptogramStep1-method

setMethod("transcriptogramStep1", "Transcriptogram",
    function(object, expression, dictionary, nCores = 1L) {
        if (is.na(object@status)) {
            stop("argument of class Transcriptogram - needs preprocessing!")
        }
        nCores <- transcriptogramerCheck("nCores", nCores)
        dictionary <- transcriptogramerCheck("dictionary",
            dictionary)
        expression <- transcriptogramerCheck("expression",
            expression)
        if (!(any(rownames(expression) %in%
            unique(dictionary$identifier)))) {
            stop("arguments expression and dictionary - does not match!")
        }
        dictionary <- dictionary[which(dictionary$protein %in%
            object@ordering$Protein), ]
        singletons <- names(which(table(dictionary$identifier) ==
            1))
        dictionary <- dictionary[which(dictionary$identifier %in%
            singletons), ]
        rm(singletons)
        rownames(dictionary) <- NULL
        nsamples <- ncol(expression)
        samples <- expression
        samples$identifier <- rownames(samples)
        message("mapping identifiers to ENSEMBL Peptide ID... step 1 of 2")
        map <- merge(object@ordering, dictionary,
            by.x = "Protein", by.y = "protein")
        map <- merge(samples, map, by.x = "identifier",
            by.y = "identifier")
        map <- map[order(map$Position), ]
        e <- unique(map[, "Protein"])
        result <- data.frame()
        map <- map[, -1]
        col <- ncol(map)
        col <- c(col - 1, col)
        ntasks <- length(e)
        pb <- progress::progress_bar$new(format = "running [:bar] :percent elapsed time :elapsed",
            total = ntasks, clear = FALSE,
            width = 60)
        cl <- snow::makeSOCKcluster(nCores)
        on.exit(snow::stopCluster(cl))
        doSNOW::registerDoSNOW(cl)
        progress <- function() {
            pb$tick()
        }
        opts <- list(progress = progress)
        message("calculating average over all identifiers ",
            "related to the same protein... step 2 of 2")
        i <- NULL
        result <- foreach::foreach(i = seq.int(1, ntasks),
            .combine = "rbind", .options.snow = opts) %dopar%
            {
                temp <- map[which(map$Protein ==
                  e[i]), ]
                if (nrow(temp) > 1) {
                  temp[1, -col] <- t(apply(temp[,
                    -col], 2, mean))
                }
                return(temp[1, ])
            }
        rownames(result) <- NULL
        result <- result[, c(col, seq.int(1, nsamples))]
        message("done!")
        object@transcriptogramS1 = result
        object@status = 1L
        return(object)
    })

# transcriptogramS2 ####

#' @rdname transcriptogramStep2-method

setMethod("transcriptogramStep2", "Transcriptogram",
    function(object, nCores = 1L) {
        if (object@status < 1L) {
            stop("argument of class Transcriptogram - be sure ",
                "to call the method transcriptogramStep1() before this one!")
        }
        nCores <- transcriptogramerCheck("nCores", nCores)
        object@transcriptogramS1 = object@transcriptogramS1[order(object@transcriptogramS1$Position),
            ]
        min <- min(object@transcriptogramS1$Position)
        max <- max(object@transcriptogramS1$Position)
        ntasks <- nrow(object@transcriptogramS1)
        pb <- progress::progress_bar$new(format = "running [:bar] :percent elapsed time :elapsed",
            total = ntasks, clear = FALSE,
            width = 60)
        result <- data.frame()
        cl <- snow::makeSOCKcluster(nCores)
        on.exit(snow::stopCluster(cl))
        doSNOW::registerDoSNOW(cl)
        progress <- function() {
            pb$tick()
        }
        opts <- list(progress = progress)
        message("applying sliding window and mounting resulting ",
            "data... step 1 of 1")
        col <- c(1, 2)
        i <- NULL
        result <- foreach::foreach(i = seq.int(1, ntasks),
            .combine = "rbind", .options.snow = opts) %dopar%
            {
                pos <- object@transcriptogramS1[i,
                  "Position"]
                l1 <- pos - object@radius
                l2 <- pos + object@radius
                if (l1 >= min && l2 <= max) {
                  temp <- object@transcriptogramS1[which(object@transcriptogramS1$Position >=
                    l1 & object@transcriptogramS1$Position <=
                    l2), -col]
                } else if (l1 < min) {
                  temp <- object@transcriptogramS1[which(object@transcriptogramS1$Position <=
                    l2), -col]
                  temp <- rbind(temp, object@transcriptogramS1[which(object@transcriptogramS1$Position >=
                    (max + 1 + (l1 - min))),
                    -col])
                } else if (l2 > max) {
                  temp <- object@transcriptogramS1[which(object@transcriptogramS1$Position >=
                    l1), -col]
                  temp <- rbind(temp, object@transcriptogramS1[which(object@transcriptogramS1$Position <=
                    (l2%%max + min - 1)),
                    -col])
                }
                temp[1, ] <- t(apply(temp,
                  2, mean))
                return(temp[1, ])
            }
        rownames(result) <- NULL
        nsamples <- ncol(result)
        result$Protein <- object@transcriptogramS1$Protein
        result$Position <- object@transcriptogramS1$Position
        result <- result[, c(nsamples + 1,
            nsamples + 2, seq.int(1, nsamples))]
        message("done!")
        object@transcriptogramS2 = result
        object@status = 2L
        return(object)
    })

# differentiallyExpressed ####

#' @rdname differentiallyExpressed-method

setMethod("differentiallyExpressed", "Transcriptogram", function(object,
    levels, pValue = 0.05, species = NULL, adjustMethod = "BH") {
    if (object@status < 2L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the methods transcriptogramStep1() and ",
            "transcriptogramStep2() before this one!")
    }
    transcriptogramerCheck("pValue", pValue)
    aux <- species
    if (is.data.frame(aux)) {
        species <- transcriptogramerCheck("species1", species)
    } else {
        transcriptogramerCheck("species1", species)
    }
    rm(aux)
    transcriptogramerCheck("adjustMethod1", adjustMethod)
    transcriptogramerCheck("levels", levels)
    if (length(levels) != (ncol(object@transcriptogramS2) -
        2)) {
        stop("argument levels - does not have a valid length!")
    }
    levels <- as.factor(levels)
    design <- stats::model.matrix(~0 + levels)
    contrasts <- "levelsFALSE-levelsTRUE"
    fit <- limma::lmFit(as.matrix(object@transcriptogramS2[,
        -c(1, 2)]), design)
    fit$Protein <- object@transcriptogramS2[, 1]
    fit$Position <- object@transcriptogramS2[, 2]
    contrasts <- limma::makeContrasts(contrasts = contrasts,
        levels = design)
    rm(design)
    message("calculating statistics... step 1 of 3")
    ct.fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts))
    res.fit <- limma::decideTests(ct.fit, method = "global",
        adjust.method = adjustMethod, p.value = pValue)
    temp <- data.frame(Protein = ct.fit$Protein, Position = ct.fit$Position,
        logFC = ct.fit$coef, pValue = ct.fit$p.value,
        degenes = unclass(res.fit), stringsAsFactors = FALSE)
    rm(contrasts)
    features <- rowSums(res.fit != 0) > 0
    DElimma <- temp[features, ]
    if (nrow(DElimma) == 0) {
        stop("no differentially expressed protein, ",
            "meeting the p-value requirement, was detected!")
    }
    rm(temp)
    colnames(DElimma)[c(3, 4, 5)] <- c("logFC", "pValue", "degenes")
    rownames(DElimma) <- NULL
    message("identifying clusters... step 2 of 3")
    pBreaks <- list()
    positions <- DElimma$Position
    clusterStartIndex <- clusterNumber <- 1
    nextIndex <- NULL
    invisible(sapply(seq.int(1, (length(positions) - 1)), function(i) {
        nextIndex <<- i + 1
        if ((positions[nextIndex] - positions[i]) > object@radius) {
            pBreaks[[clusterNumber]] <<- c(positions[clusterStartIndex],
                positions[i])
            clusterStartIndex <<- nextIndex
            clusterNumber <<- clusterNumber + 1
        }
        return(NULL)
    }))
    pBreaks[[clusterNumber]] <- c(positions[clusterStartIndex],
        positions[nextIndex])
    rm(nextIndex, clusterNumber, clusterStartIndex, positions)
    DElimma$ClusterNumber <- NA
    invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
        DElimma[which(DElimma$Position >= pBreaks[[i]][1] & DElimma$Position <=
            pBreaks[[i]][2]), "ClusterNumber"] <<- i
        return(NULL)
    }))
    DElimma <- DElimma[, c(1, 2, 6, 3, 4, 5)]
    message("generating plot... step 3 of 3")
    case <- object@transcriptogramS2[, -c(1, 2)]
    control <- case[, which(levels == TRUE)]
    case <- case[, which(levels == FALSE)]
    n <- nrow(control)
    caseValues <- vapply(seq.int(1, n), function(i) {
        result <- mean(unlist(case[i, ])) - mean(unlist(control[i,
            ]))
        return(result)
    }, numeric(1))
    smoothedLine <- stats::smooth.spline(object@transcriptogramS2$Position,
        caseValues, spar = 0.35)
    lim <- max(abs(min(caseValues)), abs(max(caseValues)))
    rm(case, control, n, caseValues)
    graphics::plot(smoothedLine, type = "l",
        ylab = "Difference of means (case - control)",
        xlab = "Gene position", main = "Differential expression",
        col = "black", lwd = 2, ylim = c(-lim, lim))
    graphics::grid(NULL, NULL, lwd = 1, lty = 1, col = "gray")
    graphics::abline(h = 0, col = "blue", lwd = 2)
    myColors <- grDevices::rainbow(length(pBreaks))
    invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
        idx <- which(smoothedLine$x >= pBreaks[[i]][1] & smoothedLine$x <=
            pBreaks[[i]][2])
        graphics::lines(x = smoothedLine$x[idx], y = smoothedLine$y[idx],
            lwd = 4, col = myColors[i])
        return(NULL)
    }))
    graphics::legend(x = "topright", legend = c("Case", "Control"),
        bty = "n", col = c("black", "blue"), lwd = 2, xpd = TRUE,
        inset = c(0, -0.15))
    if (!is.null(species)) {
        symbols <- NULL
        message("translating ENSEMBL Peptide ID to SYMBOL... extra step")
        taxonomyID <- NULL
        if (grepl("\\.", DElimma[1, 1])) {
            taxonomyID <- strsplit(DElimma[1, 1], "\\.")[[1]][1]
            taxonomyID <- paste0(taxonomyID, ".")
        }
        if (is.character(species)) {
            message("** this may take some time...")
            species <- tolower(species)
            species <- gsub("^([[:alpha:]]).* ", "\\1", species)
            species <- paste0(species, "_gene_ensembl")
            ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                dataset = species)
            proteins <- NULL
            if (grepl("\\.", DElimma[1, 1])) {
                proteins <- sapply(strsplit(DElimma[, 1], "\\."),
                    "[", 2)
            } else {
                proteins <- DElimma[, 1]
            }
            symbols <- biomaRt::getBM(filters = "ensembl_peptide_id",
                attributes = c("ensembl_peptide_id", "external_gene_name"),
                values = proteins, mart = ensembl)
        } else if (is.data.frame(species)) {
            symbols <- species
            if (grepl("\\.", symbols[1, 1])) {
                symbols$ensembl_peptide_id <- sapply(strsplit(symbols[,
                  1], "\\."), "[", 2)
            }
        }
        symbols[symbols == ""] <- NA
        symbols <- stats::na.omit(symbols)
        DElimma$Symbol <- DElimma$Protein
        invisible(sapply(seq.int(1, nrow(symbols)), function(i) {
            DElimma[which(DElimma$Protein == paste0(taxonomyID,
                symbols[i, "ensembl_peptide_id"])), "Symbol"] <<- symbols[i,
                "external_gene_name"]
            return(NULL)
        }))
    }
    object@status = 3L
    object@DE = DElimma
    message("done!")
    return(object)
})

# clusterVisualization ####

#' @rdname clusterVisualization-method

setMethod("clusterVisualization", "Transcriptogram",
    function(object,
    maincomp = FALSE, connected = FALSE,
    host = "127.0.0.1", port = 9091) {
    if (object@status < 3L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the method differentiallyExpressed() ",
            "before this one!")
    }
    symbolAsNodeAlias <- FALSE
    transcriptogramerCheck("maincomp",
        maincomp)
    transcriptogramerCheck("connected",
        connected)
    transcriptogramerCheck("host", host)
    transcriptogramerCheck("port", port)
    if ("Symbol" %in% colnames(object@DE)) {
        symbolAsNodeAlias <- TRUE
    }
    message("invoking RedeR... step 1 of 4")
    message("** this may take some time...")
    rdp <- RedeR::RedPort()
    RedeR::calld(rdp)
    message("generating the graphs... step 2 of 4")
    g <- igraph::graph.data.frame(d = object@association,
        directed = FALSE)
    n <- length(unique(object@DE$ClusterNumber))
    myColors <- grDevices::rainbow(n)
    sgList <- lapply(seq.int(1, n), function(i) {
        RedeR::subg(g = g, dat = object@DE[
            which(object@DE$ClusterNumber ==
            i), ], refcol = 1, maincomp = maincomp,
            connected = connected, transdat = TRUE)
    })
    rm(g)
    dim <- ceiling(sqrt(n))
    slice <- 100/dim
    myTheme <- list(isNest = TRUE, theme = 3,
        gscale = slice, nestFontSize = 50,
        zoom = 40)
    x <- y <- 0
    message("adding graphs into RedeR... step 3 of 4")
    if (n > 3) {
        message("** this may take some time...")
    }
    invisible(sapply(seq.int(1, n), function(i) {
        sgList[[i]] <<- RedeR::att.setv(g = sgList[[i]],
            cols = myColors[i])
        if (symbolAsNodeAlias) {
            sgList[[i]] <<- RedeR::att.setv(g = sgList[[i]],
              from = "Symbol", to = "nodeAlias")
        }
        message("** adding cluster ", i, " of ", n, "...")
        suppressMessages(RedeR::addGraph(rdp,
            sgList[[i]], theme = c(myTheme,
              nestAlias = paste0("C",
                i)), gcoord = c(x * slice,
              y * slice)))
        x <<- x + 1
        if (x == dim) {
            x <<- 0
            y <<- y + 1
        }
        return(NULL)
    }))
    RedeR::selectNodes(rdp, NULL)
    message("relaxing nodes... step 4 of 4")
    RedeR::relax(rdp)
    message("done!")
    return(rdp)
})

# clusterEnrichment ####

#' @rdname clusterEnrichment-method

setMethod("clusterEnrichment", "Transcriptogram", function(object,
    universe = NULL, species, ontology = "biological process",
    algorithm = "classic", statistic = "fisher", pValue = 0.05,
    adjustMethod = "BH", nCores = 1L) {
    if (object@status < 3L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the method differentiallyExpressed() before this one!")
    }
    nCores <- transcriptogramerCheck("nCores", nCores)
    transcriptogramerCheck("universe", universe)
    transcriptogramerCheck("pValue", pValue)
    aux <- species
    if (is.data.frame(species)) {
        species <- transcriptogramerCheck("species2", species)
    } else {
        transcriptogramerCheck("species2", species)
    }
    rm(aux)
    transcriptogramerCheck("statistic", statistic)
    transcriptogramerCheck("algorithm", algorithm)
    transcriptogramerCheck("ontology", ontology)
    transcriptogramerCheck("adjustMethod2", adjustMethod)
    if (is.null(universe)) {
        universe <- object@transcriptogramS2$Protein
    }
    if (!any(object@DE[, "Protein"] %in% universe)) {
        stop("argument of class Transcriptogram - none of ",
            "the Proteins of the DE slot are present in the argument universe!")
    }
    message("getting the terms... step 1 of 2")
    ontology <- tolower(ontology)
    ontology <- gsub(" ", "_", ontology)
    GO <- NULL
    if (grepl("\\.", universe[1])) {
        universe <- sapply(strsplit(universe, "\\."), "[", 2)
    }
    if (is.character(species)) {
        message("** this may take some time...")
        species <- tolower(species)
        species <- gsub("^([[:alpha:]]).* ", "\\1", species)
        species <- paste0(species, "_gene_ensembl")
        ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = species)
        GO <- biomaRt::getBM(filters = "ensembl_peptide_id",
            attributes = c("ensembl_peptide_id", "go_id", "namespace_1003"),
            values = universe, mart = ensembl)
        GO[GO == ""] <- NA
        GO <- stats::na.omit(GO)
        GO <- GO[which(GO$namespace_1003 == ontology), c(1, 2)]
        rm(ensembl)
    } else if (is.data.frame(species)) {
        GO <- species
        if (grepl("\\.", GO[1, 1])) {
            GO$ensembl_peptide_id <- sapply(strsplit(GO[, 1],
                "\\."), "[", 2)
        }
    }
    gene2GO <- split(GO$go_id, GO$ensembl_peptide_id)
    gene2GO <- lapply(gene2GO, unique)
    rm(species, GO)
    n <- length(unique(object@DE$ClusterNumber))
    message("running topGO enrichment for each cluster... step 2 of 2")
    message("** this may take some time...")
    ontology <- toupper(gsub("^([[:alpha:]]).*\\_([[:alpha:]]).*$",
        "\\1\\2", ontology))
    cl <- snow::makeSOCKcluster(nCores)
    on.exit(snow::stopCluster(cl))
    enrichment <- snow::parLapply(cl, seq.int(1, n), function(i){
        e <- environment()
        suppressMessages(topGO::groupGOTerms(e))
        attach(e)
        on.exit(detach(e))
        genesOfInterest <- object@DE[which(object@DE$ClusterNumber ==
            i), 1]
        if (grepl("\\.", genesOfInterest[1])) {
            genesOfInterest <- sapply(strsplit(genesOfInterest,
                "\\."), "[", 2)
        }
        geneList <- factor(as.integer(universe %in% genesOfInterest))
        names(geneList) <- universe
        myGOdata <- suppressMessages(methods::new("topGOdata",
            ontology = ontology, allGenes = geneList,
            annot = topGO::annFUN.gene2GO,
            gene2GO = gene2GO))
        result <- topGO::runTest(myGOdata, algorithm = algorithm,
            statistic = statistic)
        result <- topGO::GenTable(myGOdata, result,
            topNodes = length(result@score))
        colnames(result)[6] <- "pValue"
        result <- result[which(suppressWarnings(stats::p.adjust(result[,
            "pValue"], method = adjustMethod)) <= pValue),
            ]
        if (nrow(result) == 0) {
            return(NULL)
        }
        rownames(result) <- NULL
        return(result)
    })
    df <- data.frame()
    invisible(sapply(seq.int(1, n), function(i){
        if(!is.null(enrichment[[i]])){
            enrichment[[i]]$ClusterNumber <- i
            df <<- rbind(df, enrichment[[i]])
        }
        return (NULL)
    }))
    rm(enrichment)
    message("done!")
    return(df)
})

# getSlot ####

#' @rdname getSlot-method

setMethod("getSlot", "Transcriptogram", function(object,
    slot) {
    opts <- c("DE", "association", "ordering",
        "transcriptogramS1", "transcriptogramS2",
        "radius", "status")
    if (!is.character(slot) || length(slot) !=
        1 || !(slot %in% opts)) {
        stop("argument slot - should be any one of the options:\n",
            paste0(opts, collapse = ", "), "!")
    }
    if (slot == "DE") {
        return(object@DE)
    } else if (slot == "association") {
        return(object@association)
    } else if (slot == "ordering") {
        return(object@ordering)
    } else if (slot == "radius") {
        return(object@radius)
    } else if (slot == "status") {
        return(object@status)
    } else if (slot == "transcriptogramS1") {
        return(object@transcriptogramS1)
    } else if (slot == "transcriptogramS2") {
        return(object@transcriptogramS2)
    }
})
