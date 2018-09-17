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

# radius<- ####

#' @rdname radius-method
#' @aliases radius-method

setReplaceMethod("radius", "Transcriptogram",
    function(object, value) {
        value <- check_radius(value)
        object@radius <- value
        object
    }
)

# orderingProperties ####

#' @rdname orderingProperties-method
#' @aliases orderingProperties-method

setMethod("orderingProperties", "Transcriptogram",
    function(object, nCores = 1L) {
        if (is.na(object@status)) {
            stop("argument of class Transcriptogram - needs preprocessing!")
        }
        nCores <- check_nCores(nCores)
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
                    (max + 1 + l1 - min)),
                    -c(2, 4)])
                } else if (l2 > max) {
                  temp <- ord[which(ord$Position >=
                    l1), -c(2, 4)]
                  temp <- rbind(temp, ord[which(ord$Position <=
                    (l2 %% max + min - 1)),
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
#' @aliases connectivityProperties-method

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
#' @aliases transcriptogramStep1-method

setMethod("transcriptogramStep1", "Transcriptogram",
    function(object, expression, dictionary, nCores = 1L) {
        if (is.na(object@status)) {
            stop("argument of class Transcriptogram - needs preprocessing!")
        }
        nCores <- check_nCores(nCores)
        dictionary <- check_dictionary(dictionary)
        expression <- check_expression(expression)
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
        message("averaging over all identifiers ",
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
        object@Protein2Symbol = data.frame(ensembl_peptide_id = object@ordering$Protein,
                                           external_gene_name = object@ordering$Protein,
                                           stringsAsFactors = FALSE)
        return(object)
    })

# transcriptogramS2 ####

#' @rdname transcriptogramStep2-method
#' @aliases transcriptogramStep2-method

setMethod("transcriptogramStep2", "Transcriptogram",
    function(object, nCores = 1L) {
        if (object@status < 1L) {
            stop("argument of class Transcriptogram - be sure ",
                "to call the method transcriptogramStep1() before this one!")
        }
        nCores <- check_nCores(nCores)
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
                    (max + 1 + l1 - min)),
                    -col])
                } else if (l2 > max) {
                  temp <- object@transcriptogramS1[which(object@transcriptogramS1$Position >=
                    l1), -col]
                  temp <- rbind(temp, object@transcriptogramS1[which(object@transcriptogramS1$Position <=
                    (l2 %% max + min - 1)),
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
#' @aliases differentiallyExpressed-method

setMethod("differentiallyExpressed", "Transcriptogram", function(object,
    levels, pValue = 0.05, species = object@Protein2Symbol, adjustMethod = "BH",
    trend = FALSE, title = "Differential expression",
    boundaryConditions = FALSE) {
    if (object@status < 2L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the methods transcriptogramStep1() and ",
            "transcriptogramStep2() before this one!")
    }
    check_pValue(pValue)
    check_title(title)
    aux <- species
    if (is.data.frame(aux)) {
      species <- check_species1(species)
    } else {
      check_species1(species)
    }
    rm(aux)
    check_adjustMethod1(adjustMethod)
    check_trend(trend)
    check_boundaryConditions(boundaryConditions)
    check_levels(levels)
    if (length(levels) != (ncol(object@transcriptogramS2) -
        2)) {
        stop("argument levels - does not have a valid length!")
    }
    object@pbc = FALSE
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
    message("calculating statistics... step 1 of 4")
    ct.fit <- limma::eBayes(limma::contrasts.fit(fit, contrasts), trend = trend)
    res.fit <- limma::decideTests(ct.fit, method = "global",
        adjust.method = adjustMethod, p.value = pValue)
    temp <- data.frame(Protein = ct.fit$Protein, Position = ct.fit$Position,
        logFC = ct.fit$coef, pValue = ct.fit$p.value,
        degenes = as.integer(unclass(res.fit)), stringsAsFactors = FALSE)
    rm(contrasts)
    features <- rowSums(res.fit != 0) > 0
    DElimma <- temp[features, ]
    if (nrow(DElimma) == 0) {
        stop("no differentially expressed protein, ",
            "meeting the p-value requirement, was detected!")
    }
    rm(temp)
    colnames(DElimma)[c(3, 4, 5)] <- c("logFC", "pValue", "DEgenes")
    rownames(DElimma) <- NULL
    message("identifying clusters... step 2 of 4")
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
    if(boundaryConditions){
      aux <- c()
      min <- object@ordering$Position[1]
      max <- object@ordering$Position[nrow(object@ordering)]
      aux <- invisible(lapply(seq.int(1, length(pBreaks)), function(i) {
        l1 <- pBreaks[[i]][1] - object@radius
        l2 <- pBreaks[[i]][2] + object@radius
        return(c(l1,l2))
      }))
      elim <- list()
      invisible(lapply(seq.int(1, length(aux)), function(i) {
        if(i == length(aux)){
          elim <<- append(elim, list(c(aux[[i]][1], aux[[i]][2])))
        }else if(aux[[i]][2] >= aux[[i + 1]][1]){
          aux[[i + 1]] <<- c(aux[[i]][1], aux[[i + 1]][2])
        }else{
          elim <<- append(elim, list(c(aux[[i]][1], aux[[i]][2])))
        }
        return(NULL)
      }))
      aux <- elim
      if(aux[[1]][1] < min){
        object@pbc = TRUE
        x <- max + 1 + aux[[1]][1] - min
        if(x <= aux[[length(aux)]][2]){
          aux[[1]][1] <- min
          aux[[length(aux)]][2] <- max
        }else{
          aux[[1]][1] <- min
          aux <- append(aux, list(c(x, max)))
        }
      }else if(aux[[length(aux)]][2] > max){
        object@pbc = TRUE
        x <- aux[[length(aux)]][2] %% max + min - 1
        if(x >= aux[[1]][1]){
          aux[[1]][1] <- min
          aux[[length(aux)]][2] <- max
        }else{
          aux[[length(aux)]][2] <- max
          aux <- append(list(c(min, x)), aux)
        }
      }
      pBreaks <- aux
    }
    DElimma$ClusterNumber <- NA
    invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
      if(i == length(pBreaks) && object@pbc){
        DElimma[which(DElimma$Position >= pBreaks[[i]][1] & DElimma$Position <=
                        pBreaks[[i]][2]), "ClusterNumber"] <<- 1
      }else{
        DElimma[which(DElimma$Position >= pBreaks[[i]][1] & DElimma$Position <=
                        pBreaks[[i]][2]), "ClusterNumber"] <<- i
      }
        return(NULL)
    }))
    DElimma <- DElimma[, c(1, 2, 6, 3, 4, 5)]
    message("generating plot... step 3 of 4")
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
    lim <- max(abs(caseValues))
    rm(case, control, n, caseValues)
    lim <- round(lim, digits = 1)
    if(object@pbc){
      myColors <- grDevices::rainbow(length(pBreaks) - 1)
      myColors <- c(myColors, myColors[1])
    }else{
      myColors <- grDevices::rainbow(length(pBreaks))
    }
    df <- data.frame(x = smoothedLine$x, y = smoothedLine$y)
    rm(smoothedLine)
    p <- ggplot2::ggplot(df, ggplot2::aes_string("x", "y")) +
      ggplot2::geom_line(lwd = 1, ggplot2::aes_string(y = "0", colour = '"c1"')) +
      ggplot2::geom_line(lwd = 1, ggplot2::aes_string(y = "y", colour = '"c2"')) +
      ggplot2::scale_y_continuous(limits = c(-lim, lim), breaks = round(seq(-lim, lim, 0.1), digits = 1)) +
      ggplot2::scale_x_continuous(limits = c(0, length(object@ordering$Position) - 1),
                                  breaks = seq.int(0, length(object@ordering$Position) - 1, 1000)) +
      ggplot2::scale_colour_manual(values = c("black", "grey80"), name = "Conditions",
                                   labels =  c("Control", "Case")) +
      ggplot2::scale_linetype_manual(values = "blank", name = "Number of clusters",
                                     labels = ifelse(object@pbc, length(myColors) - 1, length(myColors))) +
      ggplot2::labs(x = "Gene position", y = "Difference of means (case - control)", title = title) +
      ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    invisible(sapply(seq.int(1, length(pBreaks)), function(i) {
      idx <- which(df$x >= pBreaks[[i]][1] & df$x <= pBreaks[[i]][2])
      aux <- data.frame(x = df$x[idx], y = df$y[idx])
      p <<- p + ggplot2::geom_line(data = aux, lwd = 1.4, col = myColors[i],
                                   ggplot2::aes_string(x = "x", y = "y")) +
        ggplot2::geom_line(ggplot2::aes_string(linetype = '"lines"'))
      return(NULL)
    }))
    suppressMessages(graphics::plot(p))
    symbols <- NULL
    message("translating ENSEMBL Peptide ID to SYMBOL... step 3 of 4")
    taxonomyID <- NULL
    if (grepl("\\.", object@ordering[1, 1])) {
      taxonomyID <- strsplit(object@ordering[1, 1], "\\.")[[1]][1]
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
      if (grepl("\\.", object@ordering[1, 1])) {
        proteins <- sapply(strsplit(object@ordering[, 1], "\\."),
                           "[", 2)
      } else {
        proteins <- object@ordering[, 1]
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
    symbols$ensembl_peptide_id <- paste0(taxonomyID,
                                         symbols$ensembl_peptide_id)
    DElimma$Symbol <- NA_character_
    object@Protein2Symbol = symbols
    invisible(sapply(seq.int(1, nrow(DElimma)), function(i) {
      DElimma$Symbol <<- symbols[match(DElimma[, "Protein"],
                                       symbols[,"ensembl_peptide_id"]),
                                 "external_gene_name"]
      return(NULL)
    }))
    if(any(is.na(DElimma$Symbol))){
      idx <- which(is.na(DElimma$Symbol))
      DElimma[idx, "Symbol"] <- DElimma[idx, "Protein"]
    }
    object@status = 3L
    object@DE = DElimma
    object@clusters = pBreaks
    message("done!")
    return(object)
})

# clusterVisualization ####

#' @rdname clusterVisualization-method
#' @aliases clusterVisualization-method

setMethod("clusterVisualization", "Transcriptogram",
    function(object,
    maincomp = FALSE, connected = FALSE,
    host = "127.0.0.1", port = 9091, clusters = NULL, onlyGenesInDE = TRUE) {
    if (object@status < 3L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the method differentiallyExpressed() before this one!")
    }
    check_maincomp(maincomp)
    check_onlyGenesInDE(onlyGenesInDE)
    check_connected(connected)
    check_host(host)
    check_port(port)
    if(is.null(clusters)){
        clusters  <- unique(object@DE$ClusterNumber)
    }else{
        if(is.numeric(clusters)){
            clusters <- as.integer(unique(clusters))
        }
        if(!is.integer(clusters) || !all(clusters %in%
            unique(object@DE$ClusterNumber))){
            error("clusters")
        }
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
    sgList <- list()
    if(onlyGenesInDE){
      sgList <- lapply(seq.int(1, n), function(i) {
        RedeR::subg(g = g, dat = object@DE[
          which(object@DE$ClusterNumber ==
                  i), ], refcol = 1, maincomp = maincomp,
          connected = connected, transdat = TRUE)
      })
    }else{
      sgList <- lapply(seq.int(1, n), function(i) {
        positions <- c()
        if(object@pbc && (i == 1)){
          positions <- c(positions,
                         object@clusters[[1]][1]:object@clusters[[1]][2])
          positions <- c(positions,
                         object@clusters[[length(object@clusters)]][1]:
                           object@clusters[[length(object@clusters)]][2])

        }else{
          positions <- c(positions,
                         object@clusters[[i]][1]:object@clusters[[i]][2])
        }
        df <- object@ordering[object@ordering$Position %in% positions, 1]
        df <- as.data.frame(cbind("Protein" = df, "Symbol" = object@Protein2Symbol[match(df, object@Protein2Symbol$ensembl_peptide_id), 2]), stringsAsFactors = FALSE)
        RedeR::subg(g = g, dat = df,
          refcol = 1, maincomp = maincomp,
          connected = connected, transdat = TRUE)
      })
    }
    rm(g)
    message("adding graphs into RedeR... step 3 of 4")
    n <- length(clusters)
    if (n > 3) {
        message("** this may take some time...")
    }
    dim <- ceiling(sqrt(n))
    slice <- 100/dim
    myTheme <- list(isNest = TRUE, theme = 3, gscale = slice,
        nestFontSize = 50, zoom = 40)
    x <- y <- 0
    invisible(sapply(clusters, function(i) {
      sgList[[i]] <<- RedeR::att.setv(g = sgList[[i]],
                                      cols = myColors[i])
      igraph::E(sgList[[i]])$edgeColor <<- "grey80"
      igraph::V(sgList[[i]])$nodeLineColor <<- "grey80"
      sgList[[i]] <<- RedeR::att.setv(g = sgList[[i]],
                                      from = "Symbol", to = "nodeAlias")
      message("** adding cluster ", i, "...")
      suppressMessages(RedeR::addGraph(rdp, sgList[[i]],
                                       theme = c(myTheme, nestAlias = paste0("C", i)),
                                       gcoord = c(x * slice, y * slice)))
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
#' @aliases clusterEnrichment-method

setMethod("clusterEnrichment", "Transcriptogram", function(object,
    universe = NULL, species, ontology = "biological process",
    algorithm = "classic", statistic = "fisher", pValue = 0.05,
    adjustMethod = "BH", nCores = 1L, onlyGenesInDE = TRUE) {
    if (object@status < 3L) {
        stop("argument of class Transcriptogram - be sure to ",
            "call the method differentiallyExpressed() before this one!")
    }
    nCores <- check_nCores(nCores)
    check_universe(universe)
    check_pValue(pValue)
    check_onlyGenesInDE(onlyGenesInDE)
    aux <- species
    if (is.data.frame(species)) {
        species <- check_species2(species)
    } else {
        check_species2(species)
    }
    rm(aux)
    check_statistic(statistic)
    check_algorithm(algorithm)
    check_ontology(ontology)
    check_adjustMethod2(adjustMethod)
    if (is.null(universe)) {
        universe <- object@ordering$Protein
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
    object@Protein2GO = GO
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
    e <- environment()
    suppressMessages(topGO::groupGOTerms(e))
    enrichment <- snow::parLapply(cl, seq.int(1, n), function(i,
                                                              onlyGenesInDE,
                                                              object, universe, e){
      attach(e)
      on.exit(detach(e))
      genesOfInterest <- c()
      if(onlyGenesInDE){
        genesOfInterest <- object@DE[which(object@DE$ClusterNumber == i), 1]
      }else{
        positions <- c()
        if(object@pbc && (i == 1)){
          positions <- c(positions,
                         object@clusters[[1]][1]:object@clusters[[1]][2])
          positions <- c(positions,
                         object@clusters[[length(object@clusters)]][1]:
                           object@clusters[[length(object@clusters)]][2])

        }else{
          positions <- c(positions,
                         object@clusters[[i]][1]:object@clusters[[i]][2])
        }
        genesOfInterest <- object@ordering[object@ordering$Position %in% positions, 1]
      }
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
      if(any(TRUE %in% grepl("^<", result[, "pValue"]))){
        result$pValue <- gsub("^<", "", result$pValue)
      }
      result$pValue <- as.numeric(result$pValue)
      result$pValue <- stats::p.adjust(result[, "pValue"], method = adjustMethod)
      result <- result[result$pValue <= pValue, ]
      if (nrow(result) == 0) {
        return(NULL)
      }
      result$ClusterNumber <- i
      rownames(result) <- NULL
      if(i == 1){
        temporary <- list()
        temporary[[1]] <- result
        temporary[2] <- myGOdata
        return(temporary)
      }
      return(result)
    }, onlyGenesInDE, object, universe, e)
    temporary <- enrichment[[1]][[2]]
    enrichment[[1]] <- enrichment[[1]][[1]]
    enrichment <- do.call("rbind", enrichment)
    object@genesInTerm = topGO::genesInTerm(temporary, unique(enrichment$GO.ID))
    rm(temporary)
    object@Terms = enrichment
    object@status = 4L
    message("done!")
    return(object)
})

# enrichmentPlot ####

#' @rdname enrichmentPlot-method
#' @aliases enrichmentPlot-method

setMethod("enrichmentPlot", "Transcriptogram",
          function(object, nCores = 1L, nTerms = 1L,
                   GOIDs = NULL, title = "Enrichment") {
            nCores <- check_nCores(nCores)
            nTerms <- check_nTerms(nTerms)
            check_title(title)
            if (object@status < 4L) {
              stop("argument of class Transcriptogram - be sure to ",
                   "call the method clusterEnrichment() before this one!")
            }
            terms <- object@Terms[order(object@Terms$pValue),]
            if(is.null(GOIDs)){
              v <- sort(unique(terms$ClusterNumber))
              GOIDs <- lapply(v, function(i){
                stats::na.omit(utils::head(terms[terms$ClusterNumber==i, c(1, 2)], nTerms))
              })
              GOIDs <- do.call("rbind", GOIDs)
            }else{
              if(is.character(GOIDs)){
                GOIDs <- unique(GOIDs)
                GOIDs <- data.frame(GO.ID = GOIDs,
                                    Term = terms[match(GOIDs, terms$GO.ID), 2],
                                    stringsAsFactors = F)
              }else{
                stop("argument GOIDs - does not have a valid value!")
              }
            }
            v <- unique(GOIDs$GO.ID)
            ord <- object@ordering
            min <- ord$Position[1]
            max <- ord$Position[nrow(object@ordering)]
            ntasks <- nrow(ord)
            pb <- progress::progress_bar$new(format = "running [:bar] :percent elapsed time :elapsed",
                                             total = ntasks, clear = FALSE,
                                             width = 60)
            message("applying sliding window and mounting resulting ",
                    "data... step 1 of 2")
            message("** this may take some time...")
            cl <- snow::makeSOCKcluster(nCores)
            on.exit(snow::stopCluster(cl))
            doSNOW::registerDoSNOW(cl)
            progress <- function() {
              pb$tick()
            }
            GOmapping <- object@Protein2GO
            opts <- list(progress = progress)
            i <- NULL
            data <- foreach::foreach(i = seq.int(1, ntasks),
                                     .combine = "rbind", .options.snow = opts) %dopar%
                                     {
                                       pos <- i-1
                                       l1 <- pos - object@radius
                                       l2 <- pos + object@radius
                                       temp <- data.frame()
                                       if (l1 >= min && l2 <= max) {
                                         temp <- ord[which(ord$Position >=
                                                             l1 & ord$Position <=
                                                             l2),]
                                       } else if (l1 < min) {
                                         temp <- ord[which(ord$Position <=
                                                             l2),]
                                         temp <- rbind(temp, ord[which(ord$Position >=
                                                                         (max + 1 + l1 - min)),])
                                       } else if (l2 > max) {
                                         temp <- ord[which(ord$Position >=
                                                             l1),]
                                         temp <- rbind(temp, ord[which(ord$Position <=
                                                                         (l2 %% max + min - 1)),])
                                       }
                                       temp <- temp[, 1]
                                       if(grepl("\\.", temp[1])){
                                         taxonomyID <- strsplit(temp[1], "\\.")[[1]][1]
                                         temp <- gsub(paste0(taxonomyID, "."), "",
                                                      temp, fixed = TRUE)
                                       }
                                       n <- length(temp)
                                       aux <- GOmapping[GOmapping$ensembl_peptide_id %in% temp,]
                                       rates <- vapply(v, function(x) {
                                         nrow(aux[aux[, 2] == x,])/n
                                       }, numeric(1))
                                       return(data.frame(Position = i - 1, t(rates)))}
            invisible(sapply(seq.int(2, ncol(data)), function(i) {
              smoothedLine <- stats::smooth.spline(data[, 1], data[, i], spar = 0.35)
              data[, i] <<- smoothedLine$y
              return(NULL)
            }))
            colnames(data) <- c("Position", paste0(GOIDs[match(v, GOIDs$GO.ID), 2], " (", v, ")"))
            data <- data[, colSums(data) != 0]
            data <- tidyr::gather(data, "key", "value", -Position)
            colnames(data) <- c("x", "Terms", "y")
            message("generating plot... step 2 of 2")
            p <- ggplot2::ggplot(data, ggplot2::aes_string(x = "x", y = "y", colour = "Terms")) +
              ggplot2::geom_line() +
              ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
              ggplot2::scale_x_continuous(limits = c(0, length(ord$Position) - 1),
                                          breaks = seq.int(0, length(ord$Position) - 1, 1000)) +
              ggplot2::labs(x = "Gene position", y = "GO term rate by window", title = title) +
              ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
            message("done!")
            suppressWarnings(graphics::plot(p))
            return(p)
          })

# radius ####

#' @rdname radius-method
#' @aliases radius-method

setMethod("radius", "Transcriptogram",
    function(object) {
        object@radius
    })

# DE ####

#' @rdname DE-method
#' @aliases DE-method

setMethod("DE", "Transcriptogram",
    function(object) {
        object@DE
    })

# Terms ####

#' @rdname Terms-method
#' @aliases Terms-method

setMethod("Terms", "Transcriptogram",
          function(object) {
            object@Terms
          })

# show ####

setMethod("show", "Transcriptogram",
    function(object) {
        cat('Slot "association":\n')
        print(object@association)
        cat('\nSlot "ordering":\n')
        print(object@ordering)
        cat('\nSlot "transcriptogramS1":\n')
        print(object@transcriptogramS1)
        cat('\nSlot "transcriptogramS2":\n')
        print(object@transcriptogramS2)
        cat('\nSlot "DE":\n')
        print(object@DE)
        cat('\nSlot "radius":\n', object@radius)
        cat('\n\nSlot "status":\n', object@status, "\n")
    })
