check_ordering <- function(argument){
    name <- "ordering"
    if (is.character(argument) && length(argument) ==
        1 && file.exists(argument)) {
        argument <- data.table::fread(argument, header = TRUE,
            stringsAsFactors = FALSE)
        argument <- as.data.frame(argument)
        colnames(argument) <- c("Protein", "Position")
        if (orderingCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.character(argument) &&
        length(argument) > 1) {
        argument <- data.frame(Protein = argument,
            Position = seq.int(0, length(argument) - 1),
            stringsAsFactors = FALSE)
        if (orderingCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.data.frame(argument)) {
        if (data.table::is.data.table(argument)) {
            argument <- as.data.frame(argument)
        }
        colnames(argument) <- c("Protein", "Position")
        if (orderingCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else {
        error(name)
    }
}

check_association <- function(argument){
    name <- "association"
    if (is.character(argument) && length(argument) ==
        1 && file.exists(argument)) {
        argument <- data.table::fread(argument, header = FALSE,
            stringsAsFactors = FALSE)
        argument <- as.data.frame(argument)
        colnames(argument) <- c("p1", "p2")
        if (associationCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.matrix(argument)) {
        argument <- as.data.frame(argument,
            stringsAsFactors = FALSE)
        colnames(argument) <- c("p1", "p2")
        if (associationCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.data.frame(argument)) {
        if (data.table::is.data.table(argument)) {
            argument <- as.data.frame(argument)
        }
        colnames(argument) <- c("p1", "p2")
        if (associationCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else {
        error(name)
    }
}

check_expression <- function(argument){
    name <- "expression"
    if (is.matrix(argument)) {
        argument <- as.data.frame(argument)
        if (is.numeric(argument[1, 1])) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.data.frame(argument) &&
        is.numeric(argument[1, 1])) {
        if (data.table::is.data.table(argument)) {
            argument <- as.data.frame(argument)
        }
        return(argument)
    } else {
        error(name)
    }
}

check_dictionary <- function(argument){
    name <- "dictionary"
    if (is.matrix(argument)) {
        argument <- as.data.frame(argument,
            stringsAsFactors = FALSE)
        colnames(argument) <- c("protein", "identifier")
        if (dictionaryCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else if (is.data.frame(argument)) {
        if (data.table::is.data.table(argument)) {
            argument <- as.data.frame(argument)
        }
        colnames(argument) <- c("protein", "identifier")
        if (dictionaryCheck(argument)) {
            return(argument)
        } else {
            error(name)
        }
    } else {
        error(name)
    }
}

check_radius <- function(argument){
    name <- "radius"
    if (is.numeric(argument)) {
        argument <- as.integer(argument)
    }
    if (!is.integer(argument) || (argument < 0) || length(argument) != 1) {
        error(name)
    } else {
        return(argument)
    }
}

check_pValue <- function(argument){
    name <- "pValue"
    if (!is.numeric(argument) || (argument < 0 || argument > 1) ||
        length(argument) != 1) {
            error(name)
    }
}

check_species1 <- function(argument){
    if(is.data.frame(argument)){
        colnames(argument) <- c("ensembl_peptide_id", "external_gene_name")
        if(species1Check(argument)){
            return (argument)
        }else{
            stop("argument species - does not have a valid value!",
                call. = FALSE)
        }
    } else if (!(is.null(argument) || (is.character(argument) &&
        length(strsplit(argument, " ")[[1]]) ==
        2 && length(argument) ==
        1))) {
            stop("argument species - does not have a valid value!",
                call. = FALSE)
    }
}

check_adjustMethod1 <- function(argument){
     opts <- c("none", "BH", "fdr", "BY", "holm")
    if (!is.character(argument) || length(argument) != 1 ||
        !(argument %in% opts)) {
            stop("argument adjustMethod - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
    }
}

check_nCores <- function(argument){
    name <- "nCores"
    nc <- parallel::detectCores()
    if(is.logical(argument)){
        if(TRUE){
            return (nc)
        }else{
            return (1L)
        }
    }else if(is.numeric(argument)){
        argument <- as.integer(argument)
        if(argument < 1 || argument > nc){
            stop("argument nCores - should be greater than 0 and ",
                "less than or equal to ", nc)
        }else{
            return (argument)
        }
    }else{
        error(name)
    }
}

check_levels <- function(argument){
    name <- "levels"
    if (!is.logical(argument)) {
        error(name)
    }
}

check_maincomp <- function(argument){
    name <- "maincomp"
    if (!is.logical(argument) || length(argument) != 1) {
        error(name)
    }
}

check_connected <- function(argument){
    name <- "connected"
    if (!is.logical(argument) || length(argument) != 1) {
        error(name)
    }
}

check_trend <- function(argument){
  name <- "trend"
  if (!is.logical(argument) || length(argument) != 1) {
    error(name)
  }
}

check_hideLegend <- function(argument){
  name <- "hideLegend"
  if (!is.logical(argument) || length(argument) != 1) {
    error(name)
  }
}

check_port <- function(argument){
    name <- "port"
    if (is.numeric(argument)) {
        argument <- as.integer(argument)
    }
    if (!is.integer(argument) || (argument <
        0 || argument > 65535) || length(argument) != 1) {
            error(name)
    }
}

check_host <- function(argument){
    name <- "host"
    if (!is.character(argument) || length(argument) != 1) {
        error(name)
    }
}

check_adjustMethod2 <- function(argument){
    opts <- c("none", "BH", "fdr", "BY", "hochberg", "hommel", "bonferroni",
        "holm")
    if (!is.character(argument) || length(argument) !=
        1 || !(argument %in% opts)) {
        stop("argument adjustMethod: should be any one of the options:\n",
            paste0(opts, collapse = ", "), "!", call. = FALSE)
    }
}

check_ontology <- function(argument){
    opts <- c("biological process", "cellular component", "molecular function")
    if (!is.character(argument) || length(argument) != 1 ||
        !(argument %in% opts) || length(strsplit(argument, " ")[[1]]) != 2) {
            stop("argument ontology - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
}

check_algorithm <- function(argument){
    opts <- c("classic", "weight01", "lea", "parentchild", "elim", "weight")
    if (!is.character(argument) || length(argument) !=
        1 || !(argument %in% opts)) {
        stop("argument algorithm - should be any one of the options:\n",
            paste0(opts, collapse = ", "), "!", call. = FALSE)
    }
}

check_statistic <- function(argument){
    opts <- c("fisher", "ks", "t", "sum", "globaltest")
    if (!is.character(argument) || length(argument) != 1 ||
        !(argument %in% opts)) {
            stop("argument statistic - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
    }
}

check_species2 <- function(argument){
    if(is.data.frame(argument)){
        colnames(argument) <- c("ensembl_peptide_id", "go_id")
        if(species2Check(argument)){
            return (argument)
        }else{
            stop("argument species - does not have a valid value!",
                call. = FALSE)
        }
    } else if (!(is.character(argument) && length(strsplit(argument,
        " ")[[1]]) == 2) || length(argument) != 1) {
            stop("argument species - does not have a valid value!",
                call. = FALSE)
    }
}

check_universe <- function(argument){
    name <- "universe"
    if (!(is.character(argument) || is.null(argument))) {
        error(name)
    }
}

orderingCheck <- function(argument) {
    if (ncol(argument) == 2 && is.character(argument$Protein) &&
        is.integer(argument$Position)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

associationCheck <- function(argument) {
    if (ncol(argument) == 2 && is.character(argument$p1) &&
        is.character(argument$p2)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

dictionaryCheck <- function(argument) {
    if (ncol(argument) == 2 && is.character(argument$protein) &&
        is.character(argument$identifier)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

species1Check <- function(argument) {
    if (ncol(argument) == 2 && is.character(argument$ensembl_peptide_id) &&
        is.character(argument$external_gene_name)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

species2Check <- function(argument) {
    if (ncol(argument) == 2 && is.character(argument$ensembl_peptide_id) &&
        is.character(argument$go_id)) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

error <- function(name) {
    stop("argument ", name, " - does not have a valid value!", call. = FALSE)
}
