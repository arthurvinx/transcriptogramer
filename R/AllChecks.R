transcriptogramer.check <- function(name,
    argument) {
    if (name == "ordering") {
        if (is.character(argument) && length(argument) ==
            1 && file.exists(argument)) {
            argument <- data.table::fread(argument)
            argument <- as.data.frame(argument)
            colnames(argument) <- c("Protein",
                "Position")
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
            colnames(argument) <- c("Protein",
                "Position")
            if (orderingCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else {
            error(name)
        }
    } else if (name == "association") {
        if (is.character(argument) && length(argument) ==
            1 && file.exists(argument)) {
            argument <- data.table::fread(argument,
                header = FALSE)
            argument <- as.data.frame(argument)
            colnames(argument) <- c("p1",
                "p2")
            if (associationCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else if (is.matrix(argument)) {
            argument <- as.data.frame(argument,
                stringsAsFactors = FALSE)
            colnames(argument) <- c("p1",
                "p2")
            if (associationCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else if (is.data.frame(argument)) {
            if (data.table::is.data.table(argument)) {
                argument <- as.data.frame(argument)
            }
            colnames(argument) <- c("p1",
                "p2")
            if (associationCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else {
            error(name)
        }
    } else if (name == "expression") {
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
    } else if (name == "dictionary") {
        if (is.matrix(argument)) {
            argument <- as.data.frame(argument,
                stringsAsFactors = FALSE)
            colnames(argument) <- c("protein",
                "identifier")
            if (dictionaryCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else if (is.data.frame(argument)) {
            if (data.table::is.data.table(argument)) {
                argument <- as.data.frame(argument)
            }
            colnames(argument) <- c("protein",
                "identifier")
            if (dictionaryCheck(argument)) {
                return(argument)
            } else {
                error(name)
            }
        } else {
            error(name)
        }
    } else if (name == "radius") {
        if (is.numeric(argument)) {
            argument <- as.integer(argument)
        }
        if (!is.integer(argument) || (argument <
            0) || length(argument) != 1) {
            error(name)
        } else {
            return(argument)
        }
    } else if (name == "pValue") {
        if (!is.numeric(argument) || (argument <
            0 || argument > 1) || length(argument) !=
            1) {
            error(name)
        }
    } else if (name == "species1") {
      if(is.data.frame(argument)){
        colnames(argument) <- c("ensembl_peptide_id",
                "external_gene_name")
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
    } else if (name == "adjustMethod1") {
        opts <- c("none", "BH", "fdr", "BY",
            "holm")
        if (!is.character(argument) || length(argument) !=
            1 || !(argument %in% opts)) {
            stop("argument adjustMethod - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
    } else if (name == "levels") {
        if (!is.logical(argument)) {
            error(name)
        }
    } else if (name %in% c("maincomp", "connected")) {
        if (!is.logical(argument) || length(argument) !=
            1) {
            error(name)
        }
    } else if (name == "port") {
        if (is.numeric(argument)) {
            argument <- as.integer(argument)
        }
        if (!is.integer(argument) || (argument <
            0 || argument > 65535) || length(argument) !=
            1) {
            error(name)
        }
    } else if (name == "host") {
        if (!is.character(argument) || length(argument) !=
            1) {
            error(name)
        }
    } else if (name == "adjustMethod2") {
        opts <- c("none", "BH", "fdr", "BY",
            "hochberg", "hommel", "bonferroni",
            "holm")
        if (!is.character(argument) || length(argument) !=
            1 || !(argument %in% opts)) {
            stop("argument adjustMethod: should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
    } else if (name == "ontology") {
        opts <- c("biological process", "cellular component",
            "molecular function")
        if (!is.character(argument) || length(argument) !=
            1 || !(argument %in% opts) ||
            length(strsplit(argument, " ")[[1]]) !=
                2) {
            stop("argument ", name, " - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
    } else if (name == "algorithm") {
        opts <- c("classic", "weight01",
            "lea", "parentchild", "elim",
            "weight")
        if (!is.character(argument) || length(argument) !=
            1 || !(argument %in% opts)) {
            stop("argument ", name, " - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
    } else if (name == "statistic") {
        opts <- c("fisher", "ks", "t", "sum",
            "globaltest")
        if (!is.character(argument) || length(argument) !=
            1 || !(argument %in% opts)) {
            stop("argument ", name, " - should be any one of the options:\n",
                paste0(opts, collapse = ", "), "!", call. = FALSE)
        }
    } else if (name == "species2") {
      if(is.data.frame(argument)){
        colnames(argument) <- c("ensembl_peptide_id",
                "go_id")
        if(species2Check(argument)){
          return (argument)
        }else{
          stop("argument species - does not have a valid value!",
                call. = FALSE)
        }
      } else if (!(is.character(argument) && length(strsplit(argument,
            " ")[[1]]) == 2) || length(argument) !=
            1) {
            stop("argument species - does not have a valid value!",
                call. = FALSE)
        }
    } else if (name == "universe") {
        if (!(is.character(argument) || is.null(argument))) {
            error(name)
        }
    } else if (name == "nCores") {
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
    return(NULL)
}
