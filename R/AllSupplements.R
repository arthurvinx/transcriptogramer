error <- function(name) {
    stop(paste0("argument ", name, " - does not have a valid value!"),
        call. = FALSE)
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
