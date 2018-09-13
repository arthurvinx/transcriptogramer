#' Class Transcriptogram
#'
#' This S4 class includes methods to use expression values with ordered
#' proteins.
#'
#' @slot association A data.frame containing two columns, with rows containing
#' ENSEMBL Peptide IDs that are connected.
#'
#' @slot ordering A data.frame containing two columns, the first one with
#' ENSEMBL Peptide IDs, and the second containing its respective position.
#'
#' @slot transcriptogramS1 A data.frame produced as the result of averaging
#' over all identifiers related to the same protein.
#'
#' @slot transcriptogramS2 A data.frame produced as the result of averaging
#' over the window.
#'
#' @slot radius An non-negative integer referring to the window radius.
#'
#' @slot status An integer used internally to check the status of the object.
#'
#' @slot DE A data.frame of differentially expressed proteins.
#'
#' @slot clusters A list indicating the first and the last position belonging to
#' each cluster.
#'
#' @slot pbc Logical value used internally to indicate the overlapping of the
#' first and the last cluster.
#'
#' @slot Protein2Symbol A data.frame containing two columns, the first one with
#' ENSEMBL Peptide IDs, and the second containing its respective Symbol.
#'
#' @slot Protein2GO A data.frame containing two columns, the first one with
#' ENSEMBL Peptide IDs, and the second containing its respective Gene Ontology accession.
#'
#' @slot Terms A data.frame containing the enriched Gene Ontology terms.
#'
#' @slot topGOdata A list of GO terms and its respective ENSEMBL Peptide IDs, feeded by
#' the clusterEnrichment() method.
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer:DE-method]{DE},
#' \link[transcriptogramer:radius-method]{radius},
#' \link[transcriptogramer:orderingProperties-method]{orderingProperties},
#' \link[transcriptogramer:connectivityProperties-method]{connectivityProperties},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2},
#' \link[transcriptogramer:differentiallyExpressed-method]{differentiallyExpressed},
#' \link[transcriptogramer:clusterVisualization-method]{clusterVisualization},
#' \link[transcriptogramer:clusterEnrichment-method]{clusterEnrichment},
#' \link[transcriptogramer:enrichmentPlot-method]{enrichmentPlot}
#'
#'
#' @author
#' Diego Morais
#'
#' @export

setClass("Transcriptogram", representation(association = "data.frame",
    ordering = "data.frame", transcriptogramS1 = "data.frame",
    transcriptogramS2 = "data.frame", DE = "data.frame",
    radius = "integer", status = "integer", clusters = "list", pbc = "logical",
    Protein2Symbol = "data.frame", Protein2GO = "data.frame",
    Terms = "data.frame", topGOdata = "list"),
    prototype = list(association = data.frame(),
        ordering = data.frame(), transcriptogramS1 = data.frame(),
        transcriptogramS2 = data.frame(),
        DE = data.frame(), radius = 0L, status = NA_integer_, clusters = list(),
        pbc = FALSE, Protein2Symbol = data.frame(), Protein2GO = data.frame(),
        Terms = data.frame(), topGOdata = list()))
