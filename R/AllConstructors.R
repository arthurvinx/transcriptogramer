#' Creates an object of class Transcriptogram
#'
#' Constructor for the Transcriptogram object.
#'
#' @param association A matrix, or data.frame, containing two columns of
#' ENSEMBL Peptide IDs (character);
#' or the path for a file containing two columns, no header, with rows
#' composed by the
#' ENSEMBL Peptide IDs of two proteins that are connected.
#'
#' @param ordering A character vector containing ordered ENSEMBL Peptide IDs;
#' a data.frame containing
#' two columns, the first one with ENSEMBL Peptide IDs (character),
#' and the second containing
#' its respective position (integer non negative); or the path for a file
#' containing two columns,
#' a row for the headers, with rows composed respectively, by a ENSEMBL Peptide
#' ID and its respective position.
#'
#' @param radius An integer, non negative, number referring to the window
#' radius required for some
#' methods.
#'
#' @return A preprocessed object of class Transcriptogram.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900)
#'
#' @seealso
#' \link[transcriptogramer]{Transcriptogram-class},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer]{Hs900}
#'
#' @author
#' Diego Morais
#'
#' @importFrom methods new
#'
#' @export

transcriptogramPreprocess <- function(association,
    ordering, radius = 0L) {
    message("preprocessing the input data... step 1 of 1")
    association = check_association(association)
    ordering = check_ordering(ordering)
    radius = check_radius(radius)
    if (!(length(unique(association$p1)) ==
        nrow(ordering) && all(ordering$Protein %in%
        unique(association$p1)))) {
        stop("arguments association and ordering - make sure that the ",
            "ordering was generated from the association!")
    }
    object <- methods::new("Transcriptogram",
        association = association, ordering = ordering,
        radius = radius, status = 0L)
    return(object)
}
