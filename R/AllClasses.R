#' Class Transcriptogram
#'
#' This S4 class includes methods to use expression values with ordered
#' proteins.
#'
#' @slot association A data frame containing two columns, and rows containing
#' proteins names
#' that are connected
#'
#' @slot ordering A data frame containing two columns, the first one with
#' proteins names,
#' and the second containing its respective position
#'
#' @slot transcriptogramS1 A data frame produced as the result of averaging
#' over all identifiers
#' related to the same protein
#'
#' @slot transcriptogramS2 A data frame produced as the result of averaging
#' over the window
#'
#' @slot radius An integer, non negative, number referring to the window radius
#'
#' @slot status An integer used internally to check the status of the object
#'
#' @slot DE A data frame of differentially expressed proteins
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{getSlot-method}
#'
#' @author
#' Diego Morais
#'
#' @export

setClass("Transcriptogram", representation(association = "data.frame",
    ordering = "data.frame", transcriptogramS1 = "data.frame",
    transcriptogramS2 = "data.frame", DE = "data.frame",
    radius = "integer", status = "integer"),
    prototype = list(association = data.frame(),
        ordering = data.frame(), transcriptogramS1 = data.frame(),
        transcriptogramS2 = data.frame(),
        DE = data.frame(), radius = 0L, status = NA_integer_))
