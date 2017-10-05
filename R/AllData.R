#' Data frame containing expression values
#'
#' Expression values, obtained by microarray, of 3 cases and 3 controls
#' referring to the innate immune responses to TREM-1 activation.
#' The data frame has 6 columns, each one contains expression values of a
#' sample,
#' the first 3 columns are case samples, and the last 3 are control samples.
#' Each row
#' contain expression values obtained by the probe mentioned in its respective
#' rowname.
#' The expression values were normalized using the affy package and, to reduce
#' the
#' storage space required for the data, this data frame is a subset from the
#' original
#' samples (GSM252443, GSM252444, GSM252445, GSM252465, GSM252466, GSM252467),
#' containing
#' only the rows on which the probes are mapped by the platform GPL570
#' dictionary.
#'
#' @examples
#' GSE9988
#'
#' @seealso
#' \link[transcriptogramer]{GPL570}
#'
#' @author
#' Diego Morais
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9988}{GSE9988}

"GSE9988"

#' Dictionary Gene2Probe
#'
#' A mapping between ENSEMBL Peptide ID and probe identifier, for the Homo
#' sapiens
#' and the platform GPL570, [HG-U133_Plus_2] Affymetrix Human Genome U133
#' Plus 2.0 Array.
#'
#' @format Each row of the data frame contains two variables:
#' \describe{
#'  \item{ENSP}{The ENSEMBL Peptide ID}
#'  \item{Probe}{The probe identifier}
#' }
#'
#' @examples
#' GPL570
#'
#' @seealso
#' \link[transcriptogramer]{GSE9988}
#'
#' @author
#' Diego Morais

"GPL570"

#' Association
#'
#' A subset of the Homo sapiens protein network data from STRINGdb,
#' release 10.5. This subset
#' contains only associations of proteins of combined score greater than or
#' equal to 900.
#'
#' @format Each row of the data frame contains two variables:
#' \describe{
#'  \item{V1}{The ENSEMBL Peptide ID of the first protein}
#'  \item{V2}{The ENSEMBL Peptide ID of the second protein}
#' }
#'
#' @examples
#' association
#'
#' @seealso
#' \link[transcriptogramer]{Hs900}
#'
#' @author
#' Diego Morais

"association"

#' Ordered Homo sapiens proteins of combined score greater than or equal to 700
#'
#' A character vector containing the Homo sapiens proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 700.
#'
#' @examples
#' Hs700
#'
#' @author
#' Diego Morais

"Hs700"

#' Ordered Homo sapiens proteins of combined score greater than or equal to 800
#'
#' A character vector containing the Homo sapiens proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 800.
#'
#' @examples
#' Hs800
#'
#' @author
#' Diego Morais

"Hs800"

#' Ordered Homo sapiens proteins of combined score greater than or equal to 900
#'
#' A character vector containing the Homo sapiens proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 900.
#'
#' @examples
#' Hs900
#'
#' @author
#' Diego Morais

"Hs900"

#' Ordered Rattus norvegicus proteins of combined score greater than or equal
#' to 700
#'
#' A character vector containing the Rattus norvegicus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 700.
#'
#' @examples
#' Rn700
#'
#' @author
#' Diego Morais

"Rn700"

#' Ordered Rattus norvegicus proteins of combined score greater than or equal
#' to 800
#'
#' A character vector containing the Rattus norvegicus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 800.
#'
#' @examples
#' Rn800
#'
#' @author
#' Diego Morais

"Rn800"

#' Ordered Rattus norvegicus proteins of combined score greater than or equal
#' to 900
#'
#' A character vector containing the Rattus norvegicus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 900.
#'
#' @examples
#' Rn900
#'
#' @author
#' Diego Morais

"Rn900"

#' Ordered Mus musculus proteins of combined score greater than or equal to 700
#'
#' A character vector containing the Mus musculus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 700.
#'
#' @examples
#' Mm700
#'
#' @author
#' Diego Morais

"Mm700"

#' Ordered Mus musculus proteins of combined score greater than or equal to 800
#'
#' A character vector containing the Mus musculus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 800.
#'
#' @examples
#' Mm800
#'
#' @author
#' Diego Morais

"Mm800"

#' Ordered Mus musculus proteins of combined score greater than or equal to 900
#'
#' A character vector containing the Mus musculus proteins, from STRINGdb
#' release 10.5,
#' of combined score greater than or equal to 900.
#'
#' @examples
#' Mm700
#'
#' @author
#' Diego Morais

"Mm900"

#' Ordered Saccharomyces cerevisiae proteins of combined score greater than or
#' equal to 700
#'
#' A character vector containing the Saccharomyces cerevisiae proteins, from
#' STRINGdb release 10.5,
#' of combined score greater than or equal to 700.
#'
#' @examples
#' Sc700
#'
#' @author
#' Diego Morais

"Sc700"

#' Ordered Saccharomyces cerevisiae proteins of combined score greater than or
#' equal to 800
#'
#' A character vector containing the Saccharomyces cerevisiae proteins, from
#' STRINGdb release 10.5,
#' of combined score greater than or equal to 800.
#'
#' @examples
#' Sc800
#'
#' @author
#' Diego Morais

"Sc800"

#' Ordered Saccharomyces cerevisiae proteins of combined score greater than or
#' equal to 900
#'
#' A character vector containing the Saccharomyces cerevisiae proteins, from
#' STRINGdb release 10.5,
#' of combined score greater than or equal to 900.
#'
#' @examples
#' Sc900
#'
#' @author
#' Diego Morais

"Sc900"

#' Dictionary Gene2Symbol
#'
#' A mapping between ENSEMBL Peptide ID and Symbol (Gene Name) of a reduced set
#' of proteins.
#'
#' @examples
#' DEsymbols
#'
#' @author
#' Diego Morais

"DEsymbols"

#' Dictionary Gene2GO
#'
#' A mapping between ENSEMBL Peptide ID and Gene Ontology, biological process,
#' terms of a set of proteins.
#'
#' @examples
#' HsBPTerms
#'
#' @author
#' Diego Morais

"HsBPTerms"
