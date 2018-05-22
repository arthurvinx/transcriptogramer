# orderingProperties ####

#' Calculates graph properties projected on the ordered proteins
#'
#' Calculates protein (node) properties, such as: degree/connectivity,
#' number of triangles and
#' clustering coefficient; and properties of the window, region of n
#' (radius * 2 + 1) proteins
#' centered at a protein, such as: connectivity, clustering coefficient
#' and modularity.
#'
#' @details
#' Connectivity/degree of a node is the number of edges it presents.
#' A triangle of a node represents a
#' pair of connected neighbors, the number of triangles on the adjacency
#' list of a node is required to
#' calculate its clustering coefficient. The clustering coefficient of a
#' node measures, in the interval [0, 1],
#' the likelihood that any two of its neighbors are themselves connected,
#' this is calculated by the ratio
#' between the number of triangles that the node has, and the maximum
#' possible number of edges on its cluster
#' (nodeTriangles / (nodeDegree * (nodeDegree - 1) / 2)). The window
#' connectivity is the average connectivity
#' calculated over the window. The window clustering coefficient, a value
#' in the interval [0, 1],
#' is the average clustering coefficient calculated over the window.
#' The window modularity, a value in the
#' interval [0, 1], is defined as the ratio between the total number
#' of edges between any two nodes of the window,
#' and the sum of the degrees of the nodes presents in the window.
#' The window considers periodic boundary
#' conditions to deal with proteins near the ends of the ordering.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param nCores An integer number, referring to the number of processing cores
#' to be used; or a logical value, TRUE indicating that all processing cores
#' should be used, and FALSE indicating the use of just one processing core.
#' The default value of this argument is 1.
#'
#' @return This method returns a data.frame containing: ENSEMBL Peptide ID,
#' its position on the ordering,
#' node degree, number of triangles and clustering coefficient, and window
#' connectivity,
#' clustering coefficient and modularity.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 2)
#' \dontrun{
#' oProperties <- orderingProperties(transcriptogram)
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom igraph graph.data.frame
#' @importFrom igraph count_triangles
#' @importFrom progress progress_bar
#' @importFrom parallel detectCores
#' @importFrom snow makeSOCKcluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association}
#'
#' @references
#' da Silva, S. R. M., Perrone, G. C., Dinis, J. M., and de Almeida, R. M. C. (2014). Reproducibility enhancement and differential expression of non predefined functional gene sets in human genome. BMC Genomics.
#'
#' de Almeida, R. M. C., Clendenon, S. G., Richards, W. G., Boedigheimer, M., Damore, M., Rossetti, S., Harris, P. C., Herbert, B. S., Xu, W. M., Wandinger-Ness, A., Ward, H. H., Glazier, J. A. and Bacallao, R. L. (2016). Transcriptome analysis reveals manifold mechanisms of cyst development in ADPKD. Human Genomics, 10(1), 1–24.
#'
#' Ferrareze, P. A. G., Streit, R. S. A., Santos, P. R. dos, Santos, F. M. dos, Almeida, R. M. C. de, Schrank, A., Kmetzsch, L., Vainstein, M. H. and Staats, C. C. (2017). Transcriptional Analysis Allows Genome Reannotation and Reveals that Cryptococcus gattii VGII Undergoes Nutrient Restriction during Infection. Microorganisms.
#'
#' Rybarczyk-Filho, J. L., Castro, M. A. A., Dalmolin, R. J. S., Moreira, J. C. F., Brunnet, L. G., and de Almeida, R. M. C. (2011). Towards a genome-wide transcriptogram: the Saccharomyces cerevisiae case. Nucleic Acids Research, 39(8), 3005-3016.
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname orderingProperties-method
#' @export

setGeneric("orderingProperties", function(object, nCores = 1L)
    standardGeneric("orderingProperties"),
    package = "transcriptogramer")

# connectivityProperties ####

#' Calculates average graph properties as function of node connectivity
#'
#' Calculates network properties as function of node connectivity/degree (k),
#' such as: probability of a protein of the graph has degree k, average
#' assortativity of the nodes of degree k, and the average clustering
#' coefficient of the nodes of degree k.
#'
#' @param object An object of class Transcriptogram.
#'
#' @details
#' The assortativity of a node can be measured by the average degree of its
#' neighbors.
#'
#' @return This method returns a data.frame containing: unique degrees (k) of
#' the nodes of the graph,
#' probability (pk) of a node of the graph has degree k, average assortativity
#' (ak) of the nodes of degree
#' k, and the average clustering coefficient (ck) of the nodes of degree k.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900)
#' \dontrun{
#' cProperties <- connectivityProperties(transcriptogram)
#' }
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association}
#'
#' @author
#' Diego Morais
#'
#' @importFrom igraph graph.data.frame
#' @importFrom igraph count_triangles
#'
#' @docType methods
#' @rdname connectivityProperties-method
#' @export

setGeneric("connectivityProperties", function(object) standardGeneric("connectivityProperties"),
    package = "transcriptogramer")

# transcriptogramS1 ####

#' Calculates the average of the expression values related to the same protein
#'
#' For each transcriptome sample, this method assigns to each protein the
#' average of the expression values of all the identifiers related to
#' it. It is necessary a
#' \code{dictionary} to map the identifiers to proteins.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param expression A matrix, or data.frame, containing normalized expression
#' values from samples of microarrays or RNA-Seq (log2-counts-per-million).
#'
#' @param dictionary A matrix, or data.frame, containing two columns, the first
#' column must contains the
#' ENSEMBL Peptide ID, and the second column must contains values that appear
#' as rownames in \code{expression},
#' in order to recognize the ENSEMBL Peptide ID of the other column.
#'
#' @param nCores An integer number, referring to the number of processing cores
#' to be used; or a logical value, TRUE indicating that all processing cores
#' should be used, and FALSE indicating the use of just one processing core.
#' The default value of this argument is 1.
#'
#' @return This method creates a data.frame to feed the transcriptogramS1
#' slot of an object
#' of class Transcriptogram. Each row of the data.frame contains: an ENSEMBL
#' Peptide ID, its
#' respective position in the ordering and the mean of the expression values
#' of the identifiers
#' related to the same protein.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom parallel detectCores
#' @importFrom snow makeSOCKcluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#'
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association}
#'
#' @author
#' Diego Morais
#'
#' @references
#' da Silva, S. R. M., Perrone, G. C., Dinis, J. M., and de Almeida, R. M. C. (2014). Reproducibility enhancement and differential expression of non predefined functional gene sets in human genome. BMC Genomics.
#'
#' de Almeida, R. M. C., Clendenon, S. G., Richards, W. G., Boedigheimer, M., Damore, M., Rossetti, S., Harris, P. C., Herbert, B. S., Xu, W. M., Wandinger-Ness, A., Ward, H. H., Glazier, J. A. and Bacallao, R. L. (2016). Transcriptome analysis reveals manifold mechanisms of cyst development in ADPKD. Human Genomics, 10(1), 1–24.
#'
#' Ferrareze, P. A. G., Streit, R. S. A., Santos, P. R. dos, Santos, F. M. dos, Almeida, R. M. C. de, Schrank, A., Kmetzsch, L., Vainstein, M. H. and Staats, C. C. (2017). Transcriptional Analysis Allows Genome Reannotation and Reveals that Cryptococcus gattii VGII Undergoes Nutrient Restriction during Infection. Microorganisms.
#'
#' Rybarczyk-Filho, J. L., Castro, M. A. A., Dalmolin, R. J. S., Moreira, J. C. F., Brunnet, L. G., and de Almeida, R. M. C. (2011). Towards a genome-wide transcriptogram: the Saccharomyces cerevisiae case. Nucleic Acids Research, 39(8), 3005-3016.
#'
#' @docType methods
#' @rdname transcriptogramStep1-method
#' @export

setGeneric("transcriptogramStep1", function(object,
    expression, dictionary, nCores = 1L)
    standardGeneric("transcriptogramStep1"),
    package = "transcriptogramer")

# transcriptogramS2 ####

#' Calculates the average of the expression values using a sliding window
#'
#' To each position of the ordering, this method assigns a value
#' equal to the average of the expression values inside a window, region of n
#' (radius * 2 + 1) proteins
#' centered at a protein. The window considers periodic boundary conditions to
#' deal
#' with proteins near the ends of the ordering.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param nCores An integer number, referring to the number of processing cores
#' to be used; or a logical value, TRUE indicating that all processing cores
#' should be used, and FALSE indicating the use of just one processing core.
#' The default value of this argument is 1.
#'
#' @return This method creates a data.frame to feed the transcriptogramS2
#' slot of an object
#' of class Transcriptogram. Each row of the data.frame contains: the ENSEMBL
#' Peptide ID used as center of the
#' window, its position on the ordering, and the mean of the expression values
#' of the window.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom parallel detectCores
#' @importFrom snow makeSOCKcluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1}
#'
#' @author
#' Diego Morais
#'
#' @references
#' da Silva, S. R. M., Perrone, G. C., Dinis, J. M., and de Almeida, R. M. C. (2014). Reproducibility enhancement and differential expression of non predefined functional gene sets in human genome. BMC Genomics.
#'
#' de Almeida, R. M. C., Clendenon, S. G., Richards, W. G., Boedigheimer, M., Damore, M., Rossetti, S., Harris, P. C., Herbert, B. S., Xu, W. M., Wandinger-Ness, A., Ward, H. H., Glazier, J. A. and Bacallao, R. L. (2016). Transcriptome analysis reveals manifold mechanisms of cyst development in ADPKD. Human Genomics, 10(1), 1–24.
#'
#' Ferrareze, P. A. G., Streit, R. S. A., Santos, P. R. dos, Santos, F. M. dos, Almeida, R. M. C. de, Schrank, A., Kmetzsch, L., Vainstein, M. H. and Staats, C. C. (2017). Transcriptional Analysis Allows Genome Reannotation and Reveals that Cryptococcus gattii VGII Undergoes Nutrient Restriction during Infection. Microorganisms.
#'
#' Rybarczyk-Filho, J. L., Castro, M. A. A., Dalmolin, R. J. S., Moreira, J. C. F., Brunnet, L. G., and de Almeida, R. M. C. (2011). Towards a genome-wide transcriptogram: the Saccharomyces cerevisiae case. Nucleic Acids Research, 39(8), 3005-3016.
#'
#' @docType methods
#' @rdname transcriptogramStep2-method
#' @export

setGeneric("transcriptogramStep2", function(object, nCores = 1L)
    standardGeneric("transcriptogramStep2"),
    package = "transcriptogramer")

# radius<- ####

#' @docType methods
#' @rdname radius-method
#' @export

setGeneric("radius<-", signature = "object",
    function(object, value) standardGeneric("radius<-"),
    package = "transcriptogramer")

# differentiallyExpressed ####

#' Identify which genes are differentially expressed
#'
#' This method uses the \pkg{limma} package to identify which genes are
#' differentially expressed,
#' meeting the \code{pValue} requirement, for the contrast "case-control".
#' The \code{levels} lenght must be
#' equal to the number of samples present in the transcriptogramS2 slot of
#' the \code{object}, and its contents
#' is related to the order that the samples appear. FALSE must be used to
#' indicate case samples,
#' and TRUE to indicate control samples. If \code{species} is NULL, no
#' translation will be done, if \code{species} is a character,
#' the \pkg{biomaRt} package will be used to translate the ENSEMBL
#' Peptide ID to Symbol
#' (Gene Name), and if \code{species} is a data.frame, it will be used
#' instead.
#' If the translation fail for some protein, its ENSEMBL
#' Peptide ID will be present
#' into the Symbol column. This method also groups the proteins detected as
#' differentially expressed
#' in clusters, and plots a graphical representation of the groupings.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param levels A logical vector that classify the columns, referring to
#' samples, of the transcriptogramS2 slot of the \code{object}. FALSE must
#' be used to indicate case samples, and TRUE to indicate control samples.
#'
#' @param pValue A numeric value between 0 and 1 giving the required
#' family-wise error rate
#' or false discovery rate. The default value of this argument is 0.05.
#'
#' @param species A character string that will be used,
#' ignoring case sensitivity,
#' to translate the ENSEMBL Peptide ID to Symbol (Gene Name); or a data.frame
#' containing two columns, the first one with ENSEMBL Peptide IDs (character),
#' which may, or not, to contain the taxonomy ID of the species as prefix,
#' and the second containing its respective Symbol (character). The default
#' value of this argument is the content of the object Protein2Symbol slot.
#'
#' @param adjustMethod Character string specifying p-value adjustment method,
#' the possible values are
#' 'none', 'BH', 'fdr' (equivalent to 'BH'), 'BY' and 'holm'. The default value
#' for this argument is 'BH'.
#'
#' @param trend Logical value, set as TRUE to use the limma-trend approach for RNA-Seq.
#' The default value of this argument is FALSE.
#'
#' @param title An overall title for the plot. The default value of this argument is "Differential expression"
#'
#' @param boundaryConditions Logical value, set as TRUE to check if nearby clusters could be merged.
#' The default value of this argument is FALSE.
#'
#' @return This method creates a data.frame to feed the DE slot of an object
#' of class Transcriptogram. This data.frame of differentially expressed
#' proteins
#' contains the log2-fold-change, the p-values and an
#' integer number that indicates if the protein is downregulated or upregulated.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01)
#'
#' ## translating ENSEMBL Peptide IDs to Symbols
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01,
#' "Homo sapiens")
#'
#' ## these calls also works
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01,
#' "H sapiens")
#'
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01,
#' DEsymbols)
#' }
#'
#' @seealso
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer]{DEsymbols},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2}
#'
#' @author
#' Diego Morais
#'
#' @importFrom limma lmFit
#' @importFrom limma makeContrasts
#' @importFrom limma eBayes
#' @importFrom limma decideTests
#' @importFrom limma contrasts.fit
#' @importFrom graphics plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 scale_linetype_manual
#' @importFrom ggplot2 scale_colour_manual
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom grDevices rainbow
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#' @importFrom stats smooth.spline
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#'
#' @docType methods
#' @rdname differentiallyExpressed-method
#' @export

setGeneric("differentiallyExpressed", function(object,
    levels, pValue = 0.05, species = object@Protein2Symbol,
    adjustMethod = "BH", trend = FALSE,
    title = "Differential expression",
    boundaryConditions = FALSE) standardGeneric("differentiallyExpressed"),
    package = "transcriptogramer")

# clusterVisualization ####

#' Displays graphs of the differentially expressed clusters
#'
#' This method uses the \pkg{RedeR} package to display graphs of the
#' differentially expressed clusters.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param maincomp Logical value, set as TRUE if you want to display only the main component of
#' each cluster. The default value of this argument is FALSE.
#'
#' @param connected Logical value, set as TRUE if you want to display only connected nodes.
#' The default value of this argument is FALSE.
#'
#' @param host The domain name of the machine that is running the RedeR XML-RPC
#' server.
#'
#' @param port An integer specifying the port on which the XML-RPC server should
#' listen.
#'
#' @param clusters An integer vector specifying the clusters to be
#' displayed, if NULL, all clusters will be displayed.
#'
#' @param onlyGenesInDE Logical value, set as TRUE to use only the genes
#' in the DE slot. Set as FALSE to use all the genes referring to the positions
#' in the clusters slot. The default value of this argument is TRUE.
#'
#' @return This function returns an object of the RedPort Class.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01,
#' DEsymbols)
#' rdp <- clusterVisualization(transcriptogram)
#' }
#'
#' @seealso
#' \link[transcriptogramer:differentiallyExpressed-method]{differentiallyExpressed},
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2},
#' \link[RedeR]{RedPort}
#'
#' @author
#' Diego Morais
#'
#' @details
#' RedeR package requirements: Java Runtime Environment (>= 6).
#'
#' @importFrom RedeR RedPort
#' @importFrom RedeR calld
#' @importFrom RedeR subg
#' @importFrom RedeR att.setv
#' @importFrom RedeR addGraph
#' @importFrom RedeR selectNodes
#' @importFrom RedeR relax
#' @importFrom grDevices rainbow
#' @importFrom igraph graph.data.frame
#' @importFrom igraph E
#' @importFrom igraph V
#'
#' @docType methods
#' @rdname clusterVisualization-method
#' @export

setGeneric("clusterVisualization", function(object,
    maincomp = FALSE,
    connected = FALSE, host = "127.0.0.1",
    port = 9091, clusters = NULL,
    onlyGenesInDE = TRUE) standardGeneric("clusterVisualization"),
    package = "transcriptogramer")

# clusterEnrichment ####

#' Term enrichment
#'
#' If \code{species} is a character, this method uses the \pkg{biomaRt} package
#' to build a Protein2GO list, if \code{species} is a data.frame, it will be used
#' instead.
#' The Protein2GO list will be used with the
#' \pkg{topGO} package to detect the most significant terms of each cluster
#' present in the DE slot of the \code{object}.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param universe A character vector containing ENSEMBL Peptide IDs, or NULL,
#' if the universe
#' is composed by all the proteins present in the transcriptogramS2 slot of
#' \code{object}.
#'
#' @param species A character string specifying the species; or a data.frame
#' containing two columns, the first one with ENSEMBL Peptide IDs (character),
#' which may, or not, to contain the taxonomy ID of the species as prefix,
#' and the second containing its respective Gene Ontology term (character).
#'
#' @param ontology A character string specifying the Gene Ontology domain,
#' ignoring case sensitivity,
#' the possible values are 'biological process', 'cellular component' and
#' 'molecular function'.
#' The default value of this argument is 'biological process'.
#'
#' @param algorithm Character string specifying which algorithm to use, the
#' possible values are
#' 'classic', 'elim', 'weight', 'weight01', 'lea' and 'parentchild'.
#' The default value of this argument is 'classic'.
#'
#' @param statistic Character string specifying which test to use, the possible
#' values are
#' 'fisher', 'ks', 't', 'sum' and 'globaltest'.
#' The default value of this argument is 'fisher'.
#'
#' @param pValue A numeric value between 0 and 1 giving the required
#' family-wise error rate or false discovery rate. The default value of this argument is 0.05.
#'
#' @param adjustMethod Character string specifying p-value adjustment method,
#' the possible values are
#' 'none', 'BH', 'fdr' (equivalent to 'BH'), 'BY', 'hochberg', 'hommel',
#' 'bonferroni', and 'holm'.
#' The default value of this argument is 'BH'.
#'
#' @param nCores An integer number, referring to the number of processing cores
#' to be used; or a logical value, TRUE indicating that all processing cores
#' should be used, and FALSE indicating the use of just one processing core.
#' The default value of this argument is 1.
#'
#' @param onlyGenesInDE Logical value, set as TRUE to use only the genes
#' in the DE slot. Set as FALSE to use all the genes referring to the positions
#' in the clusters slot. The default value of this argument is TRUE.
#'
#' @return This method creates a data.frame, containing the most significant
#' terms of each cluster, to feed the Terms slot of an object of class
#' Transcriptogram.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01)
#' transcriptogram <- clusterEnrichment(transcriptogram, species = "Homo sapiens",
#' pValue = 0.005)
#'
#' ## this call also works
#' transcriptogram <- clusterEnrichment(transcriptogram, species = HsBPTerms,
#' pValue = 0.005)
#' }
#'
#' @seealso
#' \link[transcriptogramer:differentiallyExpressed-method]{differentiallyExpressed},
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{HsBPTerms},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2}
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname clusterEnrichment-method
#' @export
#'
#' @importFrom methods new
#' @importClassesFrom topGO topGOdata
#' @importFrom topGO groupGOTerms
#' @importFrom topGO annFUN.gene2GO
#' @importFrom topGO GenTable
#' @importFrom topGO runTest
#' @importFrom biomaRt useMart
#' @importFrom biomaRt getBM
#' @importFrom stats na.omit
#' @importFrom stats p.adjust
#' @importFrom snow stopCluster
#' @importFrom snow parLapply
#' @importFrom snow makeSOCKcluster
#' @importFrom parallel detectCores

setGeneric("clusterEnrichment", function(object,
    universe = NULL, species, ontology = "biological process",
    algorithm = "classic", statistic = "fisher",
    pValue = 0.05, adjustMethod = "BH", nCores = 1L, onlyGenesInDE = TRUE)
    standardGeneric("clusterEnrichment"),
    package = "transcriptogramer")

# radius ####

#' Radius
#'
#' Retrieve or set the content of the radius slot of an object of class
#' Transcriptogram.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param value An non-negative integer referring to the window
#' radius required for some methods.
#'
#' @return This method returns the content of the radius slot of an object of
#' class Transcriptogram.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' radius(transcriptogram) <- 80
#' radius(transcriptogram)
#'
#' @seealso
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2},
#' \link[transcriptogramer:orderingProperties-method]{orderingProperties}
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname radius-method
#' @export

setGeneric("radius", function(object)
    standardGeneric("radius"),
    package = "transcriptogramer")

# DE ####

#' Get DE
#'
#' Gets the content of the DE slot of an object of class Transcriptogram.
#'
#' @param object An object of class Transcriptogram.
#'
#' @return This method returns the content of the DE slot of an object of
#' class Transcriptogram.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01)
#' DE(transcriptogram)
#' }
#'
#' @seealso
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer]{transcriptogramPreprocess}
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname DE-method
#' @export

setGeneric("DE", function(object)
    standardGeneric("DE"),
    package = "transcriptogramer")

# Terms ####

#' Get terms
#'
#' Gets the content of the Terms slot of an object of class Transcriptogram.
#'
#' @param object An object of class Transcriptogram.
#'
#' @return This method returns the content of the Terms slot of an object of
#' class Transcriptogram.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01)
#' transcriptogram <- clusterEnrichment(transcriptogram, species = "Homo sapiens",
#' pValue = 0.005)
#' Terms(transcriptogram)
#' }
#'
#' @seealso
#' \link[transcriptogramer:differentiallyExpressed-method]{differentiallyExpressed},
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{HsBPTerms},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2},
#' \link[transcriptogramer:clusterEnrichment-method]{clusterEnrichment}
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname Terms-method
#' @export

setGeneric("Terms", function(object)
  standardGeneric("Terms"),
  package = "transcriptogramer")

# enrichmentPlot ####

#' Projects Gene Ontology terms on the ordering
#'
#' Plots the rate (number of genes related to a term inside the
#' window/total number of genes in the windows) of given terms.
#'
#' @param object An object of class Transcriptogram.
#'
#' @param nCores An integer number, referring to the number of processing cores
#' to be used; or a logical value, TRUE indicating that all processing cores
#' should be used, and FALSE indicating the use of just one processing core.
#' The default value of this argument is 1.
#'
#' @param nTerms An integer number referring to the number of top terms from
#' each cluster. The default value of this argument is 1.
#'
#' @param GOIDs A character vector containing the Gene Ontology
#' accessions to be plotted. If NULL, the top \code{nTerms} of each cluster
#' will be used.
#'
#' @param title An overall title for the plot. The default value of this
#' argument is "Enrichment"
#'
#' @return This method returns an ggplot2 object.
#'
#' @examples
#' transcriptogram <- transcriptogramPreprocess(association, Hs900, 50)
#' \dontrun{
#' transcriptogram <- transcriptogramStep1(transcriptogram, GSE9988, GPL570)
#' transcriptogram <- transcriptogramStep2(transcriptogram)
#' levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#' transcriptogram <- differentiallyExpressed(transcriptogram, levels, 0.01)
#' transcriptogram <- clusterEnrichment(transcriptogram, species = "Homo sapiens",
#' pValue = 0.005)
#'
#' }
#'
#' @importFrom foreach %dopar%
#' @importFrom progress progress_bar
#' @importFrom parallel detectCores
#' @importFrom snow makeSOCKcluster
#' @importFrom snow stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom tidyr gather
#' @importFrom graphics plot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 element_text
#' @importFrom stats na.omit
#' @importFrom stats smooth.spline
#'
#' @seealso
#' \link[transcriptogramer:differentiallyExpressed-method]{differentiallyExpressed},
#' \link[transcriptogramer]{transcriptogramPreprocess},
#' \link[transcriptogramer]{GSE9988},
#' \link[transcriptogramer]{GPL570},
#' \link[transcriptogramer]{Hs900},
#' \link[transcriptogramer]{HsBPTerms},
#' \link[transcriptogramer]{association},
#' \link[transcriptogramer:transcriptogramStep1-method]{transcriptogramStep1},
#' \link[transcriptogramer:transcriptogramStep2-method]{transcriptogramStep2},
#' \link[transcriptogramer:clusterEnrichment-method]{clusterEnrichment}
#'
#' @author
#' Diego Morais
#'
#' @docType methods
#' @rdname enrichmentPlot-method
#' @export

setGeneric("enrichmentPlot", function(object, nCores = 1L, nTerms = 1L,
                                      GOIDs = NULL, title = "Enrichment")
standardGeneric("enrichmentPlot"),
package = "transcriptogramer")
