\name{transcriptogramer-News}
\title{transcriptogramer News}
\encoding{UTF-8}

\section{Version 1.5.1}{
  \itemize{
    \item Change on the differentiallyExpressed() method return: output now contains raw p-values.
    \item Change on the clusterEnrichment() method return: output now contains raw p-values.
    \item A section called "Frequently asked questions" has been added in the vignette.
  }
}

\section{Version 1.3.6}{
  \itemize{
    \item Slot genesInTerm added to the Transcriptogram class.
    \item Argument boundaryConditions from differentiallyExpressed(): default changed from FALSE to TRUE.
    \item Argument onlyGenesInDE from clusterVisualization(): default changed from TRUE to FALSE.
    \item Argument onlyGenesInDE from clusterEnrichment(): default changed from TRUE to FALSE.
    \item Argument colors added to the differentiallyExpressed() method.
    \item Argument colors added to the clusterVisualization() method.
    \item Argument colors added to the enrichmentPlot() method.
    \item Argument alpha added to the enrichmentPlot() method.
  }
}

\section{Version 1.3.4}{
  \itemize{
    \item Argument universe from clusterEnrichment(): default changed from "all the proteins present in the transcriptogramS2 slot" to "all the proteins present in the ordering slot".
  }
}

\section{Version 1.3.3}{
  \itemize{
    \item Change on the clusterEnrichment() method return.
    \item Slots Protein2GO, and Terms added to the Transcriptogram class.
    \item New methods: enrichmentPlot() and Terms().
  }
}

\section{Version 1.3.1}{
  \itemize{
    \item Argument boundaryConditions added to the differentiallyExpressed() method.
    \item Slots Protein2Symbol, clusters, and pbc added to the Transcriptogram class.
    \item Argument onlyGenesInDE added to the clusterVisualization() method.
    \item Argument onlyGenesInDE added to the clusterEnrichment() method.
  }
}

\section{Version 1.1.18}{
  \itemize{
    \item The differentiallyExpressed() method plot now uses the ggplot2 package and shows the number of clusters detected.
  }
}

\section{Version 1.1.11}{
  \itemize{
    \item Argument title added to the differentiallyExpressed() method.
  }
}

\section{Version 1.1.7}{
  \itemize{
    \item Argument trend added to the differentiallyExpressed() method to support limma-trend.
  }
}

\section{Version 1.1.6}{
  \itemize{
    \item GPL570 dataset updated.
  }
}

\section{Version 1.1.3}{
  \itemize{
    \item All species datasets (Homo sapiens, Mus musculus, Saccharomyces cerevisiae and Rattus norvegicus) were updated to obtain a better clustering.
  }
}

\section{Version 1.0.0}{
  \itemize{
    \item Released at Bioconductor [2017-10-31].
  }
}
