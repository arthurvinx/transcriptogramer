## ----message = FALSE-------------------------------------------------------
library(transcriptogramer)
t <- transcriptogram.preprocess(association = association, ordering = Hs900)

## ----message = FALSE-------------------------------------------------------
## during the preprocessing

## this command creates the object and set the radius as 0
t <- transcriptogram.preprocess(association = association, ordering = Hs900)

## this command creates the object and set the radius as 50
t <- transcriptogram.preprocess(association = association, ordering = Hs900,
                                radius = 50)

## ----message = FALSE-------------------------------------------------------
## after the preprocessing

## this command modify the radius of an existing Transcriptogram object
t <- setRadius(.Object = t, radius = 25)

## this command get the radius of an existing Transcriptogram object
r <- getSlot(.Object = t, slot = "radius")

## ----message = FALSE-------------------------------------------------------
oPropertiesR25 <- orderingProperties(.Object = t)

## ----eval = FALSE----------------------------------------------------------
#  ## slight change of radius
#  t <- setRadius(.Object = t, radius = 30)
#  
#  oPropertiesR30 <- orderingProperties(.Object = t)

## ----message = FALSE-------------------------------------------------------
cProperties <- connectivityProperties(.Object = t)

## ----message = FALSE-------------------------------------------------------
t <- transcriptogramStep1(.Object = t, expression = GSE9988,
                          dictionary = GPL570)

## ----eval = FALSE----------------------------------------------------------
#  t <- transcriptogramStep2(.Object = t)

## ----eval = FALSE----------------------------------------------------------
#  t <- setRadius(.Object = t, radius = 50)
#  t <- transcriptogramStep2(.Object = t)

## ----eval = FALSE----------------------------------------------------------
#  levels <- c(rep(FALSE, 3), rep(TRUE, 3))
#  t <- differentiallyExpressed(.Object = t, levels = levels, pValue = 0.005)

## ----eval = FALSE----------------------------------------------------------
#  ## this command also translates ENSEMBL Peptide IDs to Symbols
#  ## internet connection is required
#  t <- differentiallyExpressed(.Object = t, levels = levels, pValue = 0.005,
#                               species = "Homo sapiens")
#  
#  ## this command also translates ENSEMBL Peptide IDs to Symbols
#  ## see dataset DEsymbols
#  t <- differentiallyExpressed(.Object = t, levels = levels, pValue = 0.005,
#                               species = DEsymbols)

## ----eval = FALSE----------------------------------------------------------
#  View(getSlot(.Object = t, slot = "DE"))
#  
#  ##this also works
#  View(t@DE)

## ----eval = FALSE----------------------------------------------------------
#  rdp <- clusterVisualization(.Object = t)
#  
#  ## if the ENSEMBL Peptide IDs were translated to Symbols
#  rdp <- clusterVisualization(.Object = t, symbolAsNodeAlias = TRUE)

## ----eval = FALSE----------------------------------------------------------
#  ## internet connection is required
#  terms <- clusterEnrichment(.Object = t, species = "Homo sapiens",
#                             pValue = 0.005)
#  
#  ## see HsBPTerms dataset
#  terms <- clusterEnrichment(.Object = t, species = HsBPTerms,
#                             pValue = 0.005)

## ----echo = FALSE----------------------------------------------------------
load("terms.RData")

## --------------------------------------------------------------------------
head(terms)

## --------------------------------------------------------------------------
sessionInfo()

## --------------------------------------------------------------------------
warnings()

