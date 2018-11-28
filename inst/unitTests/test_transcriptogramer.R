testData <- function(){
    checkTrue(is.data.frame(association) &&
        nrow(association) == 547620 && ncol(association) ==
        2, "association")
    checkTrue(is.data.frame(GPL570) &&
        nrow(GPL570) == 33043 && ncol(GPL570) ==
        2, "GPL570")
    checkTrue(is.data.frame(GSE9988) &&
        nrow(GSE9988) == 24804 && ncol(GSE9988) ==
        6, "GSE9988")
    checkTrue(is.character(Hs900) && length(Hs900) ==
        11030, "Hs900")
}

ord <- Hs900[seq_len(1000)]
links <- association
links <- links[which(links$protein1 %in% ord),]
t <- transcriptogramPreprocess(links, ord, 2)

testOrderingProperties <- function(){
        oProperties <- orderingProperties(t, nCores = TRUE)
        checkTrue(is.data.frame(oProperties) &&
            ncol(oProperties) == 8 && nrow(oProperties) ==
            1000 && oProperties[757, 5] ==
            (oProperties[757, 4]/(oProperties[757,
                3] * (oProperties[757, 3] -
                1)/2)), "orderingProperties")
}

testTranscriptogramStep2 <- function(){
        t <- transcriptogramStep1(t, GSE9988,
            GPL570, nCores = TRUE)
        t <- transcriptogramStep2(t, nCores = TRUE)
        checkTrue(class(t) == "Transcriptogram" &&
            is.data.frame(t@transcriptogramS2) &&
            nrow(t@transcriptogramS2) == 846 && ncol(t@transcriptogramS2) == 8)
        checkTrue(t@transcriptogramS2[3,
            3] == sum(c(t@transcriptogramS1[2, 3], t@transcriptogramS1[3, 3],
            t@transcriptogramS1[4, 3], t@transcriptogramS1[5, 3]))/(t@radius*2+1))
}
