#' Gene Expression Signature Similarity based on GSEA
#' 
#' Gene Expression Signature Similarity based on  Gene Set Enrichment Analysis
#' 
#' @param dat gene rank matrix or data.frame
#' @param top top genes selected
#' @param ncore the number of cores used
#' @param type avg or max
#' @import 
#' @examples 
#' 
#' 
gsea <- function(dat, top=250, ncore=2, type=c("avg", "max") ) {
    
    
    ES <- matrix(-2, ncol=ncol(dat), nrow=ncol(dat) )
    for (i in 1:ncol(dat) ) {
        for (j in 1:ncol(dat) ) {
            # i = 1; j=2
            ES[i,j] <- quickenrichmentscore(which(dat[,j]<=top), which(dat[,j]>=nrow(dat)-top+1), dat[,i])
        }
    }
    if (type=="avg") {
        sim <- (ES + t(ES))/2
    } else if (type=="max") {
        sim <- pmax(ES, t(ES))
    }
    rownames(sim) <- colnames(sim) <- colnames(dat)
    
    return (sim)
    
}
