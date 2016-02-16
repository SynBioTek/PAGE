#' scorePage: Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment
#' 
#' Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment
#' 
#' @param dat logFC matrix.
#' @param top top genes selected.
#' @param center boolean - median center gene expression matrix columns prior to analysis.
#' @param type avg or max
#' @param similarity boolean
#' @return a list
#' @export 
#' @examples 
#' 
#' data(Psoriasis_Etanercept_LogFC_mat)
#' sim <- scorePage(Psoriasis_Etanercept_LogFC_mat, type="avg")$score
#' 
scorePage <- function(dat, top=250, type=c("avg", "max"), similarity=TRUE) {
    
    dat <- as.matrix(dat)

    # mu is mean of total logFC, mod is related with standard deviation of total logFC.
    mu <- apply(dat, 2, mean, na.rm=TRUE)
    mod <- top^(1/2) / apply(dat, 2, sd, na.rm=TRUE)
    
    Zscore_up <- matrix(NA, nrow=ncol(dat), ncol=ncol(dat), dimnames = list(colnames(dat), colnames(dat)) )
    Zscore_dn <- matrix(0, nrow=ncol(dat), ncol=ncol(dat), dimnames = list(colnames(dat), colnames(dat)) )
    Zscore <- matrix(0, nrow=ncol(dat), ncol=ncol(dat), dimnames = list(colnames(dat), colnames(dat)) )
    
    for (i in colnames(dat) ) {
        
        up_gene <- names(sort(dat[,i], decreasing=TRUE)[1:top])
        sm_up <- apply(dat[up_gene,], 2, mean, na.rm=TRUE)
        Zscore_up[,i] <- (sm_up -mu)*mod
        
        dn_gene <- names(sort(dat[,i], decreasing=FALSE)[1:top])
        sm_dn <- apply(dat[dn_gene,], 2, mean, na.rm=TRUE)
        Zscore_dn[,i] <- (sm_dn -mu)*mod
    }
    Zscore <- (Zscore_up - Zscore_dn)/2
    
    if (type=="avg") {
        page_score <- (Zscore + t(Zscore) )/2
    } else if (type=="max") {
        page_score <- pmax(Zscore, t(Zscore) )
    }
    
    p.val <- 2*pnorm(-abs(page_score))
    
    if (similarity) {
        page_score <- page_score / apply( abs(page_score), 2, max)
    }
    
    return (list(score=page_score, pval=p.val, Zscore_up=Zscore_up, Zscore_dn=Zscore_dn))
}



