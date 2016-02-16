#' scorePage: Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment
#' 
#' Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment
#' 
#' @param dat logFC matrix.
#' @param top top genes selected.
#' @param type avg or max
#' @param similarity boolean
#' @param isRank boolean
#' @return a list
#' @export 
#' @references 
#' 1. Kim, Seon-Young, and David J. Volsky. "PAGE: parametric analysis of gene set 
#' enrichment." BMC bioinformatics 6.1 (2005): 1.
#' 
#' 2. Iorio, Francesco, et al. "Discovery of drug mode of action and drug 
#' repositioning from transcriptional responses." Proceedings of the National 
#' Academy of Sciences 107.33 (2010): 14621-14626.
#' @examples 
#' 
#' data(Psoriasis_Etanercept_LogFC)
#' 
#' sim <- scorePage(Psoriasis_Etanercept_LogFC_mat, type="avg")$score
#' sim1 <- scorePage(Psoriasis_Etanercept_LogFC_rank, type="avg", isRank=TRUE)$score
#' 
#'  
scorePage <- function(dat, top=250, type=c("avg", "max"), similarity=TRUE, isRank=FALSE) {
    
    dat <- as.matrix(dat)
    if (isRank) { dat <- scale(dat, scale=FALSE) }

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



