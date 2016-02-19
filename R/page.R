#' Gene Expression Signature Similarity
#' 
#' 
#' Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment.
#' 
#' @param datM A Gene Expression Signature matrix including logFC or rank
#' @param geneV a gene symbol name vector
#' @param gene_up a up-regulated gene symbol name vector
#' @param gene_dn a down-regulated gene symbol name vector
#' @param isRank datM is rank or not
#' @param adjust p.adjust.methods. see \code{\link[stats]{p.adjust}}
#' 
#' @export 
#' 
#' @references 
#' 1. Kim, Seon-Young, and David J. Volsky. "PAGE: parametric analysis of gene set 
#' enrichment." BMC bioinformatics 6.1 (2005): 1.
#' 
#' 2. Iorio, Francesco, et al. "Discovery of drug mode of action and drug 
#' repositioning from transcriptional responses." Proceedings of the National 
#' Academy of Sciences 107.33 (2010): 14621-14626.
#' 
#' @examples
#' data(Psoriasis_Etanercept_LogFC)
#' datM <- apply(Psoriasis_Etanercept_LogFC, 2, rank)
#' gene_up <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=TRUE), 100) )
#' gene_dn <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=FALSE), 100) )
#' page1(datM, gene_up)
#' page1(datM, gene_dn)
#' page2(datM, gene_up, gene_dn)
#' 
#' @rdname page
page1 <- function(datM, geneV, isRank=TRUE, adjust="fdr"){
    
    if (isRank) {
        datM <- scale(datM, scale=FALSE)
    }
    mu <- apply(datM, 2, mean, na.rm=TRUE)
    sigma <- apply(datM, 2, sd, na.rm=TRUE)
    
    Sm <- apply(datM[geneV,], 2, mean, na.rm=TRUE)
    zscore <- (Sm - mu)* length(geneV)^(1/2)/sigma
    
    # similarity
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(head(sort(x, decreasing=TRUE), length(geneV) ), na.rm=TRUE) })
    z_max <- (top_mean - mu)* length(geneV)^(1/2)/sigma
    sim <- zscore/abs(z_max)
    
    pval <- 2*pnorm(-abs(zscore))
    fdr <- p.adjust(pval, method = adjust)
    cbind(zscore, sim, pval, fdr)
}

#' @rdname page
#' @export 
page2 <- function(datM, gene_up, gene_dn, isRank=TRUE, adjust="fdr"){
    
    if (isRank) {
        datM <- scale(datM, scale=FALSE)
    }
    mu <- apply(datM, 2, mean, na.rm=TRUE)
    sigma <- apply(datM, 2, sd, na.rm=TRUE)
    
    #zscore
    Sm_up <- apply(datM[gene_up,], 2, mean, na.rm=TRUE)
    zscore_up <- (Sm_up - mu)* length(gene_up)^(1/2)/sigma
    
    Sm_dn <- apply(datM[gene_dn,], 2, mean, na.rm=TRUE)
    zscore_dn <- (Sm_dn - mu)* length(gene_dn)^(1/2)/sigma
    
    # similarity
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(head(sort(x, decreasing=TRUE), length(gene_up) ), na.rm=TRUE) })
    max_up <- (top_mean - mu)* length(gene_up)^(1/2)/sigma
    sim_up <- zscore_up/max_up
    
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(tail(sort(x, decreasing=TRUE), length(gene_dn) ), na.rm=TRUE) })
    max_dn <- (top_mean - mu)* length(gene_dn)^(1/2)/sigma
    sim_dn <- zscore_dn/max_dn
    sim <- (sim_up + sim_dn)/2
    
    score <- (zscore_up - zscore_dn)/2
    pval <- 2*pnorm(-abs(score))
    fdr <- p.adjust(pval, method = adjust)
    cbind(score, sim, pval, fdr)
}
