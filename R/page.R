#' Gene Expression Signature Similarity
#' 
#' Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment.
#' 
#' @param datM A Gene Expression Signature matrix including logFC or rank
#' @param geneV a gene symbol name vector
#' @param gene_up a up-regulated gene symbol name vector
#' @param gene_dn a down-regulated gene symbol name vector
#' @param center boolean - median center gene expression matrix columns prior to analysis.
#' @param isRank datM is rank or not (the top one up-regualted gene has rank 1).
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
#' datM <- apply(-Psoriasis_Etanercept_LogFC, 2, rank)
#' gene_up <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=TRUE), 100) )
#' gene_dn <- names(head(sort(Psoriasis_Etanercept_LogFC[,1], decreasing=FALSE), 100) )
#' page1(datM, gene_up)
#' page1(datM, gene_dn)
#' page2(datM, gene_up, gene_dn)
#' 
#' 
#' @rdname page
page1 <- function(datM, geneV, center=TRUE, isRank=TRUE, adjust="fdr"){
    
    if (center) { datM <- scale(datM, scale=FALSE) }
    if (isRank) {
        datM <- -datM
    }
    
    mu <- apply(datM, 2, mean, na.rm=TRUE)
    sigma <- apply(datM, 2, sd, na.rm=TRUE)
    
    Sm <- apply(datM[geneV,], 2, mean, na.rm=TRUE)
    zscore <- (Sm - mu)* length(geneV)^(1/2)/sigma
    
    # similarity
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(head(sort(x, decreasing=TRUE), length(geneV) ), na.rm=TRUE) })
    z_max <- (top_mean - mu)* length(geneV)^(1/2)/sigma
    similarity <- zscore/abs(z_max)
    
    P_Value <- 2*pnorm(-abs(zscore))
    FDR <- p.adjust(P_Value, method = adjust)
    res <- cbind(zscore, similarity, P_Value, FDR)
    res[order(res[,"P_Value"]),]
}

#' @rdname page
#' @export 
page2 <- function(datM, gene_up, gene_dn, center=TRUE, isRank=TRUE, adjust="fdr"){
    
    
    if (center) { datM <- scale(datM, scale=FALSE) }
    if (isRank) {
        datM <- -datM
    }
    mu <- apply(datM, 2, mean, na.rm=TRUE)
    sigma <- apply(datM, 2, sd, na.rm=TRUE)
    
    gene_up <- intersect(gene_up, rownames(datM))
    gene_dn <- intersect(gene_dn, rownames(datM))
    if (length(gene_up)==0 | length(gene_dn)==0 ) {
        stop("Improper Input. Only HGCN gene symbols are accepted.")
    }
    
    # zscore
    Sm_up <- apply(datM[gene_up,], 2, mean, na.rm=TRUE)
    zscore_up <- (Sm_up - mu)* length(gene_up)^(1/2)/sigma
    
    Sm_dn <- apply(datM[gene_dn,], 2, mean, na.rm=TRUE)
    zscore_down <- (Sm_dn - mu)* length(gene_dn)^(1/2)/sigma
    
    # similarity
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(head(sort(x, decreasing=TRUE), length(gene_up) ), na.rm=TRUE) })
    max_up <- (top_mean - mu)* length(gene_up)^(1/2)/sigma
    sim_up <- zscore_up/max_up
    
    top_mean <- apply(datM, 2, function(x){names(x) <- rownames(datM);  mean(tail(sort(x, decreasing=TRUE), length(gene_dn) ), na.rm=TRUE) })
    max_dn <- (top_mean - mu)* length(gene_dn)^(1/2)/sigma
    sim_dn <- zscore_down/max_dn
    similarity <- round((sim_up + sim_dn)/2, 4)
    
    score <- round((zscore_up - zscore_down)/2, 4)
    zscore_up <- round(zscore_up, 4)
    zscore_down <- round(zscore_down, 4)
    
    P_Value <- as.numeric( format(2*pnorm(-abs(score)), digits=4) )
    FDR <- as.numeric( format(p.adjust(P_Value, method = adjust), digits=4) )
    res <- cbind(zscore_up, zscore_down, score, similarity, P_Value, FDR)
    res[order(res[,"P_Value"]),]
}
