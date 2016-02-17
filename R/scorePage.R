#' scorePage
#' 
#' Gene Expression Signature Similarity based on Parametric Analysis of Gene Set Enrichment
#' 
#' @param dat logFC matrix or list
#' @param top top genes selected
#' @param ncore the number of cores used
#' @param type avg or max
#' @param similarity boolean
#' @param adjust p.adjust.methods. see \code{\link[stats]{p.adjust}}
#' @param center boolean
#' @param verbose verbose
#' @return a list with score, pval, fdr
#' @import methods
#' @import parallel
#' @import foreach
#' @import doParallel
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
#' 
#' data(Psoriasis_Etanercept_LogFC)
#' 
#' sim <- scorePage(Psoriasis_Etanercept_LogFC, type="avg")$score
#' 
#' 
setGeneric("scorePage", 
           function(dat, top=250, ncore=2, type=c("avg", "max"), similarity=TRUE, adjust="fdr", center=TRUE, verbose=FALSE)
               standardGeneric("scorePage"))

#' @rdname scorePage
#' @aliases scorePage
setMethod("scorePage", signature(dat="matrix"),
    function(dat, top=250, ncore=2, type=c("avg", "max"), similarity=TRUE, adjust="fdr", center=TRUE, verbose=FALSE) {
    
    if (center) { dat <- scale(dat, scale=FALSE) }

    # mu is mean of total logFC, mod is related with standard deviation of total logFC.
    mu <- apply(dat, 2, mean, na.rm=TRUE)
    mod <- top^(1/2) / apply(dat, 2, sd, na.rm=TRUE)
    
    cl <- parallel::makeCluster(ncore)
    doParallel::registerDoParallel(cl)
    if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
    nc=NULL
    
    Zscore <- foreach(nc=1:ncol(dat), .combine=cbind, .verbose=verbose) %dopar% {
        i <- colnames(dat)[nc]
        up_gene <- names(sort(dat[,i], decreasing=TRUE)[1:top])
        sm_up <- apply(dat[up_gene,], 2, mean, na.rm=TRUE)
        
        dn_gene <- names(sort(dat[,i], decreasing=FALSE)[1:top])
        sm_dn <- apply(dat[dn_gene,], 2, mean, na.rm=TRUE)
        ((sm_up -mu)*mod - (sm_dn -mu)*mod)/2
    }
    parallel::stopCluster(cl)
    
    return ( zscore_sim(Zscore, type, similarity, adjust) )
    }
)

#' @rdname scorePage
#' @aliases scorePage
setMethod("scorePage", signature(dat="list"),
    function(dat, top=250, ncore=2, type=c("avg", "max"), similarity=TRUE,  adjust="fdr", center=TRUE, verbose=FALSE) {
        
        type <- match.arg(type, c("avg", "max"))
        if (center) { dat <- lapply(dat, function(x) {x-mean(x, na.rm=TRUE)}) }
        
        # mu is mean of total logFC, mod is related with standard deviation of total logFC.
        mu <- sapply(dat, mean, na.rm=TRUE)
        mod <- top^(1/2) / sapply(dat, sd, na.rm=TRUE)
        
        cl <- parallel::makeCluster(ncore)
        doParallel::registerDoParallel(cl)
        if (verbose) {print(paste("getDoParWorkers:", foreach::getDoParWorkers()))}
        nc=NULL
        
        Zscore <- foreach(nc=1:length(dat), .combine=cbind, .verbose=verbose) %dopar% {
            i <- names(dat)[nc]
            up_gene <- names(sort(dat[[i]], decreasing=TRUE)[1:top])
            sm_up <- sapply(dat, function(x) {mean(x[up_gene], na.rm=TRUE) })
            
            dn_gene <- names(sort(dat[[i]], decreasing=FALSE)[1:top])
            sm_dn <- sapply(dat, function(x) {mean(x[dn_gene], na.rm=TRUE) })
            ((sm_up -mu)*mod - (sm_dn -mu)*mod)/2
        }
        parallel::stopCluster(cl)
        
        return ( zscore_sim(Zscore, type, similarity, adjust) )
    }
)

#' zscore_sim
#' convert zscore to similarity
#' 
#' @keywords internal
zscore_sim <- function(Zscore, type="avg", similarity=TRUE, adjust="fdr") {
    colnames(Zscore) <- rownames(Zscore)
    
    if (type=="avg") {
        page_score <- (Zscore + t(Zscore) )/2
    } else if (type=="max") {
        page_score <- pmax(Zscore, t(Zscore) )
    }
    
    # ref: PGSEA
    pval <- 2*pnorm(-abs(page_score))
    fdr <- p.adjust(pval, method = adjust)
    fdr <- matrix(fdr, nrow=nrow(pval), ncol=ncol(pval), dimnames = dimnames(pval))
    
    if (similarity) {
        page_score <- page_score / apply( abs(page_score), 2, max)
    }
    return (list(score=page_score, pval=pval, fdr=fdr))
}


