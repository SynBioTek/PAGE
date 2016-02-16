load("../data/test.RData")
logFCmat <- Psoriasis_Etanercept_LogFC_mat
top=250
i <- names(logFCmat)[1]

scorePage <- function(logFCmat, top=250, center=TRUE) {
    
    if (center) {logFCmat <- scale(logFCmat, scale=FALSE) }
    
    mu <- apply(logFCmat, 2, mean, na.rm=TRUE)
    mod <- top^(1/2) / apply(logFCmat, 2, sd, na.rm=TRUE)
    
    Zscore_up <- matrix(NA, nrow=ncol(logFCmat), ncol=ncol(logFCmat), dimnames = list(colnames(logFCmat), colnames(logFCmat)) )
    Zscore_dn <- matrix(0, nrow=ncol(logFCmat), ncol=ncol(logFCmat), dimnames = list(colnames(logFCmat), colnames(logFCmat)) )
    Zscore <- matrix(0, nrow=ncol(logFCmat), ncol=ncol(logFCmat), dimnames = list(colnames(logFCmat), colnames(logFCmat)) )
    
    for (i in colnames(logFCmat) ) {
        
        up_gene <- names(sort(logFCmat[,i], decreasing=TRUE)[1:top])
        dn_gene <- names(sort(logFCmat[,i], decreasing=FALSE)[1:top])
        
        sm_up <- apply(logFCmat[up_gene,], 2, mean, na.rm=TRUE)
        Zscore_up[,i] <- (sm_up -mu)*mod
        
        sm_dn <- apply(logFCmat[dn_gene,], 2, mean, na.rm=TRUE)
        Zscore_dn[,i] <- (sm_dn -mu)*mod
        
    }
    Zscore <- (Zscore_up - Zscore_dn)/2
}
