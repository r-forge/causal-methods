.reshapeSelKnots <- function(model, w, mnk, group)
{
    w <- matrix(w, nrow=mnk)
    w <- lapply(1:ncol(w) ,function(i) {
        if (w[1,i]==0)
            return(NULL)
        sort(w[w[,i]!=0,i])})
    names(w) <- names(model$knots[[group]])
    w
}

.reshapePval <- function(model, pval, group)
{
    knots <- model$knots[[group]]
    nk <- sapply(knots, length)
    p <- length(nk)
    mnk <- max(nk)
    pval <- matrix(pval, mnk, p)
    pval[pval>1]  <- NA
    pval <- lapply(1:p, function(i) {
        if (nk[i] == 0)
            return(NA)
        k <- numeric(nk[i])
        names(k) <- names(knots[[i]])
        k[] <- pval[1:nk[i],i]
        k})
    names(pval) <- names(knots)
    list(pval = pval, knots = knots)
}

.modelPrepF <- function(model, w0, w1, pvalT=function(p) 1/log(p))
{
    z <- model$data[, model$treat]
    y <- model$data[, model$nameY]
    x <- model.matrix.cslseModel(model)[,,drop=FALSE]
    n <- length(z)
    id1 <- z==1
    n1 <- sum(id1)
    n0 <- n-n1
    p <- ncol(x)
    nk0 <- sapply(model$knots$nontreated, length)
    nk1 <- sapply(model$knots$treated, length)
    mnk0 <- max(nk0)
    tnk0 <- sum(nk0)
    mnk1 <- max(nk1)
    tnk1 <- sum(nk1)
    pvt0 <- min(pvalT((tnk0+p)/p), 1)
    pvt1 <- min(pvalT((tnk1+p)/p), 1)    
    k0 <- sapply(model$knots$nontreated, function(ki) {
        k <- numeric(mnk0)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    k1 <- sapply(model$knots$treated, function(ki) {
        k <- numeric(mnk1)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    if (missing(w0))
    {
        w0 <- sapply(1:p, function(i) c(mnk0+1, rep(0, mnk0-1)))
    } else if (is.null(w0)) {
        w0 <- matrix(0, mnk0, p)
    } else {
        w0 <- sapply(w0, function(wi) {
            w <- numeric(mnk0)
            if (!is.null(wi))
                    w[1:length(wi)] <- wi
            w})
    }
    if (missing(w1))
    {
        w1 <- sapply(1:p, function(i) c(mnk1+1, rep(0, mnk1-1)))
    } else if (is.null(w1)) {
        w1 <- matrix(0, mnk1, p)
    } else {
        w1 <- sapply(w1, function(wi) {
            w <- numeric(mnk1)
            if (!is.null(wi))
                    w[1:length(wi)] <- wi
            w})
    }    
    list(y0=y[!id1], y1=y[id1], x0=x[!id1,,drop=FALSE], x1=x[id1,,drop=FALSE],
         p=p, n1=n1, n0=n0, k0=k0, nk0=nk0, k1=k1, nk1=nk1, mnk1=mnk1, mnk0=mnk0,
         tnk0=tnk0, tnk1=tnk1, w0=w0, w1=w1, pvt0=pvt0, pvt1=pvt1)    
}

## The default HCCM is HC0 because we only want to sort the
## p-values. It avoids having to compute the hat values, which slows
## down the procedure, especially for FLSE. Being consistent is good
## enough for that.

selMod_F <- function(model, selType=c("BLSE","FLSE"),
                     selCrit = c("AIC", "BIC", "PVT"),
                     pvalT=function(p) 1/log(p),
                     vT=c("HC0", "vcov", "HC1", "HC2", "HC3"))
{
    selType <- match.arg(selType)
    met <- ifelse(selType=="BLSE", 1, 2)
    if (is.null(model$selections))
    {
        model$selections <- list()
        model$selections$knots <- model$knots
    } else {
        model$knots <- model$selections$knots
    }
    if (!is.null(model$selection[[selType]]))
    {
        warning(paste("Selection by ", selType, " has already been computed.",
                      " The selection will be replaced by the new one.", sep=""))
    }
    model$selections[[selType]] <- list()
    vT <- match.arg(vT)
    vT <- switch(vT,
                 vcov = -1,
                 HC0 = 0,
                 HC1 = 1,
                 HC2 = 2,
                 HC3 = 3)
    spec <- .modelPrepF(model, pvalT=pvalT)
    res <- .Fortran(F_selmodel, as.double(spec$y0), as.double(spec$y1), as.double(spec$x0),
                    as.double(spec$x1), as.integer(spec$n0), as.integer(spec$n1),
                    as.integer(spec$p), as.double(1e-7), as.double(spec$pvt0),
                    as.double(spec$pvt1), as.integer(met), as.integer(vT), as.integer(2),
                    as.double(spec$k0), as.integer(spec$nk0), as.integer(spec$mnk0),
                    as.integer(spec$tnk0),
                    as.double(spec$k1), as.integer(spec$nk1), as.integer(spec$mnk1),
                    as.integer(spec$tnk1),
                    pval0=double(spec$mnk0*spec$p), pval1=double(spec$mnk1*spec$p),
                    bic=double(spec$tnk0+spec$tnk1+1), aic=double(spec$tnk0+spec$tnk1+1),
                    w0bic=integer(spec$mnk0*spec$p), w0aic=integer(spec$mnk0*spec$p),
                    w0pvt=integer(spec$mnk0*spec$p),
                    w1bic=integer(spec$mnk1*spec$p), w1aic=integer(spec$mnk1*spec$p),
                    w1pvt=integer(spec$mnk1*spec$p), npval=integer(1))
    pval <- list(treated=.reshapePval(model, res$pval1, "treated"),
                 nontreated=.reshapePval(model, res$pval0, "nontreated"))
    class(pval) <- "slsePval"
    model$selections[[selType]]$pval <- pval
    model$selections[[selType]]$PVT <-
        list(treated=.reshapeSelKnots(model, res$w1pvt, spec$mnk1, "treated"),
             nontreated=.reshapeSelKnots(model, res$w0pvt, spec$mnk0, "nontreated"))
    model$selections[[selType]]$Threshold <- c(treated=spec$pvt1, nontreated=spec$pvt0)
    if (selType !="PVT")
    {
        model$selections[[selType]]$AIC <-
            list(treated=.reshapeSelKnots(model, res$w1aic, spec$mnk1, "treated"),
                 nontreated=.reshapeSelKnots(model, res$w0aic, spec$mnk0, "nontreated"))
        model$selections[[selType]]$BIC <-
            list(treated=.reshapeSelKnots(model, res$w1bic, spec$mnk1, "treated"),
                 nontreated=.reshapeSelKnots(model, res$w0bic, spec$mnk0, "nontreated"))
        model$selections[[selType]]$IC <- cbind(AIC=res$aic[1:(res$npval+1)],
                                                BIC=res$bic[1:(res$npval+1)])
    }
    model$knots <- update(model$knots, model$selections[[selType]][[selCrit]])
    attr(model$knots$treated, "pval") <- pval$treated
    attr(model$knots$nontreated, "pval") <- pval$nontreated
        attr(model$knots$treated, "curSel") <- attr(model$knots$nontreated, "curSel") <-
            attr(model$knots, "curSel") <-  list(select=selType, crit=selCrit)
    model
}


