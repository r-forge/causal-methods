.reshapeSelKnots <- function(knots, w, mnk)
{
    w <- matrix(w, nrow=mnk)
    w <- lapply(1:ncol(w) ,function(i) {
        if (w[1,i]==0)
            return(NULL)
        sort(w[w[,i]!=0,i])})
    names(w) <- names(knots)
    w
}

.reshapePval <- function(knots, pval)
{
    nk <- sapply(knots, length)
    p <- length(nk)
    mnk <- max(max(nk),1)
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
    pval <- list(pval = pval, knots = knots)
    class(pval) <- "slsePval"
    pval
}

## creates all necessary elements to call the fortran subroutine
## it is called for each group separately

.modelPrepF <- function(model, pvalT=function(p) 1/log(p))
{
    y <- model$data[, model$nameY]
    f <- reformulate(model$formX)
    x <- model.matrix(f, model$data)[,,drop=FALSE]
    if (attr(terms(f), "intercept") == 1) 
        x <- x[, -1, drop = FALSE]
    n <- length(y)
    p <- ncol(x)
    nk <- sapply(model$knots, length)
    mnk <- max(max(nk), 1)
    tnk <- sum(nk)
    pvt <- min(pvalT((tnk+p)/p), 1)
    tnk <- max(tnk,1)
    k <- sapply(model$knots, function(ki) {
        k <- numeric(mnk)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    list(y=y, x=x, p=p, n=n, k=k, nk=nk, mnk=mnk, tnk=tnk, pvt=pvt)    
}

## The default HCCM is HC0 because we only want to sort the
## p-values. It avoids having to compute the hat values, which slows
## down the procedure, especially for FLSE. Being consistent is good
## enough for that.

.selCMod <- function(model, selType=c("BLSE","FLSE"),
                     selCrit = c("AIC", "BIC", "PVT"),
                     pvalT=function(p) 1/log(p),
                     vT=c("HC0", "Classical", "HC1", "HC2", "HC3"))
{
    selType <- match.arg(selType)
    selCrit <- match.arg(selCrit)
    met <- ifelse(selType=="BLSE", 1, 2)
    vT <- match.arg(vT)
    critN <- ifelse(selCrit=="PVT", "PVT", paste("J", selCrit, sep=""))    
    for (gi in names(model))
    {
        if (is.null(model[[gi]]$selections))
        {
            model[[gi]]$selections <- list()
            model[[gi]]$selections$originalKnots <- model[[gi]]$knots
        } else {
            model[[gi]]$knots <- model[[gi]]$selections$originalKnots
        }
    }
    vT <- switch(vT,
                 Classical = -1,
                 HC0 = 0,
                 HC1 = 1,
                 HC2 = 2,
                 HC3 = 3)
    spec <- lapply(model, function(mi) .modelPrepF(mi, pvalT))
    selm <- ifelse(selCrit == "PVT", 1, 2)                  
    res <- .Fortran(F_selcmodel, as.double(spec$nontreated$y), as.double(spec$treated$y),
                    as.double(spec$nontreated$x), as.double(spec$treated$x),
                    as.integer(spec$nontreated$n), as.integer(spec$treated$n),
                    as.integer(spec$treated$p), as.double(1e-7),
                    as.double(spec$nontreated$pvt), as.double(spec$treated$pvt),
                    as.integer(met), as.integer(vT), as.integer(selm),
                    as.double(spec$nontreated$k), as.integer(spec$nontreated$nk),
                    as.integer(spec$nontreated$mnk), as.integer(spec$nontreated$tnk),
                    as.double(spec$treated$k), as.integer(spec$treated$nk),
                    as.integer(spec$treated$mnk), as.integer(spec$treated$tnk),
                    pval0=double(spec$nontreated$mnk*spec$nontreated$p),
                    pval1=double(spec$treated$mnk*spec$treated$p),
                    bic=double(spec$nontreated$tnk+spec$treated$tnk+1),
                    aic=double(spec$treated$tnk+spec$nontreated$tnk+1),
                    w0bic=integer(spec$nontreated$mnk*spec$nontreated$p),
                    w0aic=integer(spec$nontreated$mnk*spec$nontreated$p),
                    w0pvt=integer(spec$nontreated$mnk*spec$nontreated$p),
                    w1bic=integer(spec$treated$mnk*spec$treated$p),
                    w1aic=integer(spec$treated$mnk*spec$treated$p),
                    w1pvt=integer(spec$treated$mnk*spec$treated$p), npval=integer(1))
    ## For the fortran code, it is currently only treated=1 amnd nontreated=0
    giv <- c(treated=1, nontreated=0)
    for (gi in names(model))
    {
        resi <- res[grepl(giv[gi], names(res))]
        nresi <- names(resi)
        pval <- .reshapePval(model[[gi]]$knots, resi$pval)
        model[[gi]]$selections[[selType]]$pval <- pval
        model[[gi]]$selections[[selType]]$PVT <-
            .reshapeSelKnots(model[[gi]]$knots, resi[[grep("pvt", nresi)]],
                             spec[[gi]]$mnk)
        model[[gi]]$selections[[selType]]$Threshold <- spec[[gi]]$pvt        
        if (selType !="PVT")
        {
            model[[gi]]$selections[[selType]]$JAIC <-
                .reshapeSelKnots(model[[gi]]$knots, resi[[grep("aic", nresi)]],
                                 spec[[gi]]$mnk)
            model[[gi]]$selections[[selType]]$JBIC <-
                .reshapeSelKnots(model[[gi]]$knots, resi[[grep("bic", nresi)]],
                                 spec[[gi]]$mnk)
            model[[gi]]$selections[[selType]]$JIC <- cbind(AIC=res$aic[1:(res$npval+1)],
                                                           BIC=res$bic[1:(res$npval+1)])
        }
        model[[gi]]$knots <- update(model[[gi]]$knots,
                                    model[[gi]]$selections[[selType]][[critN]])
        attr(model[[gi]]$knots, "curSel") <-  list(select=selType, crit=critN)
    }
    model
}

.selMod <- function(model, selType=c("BLSE","FLSE"),
                    selCrit = c("AIC", "BIC", "PVT"),
                    pvalT=function(p) 1/log(p),
                    vT=c("HC0", "Classical", "HC1", "HC2", "HC3"))
{
    selType <- match.arg(selType)
    selCrit <- match.arg(selCrit)
    met <- ifelse(selType=="BLSE", 1, 2)
    if (is.null(model$selections))
    {
        model$selections <- list()
        model$selections$originalKnots <- model$knots
    } else {
        model$knots <- model$selections$originalKnots
    }
    vT <- match.arg(vT)
    vT <- switch(vT,
                 Classical = -1,
                 HC0 = 0,
                 HC1 = 1,
                 HC2 = 2,
                 HC3 = 3)
    spec <- .modelPrepF(model, pvalT=pvalT)
    if (sum(spec$nk) == 0)
        return(model)
    selm <- ifelse(selCrit == "PVT", 1, 2)
    res <- .Fortran(F_selmodel, as.double(spec$y), as.double(spec$x),
                    as.integer(spec$n), as.integer(spec$p), as.double(1e-7),
                    as.double(spec$pvt), as.integer(met), as.integer(vT), as.integer(selm),
                    as.double(spec$k), as.integer(spec$nk), as.integer(spec$mnk),
                    as.integer(spec$tnk), pval=double(spec$mnk*spec$p),
                    bic=double(spec$tnk+1), aic=double(spec$tnk+1),
                    wbic=integer(spec$mnk*spec$p), waic=integer(spec$mnk*spec$p),
                    wpvt=integer(spec$mnk*spec$p), npval=integer(1))    
    pval <- .reshapePval(model$knots, res$pval)
    model$selections[[selType]]$pval <- pval
    model$selections[[selType]]$PVT <- .reshapeSelKnots(model$knots, res$wpvt, spec$mnk)
    model$selections[[selType]]$Threshold <- spec$pvt
    if (selType !="PVT")
    {
        model$selections[[selType]]$AIC <- .reshapeSelKnots(model$knots, res$waic, spec$mnk)
        model$selections[[selType]]$BIC <- .reshapeSelKnots(model$knots, res$wbic, spec$mnk)
        model$selections[[selType]]$IC <- cbind(AIC=res$aic[1:(res$npval+1)],
                                                BIC=res$bic[1:(res$npval+1)])
    }
    model$knots <- update(model$knots, model$selections[[selType]][[selCrit]])
    attr(model$knots, "curSel") <-  list(select=selType, crit=selCrit)
    model
}


