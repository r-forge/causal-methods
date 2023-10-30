
.reshapeKnots <- function(model, w, mnk, group)
{
    w <- matrix(w, nrow=mnk)
    w <- lapply(1:ncol(w) ,function(i) {
        if (w[1,i]==0)
            return(NULL)
        sort(w[w[,i]!=0,i])})
    names(w) <- names(model$knots[[group]])
    update(model$knots[[group]], w)
}



.selICF <- function (model, pvalRes, pvalT = NULL, crit) 
{
    x <- model.matrix(model)
    y <- model$data[,model$nameY]
    z <- model$data[,model$treated]
    nk0 <- sapply(model$knots$nontreated, length)
    nk1 <- sapply(model$knots$treated, length)
    mnk0 <- max(nk0)
    tnk0 <- sum(nk0)
    pval <- c(do.call("c", pvalRes$nontreated$pval),
              do.call("c", pvalRes$treated$pval))
    spval <- sort(pval)
    npval <- length(spval)
    pval0 <-  sapply(pvalRes$nontreated$pval, function(pvi) {
        pv <- numeric(mnk0)
        if (!is.na(pvi[1]))
            pv[1:length(pvi)] <- pvi
        pv}) 
    knots0 <- sapply(model$knots$nontreated, function(ki) {
        k <- numeric(mnk0)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    mnk1 <- max(nk1)
    tnk1 <- sum(nk1)
    pval1 <-  sapply(pvalRes$treated$pval, function(pvi) {
        pv <- numeric(mnk1)
        if (!is.na(pvi[1]))
            pv[1:length(pvi)] <- pvi
        pv})    
    knots1 <- sapply(model$knots$treated, function(ki) {
        k <- numeric(mnk1)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    p <- ncol(x)
    nb0 <- tnk0+p+1
    nb1 <- tnk1+p+1
    n <- nrow(x)
    id0 <- z==0
    n0 <- sum(id0)
    sp <- .Fortran(F_selic, as.numeric(y[id0]), as.numeric(y[!id0]),
                   as.numeric(x[id0,]), as.numeric(x[!id0,]),
                   as.integer(n0), as.integer(n-n0),
                   as.integer(p), as.numeric(1e-7),
                   as.numeric(knots0), as.integer(nk0), as.integer(mnk0),
                   as.integer(tnk0),  
                   as.numeric(knots1), as.integer(nk1), as.integer(mnk1),                   
                   as.integer(tnk1),
                   as.numeric(pval0), as.numeric(pval1), as.numeric(spval),
                   as.integer(npval),
                   bic=numeric(npval+1), aic=numeric(npval+1),
                   w0BIC=integer(mnk0*p), w1BIC=integer(mnk1*p),
                   w0AIC=integer(mnk0*p), w1AIC=integer(mnk1*p))
    modelAIC <- model
    modelBIC <- model
    modelAIC$knots$nontreated <-.reshapeKnots(model, sp$w0AIC, mnk0, "nontreated")
    modelAIC$knots$treated <-.reshapeKnots(model, sp$w1AIC, mnk1, "treated") 
    modelBIC$knots$nontreated <-.reshapeKnots(model, sp$w0BIC, mnk0, "nontreated")
    modelBIC$knots$treated <-.reshapeKnots(model, sp$w1BIC, mnk1, "treated") 
    list(AIC=modelAIC, BIC=modelBIC)
}


selSLSE.F <- function(model, selType=c("BLSE", "FLSE"),
                      selCrit = c("AIC", "BIC", "PVT"), 
                      pvalT = function(p) 1/log(p), vcov.=vcovHC, ...)
{
    selCrit <- match.arg(selCrit)
    selType <- match.arg(selType)
    critFct <- if (selCrit == "PVT") {
                   .selPVT
               } else {
                   .selICF
               }
    if (all(sapply(model$knots$treated, function(i) is.null(i))) &
        all(sapply(model$knots$nontreated, function(i) is.null(i))))
    {
        warning("No selection needed: the number of knots is 0 for all confounders")
    } else {
        pval <- pvalSLSE(model, selType, vcov., ...)
        model <- critFct(model, pval, pvalT, selCrit)
        if (selCrit != "PVT")
            model <- model[[selCrit]]            
        attr(model$knots$treated, "pval") <- pval$treated
        attr(model$knots$nontreated, "pval") <- pval$nontreated
        attr(model$knots$treated, "curSel") <- attr(model$knots$nontreated, "curSel") <-
            attr(model$knots, "curSel") <- list(select=selType, crit=selCrit)        
    }
    model
}

mylm <- function(Y, X, type=c("vcov", "HC0", "HC1", "HC2", "HC3"))
{
    type <- match.arg(type)
    type <- switch(type,
                   vcov = -1,
                   HC0 = 0,
                   HC1 = 1,
                   HC2 = 2,
                   HC3 = 3)
    res <- .Fortran(F_myls, as.double(Y), x = as.double(X),
                    as.integer(length(Y)), as.integer(ncol(X)),
                    as.double(1e-7), as.integer(type), rank=integer(1),
                    pvot = as.integer(1:ncol(X)), e = double(length(Y)),
                    b = double(ncol(X)), vcov = double(ncol(X)^2))
    v <- matrix(res$vcov, ncol(X), ncol(X))
    list(x=X, rank=res$rank, pv=res$pvot, b=res$b,
         vcov=v)
}

mypnorm <- function(x, mu=0, sig=1)
{
    res <- .Fortran(F_mypnorm, as.double(x), as.integer(length(x)),
                    as.double(mu), as.double(sig),
                    p = double(length(x)))
    res$p
}


