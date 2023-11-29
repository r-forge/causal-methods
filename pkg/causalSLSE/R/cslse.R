## This documents include all R functions for causal-SLSE models
## Many are based on functions for SLSE models.
################################################################

 
## causal-SLSE knots
####################

##  Splines builders

llSplines <- function(object, ...)
    UseMethod("llSplines")

llSplines.cslseModel <- function (object, ...)
{
    X <- lapply(names(object), function(gi)
    {
        X2 <- lapply(names(object), function(gj)
        {
            if (gi != gj)
                object[[gi]]$data <- object[[gj]]$data
            llSplines(object[[gi]])
        })
        names(X2) <- names(object)
        X2
    })
    names(X) <- names(object)
    X
}

## CSLSE Model
###############

cslseModel <- function (form, data, nbasis = function(n) n^0.3, 
                        knots,  groupInd = c(treated=1, nontreated=0))
{
    tmp <- as.character(form)
    tmp2 <- strsplit(tmp[3], "\\|")[[1]]    
    if (!grepl("\\|", tmp[3])) 
        stop("form must be of the type y~z|~x")
    if (is.null(names(groupInd)))
        stop("groupInd must be a named vector")
    if (length(groupInd)>2)
        stop("Only one treatment is possible in this current version of the package.")
    if (length(groupInd)<2)
        stop("You need two groups in groupInd.")
    if (!all(names(groupInd) %in% c("treated","nontreated")))
        stop("The only allowed names for groupInd are currently treated and nontreated.")  
    treat <- all.vars(form)[2]  
    Z <- data[,treat]
    form <-  as.formula(paste(tmp[2], tmp2[2], sep = ""))
    if (any(naZ <- is.na(Z)))
    {
        warning("Missing values in the treatment indicator are not allowed. The observations have been removed")
        Z <- Z[!naZ,,drop=FALSE]
        data <- data[!naZ,,drop=FALSE]
    }
    if (!(all(Z %in% groupInd)))
        stop("The treatment indicator can only contain the values included in groupInd")
    data <- lapply(groupInd, function(zi) data[Z==zi,,drop=FALSE])
    names(data) <- names(groupInd)
    if (missing(knots))
    {
        obj <- lapply(names(groupInd), function(gi) slseModel(form, data[[gi]], nbasis))
    } else if (is.null(knots)) {
        obj <- lapply(names(groupInd), function(gi) slseModel(form, data[[gi]], nbasis, NULL))
    } else {
        if (!is.list(knots))
            stop("knots must be a list")
        if (is.null(names(knots)))
            stop("knots must be a named list")
        if (!all(names(knots) %in% names(groupInd)))
            stop(paste("The names of knots must be in: ",
                       paste(names(groupInd), collapse=", ", sep=""), sep=""))
        obj <- lapply(names(groupInd), function(gi)
        {
            if (gi %in% names(knots))
                slseModel(form, data[[gi]], nbasis, knots[[gi]])
            else
                slseModel(form, data[[gi]], nbasis)
        })
    }
    attr(obj, "treatedVar") <- treat
    attr(obj, "groupInd") <- groupInd
    names(obj) <- names(groupInd)
    class(obj) <- "cslseModel"
    obj
}

## cslseModel methods

.prtSel <- function(x)
{
    groups <- names(x)
    curSel <- sapply(x, function(gi) c(attr(gi$knots, "curSel")$select,
                                       attr(gi$knots, "curSel")$crit))
    initSel <- sapply(x, function(gi) c(attr(gi$knots, "initSel")$select,
                                        attr(gi$knots, "initSel")$crit))
    sameInit <- all(initSel[1, 1] == initSel[1,-1]) &  all(initSel[2,1] == initSel[2,-1])
    sameCur <- all(curSel[1, 1] == curSel[1,-1]) &  all(curSel[2,1] == curSel[2,-1])
    if (sameCur)
    {
        crit <- attr(x[[1]]$knots,"curSel")$crit
        cat("Selection method: ", attr(x[[1]]$knots,"curSel")$select, sep="")
        if (crit != "")
        {
            if (!grepl("J", crit) & crit!="PVT") {
                crit <- paste(crit, "(Sep.)")
            } else {
                crit <- gsub("J", "", crit)}
            cat("-",  crit, sep = "")
        }
        cat("\n\n")            
    } else {
        for (gi in groups)
        {
            cat("Selection method for the ", gi, ": ", sep="")
            cat(attr(x[[gi]]$knots,"curSel")$select, sep="")
            crit <- attr(x[[gi]]$knots,"curSel")$crit
            if (crit != "")
                {
                    if (!grepl("J", crit) & crit!="PVT")
                        crit <- paste(crit, "(Sep.)") 
                    cat("-",  gsub("J", "", crit), sep = "")
                }
            cat("\n")
        }
        cat("\n")
    }
    invisible()
}

print.cslseModel <- function(x, which=c("Model", "selKnots", "Pvalues"),
                              digits = max(3L, getOption("digits") - 3L), ...)
{
    groups <- names(x)
    n <- sapply(x, function(xi) nrow(xi$data))
    which <- match.arg(which)
    if (which == "Model")
    {
        cat("Causal Semiparametric LSE Model\n")
        cat("*******************************\n\n")
        for (gi in groups)
        {
            cat("Number of ", gi, ": ", n[gi], "\n")
            if (length(x[[gi]]$na))
                cat("Number of NA (", gi,"): ", length(x[[gi]]$na), "\n")
        }
        .prtSel(x)
        cat("Confounders approximated by SLSE:\n")
        for (gi in names(x))
        {
            w <- sapply(x[[gi]]$knots, is.null)
            selPW <- x[[1]]$nameX[!w]
            isApp <- if (length(selPW)) paste(selPW, collapse=", ", sep="") else "None"
            cat("\t", gi, ": ", isApp, "\n", sep="")
        }
        cat("Confounders not approximated by SLSE:\n") 
        for (gi in names(x))
        {
            w <- sapply(x[[gi]]$knots, is.null)
            nonselPW <- x[[1]]$nameX[w]
            notApp <- if (length(nonselPW)) paste(nonselPW,collapse=", ",sep="") else "None"
            cat("\t", gi, ": ", notApp, "\n", sep="")            
        }
    } else if (which == "selKnots") {
        for (gi in names(x))
        {
            cat(gi, "\n")
            cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
            print(x[[gi]]$knots, header="Select", digits=digits, ...)
        }
    } else {
        if (is.null(x$selections))
        {
            cat("No p-values are available. You must apply a selection methods first.\n")
        } else {
            selType <- attr(x,"curSel")$select
            if (is.null(x$selections[[selType]]$pval))
            {
                cat("No p-values are available in the model object\n")
            } else {
                print(x$selections[[selType]]$pval, digits=digits, ...)
            }
        }
    }
    invisible()
}

update.cslseModel <- function(object, selKnots, selType, selCrit="AIC",
                              pvalT = function(p) 1/log(p), vcov.=vcovHC, ...)
{
    treat <- attr(object, "treatedVar") 
    group <- attr(object, "groupInd") 
    if (!missing(selKnots))
    {
        if (is.null(selKnots))
        {
            object <- lapply(object, function(mi) update(mi, selKnots=NULL))
            class(object) <- "cslseModel"
            attr(object, "treatedVar") <- treat
            attr(object, "groupInd") <- group
        } else {
            if (!is.list(selKnots))
                stop("selKnots must be a list")
            if (is.null(names(selKnots)))
                stop("selKnots must be named")
            if (!all(names(selKnots) %in% names(object)))
                stop("The names of selKnots must be included in the group names of the models")
            for (gi in names(selKnots))
                object[[gi]] <- update(object[[gi]], selKnots=selKnots[[gi]])
        }
    } else if (!missing(selType)) {
        selCrit <- ifelse(selCrit=="PVT", "PVT", paste("J", selCrit, sep=""))
        object <- lapply(object, function(mi) update(object=mi, selType=selType, selCrit=selCrit,
                                                     pvalT=pvalT, vcov.=vcov., ...))
        class(object) <- "cslseModel"
        attr(object, "treatedVar") <- treat
        attr(object, "groupInd") <- group
    }
    object
}


## convert objects to models (slse or cslse)

as.model <- function(x, ...)
    UseMethod("as.model")

as.model.slseFit <- function(x, ...)
    x$model

as.model.cslseFit <- function(x, ...)
{
    mod <- lapply(x, as.model)
    class(mod) <- "cslseModel"
    attr(mod, "treatedVar") <- attr(x, "treatedVar")
    attr(mod, "groupInd") <- attr(x, "groupInd")
    mod
}

as.model.cslse <- function(x, ...)
{
    treat <- attr(x, "treatedVar")
    group <- attr(x, "groupInd")
    mod <- lapply(names(group), function(gi) x[[gi]]$model)
    names(mod) <- names(group) 
    class(mod) <- "cslseModel"
    attr(mod, "treatedVar") <- treat
    attr(mod, "groupInd") <- group
    mod
}


## The estSLSE generic function
## It estimates the model by least squares
#############################################

estSLSE <- function(model, ...)
    UseMethod("estSLSE")

estSLSE.cslseModel <- function(model, selKnots, ...)
{
    model <- update(model, selKnots)
    fit <- lapply(model, function(mi) estSLSE(mi))
    attr(fit, "treatedVar") <- attr(model, "treatedVar")
    attr(fit, "groupInd") <- attr(model, "groupInd")
    class(fit) <- "cslseFit"
    fit
}

## cslseFit methods

print.cslseFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Causal Semiparametric LSE\n")
    cat("**************************\n")
    .prtSel(lapply(x, function(xi) xi$model))
    for (gi in names(x))
    {
        cat(gi, "\n")
        cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
        print.default(coef(x[[gi]]$LSE), print.gap = 2L, digits=digits,
                      quote = FALSE)
        cat("\n")
    }
    invisible()
}

summary.cslseFit <- function (object, vcov.=vcovHC, ...) 
{
    ans <- lapply(object, function(obj) summary(obj, vcov., ...))
    class(ans) <- "summary.cslseFit"
    ans
}

print.summary.cslseFit <- function(x, groups,
                                    digits = max(3L, getOption("digits") - 3L),
                                    signif.stars = getOption("show.signif.stars"),
                                    ...)
{
    cat("Causal Semiparametric LSE\n")
    cat("**************************\n")  
    .prtSel(lapply(x, function(xi) xi$model))
    if (missing(groups))
        groups <- names(x)
    for (gi in groups)
    {
        cat(gi, "\n")
        cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
        q <- capture.output(print(x[[gi]]$lseSum, digits=digits,
                                  signif.stars=signif.stars))
        q <- q[-(1:(grep("Residuals:", q)-1))]
        q <- q[1:grep("Multiple R-squared", q)]
        q <- paste(q, collapse="\n")
        cat(q)
        cat("\n\n")
    }
    invisible()
}


## The selSLSE generic function
## It selects the knots using BLSE or FLSE

selSLSE <- function(model, ...)
    UseMethod("selSLSE")

selSLSE.cslseModel <- function(model, selType=c("BLSE", "FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"), 
                               pvalT = function(p) 1/log(p),
                               vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
                               reSelect=FALSE, ...)
{
    selCrit <- match.arg(selCrit)
    selType <- match.arg(selType)
    vcovType <- match.arg(vcovType)
    tmp <- ifelse(selCrit=="PVT", "PVT", paste("J", selCrit, sep=""))
    chk <- sapply(model, function(mi)
        !is.null(mi$selections[[selType]][[tmp]]))
    if (all(chk) &  !reSelect)
        return(update(model, selType = selType, selCrit = selCrit))
    mod <- .selCMod(model, selType, selCrit, pvalT, vcovType)
    mod
}

## The cslsePval methods
## Currently not exported and only used internally

print.cslsePval <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    groups <- names(x)
    for (gi in groups)
    {
        cat(gi, "\n")
        cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
        print(x[[gi]], digits=digits, ...)
    }
    invisible()

}


## Causal effects functions
#############################

## Hidden  functions to compute the causal effect and their SE's

.causali <- function(Z, vcov, U, beta, causal, e)
{
    id <- switch(causal,
                 ACE = rep(TRUE, length(Z)),
                 ACT = Z==1,
                 ACN = Z==0)
    n <- sum(id)
    X <- cbind(1-Z, Z, U$nontreated*(1-Z), U$treated*Z)
    SigmaX <- crossprod(X)
    X0 <- U$nontreated[id,, drop = FALSE]
    X1 <- U$treated[id,, drop = FALSE]
    Xbar0 <- colMeans(X0)
    Xbar1 <- colMeans(X1)
    beta0 <- beta$nontreated[-1]
    beta1 <- beta$treated[-1]
    est <- beta$treated[1] - beta$nontreated[1] +
        sum(beta1*Xbar1) - sum(beta0*Xbar0)
    vcovXf0 <- cov(X0)
    vcovXf1 <- cov(X1)
    vcovXf01 <- cov(X0, X1)
    Dvec0 <- c(1, Xbar0)
    Dvec1 <- c(1, Xbar1)
    Dvec <- c(-1, 1, -Xbar0, Xbar1)
    e <- do.call("c", e)[id]
    X <- X[id,,drop=FALSE]
    X1 <- c(scale(X1, scale=FALSE)%*%beta1)
    X0 <- c(scale(X0, scale=FALSE)%*%beta0)
    tmp2 <- c(e*X%*%solve(SigmaX, Dvec))
    addT <- 2*mean(X1*tmp2) - 2*mean(X0*tmp2)
    se <- (sum(Dvec0*c(crossprod(vcov$nontreated, Dvec0))) +
           sum(Dvec1*c(crossprod(vcov$treated, Dvec1))) +
           sum(beta0*c(crossprod(vcovXf0, beta0)))/n +
           sum(beta1*c(crossprod(vcovXf1, beta1)))/n -
           2 * c(beta0 %*% vcovXf01 %*% beta1)/n +
           2*mean(X1*tmp2) - 2*mean(X0*tmp2))^0.5
    ans <- c(est, se)
    names(ans) <- c("est","se")
    ans
}

.causal <- function(model, fit, vcov, U, 
                     causal=c("ACE", "ACT", "ACN", "ALL"))
{
    causal <- match.arg(causal)
    notNA <- lapply(fit, function(fi) !is.na(coef(fi$LSE))[-1])
    beta <- lapply(fit, function(fi)  na.omit(coef(fi$LSE)))
    e <- lapply(fit, function(fi) residuals(fi$LSE))
    n <- sapply(e, length)
    Z <- do.call("c", lapply(names(e), function(gi)
        if (gi == "nontreated") rep(0, length(e[[gi]])) else rep(1, length(e[[gi]]))))
    U <- lapply(names(U), function(gi)
    {
        U.gi <- do.call(rbind, U[[gi]])
        if (length(notNA[[gi]]))
            U.gi <- U.gi[, notNA[[gi]], drop=FALSE]
        U.gi
    })
    names(U) <- names(model)
    if (causal == "ALL") causal <- c("ACE","ACT","ACN")
    ans <- lapply(causal, function(ci)
        .causali(Z, vcov, U, beta, ci, e))
    names(ans) <- causal
    ans
}

## The generic function causalSLSE

causalSLSE <- function(object, ...)
{
    UseMethod("causalSLSE")
}

causalSLSE.cslseModel <- function(object,
                                  selType=c("SLSE","BLSE","FLSE"),
                                  selCrit = c("AIC", "BIC", "PVT"),
                                  selVcov = c("HC0", "Classical", "HC1", "HC2", "HC3"),
                                  causal = c("ALL","ACT","ACE","ACN"),
                                  pvalT = function(p) 1/log(p),
                                  vcov.=vcovHC, reSelect=FALSE, ...)
{
    selType <- match.arg(selType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    selVcov <- match.arg(selVcov)    
    if (selType != "SLSE")
        object <- selSLSE(object, selType, selCrit, pvalT, selVcov, reSelect)
    res <- estSLSE(object)    
    beta <- lapply(res, function(fi)  coef(fi$LSE))
    v <- lapply(res, function(fi)  vcov.(fi$LSE, ...))
    for (gi in names(object))
    {
        if (any(is.na(beta[[gi]]))) 
            warning(paste("\nThe regression is multicollinear for the ", gi, ".",
                          " The result may not be valid:", 
                          "\nThe following variables produced NA's\n",
                          paste(names(beta[[gi]])[is.na(beta[[gi]])], 
                                collapse = ", "), "\n", sep = ""))
    }
    U <- llSplines(object)
    ans <- .causal(object, res, v, U, causal)
    ans <- c(ans, res)
    class(ans) <- c("cslse", "cslseFit")
    attr(ans, "treatedVar") <- attr(res, "treatedVar")
    attr(ans, "groupInd") <- attr(res, "groupInd")    
    ans    
}

causalSLSE.cslseFit <- function(object, causal = c("ALL","ACT","ACE","ACN"),
                               vcov.=vcovHC, ...)
{
    causal <- match.arg(causal)
    beta <- lapply(object, function(fi)  coef(fi$LSE))
    v <- lapply(object, function(fi)  vcov.(fi$LSE, ...))
    for (gi in names(object))
    {
        if (any(is.na(beta[[gi]]))) 
            warning(paste("\nThe regression is multicollinear for the ", gi, ".",
                          " The result may not be valid:", 
                          "\nThe following variables produced NA's\n",
                          paste(names(beta[[gi]])[is.na(beta[[gi]])], 
                                collapse = ", "), "\n", sep = ""))
    }
    U <- llSplines(as.model(object))
    ans <- .causal(as.model(object), object, v, U, causal)
    ans <- c(ans, object)
    class(ans) <- c("cslse", "cslseFit")
    attr(ans, "treatedVar") <- attr(object, "treatedVar")
    attr(ans, "groupInd") <- attr(object, "groupInd")    
    ans        
}

causalSLSE.formula <- function(object, data, nbasis=function(n) n^0.3,
                               knots, 
                               selType=c("SLSE","BLSE","FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               selVcov = c("HC0", "Classical", "HC1", "HC2", "HC3"),                               
                               causal = c("ALL","ACT","ACE","ACN"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, reSelect=FALSE, ...)
{
    model <- cslseModel(object, data, nbasis,  knots)
    selType <- match.arg(selType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    selVcov <- match.arg(selVcov)        
    causalSLSE(object=model, selType = selType, selCrit = selCrit,
               selVcov = selVcov, causal = causal, pvalT =  pvalT, 
               vcov. = vcov., reSelect, ...)    
}


## Function to recompute the standard error
## Not currently exported. Only used in simulations
## to test the method.


cslseSE <- function(object, vcov.=vcovHC, ...)    
{
    if(!inherits(object, "cslse"))
        stop("object must be of class cslse")
    model <- as.model(object)
    group <- attr(model, "groupInd")
    treat <- attr(model, "treatedVar")
    causal <- c("ACE", "ACT", "ACN")
    causal <- causal[which(causal %in% names(object))]
    if (length(causal) > 1) 
        causal <- "ALL"
    res <- lapply(names(model), function(gi) object[[gi]])
    names(res) <- names(model)
    class(res) <- "cslseFit2"
    attr(res, "treatedVar") <- treat
    attr(res, "groupInd") <- group
    v <- lapply(res, function(fi)  vcov.(fi$LSE, ...))
    U <-  llSplines(model)
    notNA <- lapply(res, function(fi) !is.na(coef(fi$LSE))[-1])
    beta <- lapply(res, function(fi)  na.omit(coef(fi$LSE)))
    e <- lapply(res, function(fi) residuals(fi$LSE))
    n <- sapply(e, length)
    Z <- do.call("c", lapply(names(e), function(gi)
        if (gi == "nontreated") rep(0, length(e[[gi]])) else rep(1, length(e[[gi]]))))    
    U <- lapply(names(U), function(gi)
    {
        U.gi <- do.call(rbind, U[[gi]])
        if (length(notNA[[gi]]))
            U.gi <- U.gi[, notNA[[gi]], drop=FALSE]
        U.gi
    })
    names(U) <- names(model)
    if (causal == "ALL") causal <- c("ACE","ACT","ACN")
    se <- lapply(causal, function(ci) .causali(Z, v, U, beta, ci, e)["se"])
    names(se) <- causal   
    se
}

## cslse methods

print.cslse <- function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    gi <- names(x)[!(names(x)%in%c("ACE","ACT","ACN"))]
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    .prtSel(lapply(x[gi], function(xi) xi$model))
    cat("\n")
    for (causal in c("ACE","ACT","ACN"))
    {
        if (!is.null(x[[causal]]))
            cat(causal, " = ", format(x[[causal]]["est"], digits=digits),
                "\n", sep="")
    }
}


summary.cslse <- function (object, ...) 
{
    w <- sapply(c("ACE","ACT","ACN"), function(ci) !is.null(object[[ci]]))
    est <- sapply(c("ACE","ACT","ACN")[w], function(ci) {
        est <- object[[ci]]["est"]
        se <-  object[[ci]]["se"]
        t <- est/se
        pv <- 2*pnorm(-abs(t))
        c(est, se, t, pv)})
    est <- t(est)
    dimnames(est) <- list(c("ACE","ACT","ACN")[w],
                          c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    lse <- lapply(c("nontreated", "treated"), function(gi)
        summary(object[[gi]]))
    names(lse) <- c("nontreated", "treated")
    ans <- list(causal = est, lse=lse)
    class(ans) <- "summary.cslse"
    ans
}

print.summary.cslse <- function (x, digits = max(3L, getOption("digits") - 3L), 
                                signif.stars = getOption("show.signif.stars"), 
                                beta = FALSE, knots = FALSE, ...) 
{
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    .prtSel(lapply(x$lse, function(xi) xi$model))
    cat("\n")
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    if (beta) {
        cat("\nBasis function coefficients\n")
        cat("*****************************\n")
        for (gi in names(x$lse))
        {
            cat(gi, "\n")
            cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
            print(x$lse[[gi]], digits = digits, signif.stars = signif.stars)
        }
    }
    if (knots) {
        cat("\nNumber of selected knots per confounder\n")
        cat("****************************************\n")
        for (gi in names(x$lse))
        {
            cat(gi, "\n")
            cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
            print.default(format(sapply(x$lse[[gi]]$model$knots, length)), print.gap = 2L, 
                          quote = FALSE)
            cat("\n")
        }
    }
    invisible()
}

## extract methods for texreg tables
#######################################

extract.cslse <- function (model, include.nobs = TRUE,
                           include.nknots = TRUE,
                           include.numcov = TRUE, include.rsquared = TRUE,
                           include.adjrs = TRUE,
                           separated.rsquared = FALSE,
                           which=c("ALL","ACE","ACT","ACN","ACE-ACT",
                                   "ACE-ACN","ACT-ACN"), ...) 
{
    which <- match.arg(which)
    type <- c("ACE","ACT","ACN")
    if (isTRUE(include.adjrs) | isTRUE(include.rsquared))
    {
        if (!isTRUE(separated.rsquared))
        {
            groupInd <- attr(model, "groupInd")
            Y <- do.call("c", lapply(names(groupInd), function(gi) fitted(model[[gi]]$LSE)))
            e <- do.call("c", lapply(names(groupInd), function(gi) residuals(model[[gi]]$LSE)))            
            SSE <- sum((Y-mean(Y))^2)
            SSR <- sum(e^2)
            dfY <- sum(sapply(names(groupInd), function(gi) nobs(model[[gi]]$LSE)))-1
            dfe <- sum(sapply(names(groupInd), function(gi) model[[gi]]$LSE$df.residual))
            R2 <- SSE/(SSE+SSR)
            R2adj <- 1 - (1 - R2)*dfY/dfe
        }

    }
    w <- if (which == "ALL") type
         else type[sapply(type, function(ti) grepl(ti, which))]
    s <- summary(model)
    co <- s$causal
    co <- co[rownames(co) %in% w,,drop=FALSE]
    se <- co[,2]
    pval <- co[,4]
    co <- co[,1]
    names(pval) <- names(se) <- names(co) <- w
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (isTRUE(include.nknots)) {
        rs0 <- length(unlist(model$nontreated$model$knots))
        rs1 <- length(unlist(model$treated$model$knots))        
        gof <- c(gof, rs0, rs1)
        gof.names <- c(gof.names, "Num. knots (Nontreated)",
                       "Num. knots (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
   }
    if (isTRUE(include.numcov)) {
        rs3 <- length(model$treated$model$nameX)
        gof <- c(gof, rs3)
        gof.names <- c(gof.names, "Num. confounders")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    if (isTRUE(include.nobs)) {
        n1 <- nrow(model$treated$model$data)
        n0 <- nrow(model$nontreated$model$data)
        gof <- c(gof, n0, n1)
        gof.names <- c(gof.names, "Num. obs. (Nontreated)", "Num. obs. (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
    }
    if (isTRUE(include.rsquared)) {
        if (!isTRUE(separated.rsquared))
        {
            gof <- c(gof, R2)
            gof.names <- c(gof.names, "R$^2$")
            gof.decimal <- c(gof.decimal, TRUE)            
        } else {
            s0 <- summary(model$nontreated$LSE)
            s1 <- summary(model$treated$LSE)
            R2 <- c(s0$r.squared, s1$r.squared)
            gof <- c(gof, R2)
            gof.names <- c(gof.names, "R$^2$ (nontreated)", "R$^2$ (treated)")
            gof.decimal <- c(gof.decimal, TRUE, TRUE)
        }
    }
    if (isTRUE(include.adjrs)) {
        if (!isTRUE(separated.rsquared))
        {
            gof <- c(gof, R2adj)
            gof.names <- c(gof.names, "Adj. R$^2$")
            gof.decimal <- c(gof.decimal, TRUE)            
        } else {
            s0 <- summary(model$nontreated$LSE)
            s1 <- summary(model$treated$LSE)
            R2 <- c(s0$adj.r.squared, s1$adj.r.squared)
            gof <- c(gof, R2)
            gof.names <- c(gof.names, "Adj. R$^2$ (nontreated)",
                           "Adj. R$^2$ (treated)")
            gof.decimal <- c(gof.decimal, TRUE, TRUE)
        }
    }   
    tr <- createTexreg(coef.names = names(co), coef = co, se = se, 
                       pvalues = pval, gof.names = gof.names, gof = gof,
                       gof.decimal = gof.decimal)
    return(tr)
}

setMethod("extract", signature = className("cslse", "causalSLSE"),
          definition = extract.cslse)

## Other utilities for printing slseFit results

extract.slseFit <- function (model, include.rsquared = TRUE, include.adjrs = TRUE, 
                             include.nobs = TRUE, include.fstatistic = FALSE, include.rmse = FALSE,
                             ...)
{ 
    extract(model$LSE, include.rsquared, include.adjrs, 
            include.nobs, include.fstatistic, include.rmse)
}

setMethod("extract", signature = className("slseFit", "causalSLSE"),
          definition = extract.slseFit)



as.list.cslseFit <- function(x, ...)
{
    class(x) <- "list"
    x
}


## The predict and plot method for cslseFit
############################################

## utility functions

.findrep <- function(lst1, lst2 = NULL)
{
    if (length(lst2) == 0)
        return(lst1)
    if (!is.list(lst2))
        stop("Additional graphical parameters must be provided in a list")
    if (is.null(names(lst2)))
        stop("You must name your list of graphical parameters")
    if (any(duplicated(names(lst2))))
        stop("Some graphical parameters are duplicated")
    if (any(names(lst2) %in% names(lst1)))
    {
        w <- which(names(lst2) %in% names(lst1))
        lst1[names(lst2)[w]] <- lst2[w]
        lst2 <- lst2[-w]
    }
    if (length(lst2))
        lst1 <- c(lst1, lst2)
    lst1
}


.initParCSLSE <- function()
{
    treated <- list(points=list(pch = 21, col = 2),
                    lines=list(col = 2, lty = c(2, 3, 3), lwd = 2, type='l'))
    nontreated <- list(points=list(pch = 22, col = 1),
                       lines=list(col = 1, lty = c(1, 3, 3), lwd = 2, type='l'))    
    common <- list()
    legend <- list(x="topright", bty='n')
    list(treated=treated, nontreated=nontreated,
         common=common, legend=legend)
}

.cslsePar <- function(addPar, startPar=.initParCSLSE())
{
    startPar$treated$points <- .findrep(startPar$treated$points, addPar$treated$points)
    startPar$treated$lines <- .findrep(startPar$treated$lines, addPar$treated$lines) 
    startPar$nontreated$points <- .findrep(startPar$nontreated$points,
                                           addPar$nontreated$points)
    startPar$nontreated$lines <- .findrep(startPar$nontreated$lines,
                                          addPar$nontreated$lines)
    startPar$common <- .findrep(startPar$common, addPar$common)
    startPar$legend <- .findrep(startPar$legend, addPar$legend)
    startPar
}

## predict

predict.cslseFit <- function (object, interval = c("none", "confidence"),
                               se.fit = FALSE, 
                               newdata = NULL, level = 0.95, vcov. = vcovHC, ...) 
{
    interval <- match.arg(interval)
    treatedVar <- attr(object, "treatedVar")
    groupInd <- attr(object, "groupInd")
    group <- names(groupInd)
    if (!is.null(newdata))
    {
        if (!(treatedVar %in% names(newdata)))
            stop("The treatment indicator must be included in newdata")
        newdata <- lapply(group, function(gi)
        {
            dati <- newdata[newdata[,treatedVar]==groupInd[gi],,drop=FALSE]
            if (nrow(dati)) dati else NULL
        })
    } else {
        newdata <- lapply(group, function(gi) NULL)
    }
    names(newdata) <- group
    pr <- lapply(group, function(gi)
        predict(object[[gi]], interval, se.fit, newdata[[gi]], level, vcov., ...))
    names(pr) <- group
    pr
}

## plot

plot.cslseFit <- function (x, y, which = y, interval = c("none", "confidence"), 
                           level = 0.95, fixedCov = list(),
                           vcov. = vcovHC, add = FALSE, addToLegend = NULL, 
                           addPoints = FALSE, FUN = mean, plot=TRUE, graphPar=list(), ...) 
{
    interval <- match.arg(interval)
    group <- c("treated", "nontreated")
    if (length(fixedCov) == 0)
    {
        fixedCov <- lapply(group,  function(gi) NULL)
        names(fixedCov) <- group
    } else  if ( any(names(fixedCov) %in% group) ) {
        if (!all(names(fixedCov) %in% group))
            stop("When fixedCov is group specific, treated and nontreated are the only allowed names")
    } else {
        fixedCov <- list(treated=fixedCov, nontreated=fixedCov)        
    }
    pDat <- lapply(group, function(gi)
        plot(x[[gi]], which, interval=interval, level=level, fixedCov=fixedCov[[gi]],
             vcov.=vcov., FUN=FUN, plot=FALSE))
    names(pDat) <- group
    nameY <- x$treated$model$nameY
    if (!plot)
        return(pDat)
    ylim <- if(addPoints)
            {
                range(do.call("c", lapply(pDat, function(pi) pi[,nameY])))
            } else {
                range(do.call("c", lapply(pDat, function(pi) pi[,-(1:2)])))
            }
    xlim <- range(do.call("c", lapply(pDat, function(pi) pi[,which])))
    common <- list(xlab=which, ylab=nameY, ylim=ylim, xlim=xlim,
                   main=paste(nameY, " vs ", which, " using SLSE", sep = ""))
    allPar <- .cslsePar(list(common=common))
    allPar <- .cslsePar(graphPar, allPar)
    lpch <- if(addPoints)
            {
                c(allPar$treated$points$pch, allPar$nontreated$points$pch)
            } else {
                c(NA,NA)
            }
    leg <- list(legend=c("Treated", "Nontreated"),
                pch = lpch,
                col = c(allPar$treated$lines$col, allPar$nontreated$lines$col),
                lty = c(allPar$treated$lines$lty[1], allPar$nontreated$lines$lty[1]))
    allPar <- .cslsePar(list(legend=leg), allPar)
    n <- sapply(pDat, nrow)
    if (addPoints) {
        add <- TRUE
        pcol <- c(rep(allPar$treated$points$col, n["treated"]),
                  rep(allPar$nontreated$points$col, n["nontreated"]))
        pch. <- c(rep(allPar$treated$points$pch, n["treated"]),
                  rep(allPar$nontreated$points$pch, n["nontreated"]))
        pargs <- allPar$common
        pargs$pch <- pch.
        pargs$col <- pcol
        pargs$x <- do.call("c", lapply(pDat, function(pi) pi[,which]))
        pargs$y <- do.call("c", lapply(pDat, function(pi) pi[,nameY]))
        do.call("plot", pargs)
    }
    pargs1 <- c(allPar$treated$lines, allPar$common, list(add=add))
    pargs1$x <- pDat$treated[, which]
    pargs1$y <- pDat$treated[,-(1:2)]
    do.call("matplot", pargs1)
    pargs0 <- c(allPar$nontreated$lines, list(add=TRUE))
    pargs0$x <- pDat$nontreated[, which]
    pargs0$y <- pDat$nontreated[,-(1:2)]
    do.call("matplot", pargs0)
    grid()
    if (!is.null(addToLegend)) 
        allPar$legend$legend = paste(allPar$legend$legend,
                                     " (", addToLegend[1], ")", sep = "")
    do.call("legend", allPar$legend)
    invisible()
}

