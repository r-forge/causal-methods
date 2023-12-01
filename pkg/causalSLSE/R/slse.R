## This documents include all R functions for SLSE models
##########################################################


## SLSE knots 
##############

## Splines builders

llSplines.slseModel <- function(object, ...)
{
    knots <- object$knots
    all <- lapply(1:length(object$nameX), function(i) {
        Uf <- .llSpline(object, i)
        nk <- length(knots[[i]]) + 1
        colnames(Uf) <- if (nk == 1)
                        {
                            object$nameX[i]
                        } else {
                            paste(object$nameX[i], "_", 1:nk, sep = "")
                        }
        Uf})
    names(all) <- object$nameX
    all <- do.call(cbind, all)
    all
}


.llSpline <- function (model, which) 
{
    X <- model.matrix(model)[,which]
    knots <- model$knots[[which]] 
    if (is.null(knots)) 
        return(as.matrix(X))
    nz <- X!=0
    n <- length(X)
    p <- length(knots) + 1
    X <- X[nz]
    Ui <- matrix(0, nrow = n, ncol = p)
    Ui[nz, 1] <- X * (X <= knots[1]) + knots[1] * (X > knots[1])
    Ui[nz, p] <- (X - knots[p - 1]) * (X > knots[p - 1])
    if (p >= 3)
    {
         for (j in 2:(p - 1))
        {
            Ui[nz, j] <- (X - knots[j-1]) * (X >= knots[j-1]) *
                (X <= knots[j]) + (knots[j] - knots[j-1]) * (X > knots[j])
        }
    }
    Ui
}


## Utilities to validate the knots and create knots

.firstknots <- function (knames, knots) 
{
    nX <- length(knames)
    if (is.null(knots))
    {
        k <- lapply(1:nX, function(i) NULL)
        names(k) <- knames
    } else {
        k <- knots
        if (!is.list(k)) 
            stop("The knots must be a list")
        if (length(k) == 0) 
            stop("The knots cannot be an empty list")
        uknames <- names(k)
        if (length(k) > nX)
        {
            stop("You provided too many sets of knots")
        } else if (length(k) == nX) {
            if (is.null(uknames))
            {
                names(k) <- knames
            } else {
                m <- match(knames, uknames)
                if (any(is.na(m))) 
                    stop("The names of the knots must match the name of the covariates")
                k <- k[m]
            }
        } else {
            if (is.null(uknames))
                stop("To manually define a subset of knots, the list must be named")
            m <- match(uknames, knames)
            if (any(is.na(m))) 
                stop("The names of the knots must match the name of the covariates")
            tmp <- as.list(rep(NA, nX))
            names(tmp) <- knames
            tmp[m] <- k
            k <- tmp
        }        
    }       
    k
}

.chkKnots <- function(x, knots)
{
    if (is.null(knots))
        return(knots)
    noNA <- !any(is.na(knots))
    isNum <- is.numeric(knots)
    noDup <- !any(duplicated(knots))
    inRange <- (min(knots, na.rm=TRUE)>min(x))&
        (max(knots, na.rm=TRUE)<max(x))
    chk <- noNA&isNum&noDup&inRange
    if (!all(chk))
        stop(paste("If you provide knots, they must be either numeric or NULL",
                   " they cannot contain NAs, duplicated ",
                   " and must be strictly in the range of X", sep=""))
    if(is.null(names(knots)))
        names(knots) <- paste("k", 1:length(knots), sep="")
    knots
}

setKnots <- function(x, sel=1:length(x), nbasis=function(n) n^0.3,
                     knots=NA, digits.names=4)
{
    x <- x[sel]
    if (is.null(knots) | is.numeric(knots))
        return(.chkKnots(x, knots))
    n <- length(x)
    p <- max(ceiling(nbasis(n)),2)
    prop.seq <- seq(from = 0, to = 1, length.out = p + 1)
    prop.seq <- prop.seq[-c(1, p + 1)]
    knots <- quantile(x, probs = prop.seq, type = 1, digits=digits.names)
    knots <- knots[!duplicated(knots)]
    b <- knots %in% range(x)
    if (any(b))
        knots <- knots[!b]
    if (length(knots)==0)
        return(NULL)
    .chkKnots(x, knots)
}

## Constructor of slseKnots objects

slseKnots <- function(form, data, X, nbasis = function(n) n^0.3, 
                      knots)
{
    if (missing(form))
    {
        if (missing(X))
            stop("without formula, the matrix of covariances must be provided")
    } else {
        if (missing(data))
            stop("When the formula is provided, the dataset must also be provided")          
        mf <- model.frame(form, data)
        mt <- attr(mf, "terms")
        X <- model.matrix(mt, mf)
        if (attr(terms(form), "intercept") == 1) 
            X <- X[, -1, drop = FALSE]
        X <- na.omit(X)
    }
    nameX <- colnames(X)
    select <- ifelse(missing(knots), "Default", "User Based")        
    if (missing(knots))
    {
        knots <- as.list(rep(NA, length(nameX)))
        names(knots) <- nameX
    }
    knots <- .firstknots(nameX, knots)
    knots <- lapply(1:ncol(X), function(i)
        setKnots(X[,i], 1:nrow(X), nbasis, knots[[i]]))
    names(knots) <- nameX
    attr(knots, "initSel") <- attr(knots, "curSel") <-
        list(select=select, crit="")
    class(knots) <- "slseKnots"
    knots
}

## slseKnots methods

update.slseKnots <- function(object, selKnots, ...)
{
    if (missing(selKnots))
        return(object)
    if (is.null(selKnots))
    {
        selKnots <- lapply(object, function(ki) NULL)
    }
    if (!is.list(selKnots))
        stop("selKnots must be a list")
    if (is.null(names(selKnots)))
        stop("selKnots must be a named list")
    nameSel <- names(selKnots)
    if (any(duplicated(nameSel)))
        stop("selKnots has duplicated names")
    if (!all(names(nameSel) %in% names(object)))
        stop("Some names in selKnots are not in the slseKnots object")
    k <- lapply(1:length(object), function(i) {
        if (names(object)[i] %in% nameSel)
        {
            wi <- selKnots[[names(object)[i]]]
            if (is.null(wi))
                return(NULL)
            if (any(is.na(wi)))
                stop("The knots selection list cannot contain NAs")
            if (!is.integer(wi))
                stop("The knots selection list can only contain integers")
            wi <- unique(wi)
            if (any(wi<1) | any(wi>length(object[[i]])))
                stop("Knot selection out of bound.")
            object[[i]][wi]
        } else {
            object[[i]]
        }})
        names(k) <- names(object)
    attr(k, "curSel") <- list(select="Manual selection", crit="")
    attr(k, "initSel") <- attr(object, "initSel")
    class(k) <- "slseKnots"
    k        
}

print.slseKnots <- function(x, header=c("None", "All", "Select"),
                            digits = max(3L, getOption("digits") - 3L), ...)
{
    header <- match.arg(header)
    if (header == "All")
    {
        cat("Semiparametric LSE Model: Selected knots\n")
        cat("****************************************\n")
    }
    if (header %in% c("All", "Select"))
    {
        cat("Selection method: ", attr(x, "curSel")$select, sep="")
        if (attr(x, "curSel")$crit != "")
            cat("-", attr(x, "curSel")$crit, sep = "")
    }
    cat("\n\n")   
    w <- sapply(x, is.null)
    selPW <- names(x)[!w]
    nonselPW <- names(x)[w]
    isApp <- if (length(selPW)) paste(selPW, collapse=", ", sep="") else "None"
    notApp <- if (length(nonselPW)) paste(nonselPW, collapse=", ", sep="") else "None"
    cat("Covariates with no knots:\n")
    cat("\t", notApp, "\n", sep = "")
    cat("\nCovariates with knots:\n")
    if (length(selPW) == 0)
    {
        cat("None\n")
    } else {
        for (ki in selPW)
        {
            cat(ki,":\n")
            pknots <- rbind(x[[ki]])
            rownames(pknots) <- c("Knots", "P-Value")[1:nrow(pknots)]
            print.default(pknots, quote=FALSE, digits=digits, ...)
            cat("\n")
        }
    }
    invisible()
}

## SLSE Model
#############

## Constructor of slseModel objects

slseModel <- function (form, data, nbasis = function(n) n^0.3, 
                       knots)
{
    mf <- model.frame(form, data)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    if (attr(terms(form), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    nameX <- colnames(X)
    nameY <- all.vars(form)[1]    
    Y <- data[,nameY]
    nameS <- "U"
    while (TRUE)
    {
        if (!(nameS %in% c(nameX, nameY)))
            break
        nameS <- paste(nameS, "U", collapse="", sep="")
    }
    nameS <- paste(nameS, ".", sep="")
    formY <- as.formula(paste(nameY, "~", nameS), env = .GlobalEnv)
    xlevels <- .getXlevels(mt,mf)
    na <- na.omit(cbind(Y,X))
    if (!is.null(attr(na, "omit")))
    {
        na <- attr(na, "omit")
        X <- X[-na,,drop=FALSE]
        data <- data[-na,,drop=FALSE]
    } else {
        na <- NULL
    }
    knots <- slseKnots(X=X, nbasis=nbasis, knots=knots)
    obj <- list(na=na, slseForm=formY, form=form, nameY=nameY,
                knots=knots, data=data, nameX=nameX, nameS=nameS, xlevels=xlevels)
    class(obj) <- "slseModel"
    obj
}

## slseModel methods

print.slseModel <- function(x, which=c("Model", "selKnots", "Pvalues"),
                            digits = max(3L, getOption("digits") - 3L), ...)
{
    which <- match.arg(which)
    if (which == "Model")
    {
        cat("Semiparametric LSE Model\n")
        cat("************************\n\n")
        cat("Number of observations: ", nrow(x$data), "\n")
        if (length(x$na))
            cat("Number of missing values: ", length(x$na), "\n")
        cat("Selection method: ", attr(x$knots,"curSel")$select, sep="")
        if (attr(x$knots,"curSel")$crit != "")
            cat("-",  attr(x$knots,"curSel")$crit, sep = "")
        cat("\n\n")
        cat("Covariates approximated by SLSE (num. of knots):\n")
        w <- sapply(x$knots, is.null)
        selPW <- x$nameX[!w]
        nK <- sapply(x$knots, length)[!w]
        isApp <- if (length(selPW))
                 {
                     selPW <- paste(selPW, "(", nK, ")", sep="")
                     paste(selPW, collapse=", ", sep="") 
                 } else {
                     "None"
                 }
        cat("\t", isApp, "\n", sep="")
        cat("Covariates not approximated by SLSE:\n")   
        w <- sapply(x$knots, is.null)
        nonselPW <- x$nameX[w]
        notApp <- if (length(nonselPW)) paste(nonselPW,collapse=", ",sep="") else "None"
        cat("\t", notApp, "\n", sep="")            
    } else if (which == "selKnots") {
        print(x$knots, header="All", digits=digits, ...)
    } else {
        if (is.null(x$selections))
        {
            cat("No p-values are available. You must apply a selection methods first.\n")
        } else {
            selType <- attr(x$knots,"curSel")$select
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

update.slseModel <- function(object, selType, selCrit="AIC",
                             selKnots, ...)
{
    knots <- if(!is.null(object$selections))
             {
                 object$selections$originalKnots
             } else {
                 object$knots
             }
    if (!missing(selKnots))
    {
        if (is.null(object$selections))
            object$selections$originalKnots <- object$knots
        object$knots <- update(knots, selKnots)
    } else if (!missing(selType)) {
        if (selType == "None")
        {
            if (!is.null(object$selections$originalKnots))
                object$knots <- object$selections$originalKnots
        } else {
            if (!is.null(object$selections[[selType]][[selCrit]]))
            {
                object$knots <- update(knots, object$selections[[selType]][[selCrit]])
                attr(object$knots, "curSel") <-  list(select=selType, crit=selCrit)
            } else {
                stop(paste("The ", selType, " method with ", selCrit, " criterion",
                           " is not available. Use selSLSE to add it to the model", sep=""))
            }
        }
    }
    object
}

model.matrix.slseModel <- function(object, ...)
{
    tt <- delete.response(terms(object$form))
    X <- model.matrix(tt, object$data, xlev=object$xlevels)    
    if (attr(terms(object$form), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    X
}

## Knots selection for SLSE models

selSLSE.slseModel <- function(model, selType=c("BLSE", "FLSE"),
                              selCrit = c("AIC", "BIC", "PVT"), 
                              pvalT = function(p) 1/log(p),
                              vcovType = c("HC0", "Classical", "HC1", "HC2", "HC3"),
                              reSelect=FALSE, ...)
{
    selCrit <- match.arg(selCrit)
    selType <- match.arg(selType)
    vcovType <- match.arg(vcovType)    
    if (!is.null(model$selections[[selType]][[selCrit]]) & !reSelect)
        return(update(model, selType=selType, selCrit=selCrit))
    model <- .selMod(model, selType, selCrit, pvalT, vcovType)
    model
}

## The slsePval methods
## Currently not exported and only used internally

print.slsePval <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    knots <- x$knots
    pval <- x$pval
    w <- sapply(knots, is.null)
    selPW <- names(knots)[!w]
    nonselPW <- names(knots)[w]
    isApp <- if (length(selPW)) 
                 paste(selPW, collapse = ", ", sep = "")
             else "None"
    notApp <- if (length(nonselPW)) 
                  paste(nonselPW, collapse = ", ", sep = "")
              else "None"
    cat("Covariates with no knots:\n")
    cat("\t", notApp, "\n", sep = "")
    cat("\nCovariates with knots:\n")
    if (length(selPW) == 0) {
        cat("None\n")
    } else {
        for (ki in selPW) {
            cat(ki, ":\n")
            pknots <- rbind(knots[[ki]], pval[[ki]])
            rownames(pknots) <- c("Knots", "P-Value")
            print.default(pknots, quote = FALSE, digits = digits, ...)
            cat("\n")
        }
    }
    invisible()
}

## Estimation method for slseModel and slseFit methods
## It constructs the slseFit object

estSLSE.slseModel <- function(model, selKnots, ...)
{
    model <- update(model, selKnots=selKnots)
    data <- model$data
    Sname <- all.vars(model$slseForm)[2]
    data[[Sname]] <- llSplines(model)
    form <- model$slseForm
    environment(form) <-  environment()    
    fit <- lm(form, data)
    obj <- list(LSE=fit, model=model)
    class(obj) <- "slseFit"
    obj
}

## slseFit methods

print.slseFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    cat("Selection method: ", attr(x$model$knots, "curSel")$select, "\n", sep="")
    if (attr(x$model$knots, "curSel")$crit != "")
    {
        cat("Criterion: ", attr(x$model$knots,"curSel")$crit, "\n\n", sep = "")
    } else {
        cat("\n")
    }
    cat("\n")
    print.default(coef(x$LSE), print.gap = 2L, digits=digits,
                  quote = FALSE)
    invisible()
}

summary.slseFit <- function (object, vcov.=vcovHC, ...) 
{
    s <- summary(object$LSE)
    v <- vcov.(object$LSE, complete=FALSE, ...)
    se <- sqrt(diag(v))
    t <- s$coefficients[,1]/se
    pv <- 2 * pnorm(-abs(t))
    s$coefficients[,2:4] <- cbind(se, t, pv)
    ans <- list(model=object$model, lseSum=s)
    class(ans) <- "summary.slseFit"
    ans
}

print.summary.slseFit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    cat("Selection method: ", attr(x$model$knots, "curSel")$select, "\n", sep="")
    if (attr(x$model$knots, "curSel")$crit != "")
    {
        cat("Criterion: ", attr(x$model$knots,"curSel")$crit, "\n\n", sep = "")
    } else {
        cat("\n")
    }
    q <- capture.output(print(x$lseSum, digits=digits,
                              signif.stars=signif.stars))
    
    q <- q[-(1:(grep("Residuals:", q)-1))]
    q <- q[1:grep("Multiple R-squared", q)]    
    q <- paste(q, collapse="\n")
    cat(q)
    cat("\n")
    invisible()
}

## The predict method

predict.slseFit <- function (object, interval = c("none", "confidence"),
                             se.fit = FALSE, 
                             newdata = NULL, level = 0.95, vcov. = vcovHC, ...) 
{
    interval <- match.arg(interval)
    model <- object$model
    if (!is.null(newdata)) 
        model$data <- newdata
    else newdata <- model$data
    nameS <- all.vars(object$model$slseForm)[2]
    newdata[[nameS]] <- llSplines(model)    
    tt <- terms(object$LSE)
    tt <- delete.response(tt)
    m <- model.frame(tt, newdata, xlev = object$LSE$xlevels)    
    X <- model.matrix(tt, m)
    b <- coef(object$LSE)
    naCoef <- !is.na(b)
    b <- na.omit(b)    
    pr <- c(X[, naCoef] %*% b)
    if (se.fit | interval == "confidence") {
        v <- vcov.(object$LSE, ...)
        se <- apply(X[, naCoef], 1, function(x) sqrt(c(t(x) %*% 
                                                       v %*% x)))
    }
    if (interval == "confidence") {
        crit <- qnorm(0.5 + level/2)
        pr <- cbind(fit = pr, lower = pr - crit * se, upper = pr + crit * se)
    }
    if (se.fit) 
        list(fit = pr, se.fit = se)
    else pr
}

# Utility function for the plot method

.prDatak <- function(object, varCov=NULL, conCov=NULL, sel=NULL, FUN = mean)
{
    model <- object$model
    if (is.null(sel))
        sel <- 1:nrow(model$data)
    data <- model$data[sel,]
    vnames <- all.vars(model$formX)
    if (!is.null(varCov))
    {
        if (!is.character(varCov)) 
            stop("varCov must be a character vector")
        if (!all(varCov  %in% vnames))
            stop("Some variables in varCov do not exist")
        varCov <- unique(varCov)
    } else {
        varCov <- character()
    }
    if (!is.null(conCov))
    {
        if (!is.list(conCov))
            stop("conCov must be a list")
        if (any(sapply(conCov, length) !=1))
            stop("The elements of conCov must have a length of 1")
        if (is.null(names(conCov)))
            stop("conCov must be a named list")
        if (!all(names(conCov) %in% vnames))
            stop("Some variables in conCov do not exist")
        conCov <- conCov[!duplicated(names(conCov))]
        if (length(varCov))
            if (any(names(conCov) %in% varCov))
                stop("You cannot have the same covariates in conCov and varCov")
        for (i in 1:length(conCov))
            data[,names(conCov)[i]] <- conCov[[i]]
    }
    
    model$data <- data
    X <- model.matrix(model)
    data[,colnames(X)] <- X
    f <- if (ncol(X)==1)
         {
             paste("~", colnames(X))
         } else {
             paste("~", paste(colnames(X), collapse="+"))
         }   
    formX <- formula(f, .GlobalEnv)
    xlevels <- list()
    funVar <- colnames(X)[!(colnames(X) %in% c(varCov, names(conCov)))]
    for (fi in funVar)
    {
        formi <- formula(paste("~",fi,"-1",sep=""), .GlobalEnv)
        vari <- all.vars(formi)
        dati <- data[,vari,drop=FALSE]
        chk <- vari %in% varCov
        if (any(!chk))
        {
            for (vi in vari[!chk])
                dati[,vi] <- FUN(dati[,vi])
        }
        data[,fi] <- model.matrix(formi, dati)
    }
    list(data=data, formX=formX, xlevels=xlevels)
}

## The plot method


plot.slseFit <- function (x, y, which = y, interval = c("none", "confidence"), 
                          level = 0.95, fixedCov = NULL,
                          vcov. = vcovHC, add = FALSE, 
                          addPoints = FALSE, FUN = mean, plot=TRUE, graphPar=list(), ...) 
{
    interval <- match.arg(interval)
    vnames <- all.vars(x$model$form)[-1]
    if (is.numeric(which)) 
        which <- vnames[which]
    if (!is.character(which) & length(which) != 1) 
        stop("which must be a character type")
    if (!(which %in% vnames)) 
        stop("which must be one of the names of the variables")
    if (addPoints) {
        Yp <- x$model$data[, x$model$nameY]
        Xp <- x$model$data[, which]
    }
    ind <- order(x$model$data[, which])
    x$model$data <- x$model$data[ind, ]
    x$model$formX <- formula(delete.response(terms(x$model$form)))
    obj <- .prDatak(x, which, fixedCov, NULL, FUN)
    x$model$data <- obj$data
    x$model$xlevels <- obj$xlevels
    x$model$form <- obj$formX
    pr <- predict(object=x, interval = interval, level = level, newdata=x$model$data,
        se.fit = FALSE, vcov. = vcov., ...)
    if (!plot)
    {
        pr <- cbind(x$model$data[c(x$model$nameY,which)], as.data.frame(pr))
        if (interval == "none")
            names(pr)[3] <- "fit"
        return(pr)
    }
    ylim <- if(addPoints) range(Yp) else  range(c(pr))
    common <- list(xlab=which, ylab=x$model$nameY,
                   ylim=ylim, xlim=range(x$model$data[, which]),
                   main=paste(x$model$nameY, " vs ", which, " using SLSE", sep = ""))
    allPar <- .slsePar(list(common=common))
    allPar <- .slsePar(graphPar, allPar)
    lpch <- if(addPoints)
            {
                allPar$points$pch
            } else {
                NA
            }
    if (addPoints) {
        add <- TRUE
        pargs <- allPar$common
        pargs$pch <- allPar$points$pch
        pargs$col <- allPar$points$col
        pargs$x <- Xp
        pargs$y <- Yp
        do.call("plot", pargs)
    }
    pargs <- c(allPar$lines, allPar$common, list(add=add))
    pargs$x <- x$model$data[, which]
    pargs$y <- pr
    do.call("matplot", pargs)
    grid()
    invisible()
}
.slsePar <- function(addPar, startPar=.initParSLSE())
{
    startPar$points <- .findrep(startPar$points, addPar$points)
    startPar$lines <- .findrep(startPar$lines, addPar$lines) 
    startPar$common <- .findrep(startPar$common, addPar$common)
    startPar
}

.initParSLSE <- function()
{
    points <- list(pch = 21, col = 1)
    lines <- list(col = 1, lty = c(1, 3, 3), lwd = 2, type='l')
    common <- list()
    list(lines=lines, common=common)
}

## Other utilities for printing slseFit results

extract.slseFit <- function (model, include.rsquared = TRUE, include.adjrs = TRUE, 
                             include.nobs = TRUE, include.fstatistic = FALSE,
                             include.rmse = FALSE, ...)
{ 
    extract(model$LSE, include.rsquared, include.adjrs, 
            include.nobs, include.fstatistic, include.rmse)
}

setMethod("extract", signature = className("slseFit", "causalSLSE"),
          definition = extract.slseFit)

