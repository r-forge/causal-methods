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

print.slseKnots <- function(x, which=c("selKnots", "Pvalues"),
                            header=c("None", "All", "Select"),
                            digits = max(3L, getOption("digits") - 3L), ...)
{
    which <- match.arg(which)
    header <- match.arg(header)
    if (which == "Pvalues" & is.null(attr(x,"pval")))
    {
        cat("No P-values: A selection method must first be applied to the model\n\n")
        return(invisible())
    }
    if (header == "All")
    {
        if (which == "selKnots")
        {
            cat("Semiparametric LSE Model: Selected knots\n")
            cat("****************************************\n")
        } else {
            cat("Semiparametric LSE Model: Initial knots and p-values\n")
            cat("****************************************************\n")
        }
    }
    if (header %in% c("All", "Select"))
    {
        if (which == "Pvalues")
        {
            cat("Initial Selection method: ", attr(x, "initSel")$select, "\n", sep="")
            if (attr(x, "initSel")$crit != "")            
                cat("Criterion: ", attr(x, "initSel")$crit, "\n", sep = "")
        } else {
            cat("Selection method: ", attr(x, "curSel")$select, "\n", sep="")
            if (attr(x, "curSel")$crit != "")
                cat("Criterion: ", attr(x, "curSel")$crit, "\n", sep = "")
        }
        cat("\n")
    }
    knots <- if(which == "selKnots") x else attr(x, "pval")$knots
    w <- sapply(knots, is.null)
    selPW <- names(knots)[!w]
    nonselPW <- names(knots)[w]
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
            pknots <- if (which == "selKnots")
                      {
                          rbind(knots[[ki]])
                      } else {
                          rbind(knots[[ki]], attr(x,"pval")$pval[[ki]])
                      }
            rownames(pknots) <- c("Knots", "P-Value")[1:nrow(pknots)]
            print.default(pknots, quote=FALSE, digits=digits, ...)
            cat("\n")
        }
    }
    invisible()
}

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

cslseKnots <- function (form, data, X, Z, nbasis = function(n) n^0.3, 
                        knots)
{
    select <- ifelse(missing(knots), "Default", "User Based")            
    groups <- c("treated", "nontreated")
    groupsi <- c(1,0)
    names(groupsi) <- groups
    if (missing(form))
    {
        if (missing(Z) | missing(X))
            stop("without formula, the matrix of covariances and treatment indicator must be provided")
    } else {
        if (missing(data))
            stop("When the formula is provided, the dataset must also be provided")
        tmp <- as.character(form)
        if (!grepl("\\|", tmp[3])) 
            stop("form must be of the type y~z|~x")
        tmp2 <- strsplit(tmp[3], "\\|")[[1]]
        formX <- as.formula(tmp2[2], env = .GlobalEnv)
        formY <- as.formula(paste(tmp[2], "~", tmp2[1], sep = ""))
        Z <- model.matrix(formY, data)
        if (attr(terms(formY), "intercept") == 1) 
            Z <- Z[, -1, drop = FALSE]
        if (ncol(Z) > 1) 
            stop("The right hand side must be a single vector of treatment indicator")
        if (!all(Z %in% c(0, 1))) 
            stop("The right hand side must be a binary variable")
        mf <- model.frame(formX, data)
        mt <- attr(mf, "terms")
        X <- model.matrix(mt, mf)
        if (attr(terms(formX), "intercept") == 1) 
            X <- X[, -1, drop = FALSE]
    }
    na <- na.omit(cbind(X,Z))
    if (!is.null(attr(na, "omit")))
        stop("Missing values in X and/or Z are not allowed")
    nameX <- colnames(X)
    if (missing(knots))
        knots <- list()
    if (is.null(knots))
    {
        knots <- lapply(groups, function(gi) NULL)
        names(knots) <- groups
    }
    if (!is.list(knots))
        stop("knots must be a list")
    if (length(knots))
    {
        if (is.null(names(knots)))
            stop("knots must be a named list")
        if (!all(names(knots) %in% groups))
            stop(paste("The names of knots must be in: ",
                       paste(groups, collapse=", ", sep=""), sep=""))
    }
    k <- list()
    for (gi in groups)
    {
        if (!(gi %in% names(knots)))
        {
            k[[gi]] <- slseKnots(X=X[Z==groupsi[gi], , drop=FALSE],
                                 nbasis=nbasis)
        } else {
            k[[gi]] <- slseKnots(X=X[Z==groupsi[gi], , drop=FALSE],
                                 nbasis=nbasis, knots=knots[[gi]])
        }
    }
    class(k) <- "cslseKnots"
    attr(k, "initSel") <- attr(k, "curSel") <-
        list(select=select, crit="")    
    k    
}

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

update.cslseKnots <- function(object, selKnots, ...)
{
    if (missing(selKnots))
        return(object)
    if (is.null(selKnots))
        selKnots <- list(treated=NULL, nontreated=NULL)
    if (!is.list(selKnots))
        stop("selKnots must be a list")
    if (is.null(names(selKnots)))
        stop("selKnots must be named")
    if (!all(names(selKnots) %in% names(object)))
        stop("The names of selKnots must be included in the group names of the knots")
    for (gi in names(selKnots))
    {
        object[[gi]] <- update(object[[gi]], selKnots[[gi]])        
        attr(object[[gi]], "curSel") <- list(select="Manual selection", crit="")
    }
    if (any(names(selKnots) %in% names(object)))
        attr(object, "curSel") <- list(select="Manual selection", crit="")    
    object
}

print.cslseKnots <- function(x, which=c("selKnots", "Pvalues"),
                             header = c("None", "All", "Select"),
                             digits = max(3L, getOption("digits") - 3L), ...)
{
    header <- match.arg(header)
    which <- match.arg(which)
    header2 <- "None"
    curSel <- sapply(x, function(gi) c(attr(gi, "curSel")$select,
                                       attr(gi, "curSel")$crit))
    initSel <- sapply(x, function(gi) c(attr(gi, "initSel")$select,
                                       attr(gi, "initSel")$crit))
    sameInit <- all(initSel[1, 1] == initSel[1,-1]) &  all(initSel[2,1] == initSel[2,-1])
    sameCur <- all(curSel[1, 1] == curSel[1,-1]) &  all(curSel[2,1] == curSel[2,-1])
    groups <- names(x)
    if (header == "All")
    {
        if (which == "selKnots")
        {
            cat("Causal Semiparametric LSE Model: Selected knots\n")
            cat("***********************************************\n")
        } else {            
            cat("Causal Semiparametric LSE Model: Initial knots and p-values\n")
            cat("***********************************************************\n")
        }
    }
    if (header %in% c("All", "Select"))
    {
        if (which == "selKnots")
        {
            if (sameCur)
            {
                cat("Selection method: ", curSel[1,1], "\n", sep="")
                if ( curSel[2,1] != "")
                    cat("Criterion: ",  curSel[2,1], "\n", sep = "")
            } else {
                header2 <- "Select"
            }
        } else {
            if (sameInit)
            {
                cat("Initial Selection method: ", initSel[1,1], "\n", sep="")
                if (initSel[2,1] != "")
                    cat("Initial criterion: ", initSel[2,1], "\n", sep = "")
            } else {
                header2 <- "Select"
            }
        }
        cat("\n")
    }
    for (gi in groups)
    {
        cat(gi, "\n")
        cat(paste(rep("*", nchar(gi)), collapse="", sep=""), "\n")
        print(x[[gi]], which=which, header=header2, digits=digits, ...)
    }
    invisible()
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

estSLSE <- function(model, ...)
    UseMethod("estSLSE")
       
estSLSE.cslseModel <- function(model, selKnots, ...)
{
    model$knots <- update(model$knots, selKnots)
    data <- model$data
    data$Xf1 <- multiSplines(model, "treated")
    data$Xf0 <- multiSplines(model, "nontreated")
    form <- model$formY
    environment(form) <-  environment()    
    fit <- lm(form, data)
    obj <- list(lm.out=fit, model=model)
    class(obj) <- "slseFit"
    obj
}

setKnots <- function(x, sel=1:length(x), nbasis=function(n) n^0.3,
                     knots=NA, digits.names=4)
{
    x <- x[sel]
    if (is.null(knots) | is.numeric(knots))
        return(.chkKnots(x, knots))
    n <- length(x)
    p <- floor(nbasis(n))
    if (p <= 1)
        return(NULL)
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

cslseModel <- function (form, data, nbasis = function(n) n^0.3, 
                        knots)
{
    tmp <- as.character(form)
    if (!grepl("\\|", tmp[3])) 
        stop("form must be of the type y~z|~x")
    tmp2 <- strsplit(tmp[3], "\\|")[[1]]
    formX <- as.formula(tmp2[2], env = .GlobalEnv)
    formY <- as.formula(paste(tmp[2], "~", tmp2[1], sep = ""))
    Z <- model.matrix(formY, data)
    Y <- model.frame(formY, data)[[1]]
    if (attr(terms(formY), "intercept") == 1) 
        Z <- Z[, -1, drop = FALSE]
    if (ncol(Z) > 1) 
        stop("The right hand side must be a single vector of treatment indicator")
    if (!all(Z %in% c(0, 1))) 
        stop("The right hand side must be a binary variable")
    formY <- as.formula(paste(tmp[2], "~factor(", tmp2[1], ")+Xf0+Xf1-1"), 
                        env = .GlobalEnv)
    mf <- model.frame(formX, data)
    mt <- attr(mf, "terms")
    X <- model.matrix(mt, mf)
    xlevels <- .getXlevels(mt,mf)
    if (attr(terms(formX), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    na <- na.omit(cbind(Y,Z,X))
    if (!is.null(attr(na, "omit")))
    {
        na <- attr(na, "omit")
        X <- X[-na,,drop=FALSE]
        Z <- Z[-na,,drop=FALSE]
        data <- data[-na,,drop=FALSE]
    } else {
        na <- NULL
    }
    nameX <- colnames(X)
    nameY <- all.vars(formY)[1]
    nameZ <- colnames(Z)
    knots <- cslseKnots(X=X, Z=Z, nbasis=nbasis, knots=knots)
    obj <- list(na=na, formY=formY, formX=formX, treated=nameZ, nameY=nameY,
                knots=knots, data=data, nameX=nameX, xlevels=xlevels)
    class(obj) <- "cslseModel"
    obj
}

print.cslseModel <- function(x, which=c("Model", "selKnots", "Pvalues"),
                             digits = max(3L, getOption("digits") - 3L), ...)
{
    Z <- x$data[[x$treated]]
    which <- match.arg(which)
    if (which == "Model")
    {
        cat("Causal Semiparametric LSE Model\n")
        cat("*******************************\n\n")
        cat("Number of treated: ", sum(Z), "\n")
        cat("Number of nontreated: ", sum(Z==0), "\n")
        if (length(x$na))
            cat("Number of missing values: ", length(x$na), "\n")
        cat("Selection method: ", attr(x$knots,"curSel")$select, "\n", sep="")
        if (attr(x$knots,"curSel")$crit != "")
            cat("Criterion: ",  attr(x$knots,"curSel")$crit, "\n\n", sep = "")
        cat("Confounders approximated by SLSE:\n")
        for (gi in names(x$knots))
        {
            w <- sapply(x$knots[[gi]], is.null)
            selPW <- x$nameX[!w]
            isApp <- if (length(selPW)) paste(selPW, collapse=", ", sep="") else "None"
            cat("\t", gi, ": ", isApp, "\n", sep="")
        }
        cat("Confounders not approximated by SLSE:\n")   
        for (gi in names(x$knots))
        {
            w <- sapply(x$knots[[gi]], is.null)
            nonselPW <- x$nameX[w]
            notApp <- if (length(nonselPW)) paste(nonselPW,collapse=", ",sep="") else "None"
            cat("\t", gi, ": ", notApp, "\n", sep="")            
        }
    } else {
        print(x$knots, which, header="All", digits=digits, ...)
    }
    invisible()
}

.splineMatrix <- function (model, which, group="treated",
                           selObs=c("group", "all")) 
{
    selObs <- match.arg(selObs)
    Z <- model$data[[model$treated]]
    if (selObs == "group")
    {
        id <- if (group=="treated") Z==1 else Z==0
    } else {
        id <- 1:nrow(model$data)
    }
    X <- model.matrix(model)[id,which]
    knots <- model$knots[[group]][[which]] 
    if (is.null(knots)) 
        return(as.matrix(X))
    nz <- X!=0
    n <- length(X)
    p <- length(knots) + 1
    X <- X[nz]
    Xfi <- matrix(0, nrow = n, ncol = p)
    Xfi[nz, 1] <- X * (X <= knots[1]) + knots[1] * (X > knots[1])
    Xfi[nz, p] <- (X - knots[p - 1]) * (X > knots[p - 1])
    if (p >= 3)
    {
        for (j in 2:(p - 1))
        {
            Xfi[nz, j] <- (X - knots[j-1]) * (X >= knots[j-1]) *
                (X <= knots[j]) + (knots[j] - knots[j-1]) * (X > knots[j])
        }
    }
    Xfi
}

multiSplines <- function (model, group="treated", selObs=c("group", "all"))
{
    selObs <- match.arg(selObs)
    Z <- model$data[[model$treated]]
    knots <- model$knots[[group]]
    if (selObs == "group")
    {
        id <- if(group=="treated") Z==1 else Z==0
    } else {
        id <- 1:nrow(model$data)
    }
    all <- lapply(1:length(model$nameX), function(i) {
        ans <- .splineMatrix(model, i, group, selObs)
        nk <- length(knots[[i]]) + 1
        Xf <- matrix(0, nrow(model$data), nk)
        Xf[id,] <- ans
        colnames(Xf) <- if (nk == 1)
                         {
                             model$nameX[i]
                         } else {
                             paste(model$nameX[i], "_", 1:nk, sep = "")
                         }
        Xf})
    names(all) <- model$nameX
    cnames <- lapply(all, colnames)
    names(cnames) <- names(all)
    all <- do.call(cbind, all)
    attr(all, "p") <- sapply(knots, length) + 1
    attr(all, "colnames") <- cnames
    all
}

selSLSE <- function(model, ...)
    UseMethod("selSLSE")
       
selSLSE.cslseModel <- function(model, selType=c("BLSE", "FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"), 
                               pvalT = function(p) 1/log(p), vcov.=vcovHC, ...)
{
    selCrit <- match.arg(selCrit)
    selType <- match.arg(selType)
    critFct <- if (selCrit == "PVT") {
                   .selPVT
               } else {
                   .selIC
               }
    if (all(sapply(model$knots$treated, function(i) is.null(i))) &
        all(sapply(model$knots$nontreated, function(i) is.null(i))))
    {
        warning("No selection needed: the number of knots is 0 for all confounders")
    } else {
        pval <- pvalSLSE(model, selType, vcov., ...)
        model <- critFct(model, pval, pvalT, selCrit)
        attr(model$knots$treated, "pval") <- pval$treated
        attr(model$knots$nontreated, "pval") <- pval$nontreated
        attr(model$knots$treated, "curSel") <- attr(model$knots$nontreated, "curSel") <-
            attr(model$knots, "curSel") <-  list(select=selType, crit=selCrit)
    }
    model
}

.testKnots <- function(fit, model, whichK, whichX, group, vcov)
{
    wX <- ifelse(group=="treated", "Xf1", "Xf0")
    if (whichX > length(model$knots[[group]]))
        stop("whichX exceeds the number of covariates")
    if (any(whichK > length(model$knots[[group]][[whichX]])))
        stop("whichK exceeds the number of knots")   
    if (is.null(model$knots[[group]][[whichX]]))
        return(NA)
    b <- coef(fit)
    b <- na.omit(b)
    sapply(whichK, function(wi) {
        nX <- names(model$knots[[group]])[[whichX]]
        t <- c(paste(wX, nX, "_", wi, sep = ""),
               paste(wX, nX, "_", wi+1, sep = ""))
        c1 <- which(names(b) == t[1])
        c2 <- which(names(b) == t[2])
        if (length(c(c1, c2)) < 2) 
            return(NA)
        s2 <- vcov[c1, c1] + vcov[c2, c2] - 2 * vcov[c1, c2]
        ans <- 1 - pf((b[c1] - b[c2])^2/s2, 1, fit$df)
        names(ans) <- NULL
        ans})
}


pvalSLSE <- function(model, ...)
{
    UseMethod("pvalSLSE")
}

pvalSLSE.cslseModel <- function(model, method=c("BLSE", "FLSE", "SLSE"),
                                vcov.=vcovHC, ...)
{
    method <- match.arg(method)
    if (method == "SLSE")
    {
        pv <- list(pval0=lapply(length(model$knots$nontreated), function(i) NULL),
                   pval1=lapply(length(model$knots$treated), function(i) NULL))
    } else {
        pv <- if (method=="BLSE")  .getPvalB(model, vcov., ...)
              else .getPvalF(model, vcov., ...)
    }
    ans <- list(treated=list(pval=pv$pval1, knots=model$knots$treated),
                nontreated=list(pval=pv$pval0, knots=model$knots$nontreated))
    class(ans$treated) <-  class(ans$nontreated) <- "slsePval"
    ans
}

.getPvalF <- function (model, vcov.=vcovHC, ...) 
{
    selK <- list()
    selK$nontreated <- lapply(model$knots$nontreated, function(i) NULL)
    selK$treated <- lapply(model$knots$treated, function(i) NULL)
    selFct <- function(p)
    {
        sel <- if (p==1) {
                   list(1)
               } else if (p == 2) {
                   list(1:2)
               } else {
                   c(list(1:2), lapply(1:(p-2),
                                       function(j) (0:2)+j), list(c(p-1,p)))
               }
    }
    pvali <- function(i, group)
    {
        p <- length(model$knots[[group]][[i]])
        if (p==0)
            return(NA)
        sel <- selFct(p)
        c(sapply(1:length(sel), function(j) {
            selj <- as.integer(sel[[j]])
            selK[[group]][[i]] <- selj
            res <- estSLSE(model, selKnots=selK)
            v <- vcov.(res$lm.out, ...)
            if (length(selj) == 1L)
            {
                .testKnots(res$lm.out, res$model, 1L, i, group, v)
            } else if (length(sel)==1 & length(selj)==2) {
                .testKnots(res$lm.out, res$model, 1L:2L, i, group, v)
            } else {
                whichK <- ifelse(length(selj)==2 & selj[1]==1L, 1, 2)
                .testKnots(res$lm.out, res$model, whichK, i, group, v)
            }}))
    }
    pval0 <- lapply(1:length(model$knots$nontreated), function(i) pvali(i, "nontreated"))
    pval1 <- lapply(1:length(model$knots$treated), function(i) pvali(i, "treated"))
    names(pval0) <- names(pval1) <- model$nameX
    list(pval0=pval0, pval1=pval1)
}

.getPvalB <- function (model, vcov.=vcovHC, ...) 
{
    data2 <- model$data
    data2$Xf1 <- multiSplines(model, "treated")
    data2$Xf0 <- multiSplines(model, "nontreated")
    form <- model$formY
    environment(form) <- environment()
    fit <- lm(form, data2)
    v <- vcov.(fit, ...)
    pval0 <- lapply(1:length(model$knots$nontreated), function(i) {
        ki <- length(model$knots$nontreated[[i]])
        if (ki==0) NA else  .testKnots(fit, model, 1:ki, i , "nontreated", v)
    })
    pval1 <- lapply(1:length(model$knots$treated), function(i) {
        ki <- length(model$knots$treated[[i]])
        if (ki==0) NA else  .testKnots(fit, model, 1:ki, i , "treated", v)
    })
    names(pval0) <- names(pval1) <- model$nameX
    list(pval0 = pval0, pval1 = pval1)
}

.selIC <- function (model, pvalRes, pvalT = NULL, crit)
{
    pval <- c(do.call("c", pvalRes$nontreated$pval),
              do.call("c", pvalRes$treated$pval))
    pval_sort <- sort(pval)    
    q <- length(pval)
    selK <- list()
    selK$treated <- lapply(pvalRes$treated$pval, function(i) NULL)  
    selK$nontreated <- lapply(pvalRes$nontreated$pval, function(i) NULL) 
    res0 <- estSLSE(model, selK)
    icV <- ic_seq0 <- get(crit)(res0$lm.out)   
    for (i in 1:q) {
        selK$nontreated <- lapply(1:length(pvalRes$nontreated$pval), function(j)
        {
            w <- which(pvalRes$nontreated$pval[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        selK$treated <- lapply(1:length(pvalRes$treated$pval), function(j)
        {
            w <- which(pvalRes$treated$pval[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        names(selK$nontreated) <- names(pvalRes$nontreated$pval)
        names(selK$treated) <- names(pvalRes$treated$pval)
        res1 <- estSLSE(model, selK)
        ic_seq1 <- get(crit)(res1$lm.out)
        icV <- c(icV, ic_seq1)
        if (ic_seq1 < ic_seq0) {
            ic_seq0 <- ic_seq1
            res0 <- res1
        }
    }
    model <- res0$model
    class(model$knots$treated) <- class(model$knots$nontreated) <- "slseKnots"
    model
}

.selPVT <- function (model, pvalRes, pvalT = function(p) 1/log(p),
                     crit=NULL, ...)
{
    pval <- c(do.call("c", pvalRes$nontreated$pval),
              do.call("c", pvalRes$treated$pval))
    n <- nrow(model$data)
    q <- length(pval)
    p0 <-  sapply(model$knots$nontreated, function(ki) length(ki) + 1)
    p1 <-  sapply(model$knots$treated, function(ki) length(ki) + 1)
    p <- mean(c(p0, p1))
    crit <- pvalT(p)
    w0 <- lapply(1:length(model$knots$nontreated), function(i)
    {
        w <- which(pvalRes$nontreated$pval[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    w1 <- lapply(1:length(model$knots$treated), function(i)
    {
        w <- which(pvalRes$treated$pval[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    names(w0) <- names(pvalRes$nontreated$pval)
    names(w1) <- names(pvalRes$treated$pval)
    model$knots$nontreated <- update(model$knots$nontreated, w0)
    model$knots$treated <- update(model$knots$treated, w1)
    model
}

model.matrix.cslseModel <- function(object, ...)
{
    X <- model.matrix(object$formX, object$data, xlev=object$xlevels)    
    if (attr(terms(object$formX), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    X
}

.causali <- function(Z, beta, vcov, X0, X1,
                     beta0, beta1, causal, e)
{
    id <- switch(causal,
                 ACE = rep(TRUE, length(Z)),
                 ACT = Z==1,
                 ACN = Z==0)
    n <- sum(id)    
    p <- switch(causal,
                ACE = ncol(X0)+ncol(X1)+2,
                ACT = ncol(X1)+1,
                ACN = ncol(X0)+1)
    dfadj <- n/(n-p)
    X <- cbind(1-Z, Z, X0*(1-Z), X1*Z)
    SigmaX <- crossprod(X)
    X0 <- X0[id,, drop = FALSE]
    X1 <- X1[id,, drop = FALSE]
    Xbar0 <- colMeans(X0)
    Xbar1 <- colMeans(X1)
    est <- c(beta[2] - beta[1] + sum(beta1*Xbar1) - sum(beta0*Xbar0))
    vcovXf0 <- cov(X0)
    vcovXf1 <- cov(X1)
    vcovXf01 <- cov(X0, X1)
    Dvec <- c(-1, 1, -Xbar0, Xbar1)
    e <- e[id]
    X <- X[id,,drop=FALSE]
    X1 <- c(scale(X1, scale=FALSE)%*%beta1)
    X0 <- c(scale(X0, scale=FALSE)%*%beta0)
    tmp2 <- c(e*X%*%solve(SigmaX, Dvec))
    addT <- 2*mean(X1*tmp2) - 2*mean(X0*tmp2)
    se <- (sum(Dvec*c(crossprod(vcov, Dvec))) +
           sum(beta0*c(crossprod(vcovXf0, beta0)))/n +
           sum(beta1*c(crossprod(vcovXf1, beta1)))/n -
           2 * c(beta0 %*% vcovXf01 %*% beta1)/n +
           2*mean(X1*tmp2)*dfadj -
           2*mean(X0*tmp2)*dfadj)^0.5
    ans <- c(est, se)
    names(ans) <- c("est","se")
    ans
}

.causal <- function(model, fit, vcov, X0, X1,
                    causal=c("ACE", "ACT", "ACN", "ALL"))
{
    causal <- match.arg(causal)
    Z <- model$data[[model$treated]]
    p0 <- ncol(X0)
    p1 <- ncol(X1)
    beta <- coef(fit)
    beta0 <- beta[3:(p0 + 2)]
    beta1 <- beta[(p0 + 3):(p0 + p1 + 2)]
    notNA0 <- !is.na(beta0)
    notNA1 <- !is.na(beta1)
    beta0 <- na.omit(beta0)
    beta1 <- na.omit(beta1)
    if (length(notNA0))  X0 <- X0[, notNA0, drop=FALSE]
    if (length(notNA1))  X1 <- X1[, notNA1, drop=FALSE]
    if (causal == "ALL") causal <- c("ACE","ACT","ACN")                
    ans <- lapply(causal, function(ci)
        .causali(Z, beta, vcov, X0, X1, beta0, beta1, ci, residuals(fit)))
    names(ans) <- causal
    ans
}

causalSLSE <- function(object, ...)
{
    UseMethod("causalSLSE")
}

causalSLSE.cslseModel <- function(object,
                                  selType=c("SLSE","BLSE","FLSE"),
                                  selCrit = c("AIC", "BIC", "PVT"),
                                  causal = c("ALL","ACT","ACE","ACN"),
                                  pvalT = function(p) 1/log(p),
                                  vcov.=vcovHC, ...)
{
    selType <- match.arg(selType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    if (selType != "SLSE")
        object <- selSLSE(object, selType, selCrit, pvalT, vcov., ...)
    res <- estSLSE(object)
    beta <- coef(res$lm.out)
    v <- vcov.(res$lm.out, ...)
    se.beta <- sqrt(diag(v))
    if (any(is.na(beta))) 
        warning(paste("\nThe regression is multicollinear.",
                      " The result may not be valid:", 
                      "\nThe following variables produced NA's\n",
                      paste(names(beta)[is.na(beta)], 
                            collapse = ", "), "\n", sep = ""))
    X0 <- multiSplines(object, "nontreated", "all")
    X1 <- multiSplines(object, "treated", "all")
    ans <- .causal(object, res$lm.out, v, X0, X1, causal)
    ans <- c(ans, 
             list(beta = beta, se.beta = se.beta, lm.out = res$lm.out, 
                  model=object))
    class(ans) <- c("cslse", "slseFit")
    ans    
}

cslseSE <- function(object, vcov.=vcovHC, ...)    
{
    if(!inherits(object, "cslse"))
        stop("object must be of class cslse")
    Z <- object$model$data[,object$model$treate]
    causal <- c("ACE","ACT","ACN")
    causal <- causal[which(causal %in% names(object))]
    if (length(causal)>1)
        causal <- "ALL"
    X0 <- multiSplines(object$model, "nontreated", "all")
    X1 <- multiSplines(object$model, "treated", "all")
    p0 <- ncol(X0)
    p1 <- ncol(X1)
    beta <- coef(object$lm.out)
    beta0 <- beta[3:(p0 + 2)]
    beta1 <- beta[(p0 + 3):(p0 + p1 + 2)]
    notNA0 <- !is.na(beta0)
    notNA1 <- !is.na(beta1)
    beta0 <- na.omit(beta0)
    beta1 <- na.omit(beta1)
    if (length(notNA0))  X0 <- X0[, notNA0, drop=FALSE]
    if (length(notNA1))  X1 <- X1[, notNA1, drop=FALSE]
    if ("ALL" %in% causal)
        causal <- c("ACE","ACT","ACN")
    vcov <- vcov.(object$lm.out, ...)
    se <- sapply(causal, function(ci)
        .causali(Z, beta, vcov, X0, X1,
                 beta0, beta1, ci, residuals(object$lm.out))["se"])
    names(se) <- causal   
    se
}

causalSLSE.slseFit <- function(object, causal = c("ALL","ACT","ACE","ACN"),
                               vcov.=vcovHC, ...)
{
    causal <- match.arg(causal)
    beta <- coef(object$lm.out)
    v <- vcov.(object$lm.out, ...)
    se.beta <- sqrt(diag(v))
    if (any(is.na(beta))) 
        warning(paste("\nThe regression is multicollinear.",
                      " The result may not be valid:", 
                      "\nThe following variables produced NA's\n",
                      paste(names(beta)[is.na(beta)], 
                            collapse = ", "), "\n", sep = ""))
    X0 <- multiSplines(object$model, "nontreated", "all")
    X1 <- multiSplines(object$model, "treated", "all")
    ans <- .causal(object$model, object$lm.out, v, X0, X1, causal)
    ans <- c(ans, 
             list(beta = beta, se.beta = se.beta, lm.out = object$lm.out, 
                  model=object$model))
    class(ans) <- c("cslse", "slseFit")
    ans    
}

print.cslse <- function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    curSel <- sapply(x$model$knots, function(gi) c(attr(gi, "curSel")$select, 
                                                   attr(gi, "curSel")$crit))
    sameCur <- all(curSel[1, 1] == curSel[1, -1]) & all(curSel[2, 1] == curSel[2, -1])
    if (sameCur)
    {
        cat("Selection method: ", attr(x$model$knots[[1]],"curSel")$select, "\n", sep="")
        if (attr(x$model$knots[[1]],"curSel")$crit != "")
        {
            cat("Criterion: ", attr(x$model$knots[[1]],"curSel")$crit, "\n\n", sep = "")
        } else {
            cat("\n")
        }
    } else {
        for (gi in names(x$model$knots))
        {
            cat("Selection method (", gi, "): ",
                attr(x$model$knots[[gi]],"curSel")$select, "\n", sep="")
            if (attr(x$model$knots[[gi]],"curSel")$crit != "")
            {
                cat("Criterion (", gi, "): ",
                    attr(x$model$knots[[gi]],"curSel")$crit, "\n", sep = "")
            } else {
                cat("\n")
            }
        }
    }
    cat("\n")
    for (causal in c("ACE","ACT","ACN"))
    {
        if (!is.null(x[[causal]]))
            cat(causal, " = ", format(x[[causal]]["est"], digits=digits),
                "\n", sep="")
    }
}

causalSLSE.formula <- function(object, data, nbasis=function(n) n^0.3,
                               knots, 
                               selType=c("SLSE","BLSE","FLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, ...)
{
    model <- cslseModel(object, data, nbasis,  knots)
    selType <- match.arg(selType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    causalSLSE(object=model, selType = selType, selCrit = selCrit,
               causal = causal, pvalT =  pvalT,
               vcov. = vcov., ...)    
}

print.slseFit <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    curSel <- sapply(x$model$knots, function(gi) c(attr(gi, "curSel")$select, 
                                                   attr(gi, "curSel")$crit))
    sameCur <- all(curSel[1, 1] == curSel[1, -1]) & all(curSel[2, 1] == curSel[2, -1])
    if (sameCur)
    {
        cat("Selection method: ", attr(x$model$knots[[1]],"curSel")$select, "\n", sep="")
        if (attr(x$model$knots[[1]],"curSel")$crit != "")
        {
            cat("Criterion: ", attr(x$model$knots[[1]],"curSel")$crit, "\n\n", sep = "")
        } else {
            cat("\n")
        }
    } else {
        for (gi in names(x$model$knots))
        {
            cat("Selection method (", gi, "): ",
                attr(x$model$knots[[gi]],"curSel")$select, "\n", sep="")
            if (attr(x$model$knots[[gi]],"curSel")$crit != "")
            {
                cat("Criterion (", gi, "): ",
                    attr(x$model$knots[[gi]],"curSel")$crit, "\n", sep = "")
            } else {
                cat("\n")
            }
        }
    }
    cat("\n")
    print.default(coef(x$lm.out), print.gap = 2L, digits=digits,
                  quote = FALSE)
    invisible()
}

summary.slseFit <- function (object, vcov.=vcovHC, ...) 
{
    beta <- coef(object$lm.out)
    v <- vcov.(object$lm.out, ...)
    se <- sqrt(diag(v))
    t <- beta/se
    pv <- 2 * pnorm(-abs(t))
    coef <- cbind(beta, se, t, pv)
    colnames(coef) <- c("Estimate", "Std. Error", "t value", 
                        "Pr(>|t|)")
    f <- fitted(object$lm.out)
    ess <- sum((f-mean(f))^2)
    rss <- sum(residuals(object$lm.out)^2)
    r2 <- ess/(ess+rss)
    r2adj <- 1-(1-r2)*(nobs(object$lm.out)-1)/object$lm.out$df.residual
    ans <- list(coefficients = coef, model=object$model,
                r.squared=r2, adj.r.squared=r2adj)
    class(ans) <- "summary.slseFit"
    ans
}

print.summary.slseFit <- function(x, digits = max(3L, getOption("digits") - 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    curSel <- sapply(x$model$knots, function(gi) c(attr(gi, "curSel")$select, 
                                                   attr(gi, "curSel")$crit))
    sameCur <- all(curSel[1, 1] == curSel[1, -1]) & all(curSel[2, 1] == curSel[2, -1])
    if (sameCur)
    {
        cat("Selection method: ", attr(x$model$knots[[1]],"curSel")$select, "\n", sep="")
        if (attr(x$model$knots[[1]],"curSel")$crit != "")
        {
            cat("Criterion: ", attr(x$model$knots[[1]],"curSel")$crit, "\n\n", sep = "")
        } else {
            cat("\n")
        }
    } else {
        for (gi in names(x$model$knots))
        {
            cat("Selection method (", gi, "): ",
                attr(x$model$knots[[gi]],"curSel")$select, "\n", sep="")
            if (attr(x$model$knots[[gi]],"curSel")$crit != "")
            {
                cat("Criterion (", gi, "): ",
                    attr(x$model$knots[[gi]],"curSel")$crit, "\n", sep = "")
            } else {
                cat("\n")
            }
        }
    }
    cat("\n")
    printCoefmat(x$coefficients, na.print = "NA", digits = digits,
                 signif.stars = signif.stars, ...)
    cat("\nMultiple R-squared: ", formatC(x$r.squared, digits=digits, ...))
    cat(",\tAdjusted R-squared: ", formatC(x$adj.r.squared,
                                           digits=digits, ...), "\n")
    invisible()
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
    t <- na.omit(object$beta)/object$se.beta
    pv <- 2 * pnorm(-abs(t))
    beta <- cbind(na.omit(object$beta), object$se.beta, t, pv)
    colnames(beta) <- c("Estimate", "Std. Error", "t value", 
        "Pr(>|t|)")
    ans <- list(causal = est, beta = beta, 
                knots = object$model$knots,
                covNames = names(object$model$knots$treated))
    class(ans) <- "summary.cslse"
    ans
}

print.summary.cslse <- function (x, digits = max(3L, getOption("digits") - 3L), 
                                signif.stars = getOption("show.signif.stars"), 
                                beta = FALSE, knots = FALSE, ...) 
{
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    curSel <- sapply(x$knots, function(gi) c(attr(gi, "curSel")$select, 
                                             attr(gi, "curSel")$crit))
    sameCur <- all(curSel[1, 1] == curSel[1, -1]) & all(curSel[2, 1] == curSel[2, -1])
    if (sameCur)
    {
        cat("Selection method: ", attr(x$knots[[1]],"curSel")$select, "\n", sep="")
        if (attr(x$knots[[1]],"curSel")$crit != "")
        {
            cat("Criterion: ", attr(x$knots[[1]],"curSel")$crit, "\n\n", sep = "")
        } else {
            cat("\n")
        }
    } else {
        for (gi in names(x$knots))
        {
            cat("Selection method (", gi, "): ",
                attr(x$knots[[gi]],"curSel")$select, "\n", sep="")
            if (attr(x$knots[[gi]],"curSel")$crit != "")
            {
                cat("Criterion (", gi, "): ",
                    attr(x$knots[[gi]],"curSel")$crit, "\n", sep = "")
            } else {
                cat("\n")
            }
        }
    }
    cat("\n")
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)
    if (beta) {
        cat("\nBasis function coefficients\n")
        cat("*****************************\n")
        printCoefmat(x$beta, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    if (knots) {
        cat("\nNumber of selected knots per confounder\n")
        cat("****************************************\n")
        cat("Treated:\n")
        print.default(format(sapply(x$knots$treated, length)), print.gap = 2L, 
            quote = FALSE)
        cat("Nontreated:\n")
        print.default(format(sapply(x$knots$nontreated, length)), print.gap = 2L, 
            quote = FALSE)
        cat("\n")
    }
}

extract.cslse <- function (model, include.nobs = TRUE,
                                include.nknots = TRUE,
                                include.numcov = TRUE, include.rsquared = TRUE,
                                include.adjrs = TRUE, 
                                which=c("ALL","ACE","ACT","ACN","ACE-ACT",
                                        "ACE-ACN","ACT-ACN"), ...) 
{
    which <- match.arg(which)
    type <- c("ACE","ACT","ACN")
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
        rs1 <- length(unlist(model$model$knots$nontreated))
        rs2 <- length(unlist(model$model$knots$treated))        
        gof <- c(gof, rs1, rs2)
        gof.names <- c(gof.names, "Num. knots (Nontreated)", "Num. knots (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
   }
    if (isTRUE(include.numcov)) {
        rs3 <- length(model$model$nameX)
        gof <- c(gof, rs3)
        gof.names <- c(gof.names, "Num. confounders")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    if (isTRUE(include.nobs)) {
        Z <- model$model$data[[model$model$treated]]
        n1 <- sum(Z)
        n0 <- sum(1-Z)
        gof <- c(gof, n0, n1)
        gof.names <- c(gof.names, "Num. obs. (Nontreated)", "Num. obs. (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
    }
    if (isTRUE(include.rsquared)) {
        f <- fitted(model$lm.out)
        ess <- sum((f-mean(f))^2)
        rss <- sum(residuals(model$lm.out)^2)
        R2 <- ess/(ess+rss)
        gof <- c(gof, R2)
        gof.names <- c(gof.names, "R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
    }
    if (isTRUE(include.adjrs)) {
        f <- fitted(model$lm.out)
        ess <- sum((f-mean(f))^2)
        rss <- sum(residuals(model$lm.out)^2)       
        R2 <- ess/(ess+rss)
        R2adj <- 1-(1-R2)*(nobs(model$lm.out)-1)/model$lm.out$df.residual        
        gof <- c(gof, R2adj)
        gof.names <- c(gof.names, "Adj. R$^2$")
        gof.decimal <- c(gof.decimal, TRUE)
    }   
    tr <- createTexreg(coef.names = names(co), coef = co, se = se, 
                       pvalues = pval, gof.names = gof.names, gof = gof,
                       gof.decimal = gof.decimal)
    return(tr)
}

setMethod("extract", signature = className("cslse", "causalSLSE"),
          definition = extract.cslse)


predict.slseFit <- function (object, interval = c("none", "confidence"),
                             se.fit = FALSE, 
                             newdata = NULL, level = 0.95, vcov. = vcovHC, ...) 
{
    interval <- match.arg(interval)
    model <- object$model
    if (!is.null(newdata)) 
        model$data <- newdata
    else newdata <- model$data
    Z <- model$data[[model$treated]]
    if (is.null(Z)) 
        stop("newdata must contain a treatment assignment variable")
    newdata$Xf1 <- multiSplines(model, "treated", "all")
    newdata$Xf0 <- multiSplines(model, "nontreated", "all")
    newdata$Xf1[Z==0,] <- 0
    newdata$Xf0[Z==1,] <- 0
    
    tt <- terms(object$lm.out)
    tt <- delete.response(tt)
    m <- model.frame(tt, newdata, xlev = object$lm.out$xlevels)    
    X <- model.matrix(tt, m)
    b <- coef(object$lm.out)
    naCoef <- !is.na(b)
    b <- na.omit(b)    
    pr0 <- c(X[Z == 0, naCoef] %*% b)
    pr1 <- c(X[Z == 1, naCoef] %*% b)
    if (se.fit | interval == "confidence") {
        v <- vcov.(object$lm.out, ...)
        se0 <- apply(X[Z == 0, naCoef], 1, function(x) sqrt(c(t(x) %*% 
                                                              v %*% x)))
        se1 <- apply(X[Z == 1, naCoef], 1, function(x) sqrt(c(t(x) %*% 
                                                              v %*% x)))
    }
    if (interval == "confidence") {
        crit <- qnorm(0.5 + level/2)
        pr0 <- cbind(fit = pr0, lower = pr0 - crit * se0, upper = pr0 + 
            crit * se0)
        pr1 <- cbind(fit = pr1, lower = pr1 - crit * se1, upper = pr1 + 
            crit * se1)
    }
    if (se.fit) 
        list(treated = list(fit = pr1, se.fit = se1), nontreated = list(fit = pr0, 
            se.fit = se0))
    else list(treated = pr1, nontreated = pr0)
}


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

.initPar <- function()
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

.cslsePar <- function(addPar, startPar=.initPar())
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

plot.slseFit <- function (x, y, which = y, interval = c("none", "confidence"), 
                   counterfactual = FALSE, level = 0.95, fixedCov0 = NULL,
                   fixedCov1 = fixedCov0,  vcov. = vcovHC, add = FALSE, addToLegend = NULL, 
                   addPoints = FALSE, FUN = mean, plot=TRUE, graphPar=list(), ...) 
{
    interval <- match.arg(interval)
    vnames <- all.vars(x$model$formX)
    treat <- x$model$treated
    if (is.numeric(which)) 
        which <- vnames[which]
    if (!is.character(which) & length(which) != 1) 
        stop("which must be a character type")
    if (!(which %in% vnames)) 
        stop("which must be one of the names of the variables")
    if (addPoints) {
        Yp <- x$model$data[, x$model$nameY]
        Xp <- x$model$data[, which]
        Zp <- x$model$data[, treat]
    }
    ind <- order(x$model$data[, which])
    x$model$data <- x$model$data[ind, ]
    Z <-  x$model$data[, treat]
    obj1 <- .prDatak(x, which, fixedCov1, Z==1, FUN)
    obj0 <- .prDatak(x, which, fixedCov0, Z==0, FUN)    
    x$model$data <- rbind(obj1$data, obj0$data)
    Z <-  x$model$data[, treat]
    x$model$xlevels <- obj1$xlevels
    x$model$formX <- obj1$formX
    res <- predict(object=x, interval = interval, level = level, newdata=x$model$data,
        se.fit = FALSE, vcov. = vcov., counterfactual = counterfactual, 
        ...)
    pr0 <- res$nontreated
    pr1 <- res$treated
    if (!plot)
    {
        pr0 <- cbind(x$model$data[which][Z==0,], as.data.frame(pr0))
        pr1 <- cbind(x$model$data[which][Z==1,], as.data.frame(pr1))
        if (interval == "none")
            names(pr0)[2] <- names(pr1)[2] <- "fit"
        return(list(treated=pr1, nontreated=pr0))
    }
    ylim <- if(addPoints)
            {
                range(Yp)
            } else {
                range(c(pr0, pr1))
            }
    common <- list(xlab=which, ylab=x$model$nameY,
                   ylim=ylim, xlim=range(x$model$data[, which]),
                   main=paste(x$model$nameY, " vs ", which, " using SLSE", sep = ""))
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
    if (addPoints) {
        add <- TRUE
        pcol <- rep(allPar$treated$points$col, length(Zp))
        pcol[Zp == 0] <- allPar$nontreated$points$col
        pch. <- numeric(length(Zp))
        pch.[Zp == 1] <- allPar$treated$points$pch
        pch.[Zp == 0] <- allPar$nontreated$points$pch
        pargs <- allPar$common
        pargs$pch <- pch.
        pargs$col <- pcol
        pargs$x <- Xp
        pargs$y <- Yp
        do.call("plot", pargs)
    }
    pargs1 <- c(allPar$treated$lines, allPar$common, list(add=add))
    pargs1$x <- x$model$data[Z == 1, which]
    pargs1$y <- pr1
    do.call("matplot", pargs1)
    pargs0 <- c(allPar$nontreated$lines, list(add=TRUE))
    pargs0$x <- x$model$data[Z == 0, which]
    pargs0$y <- pr0
    do.call("matplot", pargs0)
    grid()
    if (counterfactual) 
        allPar$legend$legend <- paste(allPar$legend$legend, "-counterfactual", sep = "")
    if (!is.null(addToLegend)) 
        allPar$legend$legend = paste(allPar$legend$legend,
                                     " (", addToLegend[1], ")", sep = "")
    do.call("legend", allPar$legend)
    invisible()
}

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
