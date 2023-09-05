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
                     knots=NA)
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
    knots <- quantile(x, probs = prop.seq, type = 1)
    knots <- knots[!duplicated(knots)]
    b <- knots %in% range(x)
    if (any(b))
        knots <- knots[!b]
    if (length(knots)==0)
        return(NULL)
    .chkKnots(x, knots)
}

model.matrix.cslseModel <- function(object, ...)
{
    X <- model.matrix(object$formX, object$data, xlev=object$xlevels)    
    if (attr(terms(object$formX), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    X
}

.firstknots <- function (knames, knots) 
{
    nk <- length(knames)
    if (missing(knots)) {
        knots <- as.list(rep(NA, nk))
        names(knots) <- knames
    }
    else if (is.null(knots)) {
        knots <- lapply(1:nk, function(i) NULL)
        names(knots) <- knames
    }
    else {
        if (!is.list(knots)) 
            stop("The knots must be a list")
        if (length(knots) == 0) 
            stop("The knots cannot be an empty list")
        uknames <- names(knots)
        if (length(knots) > nk) {
            stop("You provided too many sets of knots")
        }
        else if (length(knots) == nk) {
            if (is.null(uknames)) {
                names(knots) <- knames
            } else {
                m <- match(knames, uknames)
                if (any(is.na(m))) 
                    stop("The names of the knots must match the name of the covariates")
                knots <- knots[m]
            }
        }
        else {
            if (is.null(uknames))
                stop("To manually define a subset of knots, the list must be named")
            m <- match(uknames, knames)
            if (any(is.na(m))) 
                stop("The names of the knots must match the name of the covariates")
            tmp <- as.list(rep(NA, nk))
            names(tmp) <- knames
            tmp[m] <- knots
            knots <- tmp
        }
    }
    knots
}
   
cslseModel <- function (form, data, nbasis = function(n) n^0.3, 
                       knots0, knots1)
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
    select <- ifelse(missing(knots0) & missing(knots1), "Default",
                     "User Based")
    knots0 <- .firstknots(nameX, knots0)
    knots1 <- .firstknots(nameX, knots1)    
    knots0 <- lapply(1:ncol(X), function(i)
        setKnots(X[,i], Z==0, nbasis, knots0[[i]]))
    knots1 <- lapply(1:ncol(X), function(i)
        setKnots(X[,i], Z==1, nbasis, knots1[[i]]))
    names(knots0) <- names(knots1) <- nameX
    obj <- list(na=na, formY=formY, formX=formX, treated=nameZ, nameY=nameY,
                knots0=knots0, knots1=knots1, data=data, nameX=nameX,
                method=list(select=select, crit=""), xlevels=xlevels)
    class(obj) <- "cslseModel"
    obj$method$pval <- pvalSLSE(obj, "SLSE")
    obj
}

print.cslseModel <- function(x, which=c("Model", "selKnots", "Pvalues"), ...)
{
    Z <- x$data[[x$treated]]
    which <- match.arg(which)
    if (which == "Model")
    {
        cat("Semiparametric LSE Model\n")
        cat("************************\n\n")
        cat("Number of treated: ", sum(Z), "\n")
        cat("Number of nontreated: ", sum(Z==0), "\n")
        cat("Number of missing values: ", length(x$na), "\n")
        cat("Selection method: ", x$method$select, "\n", sep="")
        if (x$method$crit != "")
            cat("Criterion: ", x$method$crit, "\n\n", sep = "")
        cat("Covariates approximated by SLSE:\n")
        w0 <- sapply(x$knots0, is.null)
        w1 <- sapply(x$knots1, is.null)    
        selPW0 <- x$nameX[!w0]
        selPW1 <- x$nameX[!w1]    
        nonselPW0 <- x$nameX[w0]
        nonselPW1 <- x$nameX[w1]
        allSame <- isTRUE(all.equal(w0,w1))
        isApp0 <- if (length(selPW0)) paste(selPW0, collapse=", ", sep="") else "None"
        isApp1 <- if (length(selPW1)) paste(selPW1, collapse=", ", sep="") else "None"
        notApp0 <- if (length(nonselPW0)) paste(nonselPW0, collapse=", ", sep="") else "None"
        notApp1 <- if (length(nonselPW1)) paste(nonselPW1, collapse=", ", sep="") else "None"
        if (!allSame)
        {
            cat("\tTreated: ", isApp1, "\n", sep="")
            cat("\tNontreated: ", isApp0, "\n", sep="")
        } else {
            cat("\t", isApp1, "\n", sep="")
        }
        cat("Covariates not approximated by SLSE:\n")   
        if (!allSame)
        {
            cat("\tTreated: ", notApp1, "\n", sep="")
            cat("\tNontreated: ", notApp0, "\n", sep="")
        } else {
            cat("\t", notApp1, "\n", sep="")
        }
    } else if (which == "selKnots") {
        cat("Semiparametric LSE Model: Selected knots\n")
        cat("****************************************\n")
        cat("Selection method: ", x$method$select, "\n", sep="")
        if (x$method$crit != "")
            cat("Criterion: ", x$method$crit, "\n\n", sep = "")
        
        cat("Treated:\n")
 	cat("********\n")
 	for (sel in 1:length(x$knots1))
 	{
            cat(x$nameX[sel],":\n", sep="")
            if (is.null(x$knots1[[sel]]))
                cat("None\n")
            else
                print.default(format(x$knots1[[sel]], ...), print.gap = 2L,
                              quote = FALSE)
        }
        cat("\nNontreated:\n")
 	cat("*************\n")
 	for (sel in 1:length(x$knots0))
 	{
            cat(x$nameX[sel],":\n", sep="")
            if (is.null(x$knots0[[sel]]))
                cat("None\n")
            else
                print.default(format(x$knots0[[sel]], ...), print.gap = 2L,
                              quote = FALSE)
        }
    } else {
        cat("Semiparametric LSE Model: P-Values of original knots\n")
        cat("****************************************************\n")

        if (!(x$method$select %in% c("BTLSE","FTLSE")))
        {
            cat("The knots have not been selected, so no p-values are available\n")
        } else {        
            cat("Selection method: ", x$method$select, "\n", sep="")
            if (x$method$crit != "")
                cat("Criterion: ", x$method$crit, "\n\n", sep = "")
        
            cat("Treated\n")
            cat("******************************\n")
            print(x$method$pval$treated, ...)        
            cat("\nNontreated\n")
            cat("*********************************\n")
            print(x$method$pval$nontreated, ...)
        }
    }
    invisible()
}

.splineMatrix <- function (model, which, treated=TRUE,
                           selObs=c("group", "all")) 
{
    selObs <- match.arg(selObs)
    Z <- model$data[[model$treated]]
    if (selObs == "group")
    {
        id <- if (treated) Z==1 else Z==0
    } else {
        id <- 1:nrow(model$data)
    }
    X <- model.matrix(model)[id,which]
    knots <- if(treated) model$knots1[[which]] else model$knots0[[which]]
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

multiSplines <- function (model, treated=TRUE, selObs=c("group", "all"))
{
    selObs <- match.arg(selObs)
    Z <- model$data[[model$treated]]
    knots <- if(treated) model$knots1 else model$knots0
    if (selObs == "group")
    {
        id <- if(treated) Z==1 else Z==0
    } else {
        id <- 1:nrow(model$data)
    }
    all <- lapply(1:length(model$nameX), function(i) {
        ans <- .splineMatrix(model, i, treated, selObs)
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

.analSe <- function(model, fit, X0, X1, beta0, beta1,
                   causal=c("ACE", "ACT", "ACN", "ALL"))
{
    causal <- match.arg(causal)
    id0 <- model$data[[model$treated]]==0
    n <- nrow(model$data)
    phat0 <- sum(id0)/n
    e <- residuals(fit)
    sig <- var(e[id0])/phat0 + var(e[!id0])/(1-phat0)
    Sigma11 <- var(X1[!id0,,drop=FALSE])
    Sigma00 <- var(X0[id0,,drop=FALSE])

    if (causal %in% c("ACT","ALL"))
    {
        Ubar01 <- colMeans(X0[!id0,,drop=FALSE])
        Sigma01 <- var(X0[!id0,,drop=FALSE])
        Psi1 <- cov(X0[!id0,,drop=FALSE],X1[!id0,,drop=FALSE])  
    }
        
    if (causal %in% c("ACN","ALL"))   
    {
        Ubar10 <- colMeans(X1[id0,,drop=FALSE])
        Sigma10 <- var(X1[id0,,drop=FALSE])    
        Psi0 <- cov(X1[id0,,drop=FALSE],X0[id0,,drop=FALSE])
    }

    if (causal %in% c("ACE","ALL"))
    {
        Ubar0 <- colMeans(X0)
        Ubar1 <- colMeans(X1)    
        Sigma1 <- var(X1)
        Sigma0 <- var(X0)
        Psi <- cov(X0,X1)
    }
    
    if (causal != "ACN")       
        {
            Ubar00 <- colMeans(X0[id0,,drop=FALSE])
            Omega00 <- var(scale(X0[id0,,drop=FALSE], scale=FALSE)*e[id0])
            nu00 <- c(cov(e[id0],scale(X0[id0,,drop=FALSE],scale=FALSE)*e[id0]))
        }

    if (causal != "ACT")
    {
        Ubar11 <- colMeans(X1[!id0,,drop=FALSE])
        Omega11 <- var(scale(X1[!id0,,drop=FALSE], scale=FALSE)*e[!id0])    
        nu11 <- c(cov(e[!id0],scale(X1[!id0,,drop=FALSE],scale=FALSE)*e[!id0]))    
    }

    se <- c()
    if (causal %in% c("ACE","ALL"))
    {
    se.ace <- sig +
    (crossprod(Ubar00 - Ubar0,
               solve(Sigma00) %*% Omega00 %*% solve(Sigma00, Ubar00 - Ubar0)) - 
    2 * crossprod(Ubar00 - Ubar0, solve(Sigma00, nu00)))/phat0 + 
    (crossprod(Ubar11 - Ubar1, 
     solve(Sigma11) %*% Omega11 %*% solve(Sigma11, Ubar11 - Ubar1)) - 
    2 * crossprod(Ubar11 - Ubar1, solve(Sigma11, nu11)))/(1-phat0) + 
    crossprod(beta1, Sigma1 %*% beta1) +
    crossprod(beta0, Sigma0 %*% beta0) - 
    2 * crossprod(beta0, Psi %*% beta1)
    se <- c(se, ACE=sqrt(se.ace/n))
    }

    if (causal %in% c("ACT","ALL"))
    {  
    se.act <- sig + 
    (crossprod(Ubar00 - Ubar01, 
               solve(Sigma00) %*% Omega00 %*% solve(Sigma00, Ubar00 - Ubar01)) - 
     2 * crossprod(Ubar00 - Ubar01, solve(Sigma00, nu00)))/phat0 +
    (crossprod(beta1, Sigma11 %*% beta1) +
     crossprod(beta0, Sigma01 %*% beta0) - 
     2 * crossprod(beta0, Psi1 %*% beta1))/(1-phat0)
    se <- c(se, ACT=sqrt(se.act/n))    
    }

    if (causal %in% c("ACN","ALL"))
    {  
    se.acn <- sig + 
    (crossprod(Ubar10 - Ubar11, 
               solve(Sigma11) %*% Omega11 %*% solve(Sigma11, Ubar10 - Ubar11)) - 
     2 * crossprod(Ubar10 - Ubar11, solve(Sigma11, nu11)))/(1-phat0) +
        (crossprod(beta1, Sigma10 %*% beta1) +
         crossprod(beta0, Sigma00 %*% beta0) - 
    2 * crossprod(beta1, Psi0 %*% beta0))/phat0
    se <- c(se, ACN=sqrt(se.acn/n))        
    }
    se
}                      

.causali <- function(model, beta, vcov, X0, X1,
                     beta0, beta1, causal, islm=TRUE)
{
    Z <- model$data[[model$treated]]
    id <- switch(causal,
                 ACE=rep(TRUE, nrow(model$data)),
                 ACT=Z==1,
                 ACN=Z==0)
    n <- sum(id)
    X0 <- X0[id,, drop = FALSE]
    X1 <- X1[id,, drop = FALSE]
    Xbar0 <- colMeans(X0)
    Xbar1 <- colMeans(X1)
    est <- c(beta[2] - beta[1] + sum(beta1*Xbar1) - sum(beta0*Xbar0))
    if (!islm)
        return(est)
    vcovXf0 <- cov(X0)
    vcovXf1 <- cov(X1)
    vcovXf01 <- cov(X0, X1)
    Dvec <- c(-1, 1, -Xbar0, Xbar1)
    se <- (sum(Dvec*c(crossprod(vcov, Dvec))) +
           sum(beta0*c(crossprod(vcovXf0, beta0)))/n +
           sum(beta1*c(crossprod(vcovXf1, beta1)))/n -
           2 * c(beta0 %*% vcovXf01 %*% beta1)/n)^0.5
    ans <- c(est, se)
    names(ans) <- c("est","se")
    ans
}

.causal <- function(model, fit, vcov, X0, X1,
                    causal=c("ACE", "ACT", "ACN", "ALL"),
                    type=c("analytic","lm"))
{
    causal <- match.arg(causal)
    type <- match.arg(type)
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
    if (type == "lm")
    {
        if (causal == "ALL") causal <- c("ACE","ACT","ACN")                
        ans <- lapply(causal, function(ci)
            .causali(model, beta, vcov, X0, X1,
                     beta0, beta1, ci, TRUE))
    } else {
        se <- .analSe(model, fit, X0, X1, beta0, beta1,
                     causal)
        if (causal == "ALL") causal <- c("ACE","ACT","ACN")
        ans <- lapply(1:length(causal), function(i)
        {
            est <- .causali(model, beta, vcov, X0, X1,
                            beta0, beta1, causal[i], FALSE)
            sei  <- se[causal[i]]
            est <- c(est, sei)
            names(est) <- NULL
            names(est) <- c("est","se")
            est
        })
    }
    names(ans) <- causal
    ans
}

causalSLSE <- function(object, ...)
{
    UseMethod("causalSLSE")
}

causalSLSE.cslseModel <- function(object,
                                 selType=c("SLSE","BSLSE","FSLSE"),
                                 selCrit = c("AIC", "BIC", "PVT"),
                                 causal = c("ALL","ACT","ACE","ACN"),
                                 seType=c("analytic", "lm"),
                                 pvalT = function(p) 1/log(p),
                                 vcov.=vcovHC, ...)
{
    selType <- match.arg(selType)
    seType <- match.arg(seType)
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
    X0 <- multiSplines(object, FALSE, "all")
    X1 <- multiSplines(object, TRUE, "all")
    ans <- .causal(object, res$lm.out, v, X0, X1, causal, seType)
    ans <- c(ans, 
             list(beta = beta, se.beta = se.beta, lm.out = res$lm.out, 
                  model=object))
    class(ans) <- c("cslse", "slseFit")
    ans    
}

causalSLSE.slseFit <- function(object, seType=c("analytic", "lm"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               vcov.=vcovHC, ...)
{
    seType <- match.arg(seType)
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
    X0 <- multiSplines(object$model, FALSE, "all")
    X1 <- multiSplines(object$model, TRUE, "all")
    ans <- .causal(object$model, object$lm.out, v, X0, X1, causal, seType)
    ans <- c(ans, 
             list(beta = beta, se.beta = se.beta, lm.out = object$lm.out, 
                  model=object$model))
    class(ans) <- c("cslse", "slseFit")
    ans    
}


print.cslse <- function (x, ...) 
{
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    cat("Selection method: ", x$model$method$select, "\n", sep="")
    if (x$model$method$crit != "")
    {
        cat("Criterion: ", x$model$method$crit, "\n\n", sep = "")
    } else {
        cat("\n")
    }
    for (causal in c("ACE","ACT","ACN"))
    {
        if (!is.null(x[[causal]]))
            cat(causal, " = ", x[[causal]]["est"], "\n", sep="")
    }
}

causalSLSE.formula <- function(object, data, nbasis=function(n) n^0.3,
                               knots0, knots1, 
                               selType=c("SLSE","BSLSE","FSLSE"),
                               selCrit = c("AIC", "BIC", "PVT"),
                               causal = c("ALL","ACT","ACE","ACN"),
                               seType=c("analytic", "lm"),
                               pvalT = function(p) 1/log(p),
                               vcov.=vcovHC, ...)
{
    model <- cslseModel(object, data, nbasis,  knots0, knots1)
    selType <- match.arg(selType)
    seType <- match.arg(seType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    causalSLSE(object=model, selType=selType, selCrit = selCrit,
               causal = causal, seType=seType, pvalT =  pvalT,
               vcov.=vcov., ...)    
}

print.slseFit <- function(x, ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    cat("Selection method: ", x$model$method$select, "\n", sep="")
    if (x$model$method$crit != "")
    {
        cat("Criterion: ", x$model$method$crit, "\n\n", sep = "")
    } else {
        cat("\n")
    }
    print.default(format(coef(x$lm.out), ...), print.gap = 2L, 
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

print.summary.slseFit <- function(x, digits = 4,
                                  signif.stars = getOption("show.signif.stars"),
                                  ...)
{
    cat("Semiparametric LSE\n")
    cat("******************\n")
    cat("Selection method: ", x$model$method$select, "\n", sep="")
    if (x$model$method$crit != "")
    {
        cat("Criterion: ", x$model$method$crit, "\n\n", sep = "")
    } else {
        cat("\n")
    }
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
    nb0 <- ifelse(is.null(object$model$knots0), 1,
                  length(object$model$knots0))
    nb1 <- ifelse(is.null(object$model$knots1), 1,
                  length(object$model$knots1))
    t <- object$beta/object$se.beta
    pv <- 2 * pnorm(-abs(t))
    beta <- cbind(object$beta, object$se.beta, t, pv)
    colnames(beta) <- c("Estimate", "Std. Error", "t value", 
        "Pr(>|t|)")
    ans <- list(causal = est, beta = beta, crit = object$model$method$crit, 
                knots0 = object$model$knots0,
                knots1 = object$model$knots1,
                covNames = names(object$model$knots0),
                select=object$model$method$select)
    class(ans) <- "summary.cslse"
    ans
}

print.summary.cslse <- function (x, digits = 4,
                                signif.stars = getOption("show.signif.stars"), 
                                beta = FALSE, knots = FALSE, ...) 
{
    cat("Causal Effect using Semiparametric LSE\n")
    cat("**************************************\n")
    cat("Selection method: ", x$select, "\n", sep="")
    if (x$crit != "")
        cat("Criterion: ", x$model$method$crit, "\n\n", sep = "")
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)
    if (beta) {
        cat("\nPiecewise polynomials coefficients\n")
        cat("**********************************\n")
        printCoefmat(x$beta, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    if (knots) {
        cat("\nNumber of selected knots per covariate\n")
        cat("****************************************\n")
        cat("Treated:\n")
        print.default(format(sapply(x$knots1, length)), print.gap = 2L, 
            quote = FALSE)
        cat("Nontreated:\n")
        print.default(format(sapply(x$knots0, length)), print.gap = 2L, 
            quote = FALSE)
        cat("\n")
    }
}

## The BSLSE

pvalSLSE <- function(model, ...)
{
    UseMethod("pvalSLSE")
}

pvalSLSE.cslseModel <- function(model, method=c("BSLSE", "FSLSE", "SLSE"),
                                vcov.=vcovHC, ...)
{
    method <- match.arg(method)
    if (method == "SLSE")
    {
        pv <- list(p0=sapply(model$knots0, function(ki) length(ki)+1),
                   p1=sapply(model$knots1, function(ki) length(ki)+1),
                   pval0=lapply(length(model$knots0), function(i) NULL),
                   pval1=lapply(length(model$knots1), function(i) NULL))
    } else {
        pv <- if (method=="BSLSE")  .getPvalB(model, vcov., ...)
              else .getPvalF(model, vcov., ...)
    }
    ans <- list(treated=list(pval=pv$pval1, knots=model$knots1, p=pv$p1,
                              method=method),
                nontreated=list(pval=pv$pval0, knots=model$knots0, p=pv$p0,
                                 method=method))
    class(ans$treated) <-  class(ans$nontreated) <- "slsePval"
    ans
}

print.slsePval <- function(x, ...)
{
    w <- sapply(x$knots, is.null)
    selPW <- names(x$knots)[!w]
    nonselPW <- names(x$knots)[w]
    notApp <- if (length(nonselPW))
              {
                  paste(nonselPW, collapse = ", ", sep = "")
              } else {
                  "None"
              }
    cat("\nCovariates with no knots:\n")
    cat("\t", notApp, "\n", sep = "")
    cat("\nCovariates with knots:\n")
    if (length(selPW) == 0)
    {
        cat("None\n")
    } else {
        for (ki in selPW)
        {
            cat(ki,":\n")
            pvi <- rbind(x$knots[[ki]], x$pval[[ki]])
            rownames(pvi) <- c("Knots", "P-Value")[1:nrow(pvi)]
            print.default(pvi, quote=FALSE, ...)
            cat("\n")
        }
    }
    invisible()
}
    

.getPvalB <- function (model, vcov.=vcovHC, ...) 
{
    data2 <- model$data
    data2$Xf1 <- multiSplines(model, TRUE)
    data2$Xf0 <- multiSplines(model, FALSE)
    form <- model$formY
    environment(form) <- environment()
    fit <- lm(form, data2)
    p0 <- attr(data2$Xf0, "p")
    p1 <- attr(data2$Xf1, "p")
    v <- vcov.(fit, ...)
    pval0 <- lapply(1:length(model$knots0), function(i) {
        ki <- length(model$knots0[[i]])
        if (ki==0) NA else  .testKnots(fit, model, 1:ki, i , FALSE, v)
    })
    pval1 <- lapply(1:length(model$knots1), function(i) {
        ki <- length(model$knots1[[i]])
        if (ki==0) NA else  .testKnots(fit, model, 1:ki, i , TRUE, v)
    })
    names(pval0) <- names(pval1) <- names(p0) <- names(p1) <- model$nameX
    list(pval0 = pval0, pval1 = pval1, p0=p0, p1=p1)
}

.selIC <- function (model, pvalRes, pvalT = NULL, crit)
{
    pval <- c(do.call("c", pvalRes$nontreated$pval),
              do.call("c", pvalRes$treated$pval))
    pval_sort <- sort(pval)    
    q <- length(pval)
    p <- sum(pvalRes$treated$p) + sum(pvalRes$nontreated$p)
    w0 <- lapply(1:length(pvalRes$nontreated$pval), function(i) NULL)
    w1 <- lapply(1:length(pvalRes$treated$pval), function(i) NULL)
    res0 <- estSLSE(model, w0, w1)
    icV <- ic_seq0 <- get(crit)(res0$lm.out)   
    for (i in 1:q) {
        w0 <- lapply(1:length(pvalRes$nontreated$pval), function(j)
        {
            w <- which(pvalRes$nontreated$pval[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        w1 <- lapply(1:length(pvalRes$treated$pval), function(j)
        {
            w <- which(pvalRes$treated$pval[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        res1 <- estSLSE(model, w0, w1)
        ic_seq1 <- get(crit)(res1$lm.out)
        icV <- c(icV, ic_seq1)
        if (ic_seq1 < ic_seq0) {
            ic_seq0 <- ic_seq1
            res0 <- res1
        }
    }
    res0$model
}

.selPVT <- function (model, pvalRes, pvalT = function(p) 1/log(p),
                     crit=NULL, ...)
{
    pval <- c(do.call("c", pvalRes$nontreated$pval),
              do.call("c", pvalRes$treated$pval))
    n <- nrow(model$data)
    q <- length(pval)
    p <- mean(c(pvalRes$nontreated$p,pvalRes$treated$p))
    crit <- pvalT(p)
    w0 <- lapply(1:length(model$knots0), function(i)
    {
        w <- which(pvalRes$nontreated$pval[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    w1 <- lapply(1:length(model$knots1), function(i)
    {
        w <- which(pvalRes$treated$pval[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    model$knots0 <- .chkSelKnots(model, w0, FALSE)
    model$knots1 <- .chkSelKnots(model, w1, TRUE)
    model
}

selSLSE <- function(model, ...)
{
    UseMethod("selSLSE")
}


selSLSE.cslseModel <- function(model, selType=c("BSLSE", "FSLSE"),
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
    pval <- pvalSLSE(model, selType, vcov., ...)
    model <- critFct(model, pval, pvalT, selCrit)
    model$method <- list(select=selType, crit=selCrit, pval=pval)
    model
}

.chkSelKnots <- function(model, w, treated=TRUE)
{
    wK <- ifelse(treated, "knots1", "knots0")
    knots <- model[[wK]]
    if (is.null(w))
        return(knots)
    if (!is.list(w))
        stop("The knots selection must be included in a list")
    if (length(w) != length(knots))
        stop(paste("The length of the knots selection list does not match the length of ",
                   wK, sep=""))
    k <- lapply(1:length(w), function(i) {
        ki <- knots[[i]]
        wi <- w[[i]]
        if (is.null(ki))
            return(NULL)
        if (is.null(wi))
            return(NULL)
        if (any(is.na(wi)))
            stop("The knots selection list cannot contain NAs")
        if (!is.integer(wi))
            stop("The knots selection list can only contain integers")
        wi <- unique(wi)
        if (any(wi<1) | any(wi>length(ki)))
            stop(paste("Knot selection out of bound in ", wK, sep=""))
        ki <- ki[wi]})
    names(k) <- model$nameX
    k
}

estSLSE <- function(model, ...)
{
    UseMethod("estSLSE")
}

estSLSE.cslseModel <- function(model, w0=NULL, w1=NULL, ...)
{
    if (!inherits(model, "cslseModel"))
        stop("model must be an object of class cslseModel")
    model$knots0 <- .chkSelKnots(model, w0, treated=FALSE)
    model$knots1 <- .chkSelKnots(model, w1, treated=TRUE)
    data <- model$data
    data$Xf1 <- multiSplines(model, TRUE)
    data$Xf0 <- multiSplines(model, FALSE)
    form <- model$formY
    environment(form) <-  environment()    
    fit <- lm(form, data)
    obj <- list(lm.out=fit, model=model)
    class(obj) <- "slseFit"
    obj
}

### The function tests the difference on each side of knots whichK
### (more than one is allowed) for the whichX covariate.
### Knots for treated are tested if treated is set to TRUE
### vcov is the covariance matrix of fit.


.testKnots <- function(fit, model, whichK, whichX, treated, vcov)
{
    wK <- ifelse(treated, "knots1", "knots0")
    wX <- ifelse(treated, "Xf1", "Xf0")
    if (whichX > length(model[[wK]]))
        stop("whichX exceeds the number of covariates")
    if (any(whichK > length(model[[wK]][[whichX]])))
        stop("whichK exceeds the number of knots")   
    if (is.null(model[[wK]][[whichX]]))
        return(NA)
    b <- coef(fit)
    b <- na.omit(b)
    sapply(whichK, function(wi) {
        nX <- names(model[[wK]])[[whichX]]
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

## The FSLSE

.getPvalF <- function (model, vcov.=vcovHC, ...) 
{
    w0 <- lapply(1:length(model$knots0), function(i) NULL)
    w1 <- lapply(1:length(model$knots1), function(i) NULL)
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
    pvali <- function(i, treated=TRUE)
    {
        p <- if(treated) length(model$knots1[[i]]) else length(model$knots0[[i]])
        if (p==0)
            return(NA)
        sel <- selFct(p)
        c(sapply(1:length(sel), function(j) {
            selj <- as.integer(sel[[j]])
            if (treated)
                w1[[i]] <- selj
            else
                w0[[i]] <- selj
            res <- estSLSE(model, w0, w1)
            v <- vcov.(res$lm.out, ...)
            if (length(selj) == 1L)
            {
                .testKnots(res$lm.out, res$model, 1L, i, treated, v)
            } else if (length(sel)==1 & length(selj)==2) {
                .testKnots(res$lm.out, res$model, 1L:2L, i, treated, v)
            } else {
                whichK <- ifelse(length(selj)==2 & selj[1]==1L, 1, 2)
                .testKnots(res$lm.out, res$model, whichK, i, treated, v)
            }}))
    }
    pval0 <- lapply(1:length(model$knots0), function(i) pvali(i, FALSE))
    pval1 <- lapply(1:length(model$knots1), function(i) pvali(i))
    p0 <- sapply(model$knots0, function(ki) length(ki)+1)
    p1 <- sapply(model$knots1, function(ki) length(ki)+1)    
    names(pval0) <- names(pval1) <- names(p0) <- names(p1) <- model$nameX
    list(pval0=pval0, pval1=pval1, p0=p0, p1=p1)
}


extract.cslse <- function (model, include.nobs = TRUE,
                                include.nknots = TRUE,
                                include.numcov = TRUE, include.rsquared = TRUE,
                                include.adjrsquared=TRUE, 
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
        rs1 <- length(unlist(model$model$knots0))
        rs2 <- length(unlist(model$model$knots1))        
        gof <- c(gof, rs1, rs2)
        gof.names <- c(gof.names, "Num. knots (Nontreated)", "Num. knots (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
   }
    if (isTRUE(include.numcov)) {
        rs3 <- length(model$model$nameX)
        gof <- c(gof, rs3)
        gof.names <- c(gof.names, "Num. covariates")
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
    if (isTRUE(include.adjrsquared)) {
        f <- fitted(model$lm.out)
        ess <- sum((f-mean(f))^2)
        rss <- sum(residuals(model$lm.out)^2)       
        R2 <- ess/(ess+rss)
        R2adj <- 1-(1-R2)*(nobs(model$lm.out)-1)/model$lm.out$df.residual        
        gof <- c(gof, R2adj)
        gof.names <- c(gof.names, "R$^2_{adj}$")
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
    newdata$Xf1 <- multiSplines(model, TRUE, "all")
    newdata$Xf0 <- multiSplines(model, FALSE, "all")
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

plotOLD <- function (x, y, which = y, interval = c("none", "confidence"),
                          counterfactual=FALSE,
                          level = 0.95, newdata = NULL, legendPos = "topright",
                          vcov. = vcovHC, 
                          col0=1, col1=2, lty0=1, lty1=2,  add.=FALSE,
                          addToLegend=NULL,
                          cex=1, ylim.=NULL, xlim.=NULL, addPoints=FALSE,
                          FUN=mean, main=NULL, ...) 
{
    interval <- match.arg(interval)
    vnames <- all.vars(x$model$formX)
    treat <- x$model$treated    
    if (!is.null(newdata)) {
        if (!is.numeric(newdata)) 
            stop("newdata must be a vector of numeric values")
        if (is.null(names(newdata))) 
            stop("newdata must be a named vector")
        if (any(names(newdata) == "")) 
            stop("All elements of newdata must be named")
        nd <- names(newdata)
    }
    if (is.numeric(which)) 
        which <- vnames[which]
    if (!is.character(which) & length(which) != 1) 
        stop("which must be a character type")
    if (!(which %in% vnames)) 
        stop("which must be one of the names of the variables")
    if (addPoints)
    {
        Yp <- x$model$data[,x$model$nameY]
        Xp <- x$model$data[,which]
        Zp <- x$model$data[,treat]
    }
    ind <- order(x$model$data[, which])
    data <- x$model$data[ind, ]
    Z <- data[,treat]    
    data[Z==1, !(names(data) %in% c(which, treat))] <- sapply(which(!(names(data) %in% 
        c(which, treat))), function(i) rep(FUN(data[Z==1, i]), 
                                           sum(Z)))
    data[Z==0, !(names(data) %in% c(which, treat))] <- sapply(which(!(names(data) %in% 
        c(which, treat))), function(i) rep(FUN(data[Z==0, i]), 
                                           sum(1-Z)))    
    if (!is.null(newdata)) 
        for (ndi in nd) data[[ndi]] <- newdata[ndi]
    res <- predict(x, interval = interval, level = level,
                   newdata = data, se.fit = FALSE, vcov. = vcov.,
                   counterfactual=counterfactual, ...)
    pr0 <- res$nontreated
    pr1 <- res$treated
    if (interval == "confidence") {
        lty0 = c(lty0, 3, 3)
        lty1 = c(lty1, 3, 3)
        lwd = c(2, 1, 1)
    }
    else {
        lwd = 2
    }
    if (is.null(main))
        main <- paste("Outcome VS ", which, " using SLSE", 
                      sep = "")
    if (addPoints)
    {
        add.=TRUE
        pcol <- rep(col1, length(Zp))
        pcol[Zp==0] <- col0
        pch. <- 21*(Zp==1)+22*(Zp==0)
        plot(Xp, Yp, pch=pch., col=pcol, 
             main=main, ylab = x$model$nameY, xlab = which)
    }
    if (is.null(ylim.))
        ylim. <- range(c(pr0, pr1))
    if (is.null(xlim.))
        xlim. <- range(data[, which])
    matplot(data[Z == 1, which], pr1, col = col1, ylim = ylim., type = "l", 
        lty = lty1, lwd = lwd, main = main, ylab = x$model$nameY, 
        xlab = which, add=add., xlim=xlim.)
    matplot(data[Z == 0, which], pr0, col = col0, type = "l", lty = lty0,
            lwd = lwd, add = TRUE)
    grid()
    Leg=c("Treated", "Nontreated")
    if(counterfactual)
        Leg <- paste(Leg, "-counterfactual", sep="")
    if (!is.null(addToLegend))
        Leg=paste(Leg, " (", addToLegend[1], ")", sep="")
    if (addPoints) pch. <- c(21,22) else pch.=c(NA,NA)
    
    legend(legendPos, Leg, col = c(col1, col0), pch=pch., 
           lty = c(lty1[1], lty0[1]), lwd = 2, bty = "n", cex=cex)
    invisible()
 }

plot.slseFit <- function (x, y, which = y, interval = c("none", "confidence"), 
                   counterfactual = FALSE, level = 0.95, fixedCov0 = NULL,
                   fixedCov1 = fixedCov0, legendPos = "topright", 
                   vcov. = vcovHC, col0 = 1, col1 = 2, lty0 = 1, lty1 = 2, add. = FALSE, 
                   addToLegend = NULL, cex = 1, ylim. = NULL, xlim. = NULL, 
                   addPoints = FALSE, FUN = mean, main = NULL, plot=TRUE, ...) 
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
    if (interval == "confidence") {
        lty0 = c(lty0, 3, 3)
        lty1 = c(lty1, 3, 3)
        lwd = c(2, 1, 1)
    }
    else {
        lwd = 2
    }
    outcome <- x$model$nameY
    if (is.null(main)) 
        main <- paste(outcome, " vs ", which, " using SLSE", sep = "")
    if (addPoints) {
        add. = TRUE
        pcol <- rep(col1, length(Zp))
        pcol[Zp == 0] <- col0
        pch. <- 21 * (Zp == 1) + 22 * (Zp == 0)
        plot(Xp, Yp, pch = pch., col = pcol, main = main, ylab = x$model$nameY, 
            xlab = which)
    }
    if (is.null(ylim.)) 
        ylim. <- range(c(pr0, pr1))
    if (is.null(xlim.)) 
        xlim. <- range(x$model$data[, which])
    matplot(x$model$data[Z == 1, which], pr1, col = col1, ylim = ylim., 
        type = "l", lty = lty1, lwd = lwd, main = main, ylab = x$model$nameY, 
        xlab = which, add = add., xlim = xlim.)
    matplot(x$model$data[Z == 0, which], pr0, col = col0, type = "l", 
        lty = lty0, lwd = lwd, add = TRUE)
    grid()
    Leg = c("Treated", "Nontreated")
    if (counterfactual) 
        Leg <- paste(Leg, "-counterfactual", sep = "")
    if (!is.null(addToLegend)) 
        Leg = paste(Leg, " (", addToLegend[1], ")", sep = "")
    if (addPoints) 
        pch. <- c(21, 22)
    else pch. = c(NA, NA)
    legend(legendPos, Leg, col = c(col1, col0), pch = pch., lty = c(lty1[1], 
        lty0[1]), lwd = 2, bty = "n", cex = cex)
    invisible()
}

selSLSE2 <- function (model, selType = c("FSLSE", "BSLSE"), selCrit = c("AIC", 
    "BIC", "PVT"), pvalT = function(p) 1/log(p), vcov. = vcovHC, 
    ...) 
{
    selCrit <- match.arg(selCrit)
    selType <- match.arg(selType)
    critFct <- if (selCrit == "PVT") {
        .selPVT
    }
    else {
        .selICF
    }
    if (selType == "BSLSE") 
        pval <- .getPvalB(model, vcov., ...)
    else pval <- .getPvalF(model, vcov., ...)
    critFct(model, pval, pvalT, selCrit)
}

.reshapeKnots <- function(model, w, mnk, treated)
{
    w <- matrix(w, nrow=mnk)
    w <- lapply(1:ncol(w) ,function(i) {
        if (w[1,i]==0)
            return(NULL)
        sort(w[w[,i]!=0,i])})
    .chkSelKnots(model, w, treated)
}

.selICF <- function (model, pvalRes, pvalT = NULL, crit) 
{
    x <- model.matrix(model)
    y <- model$data[,model$nameY]
    z <- model$data[,model$treated]
    nk0 <- sapply(model$knots0, length)
    nk1 <- sapply(model$knots1, length)
    mnk0 <- max(nk0)
    tnk0 <- sum(nk0)
    pval <- c(do.call("c", pvalRes$pval0), do.call("c", pvalRes$pval1))
    spval <- sort(pval)
    npval <- length(spval)
    pval0 <-  sapply(pvalRes$pval0, function(pvi) {
        pv <- numeric(mnk0)
        if (!is.na(pvi[1]))
            pv[1:length(pvi)] <- pvi
        pv}) 
    knots0 <- sapply(model$knots0, function(ki) {
        k <- numeric(mnk0)
        if (length(ki))
            k[1:length(ki)] <- ki
        k})
    mnk1 <- max(nk1)
    tnk1 <- sum(nk1)
    pval1 <-  sapply(pvalRes$pval1, function(pvi) {
        pv <- numeric(mnk1)
        if (!is.na(pvi[1]))
            pv[1:length(pvi)] <- pvi
        pv})    
    knots1 <- sapply(model$knots1, function(ki) {
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
    modelAIC$knots0 <-.reshapeKnots(model, sp$w0AIC, mnk0, FALSE)
    modelAIC$knots1 <-.reshapeKnots(model, sp$w1AIC, mnk1, FALSE) 
    modelBIC$knots0 <-.reshapeKnots(model, sp$w0BIC, mnk0, FALSE)
    modelBIC$knots1 <-.reshapeKnots(model, sp$w1BIC, mnk1, FALSE) 
    list(AIC=modelAIC, BIC=modelBIC)
}

