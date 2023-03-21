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
    knots
}


setKnots <- function(x, sel=1:length(x), nknots=function(n) n^0.3,
                     knots=NA)
{
    x <- x[sel]
    if (!is.null(knots))
        {
            if (all(x %in% c(1,0)))
                return(NULL)  
            if (is.na(knots))
            {
                nzero <- x==0
                n <- sum(!nzero)
                p <- floor(nknots(n))
                if (p == 0)
                    return(NULL)
                prop.seq <- seq(from = 0, to = 1, length.out = p + 1)
                prop.seq <- prop.seq[-c(1, p + 1)]
                knots <- quantile(x[!nzero], probs = prop.seq, type = 1)
                if (any(duplicated(knots)))
                    knots <- unique(knots)
            }
        }
    .chkKnots(x, knots)
}

model.matrix.tlseModel <- function(object, ...)
{
    X <- model.matrix(object$formX, object$data)    
    if (attr(terms(object$formX), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    X
}

    
setModel <- function (form, data, nknots = function(n) n^0.3, 
                       knots0 = NA, knots1 = NA, userRem=NULL, ...)
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
    X <- model.matrix(formX, data)
    if (attr(terms(formX), "intercept") == 1) 
        X <- X[, -1, drop = FALSE]
    na <- na.omit(cbind(Y,Z,X))
    if (is.na(knots0))
        knots0 <- as.list(rep(NA, ncol(X)))
    if (is.na(knots1))
        knots1 <- as.list(rep(NA, ncol(X)))
    if (!all(c(is.list(knots0), is.list(knots1))))
        stop("knots0 and knots1 must be a list")
    if ( (length(knots0)!=ncol(X)) | (length(knots1)!=ncol(X)) )
        stop("The length of knots must be equal to the number of covariates")    
    if (!is.null(attr(na, "omit")))
    {
        na <- attr(na, "omit")
        X <- X[-na,,drop=FALSE]
        Z <- Z[-na,,drop=FALSE]
        data <- data[-na,,drop=FALSE]
    } else {
        na <- NULL
    }
    if (!is.null(userRem))
    {
        w <- which(colnames(X) %in% userRem)
        if (length(w))
        {
            knots0[w] <- lapply(w, function(i) NULL)
            knots1[w] <- lapply(w, function(i) NULL)
        }
    }
    nameX <- colnames(X)
    nameY <- all.vars(formY)[1]
    nameZ <- colnames(Z)
    knots0 <- lapply(1:ncol(X), function(i)
        setKnots(X[,i], Z==0, nknots, knots0[[i]]))
    knots1 <- lapply(1:ncol(X), function(i)
        setKnots(X[,i], Z==1, nknots, knots1[[i]]))
    names(knots0) <- names(knots1) <- nameX
    obj <- list(na=na, formY=formY, formX=formX, treated=nameZ, nameY=nameY,
                knots0=knots0, knots1=knots1, data=data, nameX=nameX)
    class(obj) <- "tlseModel"
    obj
}

print.tlseModel <- function(x, knots=FALSE, ...)
{
    Z <- x$data[[x$treated]]
    cat("Semiparametric Thresholding LSE Model\n")
    cat("*************************************\n\n")
    cat("Number of treated: ", sum(Z), "\n")
    cat("Number of control: ", sum(Z==0), "\n")
    cat("Number of missing values: ", length(x$na), "\n")
    cat("Covariates being approximated by a piecewise function:\n")
    w <- sapply(x$knots0, is.null)
    selPW <- x$nameX[!w]
    nonselPW <- x$nameX[w]    
    if (length(selPW))
    {
        cat("\t", paste(selPW, collapse=", ", sep=""), "\n")
    } else {
        cat("\t None\n")
    }
    cat("Covariates not being approximated by a piecewise function:\n")    
    if (length(nonselPW))
    {
        cat("\t", paste(nonselPW, collapse=", ", sep=""), "\n")
    } else {
        cat("\t None\n")
    }
    if (knots)
    {
        cat("Lists of knots for the treated group\n")
        cat("************************************\n")
        for (sel in which(!w))
        {
            cat(x$nameX[sel],":\n", sep="")
            print.default(format(x$knots1[[sel]], ...), print.gap = 2L, 
                          quote = FALSE)

        }
        cat("Lists of knots for the Control group\n")
        cat("************************************\n")
        for (sel in which(!w))
        {
            cat(x$nameX[sel],":\n", sep="")
            print.default(format(x$knots0[[sel]], ...), print.gap = 2L, 
                          quote = FALSE)

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
    Sigma11 <- var(X1[!id0,])
    Sigma00 <- var(X0[id0,])

    if (causal %in% c("ACT","ALL"))
    {
        Ubar01 <- colMeans(X0[!id0,])
        Sigma01 <- var(X0[!id0,])
        Psi1 <- cov(X0[!id0,],X1[!id0,])   
    }
        
    if (causal %in% c("ACN","ALL"))   
    {
        Ubar10 <- colMeans(X1[id0,])
        Sigma10 <- var(X1[id0,])    
        Psi0 <- cov(X1[id0,],X0[id0,])
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
            Ubar00 <- colMeans(X0[id0,])
            Omega00 <- var(scale(X0[id0,], scale=FALSE)*e[id0])
            nu00 <- c(cov(e[id0],scale(X0[id0,],scale=FALSE)*e[id0]))
        }

    if (causal != "ACT")
    {
        Ubar11 <- colMeans(X1[!id0,])
        Omega11 <- var(scale(X1[!id0,], scale=FALSE)*e[!id0])    
        nu11 <- c(cov(e[!id0],scale(X1[!id0,],scale=FALSE)*e[!id0]))    
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
                    type=c("analytical","lm"))
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

tlse <- function (model, selType=c("SLSE","BTLSE","FTLSE"),
                  selCrit = c("ASY", "AIC", "BIC"),
                  causal = c("ALL","ACT","ACE","ACN"),
                  seType=c("analytical", "lm"),
                  minPV = function(p) 1/log(p), vcov.=NULL, ...)
{
    selType <- match.arg(selType)
    seType <- match.arg(seType)
    causal <- match.arg(causal)
    selCrit <- match.arg(selCrit)
    if (selType != "SLSE")
        model <- selTLSE(model, selType, selCrit, minPV, vcov., ...)
    res <- estModel(model)
    beta <- coef(res$fit)
    if(is.null(vcov.))
        v <- vcov(res$fit)
    else v <- vcov.(res$fit, ...)
    se.beta <- sqrt(diag(v))
    if (any(is.na(beta))) 
        warning(paste("\nThe final regression is multicollinear.",
                      " The result may not be valid:", 
                      "\nThe following variables produced NA's\n",
                      paste(names(beta)[is.na(beta)], 
                            collapse = ", "), "\n", sep = ""))
    X0 <- multiSplines(model, FALSE, "all")
    X1 <- multiSplines(model, TRUE, "all")
    ans <- .causal(model, res$fit, v, X0, X1, causal, seType)
    ans <- c(ans, 
             list(beta = beta, se.beta = se.beta, lm.out = res$fit, 
                  crit = selCrit, model=model, type=selType))
    class(ans) <- "tlse"
    ans    
}

print.tlse <- function (x, ...) 
{
    cat("Causal Effect using Thresholding Least Squares\n")
    cat("**********************************************\n")
    cat("Type: ", x$type, "\n", sep="")
    if (x$type != "SLSE")
        cat("Selection method: ", x$crit, "\n\n", sep = "")
    for (causal in c("ACE","ACT","ACN"))
    {
        if (!is.null(x[[causal]]))
            cat(causal, " = ", x[[causal]]["est"], "\n", sep="")
    }
}

summary.tlse <- function (object, ...) 
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
    ans <- list(causal = est, beta = beta, crit = object$crit, 
                knots0 = object$model$knots0,
                knots1 = object$model$knots1,
                covNames = names(object$model$knots0),
                type=object$type)
    class(ans) <- "summary.tlse"
    ans
}

print.summary.tlse <- function (x, digits = 4,
                                signif.stars = getOption("show.signif.stars"), 
                                beta = FALSE, knots = FALSE, ...) 
{
    cat("Causal Effect using Thresholding Least Squares\n")
    cat("**********************************************\n")
    cat("Type: ", x$type, "\n", sep="")
    if (x$type != "SLSE")
        cat("Selection method: ", x$crit, "\n\n", sep = "")
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
        na.print = "NA", ...)
    if (beta) {
        cat("\nPiecewise polynomials coefficients\n")
        cat("**********************************\n")
        printCoefmat(x$beta, digits = digits, signif.stars = signif.stars, 
            na.print = "NA", ...)
    }
    if (knots) {
        cat("\nNumber of selected knots per variables\n")
        cat("***************************************\n")
        cat("Treated group:\n")
        print.default(format(sapply(x$knots1, length)), print.gap = 2L, 
            quote = FALSE)
        cat("Control group:\n")
        print.default(format(sapply(x$knots0, length)), print.gap = 2L, 
            quote = FALSE)
        cat("\n")
    }
}

## The BTLSE

.getPvalB <- function (model, vcov.=NULL, ...) 
{
    data2 <- model$data
    data2$Xf1 <- multiSplines(model, TRUE)
    data2$Xf0 <- multiSplines(model, FALSE)
    form <- model$formY
    environment(form) <- environment()
    fit <- lm(form, data2)
    p0 <- attr(data2$Xf0, "p")
    p1 <- attr(data2$Xf1, "p")
    if(is.null(vcov.))
        v <- vcov(fit)
    else v <- vcov.(fit, ...)
    pval0 <- lapply(1:length(model$knots0), function(i) {
        ki <- length(model$knots0[[i]])
        if (ki==0) NA else  .testKnots(fit, model, 1:ki, i , FALSE, v)
    })
    pval1 <- lapply(1:length(model$knots1), function(i) {
        ki <- length(model$knots1[[i]])
        if (is.null(ki)) NA else  .testKnots(fit, model, 1:ki, i , TRUE, v)
    })
    names(pval0) <- names(pval1) <- names(p0) <- names(p1) <- model$nameX
    list(pval0 = pval0, pval1 = pval1, p0=p0, p1=p1)
}

.selIC <- function (model, pvalRes, minPV = NULL, crit)
{
    pval <- c(do.call("c", pvalRes$pval0), do.call("c", pvalRes$pval1))
    pval_sort <- sort(pval)    
    q <- length(pval)
    p <- sum(pvalRes$p0) + sum(pvalRes$p1)
    w0 <- lapply(1:length(pvalRes$pval0), function(i) NULL)
    w1 <- lapply(1:length(pvalRes$pval1), function(i) NULL)
    res0 <- estModel(model, w0, w1)
    icV <- ic_seq0 <- get(crit)(res0$fit)   
    for (i in 1:q) {
        w0 <- lapply(1:length(pvalRes$pval0), function(j)
        {
            w <- which(pvalRes$pval0[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        w1 <- lapply(1:length(pvalRes$pval1), function(j)
        {
            w <- which(pvalRes$pval1[[j]] <= pval_sort[i])
            if (length(w)==0)
                w <- NULL
            w})
        res1 <- estModel(model, w0, w1)
        ic_seq1 <- get(crit)(res1$fit)
        icV <- c(icV, ic_seq1)
        if (ic_seq1 < ic_seq0) {
            ic_seq0 <- ic_seq1
            res0 <- res1
        }
    }
    res0$model
}


.selASY <- function (model, pvalRes, minPV = function(p) 1/log(p),
                     crit=NULL)
{
    pval <- c(do.call("c", pvalRes$pval0), do.call("c", pvalRes$pval1))
    n <- nrow(model$data)
    q <- length(pval)
    p <- mean(c(pvalRes$p0,pvalRes$p1))
    crit <- minPV(p)
    w0 <- lapply(1:length(model$knots0), function(i)
    {
        w <- which(pvalRes$pval0[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    w1 <- lapply(1:length(model$knots1), function(i)
    {
        w <- which(pvalRes$pval1[[i]]<crit)
        if (length(w) == 0)
            w <- NULL
        w})
    model$knots0 <- .chkSelKnots(model, w0, FALSE)
    model$knots1 <- .chkSelKnots(model, w1, TRUE)
    model
}

selTLSE <- function(model, method=c("FTLSE", "BTLSE"),
                    crit = c("ASY", "AIC", "BIC"), 
                    minPV = function(p) 1/(p * log(p)), vcov.=NULL, ...)
{
    crit <- match.arg(crit)
    method <- match.arg(method)
    critFct <- if (crit == "ASY") {
                   get(paste(".sel", crit, sep=""))
               } else {
                   .selIC
               }
    if (method == "BTLSE")
        pval <- .getPvalB(model, vcov., ...)
    else
        pval <- .getPvalF(model, vcov., ...)
    critFct(model, pval, minPV, crit)
}

.chkSelKnots <- function(model, w, treated=TRUE)
{
    wK <- ifelse(treated, "knots1", "knots0")
    knots <- model[[wK]]
    if (is.null(w))
        return(knots)
    if (!is.list(w))
        stop("The knot selection must be inclluded in a list")
    if (length(w) != length(knots))
        stop(paste("The length of the knot selection list does not match the length of ",
                   wK, sep=""))
    k <- lapply(1:length(w), function(i) {
        ki <- knots[[i]]
        wi <- w[[i]]
        if (is.null(ki))
            stop(paste("There are no knots in the ", i, "th covariate of ", wK, sep=""))
        if (is.null(wi))
            return(NULL)
        if (any(is.na(wi)))
            stop("The knot selection list cannot contain NAs")
        if (!is.integer(wi))
            stop("The knot selection list can only contain integers")
        wi <- unique(wi)
        if (any(wi<1) | any(wi>length(ki)))
            stop(paste("Knot selection out of bound in ", wK, sep=""))
        ki <- ki[wi]})
    names(k) <- model$nameX
    k
}

estModel <- function(model, w0=NULL, w1=NULL)
{
    model$knots0 <- .chkSelKnots(model, w0, treated=FALSE)
    model$knots1 <- .chkSelKnots(model, w1, treated=TRUE)
    data <- model$data
    data$Xf1 <- multiSplines(model, TRUE)
    data$Xf0 <- multiSplines(model, FALSE)
    form <- model$formY
    environment(form) <-  environment()    
    fit <- lm(form, data)
    list(fit=fit, model=model, data=data)
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

.getPvalF <- function (model, vcov.=NULL, ...) 
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
            res <- estModel(model, w0, w1)
            if(is.null(vcov.))
                v <- vcov(res$fit)
            else v <- vcov.(res$fit, ...)
            if (length(selj) == 1L)
            {
                .testKnots(res$fit, res$model, 1L, i, treated, v)
            } else if (length(sel)==1 & length(selj)==2) {
                .testKnots(res$fit, res$model, 1L:2L, i, treated, v)
            } else {
                whichK <- ifelse(length(selj)==2 & selj[1]==1L, 1, 2)
                .testKnots(res$fit, res$model, whichK, i, treated, v)
            }}))
    }
    pval0 <- lapply(1:length(model$knots0), function(i) pvali(i, FALSE))
    pval1 <- lapply(1:length(model$knots1), function(i) pvali(i))
    p0 <- sapply(model$knots0, function(ki) length(ki)+1)
    p1 <- sapply(model$knots1, function(ki) length(ki)+1)    
    names(pval0) <- names(pval1) <- names(p0) <- names(p1) <- model$nameX
    list(pval0=pval0, pval1=pval1, p0=p0, p1=p1)
}


extract.tlse <- function (model, include.nobs = TRUE, include.nknots = TRUE,
                           include.numcov = TRUE,
                           which=c("ALL","ACE","ACT","ACN","ACE-ACT","ACE-ACN","ACT-ACN"),
                           ...) 
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
        gof.names <- c(gof.names, "Num. knots (Control)", "Num. knots (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
   }
    if (isTRUE(include.numcov)) {
        rs3 <- length(model$model$nameX)
        gof <- c(gof, rs3)
        gof.names <- c(gof.names, "Num. covariates")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    if (isTRUE(include.nobs)) {
        n <- nrow(model$model$data)
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num. obs.")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    tr <- createTexreg(coef.names = names(co), coef = co, se = se, 
        pvalues = pval, gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
    return(tr)
}

setMethod("extract", signature = className("tlse", "causalTLSE"),
          definition = extract.tlse)


predict.tlse <- function (object, interval = c("none", "confidence"), se.fit = FALSE, 
    newdata = NULL, level = 0.95, vcov. = NULL, ...) 
{
    interval <- match.arg(interval)
    model <- object$model
    if (!is.null(newdata)) 
        model$data <- newdata
    else newdata <- model$data
    Z <- model$data[[model$treated]]
    if (is.null(Z)) 
        stop("newdata must contain a treatment assignment variable")
    newdata$Xf1 <- multiSplines(model, TRUE)
    newdata$Xf0 <- multiSplines(model, FALSE)
    X <- model.matrix(model$formY, newdata)
    b <- coef(object$lm.out)
    pr0 <- c(X[Z == 0, ] %*% b)
    pr1 <- c(X[Z == 1, ] %*% b)
    if (se.fit | interval == "confidence") {
        v <- if (is.null(vcov.)) 
            vcov(object$lm.out)
        else vcov.(object$lm.out, ...)
        se0 <- apply(X[Z == 0, ], 1, function(x) sqrt(c(t(x) %*% 
            v %*% x)))
        se1 <- apply(X[Z == 1, ], 1, function(x) sqrt(c(t(x) %*% 
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
        list(treated = list(fit = pr1, se.fit = se1), control = list(fit = pr0, 
            se.fit = se0))
    else list(treated = pr1, control = pr0)
}


plot.tlse <- function (x, y, which = y, interval = c("none", "confidence"), 
                       level = 0.95, newdata = NULL, legendPos = "topright", vcov. = NULL,
                       col0=2, col1=5, ...) 
 {
    interval <- match.arg(interval)
    vnames <- all.vars(x$model$formX)
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
    ind <- order(x$model$data[, which])
    treat <- x$model$treated
    data <- x$model$data[ind, ]
    Z <- data[, treat]
    data[, !(names(data) %in% c(which, treat))] <- sapply(which(!(names(data) %in% 
        c(which, treat))), function(i) rep(mean(data[, i], na.rm = TRUE), 
                                           nrow(data)))
    if (!is.null(newdata)) 
        for (ndi in nd) data[[ndi]] <- newdata[ndi]
    res <- predict(x, interval = interval, level = level, newdata = data, 
                    se.fit = FALSE, vcov. = vcov., ...)
    pr0 <- res$control
    pr1 <- res$treated
    if (interval == "confidence") {
        lty = c(1, 3, 3)
        lwd = c(2, 1, 1)
    }
    else {
        lty = 1
        lwd = 2
    }
    main <- paste("Outcome versus ", which, " using piecewise polynomials", 
        sep = "")
    ylim <- range(c(pr0, pr1))
    matplot(data[Z == 1, which], pr1, col = col1, ylim = ylim, type = "l", 
        lty = lty, lwd = lwd, main = main, ylab = x$model$nameY, 
        xlab = which)
    matplot(data[Z == 0, which], pr0, col = col0, type = "l", lty = lty, 
        lwd = lwd, add = TRUE)
    grid()
    legend(legendPos, c("Treated", "Control"), col = c(col1, col0), 
        lty = lty, lwd = 2, bty = "n")
    invisible()
 }

