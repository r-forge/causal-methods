###################################
########  Matching methods
##################################

## units in x2 are matched to each element of x1
## minn is the minimum number of matches
## tol: if other individuals are within tol of 
## the highest distance in the first minn matches, then they are also included.
## y2 is y from group 2
## it returns:
##             list in indices for the match
##             E(y2 | z=1, x) (n1 x 1 vector of estimated missing Y(2))

getNN <- function(x1,x2,y2, minn, tol=1e-7)
    {
        x1 <- as.matrix(x1)
        x2 <- as.matrix(x2)
        n1 <- nrow(x1)
        n2 <- nrow(x2)
        k <- ncol(x1)
        res <- .Fortran(F_findnn, as.double(x1), as.double(x2),as.double(y2),
                        as.double(tol),
                        as.integer(n1), as.integer(n2), as.integer(k),
                        as.integer(minn), K=double(n2), nn=integer(n1),
                        ind=integer(n1*n2), y2hat=double(n1))
        nn <- res$nn
        ind <- matrix(res$ind,nrow=n1)
        list(matches=lapply(1:n1, function(i) ind[i,1:nn[i]]), K=res$K,
             y2hat = res$y2hat)
    }

## Function to get the missing counterfactual in y
## if which="y0", y[z==1] are replace by an estimate Y(0)|z=1
## if which="y1", y[z==0] are replace by an estimate Y(1)|z=0
## if type="match", only match method is used
## if type="reg", only regression is used
## if type="bc_match" the bias-corrected method is used.
## x is the matrix used for matching
## regMat is the regression matrix for E(Y|X,Z=1) or  E(Y|X,Z=0)

getMissingY <- function(y, z, x, which=c("y0", "y1"),
                        type=c("match","bc_match","reg","all"), minn=4,
                        tol=1e-7, regMat=x)
    {
        which <- match.arg(which)
        type <- match.arg(type)
        ## regMat should have an intercept included
        x <- as.matrix(x)
        if (which == "y0")
            sel <- z==1
        else
            sel <- z==0
        x1 <- x[sel,,drop=FALSE]
        x2 <- x[!sel,,drop=FALSE]
        y2 <- y[!sel]
        res <- getNN(x1,x2,y2, minn, tol)
        if (type != "match" | type == "all")
            b <- lm.wfit(regMat[!sel,,drop=FALSE], y[!sel],
                         res$K)$coefficients
        if (type == "match")
            {
                y[sel] <- res$y2hat
                y <- as.matrix(y)
                colnames(y) <- "match"
            } else if (type == "reg") {
                y[sel] <- c(regMat[sel,,drop=FALSE]%*%b)
                y <- as.matrix(y)
                colnames(y) <- "reg"
            } else if (type=="bc_match") {
                fity <- c(regMat%*%b)
                bc <- fity[sel] - sapply(res$matches,
                                         function(l) mean(fity[!sel][l]))
                y[sel] <- res$y2hat + bc
                y <- as.matrix(y)
                colnames(y) <- "bc"
            } else {
                ym <- y
                ym[sel] <- res$y2hat
                yreg <- y
                yreg[sel] <- c(regMat[sel,,drop=FALSE]%*%b)
                ybc <- y
                fity <- c(regMat%*%b)
                bc <- fity[sel] - sapply(res$matches,
                                         function(l) mean(fity[!sel][l]))
                ybc[sel] <- res$y2hat + bc
                y <- cbind(ym, ybc, yreg)
                colnames(y) <- c("match","bc","reg")
            }
        attr(y, "which") <- which
        attr(y, "K") <- res$K
        y
    }

### This function is use to compute sigma(X, Z) = Var(Y|X,Z)
### It uses the Equatrion (14) of Abadie and Imbens (2006)

.sigmaX <- function(y, z, x, minn=4, tol=1e-7)
{
    x <- as.matrix(x)
    yhat <- numeric(length(y))
    yhat[z==1] <- getNN(x[z==1,],x[z==1,],y[z==1], minn, tol)$y2hat
    yhat[z==0] <- getNN(x[z==0,],x[z==0,],y[z==0], minn, tol)$y2hat
    (y-yhat)^2*minn/(minn+1)
}

## PSForm is the formula used to compute the propensity score.
## BCorForm is the bias correction regression model

matching <- function(form,  balm, data, type=c("ACE","ACT","ACC"), M=4,
                      psForm=NULL, bcForm=NULL, vJ=4)
{
    ## If bcForm does not have an intercept it will be added.
    type <- match.arg(type)
    Call <- match.call()
    matchPS <- inherits(psForm, "formula")
    adjust <- inherits(bcForm, "formula")
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    ptype <- type
    if (type == "ACC")
    {
        T <- 1-T
        type <- "ACT"
        ptype <- "ACC"
    }
    BC <- NA
    se.BC <- NA
    if (!matchPS)
    {
        balm <- attr(terms(balm), "term.labels")
        balm <- reformulate(balm, intercept=FALSE)
        X <- model.matrix(balm, data)
    } else {
        res <- glm(psForm, data=data, family=binomial())
        X <- as.matrix(fitted(res))
    }
    X <- scale(X)
    if (adjust)
    {
        bcForm <- reformulate(attr(terms(bcForm), "term.labels"))
        Z <- model.matrix(bcForm, data)
    } else {
        Z <- NULL
    }
    typeCor <- ifelse(adjust, "all", "match")
    sig <- .sigmaX(Y, T, X, vJ)
    if (type == "ACT")
    {
        N1 <- sum(T)
        Y0 <- getMissingY(Y, T, X, which="y0", type=typeCor, minn=M,
                          regMat=Z)
        NN <- mean(Y[T==1]-Y0[T==1,"match"])
        if (ptype == "ACC") NN <- -NN
        K <- attr(Y0, "K")*M
        s2.NN <- sum((Y[T==1]-Y0[T==1,"match"] - NN)^2)/N1
        s2.2 <- sum(K*(K-1)/M^2*sig[T==0])/N1
        se.NN <- sqrt(s2.NN+s2.2)/sqrt(N1)
        if (adjust)
        {
            BC <- mean(Y[T==1]-Y0[T==1,"bc"])
            if (ptype == "ACC") BC <- -BC
            s2.BC <- sum((Y[T==1]-Y0[T==1,"bc"] - BC)^2)/N1
            se.BC <- sqrt(s2.BC+s2.2)/sqrt(N1)
        }       
    } else {
        N <- length(Y)
        Y0 <- getMissingY(Y, T, X, which="y0", type=typeCor, minn=M,
                          regMat=Z)
        Y1 <- getMissingY(Y, T, X, which="y1", type=typeCor, minn=M,
                          regMat=Z)
        NN <- mean(Y1[,"match"]-Y0[,"match"])
        K <- numeric(length(Y))        
        K[T==0] <- attr(Y0, "K")
        K[T==1] <- attr(Y1, "K")
        s2.NN <- mean((Y1[,"match"]-Y0[,"match"] - NN)^2)
        s2.2 <- mean((K^2+(2*M-1)*K/M)*sig)
        se.NN <- sqrt(s2.NN+s2.2)/sqrt(N)
        if (adjust)
        {
            BC <- mean(Y1[,"bc"]-Y0[,"bc"])
            s2.BC <- mean((Y1[,"match"]-Y0[,"match"] - BC)^2)
            se.BC <- sqrt(s2.BC+s2.2)/sqrt(N)
        }
    }
    coef <- c(NN, BC)
    se <- c(se.NN, se.BC)
    names(se) <- names(coef) <- c("NN","BC.NN")
    method <- ifelse(matchPS, "Propensity Score Matching Method",
                         "Covariate Matching Method")
    details <- list(ifelse(matchPS,
                           paste("Propensity Score formula: ",
                                 deparse(psForm), sep=""), 
                           paste("Covariate formula: ",
                                 deparse(balm), sep="")),
                    paste("Mininum number of matches: ", M, sep=""))

    if (adjust)
        details <- c(details, list(paste("Bias-correction formula: ",
                                         deparse(bcForm),sep="")))
    form <- list(balForm=balm, psForm=psForm, bcForm=bcForm)
    ans <- list(estim=coef, se=se, type=type, method=method,
        form=form, details=details, info=list(), data=data, call=Call,
        coefNames = c(type, paste(type,"_BC", sep="")))
    class(ans) <- "altCausal"
    ans
}


### ACE using local linear regression matching
### if h is set, it is used, if not, the optimal h is computed
### by minimizing the CV. A grid in seq(from, to, length=nh) is
### first obtain and the brent method is used to complete.
### In the output, if info=0, it converged, if info=1, the minimum
### is at from or to, so the brent merhod was not launch, and for info = 2
### the brent reached the minimum iteration. 
### the method requires propensity score so PSForm must be provided

LLmatching <- function(form, psForm, data, type=c("ACE","ACT","ACC"),
                       kern=c("Gaussian","Epanechnikov"),tol=1e-4,
                       h=NULL, from=.00001, to=5, ngrid=10, maxit=100,
                       hMethod=c("Brent","Grid"))
{
    type <- match.arg(type)
    Call <- match.call()
    kern <- match.arg(kern)
    hMethod <- match.arg(hMethod)
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    ptype <- type
    if (type == "ACC")
    {
        T <- 1-T
        ptype <- "ACC"
        type <- "ACT"
    }
    res <- glm(psForm, data=data, family=binomial())
    X <- as.matrix(fitted(res))
    if (type == "ACT")
    {
        Y0 <- getMissingY.LL(Y, T, X, "y0", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        ACE <- mean(Y[T==1]-Y0[T==1], na.rm=TRUE)
        if (ptype == "ACC") ACE <- -ACE
        info <- attr(Y0, "h")
    } else {
        Y0 <- getMissingY.LL(Y, T, X, "y0", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        Y1 <- getMissingY.LL(Y, T, X, "y1", h, kern,
                             from, to, ngrid, tol, maxit, hMethod)
        info <- list(Y0=attr(Y0, "h"), Y1=attr(Y1, "h"))
        ACE <- mean(Y1-Y0, na.rm=TRUE)
    }
    estim <- c(LL=ACE)
    method <- "Local-Linear Regression Method"
    details <- list(paste("Kernel: ", kern, sep=""),
                    paste("PS formula: ", deparse(psForm), sep=""),
                    ifelse(is.null(h),
                           paste("Bandwidth select: ",hMethod,sep=""),
                           "Bandwidth select: Fixed"))
    if (type == "ACE")
    {
        tmp <- paste("h (Y(0)) = ", round(info[[1]][[1]], 6), sep="")
        tmp2 <- paste("h (Y(1)) = ", round(info[[2]][[1]], 6), sep="")
        details <- c(details, tmp, tmp2)                                 
    } else {
        tmp <- paste("h = ", round(info[[1]], 6), sep="")
        details <- c(details, tmp)                
    }                        
    if (is.null(h))
    {
        if (type == "ACE")
        {
            tmp <- paste("Convergence code for bandwidth (Y(0)): ",
                         info[[1]][[2]], sep="")
            tmp2 <- paste("Convergence code for bandwidth (Y(1)): ",
                          info[[2]][[2]], sep="")
            details <- c(details, tmp, tmp2)
        } else {
            details <- c(details,
                             list(paste("Convergence code for bandwidth: ",
                                        info[[2]],sep="")))
        }
    }
    if (type == "ACE")
        info <- list(meanh = mean(c(info[[1]][[1]], info[[2]][[1]])),
                     convergence = max(info[[1]][[2]], info[[2]][[2]]))
    se <- as.numeric(NA)
    ans <- list(estim=estim, se=se, type=ptype, method=method,
        form=list(psForm=psForm),data=data, coefNames=type,
        details=details, info=info, call=Call)
    class(ans) <- "altCausal"
    ans
}

## p is the fitted propensity score
## The function estimate the counterfactual Y(0) or Y(1) 

getMissingY.LL <- function(y, z, p, which=c("y0", "y1"), h=NULL,
                           kern=c("Gaussian","Epanechnikov"),
                           from, to, ngrid, tol, maxit,
                           hMethod=c("Brent","Grid"))
{
    which <- match.arg(which)
    hMethod <- match.arg(hMethod)
    kern <- match.arg(kern)
    if (which == "y0")
        sel <- z==1
    else
        sel <- z==0
    missY <- llrF(p[sel], p[!sel], y[!sel], h,
                  kern, from, to, ngrid, tol, maxit, hMethod)
    y[sel] <- missY
    attr(y, "which") <- which
    attr(y, "h") <- attr(missY, "h")
    y
}

cvF <- function(p, y, h, kern=c("Gaussian","Epanechnikov"))
{
    kern <- match.arg(kern)
    nh <- length(h)
    kernF <- tolower(strtrim(kern,1))
    n <- length(p)
    res <- .Fortran(F_cvfct, as.double(p), as.double(y),
                    as.integer(n), as.double(h), as.character(kernF),
                    as.integer(nh), cv=double(nh))
    res$cv
}

llrF <- function(p1, p0, y0, h=NULL, kern=c("Gaussian","Epanechnikov"), 
                 from, to, ngrid, tol, maxit, hMethod=c("Brent","Grid"))
{
    hMethod <- match.arg(hMethod)
    kern <- match.arg(kern)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)
    n1 <- length(p1)
    n2 <- length(p0)
    convergence <- NA
    if (is.null(h))
    {
        resh <- optCVF(p0, y0, kern, from, to, ngrid,
                       tol, maxit, hMethod)
        h <- resh$h
        convergence <- resh$info
    }
    res <- .Fortran(F_llr, as.double(p1), as.double(p0), as.double(y0),
                    as.integer(n1), as.integer(n2), as.double(h),
                    as.integer(kernF), info = integer(n1),
                    y0hat = double(n1))
    y0hat <- res$y0hat
    y0hat[res$info == 1] <- NA
    attr(y0hat, "h") <- list(h=h,convergence=convergence)
    y0hat
}


gridcvF <- function(p, y, kern=c("Gaussian","Epanechnikov"), from, to, nh)
{
    kern <- match.arg(kern)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)
    n <- length(p)
    res <- .Fortran(F_gridcv, as.double(p), as.double(y),
                    as.integer(n), as.integer(kernF),
                    as.double(from), as.double(to), as.integer(nh),
                    info=integer(1), h=double(3), cv=double(3))
    cbind(res$h,res$cv)
}

### The gricv fortran subroutine returns, if info=0,  a 3x1 vector of h and cv, 
### such that cv[1]>cv[2] and cv[3]>cv[2]. The minimum is therefore the middle one.

optCVF <- function(p, y, kern=c("Gaussian","Epanechnikov"), from, to, nh,
                   tol=1e-4, maxit=100, method=c("Brent","Grid"))
{
    kern <- match.arg(kern)
    method <- match.arg(method)
    kernF <- ifelse(kern=="Gaussian", 1L, 2L)        
    n <- length(p)
    res <- .Fortran(F_gridcv, as.double(p), as.double(y),
                    as.integer(n), as.integer(kernF),
                    as.double(from), as.double(to), as.integer(nh),
                    info=integer(1), h=double(3), cv=double(3))
    if (res$info == 1)
    {
        return(list(cv=res$cv[1], h=res$h[1], info=res$info))
    }
    if (method == "Grid")
    {
        return(list(cv=res$cv[2], h=res$h[2], info=res$info))
    }
    res2 <- .Fortran(F_brentcv, as.double(p), as.double(y),
                     as.integer(n), as.integer(kernF),
                     as.double(tol), as.integer(maxit),
                     as.double(res$h), as.double(res$cv),
                     info=integer(1), h=double(1), cv=double(1))
    info <- ifelse(res2$info == 0, 0, 2)
    list(cv=res2$cv, h=res2$h, info=info)
}


### Weighting Methods
########################

ipw <- function(form, psForm, data, type=c("ACE","ACT","ACC"),
                normalized=FALSE, ...)
{
    type <- match.arg(type)
    Call <- match.call()
    mf <- model.frame(form, data)
    Y <- c(mf[[1]])
    T <- c(mf[[2]])
    n <- length(Y)
    mnorm <- normalized
    info <- list()
    ptype <- type
    if (normalized == "GPE")
    {
        res <- getGPE.PS(form, psForm, data, ...)
        X <- res$phat
        info <- res$info
        normalized <- TRUE
    } else {                        
        res <- glm(psForm, data=data, family=binomial())
        X <- as.matrix(fitted(res))
    }
    Rx <- model.matrix(psForm, data)
    if (type == "ACC")
    {
        T <- 1-T
        X <- 1-X
        type <- "ACT"
        ptype <- "ACC"
    }
    w1 <- T/X
    w0 <- (1-T)/(1-X)
    gX <- if (type == "ACE") rep(1, n) else X
    psi <- gX*(w1*Y-w0*Y)
    if (normalized)
        ACE <- sum(Y*((gX*w1)/sum(gX*w1)-(gX*w0)/sum(gX*w0)))
    else
        ACE <- sum(Y*gX*(w1-w0))/sum(gX)
    
    if (ptype == "ACC")
            ACE <- -ACE    
    b <- qr.solve(Rx,Y*(w1^2+w0^2))
    alpha <- -gX*(T-X)*c(Rx%*%b)
    V <- mean((psi-ACE*gX+alpha)^2)*n/sum(gX)^2
    se <- c(IPW=sqrt(V))
    estim <- c(IPW=ACE)
    method <- "Inverse Probability Weight Method"
    details <- list(paste("PS formula: ", deparse(psForm), sep=""))    
    if (mnorm == "GPE")
        details <- c(details,
                     list(paste("Method: normalized with GPE method")))
    else
        details <- c(details,
                     list(paste("Method: ", ifelse(mnorm, "normalized",
                                                   "non-normalized"), sep="")))
    if (mnorm == "GPE")
        details <- c(details,
                     paste("Convergence code of optim for GPE: ", info[2],
                           sep=""))
    ans <- list(estim=estim, se=se, type=ptype, method=method,
                form=list(psForm=psForm), details=details, data=data,
                coefNames = type, info=info, call=Call)
    class(ans) <- "altCausal"
    ans
}


getGPE.PS <- function(form, PSForm=z~kX, data, ...)
    {
        mf <- model.frame(form, data)
        Y <- c(mf[[1]])
        T <- c(mf[[2]])
        X <- model.matrix(PSForm, data)
        g <- function(theta, T, X)
            {
                kX <- c(X%*%theta)
                G <- plogis(kX)
                e <- (T-G)/(1-G)
                m <- as.matrix(e*X)
                mean(colMeans(m, na.rm=TRUE)^2)
            }
        res <- glm(PSForm, data=data, family=binomial())        
        res1 <- optim(coef(res), g, method="BFGS", T=T, X=X, ...)
        b <- res1$par
        fit <- c(X%*%b)
        phat <- plogis(fit)
        phat[phat>.999] <- .999
        list(phat=phat, info=list(obj=res1$value, convergence=res1$convergence))
    }


## print

print.altCausal <- function(x, digits=5, ...) {
    cat(x$method, "\n")
    cat("**************************\n")
    cat("Type of causal effect: ", x$type, "\n", sep="")
    if (length(x$estim)==1)
    {
        cat("Causal estimate: ",
            format(x$estim, digits=digits, ...), "\n")
    } else {
        cat("Causal estimates: \n")
        for (i in 1:length(x$estim))
            cat("\t", names(x$estim)[i], " = ",
                format(x$estim[i], digits=digits, ...),
                "\n", sep="")
    }
    if (length(x$details))
    {
        cat("Details:\n")
        for (i in 1:length(x$details))
            cat("\t", x$details[[i]], "\n")
    }
    invisible(x)
}

summary.altCausal <- function(object, ...)
{
    t <- object$estim/c(object$se) 
    pv <- 2 * pnorm(-abs(t))
    ace <- cbind(object$estim, object$se, t, pv)
    dimnames(ace) <- list(object$coefNames,
                          c("Estimate", 
                            "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(coef=ace,
                type=object$type, method=object$method,
                form=object$form, details=object$details,
                info=object$info)
    class(ans) <- "summaryAltCausal"
    ans
}


print.summaryAltCausal <-  function(x, digits=5,
                                    signif.stars = getOption("show.signif.stars"), ...)
          {
              cat(x$method, "\n", sep="")
              cat(paste(rep("*", nchar(x$method)), collapse=""), "\n")
              cat("Causal Effect type: ", x$type, "\n")
              for (i in 1:length(x$details))
                  cat(x$details[[i]], "\n")
              cat("\n")
              printCoefmat(x$coef, digits = digits,
                           signif.stars = signif.stars, na.print = "NA", ...)
          }

setMethod("extract", signature = className("altCausal", package='causalSLSE'),
          definition = function (model, include.nobs = TRUE, ...) 
          {
              s <- summary(model)
              co <- s$coef[,1]
              se <- s$coef[,2]
              pval <- s$coef[,4]
              names(co) <- names(se) <- names(pval) <- rownames(s$coef)
              gof <- numeric()
              gof.names <- character()
              gof.decimal <- logical()
              if (isTRUE(include.nobs)) {
                  n <- nrow(model$data)
                  gof <- c(gof, n)
                  gof.names <- c(gof.names, "Num. obs.")
                  gof.decimal <- c(gof.decimal, FALSE)
              }
              tr <- createTexreg(coef.names = names(co), coef = co, se = se, 
                                 pvalues = pval, gof.names = gof.names,
                                 gof = gof, gof.decimal = gof.decimal)
              return(tr)
          })

