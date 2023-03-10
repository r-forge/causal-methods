.splineMatrix <- function (X, knots = NA, pFact = 0.3,
                          deg = 1, method = c("manual", "bs")) 
{
    method <- match.arg(method)
    n <- length(X)
    if (is.null(knots)) 
        return(as.matrix(X))
    if (any(is.na(knots))) {
        p <- floor(n^pFact)
        prop.seq <- seq(from = 0, to = 1, length.out = p + 1)
        prop.seq <- prop.seq[-c(1, p + 1)]
        knots <- quantile(X, probs = prop.seq, type = 1)
    }
    if (method == "bs") {
        Xfi <- bs(x = X, knots = knots, degree = deg)
    }
    else {
        p <- length(knots) + 1
        Xfi <- matrix(nrow = n, ncol = p)
        Xfi[, 1] <- X * (X <= knots[1]) + knots[1] * (X > knots[1])
        Xfi[, p] <- (X - knots[p - 1]) * (X > knots[p - 1])
        if (p >= 3) {
            for (j in 2:(p - 1)) {
                Xfi[, j] <- (X - knots[j - 1]) * (X >= knots[j - 
                  1]) * (X <= knots[j]) + (knots[j] - knots[j - 
                  1]) * (X > knots[j])
            }
        }
        attr(Xfi, "knots") <- knots
    }
    Xfi
}

ppSplines <- function(form, data, knots = NA, pFact = 0.3,
                     deg = 1, method = c("manual", "bs"), minZeroProp=0.1,
                     subGroup=NULL)
{
    X <- model.matrix(form, data)
    if (attr(terms(form), "intercept") == 1)
        X <- X[,-1,drop=FALSE]
    if (!is.null(subGroup))
    {
        if (length(subGroup) != nrow(data))
            stop(paste("The length of the subGroup vector should be ", nrow(data),
                       sep=""))
        if (!is.logical(subGroup))
        {
            if (!all(subGroup %in% c(1,0)))
                stop("subGroup must be logical or binary")
            subGroup <- subGroup == 1
        }
        X <- X[subGroup,,drop=FALSE]
    }
    
    if (!is.list(knots))
    {
        if (length(knots) > 1)
            warning("knots is either a list or a scalar. \nOnly the first element is used")
        knots <- lapply(1:ncol(X), function(i) knots[1])
    }
    isBinary <- apply(X, 2, function(x) all(x %in% c(1,0)))
    isTruncated <- sapply(1:ncol(X),
                          function(i) (!isBinary[i])&(mean(X[,i]==0)>=minZeroProp))    
    if (length(knots) != ncol(X))
        stop(paste("The length of knots must be equal to ", ncol(X), sep=""))
        
    all <- lapply(1:ncol(X), function(i){
        k <- if (isBinary[i]) NULL else knots[[i]]
        if (!isTruncated[i])
        {
            ans <- .splineMatrix(X[,i] , k, pFact, deg, method)
        } else {
            id <- X[,i] != 0
            tmp <- .splineMatrix(X[id,i] , k, pFact, deg, method)
            ans <- matrix(0, nrow(X), ncol(tmp))
            ans[id,] <- tmp
            attr(ans, "knots") <- attr(tmp, "knots")
        }
        if (!is.null(subGroup))
        {
            tmp <- ans
            ans <- matrix(0, nrow(data), ncol(ans))
            ans[subGroup,] <- tmp
            attr(ans, "knots") <- attr(tmp, "knots")
        }
        nk <- length(attr(ans, "knots"))+1
        colnames(ans) <- if (nk==1)
                         {
                             colnames(X)[i]
                         } else {
                             paste(colnames(X)[i], "_", 1:nk, sep="")
                         }
        ans})
    names(all) <- colnames(X)
    knots <- lapply(all, function(li) attr(li, "knots"))
    names(knots) <- names(all)
    cnames <- lapply(all, colnames)
    names(cnames) <- names(all)
    all <- do.call(cbind, all)
    attr(all, "knots") <- knots
    attr(all, "p") <- sapply(knots, length) + 1
    attr(all, "colnames") <- cnames
    all
}

.getPval <- function (form,  data,  ppow, splineMet, HCtype="HC",
                     mZeroProp=0.1, knots0., knots1.) 
{
    tmp <- as.character(form)
    if (!grepl("\\|", tmp[3]))
        stop("form must be of the type y~z|~x")
    tmp2 <- strsplit(tmp[3], "\\|")[[1]]
    formX <- as.formula(tmp2[2], env=.GlobalEnv)
    formY <- as.formula(paste(tmp[2], "~",tmp2[1],sep=""))
    Z <- model.matrix(formY, data)
    Y <- model.frame(formY, data)[[1]]
    if (attr(terms(formY), "intercept") == 1)
        Z <- Z[,-1,drop=FALSE]
    if (ncol(Z)>1)
        stop("The right hand side must be a single vector of treatment indicator")
    if (!all(Z %in% c(0,1)))
        stop("The right hand side must be a binary variable")
    formY <- as.formula(paste(tmp[2], "~factor(",tmp2[1],")+Xf0+Xf1-1"),
                        env=.GlobalEnv)
    n <- length(Z)    
    id0 <- Z == 0

    data$Xf0 <- ppSplines(form=formX, data=data, pFact=ppow, knots=knots0., 
                          method=splineMet, subGroup=id0, minZeroProp=mZeroProp)
    data$Xf1 <- ppSplines(form=formX, data=data, pFact=ppow, knots=knots1.,
                          method=splineMet, subGroup=!id0, minZeroProp=mZeroProp)
    fit <- lm(formY, data)
    naCoef <- is.na(coef(fit))
    if (any(naCoef))       
        warning(paste("\nThe coefficients of the following variables are NA's:\n",
                      paste(names(coef(fit)[naCoef]), collapse=", "), "\n", sep=""))
    knots0 <- attr(data$Xf0, "knots")
    knots1 <- attr(data$Xf1, "knots")
    p0 <- attr(data$Xf0, "p")
    p1 <- attr(data$Xf1, "p")
    cn0 <- attr(data$Xf0, "colnames")
    cn1 <- attr(data$Xf1, "colnames")
    v <- vcovHC(fit, type = HCtype)
    pval0 <- lapply(1:length(p0), function(i){
        if (p0[i] == 1)
            return(NA)
        sapply(1:(p0[i] - 1),  function(j){
            t <- c(paste("Xf0", cn0[[i]][j], sep=""),
                   paste("Xf0", cn0[[i]][j+1], sep=""))
            .lintest(fit, t[1], t[2], v)
        })})
    names(pval0) <- NULL
    names(pval0) <- names(p0)
    pval1 <- lapply(1:length(p1), function(i){
        if (p1[i] == 1)
            return(NA)        
        sapply(1:(p1[i] - 1),  function(j){
            t <- c(paste("Xf1", cn1[[i]][j], sep=""),
                   paste("Xf1", cn1[[i]][j+1], sep=""))
            .lintest(fit, t[1], t[2], v)            
        })})
    names(pval1) <- NULL
    names(pval1) <- names(p1)
    list(pval0 = pval0, pval1 = pval1, p0 = p0, p1 = p1, knots0 = knots0, 
         knots1 = knots1, Z=Z, Y=Y, formX=formX, formY=formY,
         treatment=colnames(Z))
}

.lintest <- function(obj, c1, c2, v)
{
    b <- coef(obj)
    b <- na.omit(b)
    c1 <- which(names(b) == c1)
    c2 <- which(names(b) == c2)
    if (length(c(c1,c2))<2)
        return(NA)
    s2 <- v[c1,c1]+v[c2,c2]-2*v[c1,c2]
    ans <- 1-pf((b[c1]-b[c2])^2/s2, 1, obj$df)
    names(ans) <- NULL
    ans
}


selASY <- function (form, data, pFact = 0.3, splineMet = c("manual", "bs"),
                    HCtype="HC", mZeroProp=0.1, knots0=NA, knots1=NA,
                    minPV=function(p) 1/(p*log(p)))
{
    splineMet <- match.arg(splineMet)
    res <- .getPval(form, data, pFact, splineMet, HCtype, mZeroProp, knots0,
                    knots1)
    pval <- c(do.call("c", res$pval0), do.call("c", res$pval1))
    n <- nrow(data)
    q <- length(pval)
    p <- sum(res$p0) + sum(res$p1)
    crit <- minPV(p)
    id0 <- res$Z == 0
    knots0 <- lapply(1:length(res$pval0), function(i)
    {      
        jhat <- res$pval0[[i]]  <= crit
        ans <- if (all(!jhat) | any(is.na(jhat)))
               {
                   NULL
               } else {
                   res$knots0[[i]][jhat]
               }
    })
    knots1 <- lapply(1:length(res$pval1), function(i)
    {
        jhat <- res$pval1[[i]]  <= crit
        ans <- if (all(!jhat) | any(is.na(jhat)))
               {
                   NULL
               } else {
                   res$knots1[[i]][jhat]
               }
    })
    Xf0 <- ppSplines(res$formX, data, knots0, pFact, 1, splineMet,
                     subGroup=id0, minZeroProp=mZeroProp)
    Xf1 <- ppSplines(res$formX, data, knots1, pFact, 1, splineMet,
                     subGroup=!id0, minZeroProp=mZeroProp)
    names(knots0) <- names(attr(Xf0, "colnames"))
    names(knots1) <- names(attr(Xf1, "colnames"))        
    list(Xf1 = Xf1, Xf0 = Xf0, knots0 = knots0, knots1 = knots1, 
         pval = pval, id0=id0, formY=res$formY, formX=res$formX,
         treatment=res$treatment)
}

.getCV <- function(Y, Z, X0, X1)
{
    n <- length(Y)
    X0 <- as.matrix(X0)
    X1 <- as.matrix(X1)
    myK <- floor(log(n))
    cv.outi <- cvFolds(n, K = myK)
    mspe_pred <- rep(NA, myK)
    X <- cbind(1-Z,Z,X0,X1)
    for(k in 1 : myK)
    {
        id.train <- cv.outi$subsets[cv.outi$which != k]
        id.valid <- cv.outi$subsets[cv.outi$which == k]
        lm.outk <- lm.fit(X[id.train,],Y[id.train])
        sel <- !is.na(lm.outk$coefficients)
        pred.outk <- c(X[id.valid,sel]%*%lm.outk$coefficients[sel])
        mspe_pred[k] <- mean((Y[id.valid]-pred.outk)^2, na.rm=TRUE)
    }    
    mean(mspe_pred, na.rm = TRUE)
}

selIC <- function(form, data, pFact = 0.3, type=c("AIC", "BIC", "CV"),
                  splineMet = c("manual", "bs"), HCtype="HC",
                  mZeroProp=0.1, knots0=NA, knots1=NA) 
{
    type <- match.arg(type)
    splineMet <- match.arg(splineMet)
    res <- .getPval(form, data, pFact, splineMet, HCtype, mZeroProp, knots0,
                    knots1)
    pval <- c(do.call("c", res$pval0), do.call("c", res$pval1))
    n <- nrow(data)
    q <- length(pval)
    p <- sum(res$p0) + sum(res$p1)
    id0 <- res$Z == 0
    pval_sort <- sort(pval)
    data$Xf0 <- ppSplines(form = res$formX, data = data, knots = NULL,
                                 pFact = pFact, deg = 1, method = splineMet,
                                 subGroup=id0, minZeroProp=mZeroProp)
    data$Xf1 <- ppSplines(form = res$formX, data = data, knots = NULL,
                                 pFact = pFact, deg = 1, method = splineMet,
                                 subGroup=!id0, minZeroProp=mZeroProp)
    if (type != "CV")
    {
        fit0 <- lm(res$formY, data)
        icV <- ic_seq0 <- get(type)(fit0)
    } else {
        icV <- ic_seq0 <- .getCV(res$Y, res$Z, data$Xf0, data$Xf1)
    }
    Xf0 <- data$Xf0
    Xf1 <- data$Xf1    
    knots0 <- lapply(1:length(res$pval0), function(i) NULL)
    knots1 <- lapply(1:length(res$pval1), function(i) NULL)
    for(i in 1 : q)
    {
        knots02 <- lapply(1:length(res$pval0), function(j) {
            jhat <- res$pval0[[j]] <= pval_sort[i]
            if(all(!jhat) | any(is.na(jhat))) NULL else res$knots0[[j]][jhat]
        })
        knots12 <- lapply(1:length(res$pval1), function(j) {
            jhat <- res$pval1[[j]] <= pval_sort[i]
            if(all(!jhat) | any(is.na(jhat))) NULL else res$knots1[[j]][jhat]
        })
        data$Xf0 <- ppSplines(form = res$formX, data = data, knots = knots02,
                              pFact = pFact, deg = 1, method = splineMet,
                              subGroup=id0, minZeroProp=mZeroProp)
        data$Xf1 <- ppSplines(form = res$formX, data = data, knots = knots12,
                              pFact = pFact, deg = 1, method = splineMet,
                              subGroup=!id0, minZeroProp=mZeroProp)
        if (type != "CV")
        {
            fit1 <- lm(res$formY, data)
            ic_seq1 <- get(type)(fit1)
        } else {
            ic_seq1 <- .getCV(res$Y, res$Z, data$Xf0, data$Xf1)
        }
        icV <- c(icV, ic_seq1)
        if (ic_seq1<ic_seq0)
        {
            ic_seq0 <- ic_seq1
            Xf0 <- data$Xf0
            Xf1 <- data$Xf1
            knots0 <- knots02
            knots1 <- knots12
        } 
    }
    names(knots0) <- names(attr(Xf0, "colnames"))
    names(knots1) <- names(attr(Xf1, "colnames"))           
    list(Xf1 = Xf1, Xf0 = Xf0, knots0 = knots0, knots1 = knots1, 
         IC=icV, pval = pval, id0=id0, formY=res$formY, formX=res$formX,
         treatment=res$treatment)
}

.selNONE <- function(form, data, pFact = 0.3, type=c("AIC", "BIC", "CV"),
                    splineMet = c("manual", "bs"), HCtype="HC",
                    mZeroProp=0.1, knots0=NA, knots1=NA) 
{
    type <- match.arg(type)
    splineMet <- match.arg(splineMet)    
    tmp <- as.character(form)
    if (!grepl("\\|", tmp[3]))
        stop("form must be of the type y~z|~x")
    tmp2 <- strsplit(tmp[3], "\\|")[[1]]
    formX <- as.formula(tmp2[2], env=.GlobalEnv)
    formY <- as.formula(paste(tmp[2], "~",tmp2[1],sep=""))
    Z <- model.matrix(formY, data)
    if (attr(terms(formY), "intercept") == 1)
        Z <- Z[,-1,drop=FALSE]
    id0 <- Z == 0   
    formY <- as.formula(paste(tmp[2], "~factor(",tmp2[1],")+Xf0+Xf1-1"),
                        env=.GlobalEnv)
    Xf0 <- ppSplines(form = formX, data = data, knots = knots0,
                     pFact = pFact, deg = 1, method = splineMet,
                     subGroup=id0, minZeroProp=mZeroProp)
    Xf1 <- ppSplines(form = formX, data = data, knots = knots1,
                     pFact = pFact, deg = 1, method = splineMet,
                     subGroup=!id0, minZeroProp=mZeroProp)
    knots0 <- attr(Xf0, "knots")
    knots1 <- attr(Xf1, "knots")
    list(Xf1 = Xf1, Xf0 = Xf0, knots0 = knots0, knots1 = knots1, 
         pval = NULL, id0=id0, formY=formY, formX=formX,
         treatment=colnames(Z))
}


otlse <- function(form, data, crit = c("ASY", "AIC", "BIC", "CV", "NONE"),
                  pFact=0.3, splineMet=c("manual","bs"), HCtype="HC",
                  mZeroProp=0.1, knots0=NA, knots1=NA, ...)
{
    crit <- match.arg(crit)
    splineMet <- match.arg(splineMet)
    optBasis <- switch(crit,
                       ASY = selASY(form, data, pFact, splineMet, HCtype, mZeroProp,
                                    knots0, knots1, ...),
                       AIC = selIC(form, data, pFact, "AIC", splineMet, HCtype, mZeroProp,
                                   knots0, knots1),
                       BIC = selIC(form, data, pFact, "BIC", splineMet, HCtype, mZeroProp,
                                   knots0, knots1),
                       CV = selIC(form, data, pFact, "CV", splineMet, HCtype, mZeroProp,
                                  knots0, knots1),
                       NONE = .selNONE(form, data, pFact, "CV", splineMet, HCtype, mZeroProp,
                                      knots0, knots1)) 
    data2 <- data
    data2$Xf0 <- optBasis$Xf0
    data2$Xf1 <- optBasis$Xf1
    lm.out <- lm(optBasis$formY, data2, na.action="na.exclude")
    n <- nrow(data2)
    id0 <- optBasis$id0
    n0 <- sum(id0)
    n1 <- n-n0
    knots0 <- optBasis$knots0
    knots1 <- optBasis$knots1
    pval <- optBasis$pval
    p00 <- ncol(optBasis$Xf0)
    p10 <- ncol(optBasis$Xf1)
    vcov <- vcovHC(lm.out, type=HCtype)
    idb0 <- 3 : (p00 + 2)
    idb1 <- (p00 + 3) : (p00 + p10 + 2)
    beta <- coef(lm.out)
    se.beta <- sqrt(diag(vcov))
    if (any(is.na(beta)))
        warning(paste("\nThe final regression is multicoliear. The result may not be valid:",
                      "\nThe following variables produced NA's\n",
                      paste(names(beta)[is.na(beta)], collapse=", ")), "\n", sep="")
    notNA0 <- !is.na(beta[idb0])
    notNA1 <- !is.na(beta[idb1])
    beta0 <- na.omit(beta[idb0])
    beta1 <- na.omit(beta[idb1])
    p00 <- length(beta0)
    p10 <- length(beta1)
    idb0 <- 3 : (p00 + 2)
    idb1 <- (p00 + 3) : (p00 + p10 + 2)
    
    ## ACE

    X0 <- ppSplines(optBasis$formX, data, knots0, pFact, 1, splineMet, mZeroProp,
                    NULL)
    X1 <- ppSplines(optBasis$formX, data, knots1, pFact, 1, splineMet, mZeroProp,
                    NULL)    

    Xbar0 <- apply(X0, 2, mean)[notNA0]
    Xbar1 <- apply(X1, 2, mean)[notNA1]
    vcovXf0 <- cov(X0[,notNA0,drop=FALSE])
    vcovXf1 <- cov(X1[,notNA1,drop=FALSE])
    vcovXf01 <- cov(X0[,notNA0,drop=FALSE], X1[,notNA1,drop=FALSE])
    ace <- c(beta[2] - beta[1] + crossprod(beta1, Xbar1) - 
             crossprod(beta0, Xbar0))
    
    Dvec <- rep(0, p00 + p10 + 2)
    Dvec[1] <- - 1
    Dvec[2] <- 1
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.ace <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta0, crossprod(vcovXf0, beta0)) / n + 
                 crossprod(beta1, crossprod(vcovXf1, beta1)) / n -
                 2 * beta0 %*% vcovXf01 %*% beta1 / n )^.5)

    ## ACT
  
    Xbar0 <- apply(X0[!id0, notNA0, drop = FALSE], 2, mean)
    Xbar1 <- apply(X1[!id0, notNA1, drop = FALSE], 2, mean)
    vcovXf0 <- cov(X0[!id0, notNA0, drop = FALSE])
    vcovXf1 <- cov(X1[!id0, notNA1, drop = FALSE])
    vcovXf01 <- cov(X0[!id0, notNA0, drop = FALSE], X1[!id0, notNA1, drop = FALSE])
    act <- c(beta[2] - beta[1] + crossprod(beta1, Xbar1) - 
             crossprod(beta0, Xbar0))
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.act <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta0, crossprod(vcovXf0, beta0)) / n1 + 
                 crossprod(beta1, crossprod(vcovXf1, beta1)) / n1 -
                 2 * beta0 %*% vcovXf01 %*% beta1 / n1)^.5)
   
    ## ACN
    
    Xbar0 <- apply(X0[id0, notNA0, drop = FALSE], 2, mean)
    Xbar1 <- apply(X1[id0, notNA1, drop = FALSE], 2, mean)
    vcovXf0 <- cov(X0[id0, notNA0, drop = FALSE])
    vcovXf1 <- cov(X1[id0, notNA1, drop = FALSE])
    vcovXf01 <- cov(X0[id0, notNA0, drop = FALSE], X1[id0, notNA1, drop = FALSE])
    acn <- c(beta[2] - beta[1] + crossprod(beta1, Xbar1) - 
             crossprod(beta0, Xbar0))
    Dvec[idb0] <- - Xbar0
    Dvec[idb1] <- Xbar1
    se.acn <- c((crossprod(Dvec, crossprod(vcov, Dvec)) + 
                 crossprod(beta0, crossprod(vcovXf0, beta0)) / n0 + 
                 crossprod(beta1, crossprod(vcovXf1, beta1)) / n0 -
                 2 * beta0 %*% vcovXf01 %*% beta1 / n0)^.5)
    ans <- list(beta = beta, se.beta = se.beta, 
                lm.out = lm.out, ace = ace, se.ace = se.ace, 
                act = act, se.act = se.act, acn = acn, se.acn = se.acn,
                knots0 = knots0, knots1 = knots1, pval = pval, crit=crit,
                data=data, formX=optBasis$formX, formY=optBasis$formY,
                pFact=pFact, splineMet=splineMet, HCtype=HCtype,
                mZeroProp=mZeroProp, treatment=optBasis$treatment)
    class(ans) <- "otlse"
    ans
}

print.otlse <- function (x, ...)
{
    cat("Causal Effect using Optimal Thresholding Least Squares\n")
    cat("******************************************************\n")
    cat("Selection method: ", x$crit, "\n\n", sep="")
    cat("ACE = ", x$ace, "\nACT = ", x$act, "\nACN = ", x$acn,"\n")
}


summary.otlse <- function(object, ...)
{
    t <- c(object$ace, object$act, object$acn)/
        c(object$se.ace, object$se.act, object$se.acn)
    pv <- 2*pnorm(-abs(t))
    ace <- cbind(c(object$ace, object$act, object$acn),
                 c(object$se.ace, object$se.act, object$se.acn),
                 t, pv)
    dimnames(ace) <- list(c("ACE","ACT","ACN"),
                           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    nb0 <- ifelse(is.null(object$knots0), 1, length(object$knots0))
    nb1 <- ifelse(is.null(object$knots1), 1, length(object$knots1))
    t <- object$beta/object$se.beta
    pv <- 2*pnorm(-abs(t))
    beta <- cbind(object$beta, object$se.beta, t, pv)
    colnames(beta) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")    
    ans <- list(causal=ace, beta=beta, crit=object$crit, knots0=object$knots0,
                knots1=object$knots1, covNames = names(object$knots0))
    class(ans) <- "summary.otlse"
    ans
}

print.summary.otlse <- function(x, digits = 4,
                                signif.stars = getOption("show.signif.stars"),
                                beta=FALSE, knots=FALSE, ...)
{
    cat("Causal Effect using Optimal Thresholding Least Squares\n")
    cat("******************************************************\n")
    cat("Selection method: ", x$crit, "\n\n", sep="")    
    printCoefmat(x$causal, digits = digits, signif.stars = signif.stars, 
                 na.print = "NA", ...)
    if (beta)
    {
        cat("\nPiecewise polynomials coefficients\n")
        cat("**********************************\n")
        printCoefmat(x$beta, digits = digits, signif.stars = signif.stars, 
                     na.print = "NA", ...)        
    }
    if (knots)
    {
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

plot.otlse <- function(x, y, which=y, addInterval=FALSE, level=0.95, 
                       leg.pos="topright", ...)
{
    fit <- x$lm.out
    vnames <- all.vars(x$formX)
    if (is.numeric(which))
        which <- vnames[which]
    if (!is.character(which) & length(which) != 1)
        stop("which must be a character type")
    if (!(which %in% vnames))
        stop("which must be one of the names of the variables")
    data <- x$data[order(x$data[,which]),]
    treat <- x$treatment
    Z <- data[,treat]
    data[,!(names(data)%in%c(which, treat))] <-
        sapply(which(!(names(data)%in%c(which, treat))), function(i)
            rep(mean(data[,i], na.rm=TRUE), nrow(data)))
    data$Xf0 <- ppSplines(x$formX, data, x$knots0, x$pFact, 1,
                          x$splineMet, x$mZeroProp, Z==0)
    data$Xf1 <- ppSplines(x$formX, data, x$knots1, x$pFact, 1,
                          x$splineMet, x$mZeroProp, Z==1)
    main <- paste("Outcome versus ", which, " using piecewise polynomials", sep="")
    if (addInterval)
    {
        int <- "confidence"
        lty=c(1,3,3)
        lwd=c(2,1,1)
    } else {
        int <- "none"
        lty=1
        lwd=2
    }
    pr.t <- predict(x$lm.out, newdata=data[Z==1,], interval=int, level=level)
    pr.c <- predict(x$lm.out, newdata=data[Z==0,], interval=int, level=level)    
    ylim <- range(c(pr.c, pr.t))
    matplot(data[Z==1,which], pr.t, col=2, ylim=ylim, type='l',
            lty=lty, lwd=lwd, main=main, ylab="Outcome", xlab=which)
    matplot(data[Z==0,which], pr.c, col=5, type='l',
            lty=lty, lwd=lwd, add=TRUE)
    grid()
    legend(leg.pos, c("Treated","Control"), col=c(2,5), lty=lty, lwd=2,
           bty='n')
    invisible()
}  


extract.otlse <- function (model, include.nobs = TRUE, include.nknots = TRUE,
                           include.numcov = TRUE,
                           which=c("ALL","ACE","ACT","ACN","ACE-ACT","ACE-ACN","ACT-ACN"),
                           ...) 
{
    which <- match.arg(which)
    type <- c("ACE","ACT","ACN")
    w <- if (which == "ALL") type else type[sapply(type, function(ti) grepl(ti, which))]
    wl <- tolower(w)
    co <- unlist(model[wl])
    names(co) <- toupper(names(co))
    se <- unlist(model[paste("se.",wl,sep="")])
    names(se) <- toupper(names(se))
    pval <- 2*pnorm(-abs(co/se))
    names(pval) <- toupper(names(pval))
    gof <- numeric()
    gof.names <- character()
    gof.decimal <- logical()
    if (isTRUE(include.nknots)) {
        rs1 <- length(unlist(model$knots0))
        rs2 <- length(unlist(model$knots1))        
        gof <- c(gof, rs1, rs2)
        gof.names <- c(gof.names, "Num. knots (Control)", "Num. knots (Treated)")
        gof.decimal <- c(gof.decimal, FALSE, FALSE)
   }
    if (isTRUE(include.numcov)) {
        rs3 <- length(model$knots0)
        gof <- c(gof, rs3)
        gof.names <- c(gof.names, "Num. covariates")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    if (isTRUE(include.nobs)) {
        n <- nrow(model$data)
        gof <- c(gof, n)
        gof.names <- c(gof.names, "Num. obs.")
        gof.decimal <- c(gof.decimal, FALSE)
    }
    tr <- createTexreg(coef.names = names(co), coef = co, se = se, 
        pvalues = pval, gof.names = gof.names, gof = gof, gof.decimal = gof.decimal)
    return(tr)
}

setMethod("extract", signature = className("otlse", "causalOTLSE"),
          definition = extract.otlse)

