
W = as.matrix(read.table("Wnew.dat"))
X = as.matrix(read.table("Xnew.dat"))
Y = scan("Ypnew.dat")
des = as.matrix(read.table("design3934.dat"))

beta = as.matrix(read.table("betadraws.dat"))
mbeta=apply(beta,2,mean)
delta = as.matrix(read.table("deltadraws.dat"))
mdelta=apply(delta,2,mean)
Z = scan("Z.dat")
curve = scan("curve.dat")
P = 100
lower = P + 1
upper = length(Y)
quartz()
plot(Y,type="l",col="gray75")
lines(lower:upper,(W%*%mdelta+X%*%mbeta+Z+curve[lower:upper]),col=4,lwd=2)
lines(lower:upper,W%*%mdelta,col=7)
lines(curve,col=5,lwd=2)
lines(lower:upper,Z,col=3,lwd=2)
lines(lower:upper,X%*%mbeta,col=2,lwd=2)
matlines(des-2,lty=1,col=rep(1:2,each=6))

plot(2*(Mod(fft(Y[lower:upper]-W%*%mdelta - X%*%mbeta-curve[lower:upper])/sqrt(upper-P))[1:((upper-P)/2)])^2,type="h")
plot(2*(Mod(fft(Y[lower:upper]-W%*%mdelta - X%*%mbeta-Z-curve[lower:upper])/sqrt(upper-P))[1:((upper-P)/2)])^2,type="h")

#plot((4001:4400)/20,Y[4101:4500],type="l")
#lines((4001:4400)/20,(W%*%mdelta+X%*%mbeta+Z+curve[101:10820])[4001:4400],type="l",col=2)

acf(Y[lower:upper]-W%*%mdelta-X%*%mbeta-Z-curve[lower:upper],200)
pacf(Y[lower:upper]-W%*%mdelta-X%*%mbeta-Z-curve[lower:upper],200)

ldelta = apply(delta,2,quantile,prob=c(0.025,0.975))
matplot(t(ldelta),type="n")
abline(h=0,col="gray75")
points(mdelta,pch=19,cex=0.5)
segments(1:length(mdelta),ldelta[1,],1:length(mdelta),ldelta[2,],col=2)

fit = (curve[(P+1):length(curve)]+W%*%mdelta+X%*%mbeta+Z)
res = (Y[(P+1):length(curve)]-(curve[(P+1):length(curve)]+W%*%mdelta+X%*%mbeta+Z))
plot(fit,res,pch=19,cex=0.25)

std.res = (res-mean(res))/sd(res)
qqnorm(std.res);abline(0,1,col=2,lwd=2)

#####

l = 5410  - 10820*0.1 + 1  # top 20% high frequencies
u =  5410

1/(2*sum((Mod(fft(Y)/sqrt(10820))^2)[l:u]*(1/10820)))

#####

my.ar.ols = function (x, aic = TRUE, order.max = NULL, na.action = na.fail, 
    demean = TRUE, intercept = demean, series = NULL, ...) 
{
    if (is.null(series)) 
        series <- deparse(substitute(x))
    rescale <- TRUE
    ists <- is.ts(x)
    x <- na.action(as.ts(x))
    if (anyNA(x)) 
        stop("NAs in 'x'")
    if (ists) 
        xtsp <- tsp(x)
    xfreq <- frequency(x)
    x <- as.matrix(x)
    if (!is.numeric(x)) 
        stop("'x' must be numeric")
    n.used <- nrow(x)
    nser <- ncol(x)
    iser <- seq_len(nser)
    if (rescale) {
        sc <- sqrt(drop(apply(x, 2L, var)))
        sc[sc == 0] <- 1
        x <- x/rep.int(sc, rep.int(n.used, nser))
    }
    else sc <- rep.int(1, nser)
    order.max <- if (is.null(order.max)) 
        min(n.used - 1L, floor(10 * log10(n.used)))
    else round(order.max)
    if (order.max < 0L) 
        stop("'order.max' must be >= 0")
    if (order.max >= n.used) 
        stop("'order.max' must be < 'n.used'")
    order.min <- if (aic) 
        0L
    else order.max
    varE <- seA <- A <- vector("list", order.max - order.min + 
        1L)
    xaic <- rep.int(Inf, order.max - order.min + 1L)
    det <- function(x) max(0, prod(diag(qr(x)$qr)) * (-1)^(ncol(x) - 
        1))
    if (demean) {
        xm <- colMeans(x)
        x <- sweep(x, 2L, xm, check.margin = FALSE)
    }
    else xm <- rep.int(0, nser)
    for (m in order.max:order.max) {
        y <- embed(x, m + 1L)
        if (intercept) {
            if (m) 
                X <- cbind(rep.int(1, nrow(y)), y[, (nser + 1L):ncol(y)])
            else X <- as.matrix(rep.int(1, nrow(y)))
        }
        else {
            if (m) 
                X <- y[, (nser + 1L):ncol(y)]
            else X <- matrix(0, nrow(y), 0)
        }
        Y <- t(y[, iser])
        N <- ncol(Y)
        XX <- t(X) %*% X
        rank <- qr(XX)$rank
        if (rank != nrow(XX)) {
            warning(paste("model order: ", m, "singularities in the computation of the projection matrix", 
                "results are only valid up to model order", m - 
                  1L), domain = NA)
            break
        }
        P <- if (ncol(XX) > 0) 
            solve(XX)
        else XX
        A[[m - order.min + 1L]] <- Y %*% X %*% P
        YH <- A[[m - order.min + 1L]] %*% t(X)
        E <- (Y - YH)
        varE[[m - order.min + 1L]] <- tcrossprod(E)/N
        varA <- P %x% (varE[[m - order.min + 1L]])
        seA[[m - order.min + 1L]] <- if (ncol(varA) > 0) 
            sqrt(diag(varA))
        else numeric()
        xaic[m - order.min + 1L] <- n.used * log(det(varE[[m - 
            order.min + 1L]])) + 2 * nser * (nser * m + intercept)
    }
    m <- if (aic) 
        which.max(xaic == min(xaic)) + order.min - 1L
    else order.max
    y <- embed(x, m + 1L)
    AA <- A[[m - order.min + 1L]]
    if (intercept) {
        xint <- AA[, 1L]
        ar <- AA[, -1L]
        X <- if (m) 
            cbind(rep.int(1, nrow(y)), y[, (nser + 1L):ncol(y)])
        else as.matrix(rep.int(1, nrow(y)))
    }
    else {
        X <- if (m) 
            y[, (nser + 1L):ncol(y)]
        else matrix(0, nrow(y), 0L)
        xint <- NULL
        ar <- AA
    }
    Y <- t(y[, iser, drop = FALSE])
    YH <- AA %*% t(X)
    E <- drop(rbind(matrix(NA, m, nser), t(Y - YH)))
    maic <- min(aic)
    xaic <- setNames(if (is.finite(maic)) 
        xaic - min(xaic)
    else ifelse(xaic == maic, 0, Inf), order.min:order.max)
    dim(ar) <- c(nser, nser, m)
    ar <- aperm(ar, c(3L, 1L, 2L))
    ses <- seA[[m - order.min + 1L]]
    if (intercept) {
        sem <- ses[iser]
        ses <- ses[-iser]
    }
    else sem <- rep.int(0, nser)
    dim(ses) <- c(nser, nser, m)
    ses <- aperm(ses, c(3L, 1L, 2L))
    var.pred <- varE[[m - order.min + 1L]]
    if (nser > 1L) {
        snames <- colnames(x)
        dimnames(ses) <- dimnames(ar) <- list(seq_len(m), snames, 
            snames)
        dimnames(var.pred) <- list(snames, snames)
        names(sem) <- colnames(E) <- snames
    }
    if (ists) {
        attr(E, "tsp") <- xtsp
        attr(E, "class") <- "ts"
    }
    if (rescale) {
        xm <- xm * sc
        if (!is.null(xint)) 
            xint <- xint * sc
        aa <- outer(sc, 1/sc)
        if (nser > 1L && m) 
            for (i in seq_len(m)) ar[i, , ] <- ar[i, , ] * aa
        var.pred <- var.pred * drop(outer(sc, sc))
        E <- E * rep.int(sc, rep.int(NROW(E), nser))
        sem <- sem * sc
        if (m) 
            for (i in seq_len(m)) ses[i, , ] <- ses[i, , ] * 
                aa
    }
    res <- list(order = m, ar = ar, var.pred = var.pred, x.mean = xm, 
        x.intercept = xint, aic = xaic, n.used = n.used, order.max = order.max, 
        partialacf = NULL, resid = E, method = "Unconstrained LS", 
        series = series, frequency = xfreq, call = match.call(), 
        asy.se.coef = list(x.mean = sem, ar = drop(ses)))
    class(res) <- "ar"
    #res
    X
}




plot(2*(Mod(fft(res)/sqrt(length(res)))[1:((length(res))/2)])^2,type="h")
pY = as.matrix(read.table("sub_precYdraws00.dat",skip=5001))
d1 = as.matrix(read.table("deltadraws00.dat",skip=5001))
d2 = as.matrix(read.table("delta2draws00.dat",skip=5001))
