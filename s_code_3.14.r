# greedy関数（コード3.2）
greedy <- function(X, y) {
    mmaxp1 <- ncol(X)
    sidx <- NULL
    cidx <- 1:mmaxp1
    for(m in 1:mmaxp1){
        errtmp <- NULL
        for(k in cidx){
            idx <- c(sidx, k)
            Xm <- X[, idx]
            bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
            fm <- Xm %*% bm
            errtmp <- c(errtmp, mean((y - fm)^2))
        }
        pos <- which.min(errtmp)     
        sidx <- c(sidx, cidx[pos])
        cidx <- cidx[-pos]
    }
    return(sidx)
}
# コード3.14
cal_cverr_greedy <- function(X, y, folds) {
# shuffling order of samples
  shuf <- sample(1:n, replace = F)
  X <- X[shuf, ]
  y <- y[shuf]
# calculation of start and end points of each validation block
  n1b <- floor(n / folds)
  nB <- rep(n1b, folds)
  if (n != n1b * folds) nB[1:(n - n1b * folds)] <- nB[1:(n - n1b * folds)] + 1
  st <- c(1, cumsum(nB)[-folds] + 1)
  ed <- c(cumsum(nB))
# loop for calculation of validation error for each candidate value
  cvej <- numeric(ncol(X))
  for (j in 1:folds) {
# constructing estimation and validation sets
    idx <- (st[j]:ed[j])
    Xe <- X[-idx, ]
    ye <- y[-idx]
    Xv <- X[idx, ]
    yv <- y[idx]
# application of greedy method for estimation set
    sidx <- greedy(Xe, ye)
# estimation of coefficients and calculation of CV error in order of sidx
    for (m in 1:length(sidx)) {
      Xvt <- matrix(Xv[, sidx[1:m]], nrow(Xv), m)
      Xet <- matrix(Xe[, sidx[1:m]], nrow(Xe), m)
      coef <- solve(t(Xet) %*% Xet) %*% t(Xet) %*% ye
      errtmp <- sum((yv - Xvt %*% coef)^2)
      cvej[m] <- cvej[m]+errtmp
    }
  }
  cve <- cvej / nrow(X)
# determine the number of variables, which minimizes CV error
  pos <- which.min(cve)
# application of greedy method for whole sample
  rlist <- greedy(X, y)
    
# estimation of coefficients under determined variables
  sidx <- sidx[1:pos]
  Xpos <- matrix(X[, sidx], nrow(X), pos)
  coef <- solve(t(Xpos) %*% Xpos) %*% t(Xpos) %*% y
  return(list(cve = cve, sidx = sidx, coef = coef))
}
# コード3.14を利用した数値例
# 数値実験の手順はテキストを参照のこと
set.seed(2017)
n<-100; mmax<-8; mmaxp1<-mmax+1; sigma<-0.4; S1<-100; S2<-100
x <- 2*(0:(n-1))/(n-1)-1
X <- rep(1, n)
for(m in 1:mmax) X <- cbind(X, x^m)
tf <- 0.5*sin(pi*x)
pse <- NULL; cve <- NULL; perr <- NULL;
folds <- 10
z <- NULL    
for(k in 1:S2) z <- cbind(z, tf + rnorm(n, 0, sigma))
for(s1 in 1:S1){
    y <- tf + rnorm(n, 0, sigma)
    bf <- solve(t(X) %*% X) %*% t(X) %*% y
    ff <- X %*% bf
    esigma2 <- sum((y - ff)^2)/(n - mmaxp1)
    perrm <- NULL; psem <- NULL
    sidx <- greedy(X, y)
    for(m in 1:mmaxp1){
        idx <- sidx[1:m]
        Xm <- matrix(X[, idx], n, m)
        bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
        fm <- Xm %*% bm
        perrtmp <- NULL
        for(k in 1:S2) perrtmp <- c(perrtmp, mean((z[, k] - fm)^2))
        perrm <- c(perrm, mean(perrtmp))
        psem <- c(psem, mean((y - fm)^2)+2 * esigma2 * m / n)
    }
    perr <- rbind(perr, perrm)
    pse <- rbind(pse, psem)
    ret <- cal_cverr_greedy(X, y, folds)
    cve <- rbind(cve, ret$cve)
}
# PSE，評価誤差および予測二乗誤差の平均の描画
PSE <- apply(pse, 2, mean)
CVE <- apply(cve, 2, mean)
PERR <- apply(perr, 2, mean)
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(1, mmaxp1)
ylim <- c(0.15, 0.25)
plot(1:mmaxp1, PERR, type = "l", xlim = xlim, ylim = ylim, xlab = "Number of Terms", ylab = "Errors", main = "")
par(new = T)
plot(1:mmaxp1, PSE, type = "p", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
par(new = T)
plot(1:mmaxp1, CVE, type = "p", pch = 16, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
labels <- c("Prediction Error", "CV", "PSE")
pchs <- c(-1, 16, 1)
ltys <- c(1, -1, -1)
legend("topright",  legend  =  labels,  pch  =  pchs,  lty  =  ltys, bty = "n", y.intersp = 2.5, pt.cex = 2)
