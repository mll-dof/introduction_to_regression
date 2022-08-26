# Soft-Thresholding 関数
st <- function(x, t) {
  return(sign(x) * pmax(abs(x) - t, 0))
}
# 座標降下法
lrcd <- function(X, y, tmax, lambda, alpha) {
  mmax <- ncol(X)
  coef <- numeric(mmax)
  for (t in 1:tmax) {
    for (j in 1:mmax) {
      ytj <- X[, -j] %*% coef[-j]
      gammaj <- mean(X[, j]^2) + lambda*(1 - alpha)
      xij <- mean(X[, j]*(y - ytj))
      coef[j] <- st(xij, lambda * alpha) / gammaj
    }
  }
  return(coef)
}
# 多項式回帰に対する座標降下法の数値例
# サンプルの生成
set.seed(2017)
n <- 100; mmax <- 8; sigma <- 0.1
x <- 2*(0:(n - 1)) / (n - 1) - 1
X <- NULL
for (m in 1:mmax) X <- cbind(X, x^m)
X <- cbind(rep(1, n), X)
tf <- x^2 - 1
y <- tf + rnorm(n, 0, sigma)
# 座標降下法によるパラメータ推定
coef <- lrcd(X, y, 100, 0.01, 1)
fhat <- X %*% coef
# サンプル，真の回帰式の出力および推定された出力の描画
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(-1, 1)
ylim <- c(-1.5, 0.5)
plot(x, fhat, type = "l", lty = "solid", xlim = xlim, ylim = ylim, xlab = "x", ylab = "Data & Output", main = "")
par(new = T)
plot(x, tf, type = "l", lty = "dashed", xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
par(new = T)
plot(x, y, type = "p", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
labels <- c("Data", "True", "Estimated")
pchs <- c(1, -1, -1)
ltys <- c(0, 2, 1)
legend("topright",  legend  =  labels,  pch = pchs, lty = ltys,  bty = "n", y.intersp = 2.0, pt.cex = 2)
