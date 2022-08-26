# サンプルの生成
set.seed(2017)
n <- 100; sigma <- 0.2
x <- 2*(0:(n - 1)) / (n - 1) - 1
tf <- sin(2 * pi * x)
y <- tf + rnorm(n, 0, sigma)
# デザイン行列の構成
K <- NULL; tau <- 0.02
for (j in 1:n) K <- cbind(K, exp(-(x - x[j])^2 / tau))
# 正則化パラメータの候補集合の生成
rpset <- seq(-8, 2, 2); rpset <- 10^rpset; rpset <- sort(rpset)
# LOOCV 誤差の計算
loocve <- NULL
for (rp in rpset) {
  H <- K %*% solve(t(K) %*% K + rp * diag(n)) %*% t(K)
  f <- H %*% y
  loocve <- c(loocve, sum((y - f)^2 / (1 - diag(H))^2) / n)
}
# LOOCV 誤差を最小とする正則化パラメータ値の特定
minpos <- which.min(loocve); rpmin <- rpset[minpos]
# 各正則化パラメータ値での出力
H <- K  %*% solve(t(K) %*% K + rpmin*diag(n)) %*% t(K)
fest <- H %*% y
H <- K %*% solve(t(K) %*% K + rpset[1]*diag(n)) %*% t(K)
fmin <- H %*% y
H <- K %*% solve(t(K) %*% K + rpset[length(rpset)]*diag(n)) %*% t(K)
fmax <- H %*% y
# 出力の描画
par(omi = c(0, 0, 0, 0))
par(mai = c(0.8, 1.0, 0.2, 0.2))
par(ps = 14)
xlim <- c(-1, 1)
ylim <- c(-1.5, 1.5)
par(mfrow = c(2, 2))
par(new = F)
plot(x, y, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y and target")
par(new = T)
plot(x, tf, type = "l", lwd = 1, col = grey(0.6), xlim = xlim, ylim = ylim, ann = F)
par(new = F)
plot(x, tf, type = "l", lwd = 1, col = grey(0.6), xlim = xlim, ylim = ylim, xlab = "x", ylab = "target and estimate")
par(new = T)
plot(x, fmin, type = "l", lwd = 2, xlim = xlim, ylim = ylim, ann = F)
par(new = F)
plot(x, tf, type = "l", lwd = 1, col = grey(0.6), xlim = xlim, ylim = ylim, xlab = "x", ylab = "target and estimate")
par(new = T)
plot(x, fmax, type = "l", lwd = 2, xlim = xlim, ylim = ylim, ann = F)
par(new = F)
plot(x, tf, type = "l", lwd = 1, col = grey(0.6), xlim = xlim, ylim = ylim, xlab = "x", ylab = "target and estimate")
par(new = T)
plot(x, fest, type = "l", lwd = 2, xlim = xlim, ylim = ylim, ann = F)
