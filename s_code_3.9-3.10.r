# コード3.9
cal_cverr_l2reg <- function(X, y, folds, rpset) {
  n <- length(y)
# サンプルの順番をシャッフルする
  shuf <- sample(1:n, n, replace = F)
  X <- X[shuf, ]
  y <- y[shuf]
# ブロックごとのサンプル数とブロックの始点と終点を計算する
  n1b <- floor(n / folds)
  nB <- rep(n1b, folds)
  if (n %% folds != 0)  nB[1:(n - n1b * folds)] <- nB[1:(n - n1b * folds)] + 1
  st <- c(1, cumsum(nB)[-folds] + 1)
  ed <- cumsum(nB)
# 検証誤差を計算する
  cve <- NULL
  for (rp in rpset) {
    cvetmp <- NULL
    for (j in 1:folds) {
# 推定用セットと検証用セットを構成する
      idx <- (st[j]:ed[j])
      Xe <- X[-idx, ]
      ye <- y[-idx]
      Xv <- X[idx, ]
      yv <- y[idx]
# 推定用セットに対するパラメータ推定
      coef <- solve(t(Xe) %*% Xe + rp * diag(ncol(Xe))) %*% t(Xe) %*% ye
# 検証用セットに対する出力の計算
      fv <- Xv%*%coef
# 検証用セットに対する誤差の計算
      cvetmp <- c(cvetmp, (yv - fv)^2)
    }
    cve <- c(cve, mean(cvetmp))
  }
# 検証誤差を最小化する正則化パラメータ値を求める
  minpos <- which.min(cve) 
  rlist <- list(cve, minpos)
  names(rlist) <- c("cve", "minpos")
  return(rlist)
}
# コード3.10
# 説明変数のサンプルの生成
set.seed(2017)
n <- 50; mmax <- 25; sigma <- 0.4; S <- 1000
x <- 2*(0:(n - 1)) / (n - 1) - 1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tf <- sin(pi * x)
# 正則化パラメータの候補集合の生成
rpset <- seq(-3, 1, 1)
rpset <- c(1*10^rpset, 5*10^rpset)
rpset <- sort(rpset)
freq <- numeric(length(rpset))
folds <- 10
cve <- NULL
# 正則化パラメータ値をクロスバリデーションによって求める数値実験
for (s in 1:S) {
  y <- tf + rnorm(n, 0, sigma) # 目的変数のサンプルの生成
  rlist <- cal_cverr_l2reg(X, y, folds, rpset)
  pos <- rlist$minpos
  cve <- rbind(cve, rlist$cve)
}
# 評価誤差の平均の描画
CVE <- apply(cve, 2, mean);
par(omi = c(0, 0, 0, 0));
par(mai = c(1.0, 1.5, 0.5, 0.5));
par(ps = 18);
xlim <- c(rpset[1], rpset[length(rpset)]);
ylim <- c(0.1, 1.0);
plot(rpset, CVE, type = "p", log = "xy", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "reguralization param.", ylab = "validation error");

