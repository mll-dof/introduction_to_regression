# 数値実験の準備
set.seed(2017)
n <- 20 
mmax <- 15 
S1 <- 100
S2 <- 100
sigma <- 0.4
x <- 2 * (0:(n - 1)) / (n - 1) - 1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tb <- rep(0, mmax + 1)
tb[1:3] <- c(2.5, -0.5, -2)
tf <- X %*% tb
# 予測用サンプルセットの生成
Z <- NULL
for (s2 in 1:S2) Z <- cbind(Z, tf + rnorm(n, 0, sigma))
# 各推定用サンプルに対して最小二乗推定を行い，予測二乗誤差などを計算する
pe1 <- pe2 <- ee <- NULL
for (s1 in 1:S1) {
# 推定用サンプルの生成
  y <- tf + rnorm(n, 0, sigma)
  PEtmp1 <- PEtmp2 <- EEtmp <- NULL
# 次数を変えながら最小二乗推定量を計算
  for (m in 1:(mmax + 1)) {
    Xm <- X[, 1:m]
    bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
    fm <- Xm %*% bm
    ptmp <- 0
    for (s2 in 1:S2) {
      ptmp <- ptmp + mean((Z[, s2] - fm)^2)
    }
# 予測二乗誤差の計算
    PEtmp1 <- append(PEtmp1, ptmp / S2)
# リスクの計算
    PEtmp2 <- append(PEtmp2, mean((tf - fm)^2))
# 残差二乗平均の計算
    EEtmp <- append(EEtmp, mean((y - fm)^2))
  }
  pe1 <- rbind(pe1, PEtmp1)
  pe2 <- rbind(pe2, PEtmp2)
  ee <- rbind(ee, EEtmp)
}
PE1 <- apply(pe1, 2, mean)
PE2 <- apply(pe2, 2, mean)
EE <- apply(ee, 2, mean)
# 残差二乗平均，予測二乗誤差およびリスクの描画
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(0, mmax)
ylim<- c(0.01, 2)
ltys <- c(1, 2, -1)
pchs <- c(-1, -1, 1)
plot(0:mmax, PE1, log = "y", type = "l", lty = ltys[1], xlim = xlim, ylim = ylim, xlab = "Variables", ylab = "Error", main = "")
par(new = T)
plot(0:mmax, PE2, log = "y", type = "l", lty = ltys[2], xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
par(new = T)
plot(0:mmax, EE, log = "y", pch = pchs[3], cex = 1.5, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
labels <- c("Prediction Error", "Risk", "Residual")
legend("topright",  legend  =  labels,  pch = pchs, lty = ltys,  bty = "n", y.intersp = 2.5, pt.cex = 2)

