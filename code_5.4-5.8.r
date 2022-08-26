# コード5.4：分解アルゴリズム
wtf_1d <- function(cj,  p) {
  K <- length(p)
  Mj <- length(cj)
  Mj1 <- Mj / 2
  cj1 <- numeric(Mj1)
  dj1 <- numeric(Mj1)
  q <- (-1)^(0:(K - 1)) * p[K:1]
  for (l in 1:Mj1) {
    for (k in 1:K) {
      pos <- ((2 * (l - 1) + (k - 1)) %% Mj) + 1
      cj1[l] <- cj1[l] + cj[pos] * p[k]
      dj1[l] <- dj1[l] + cj[pos] * q[k]
    }
  }
  return(rbind(cj1, dj1))
}
# コード5.5：再構成アルゴリズム
iwtf_1d <- function(cj1, dj1, p) {
  K <- length(p)
  Mj1 <- length(cj1)
  Mj <- 2 * Mj1
  cj <- numeric(Mj)
  q <- (-1)^(0:(K - 1)) * p[K:1]
  for (l in 1:(Mj / 2)) {
    for (k in 1:(K / 2)) {
      pos <- (((l - 1) - (k - 1) + Mj1) %% Mj1) + 1
      cj[2 * l - 1] <- cj[2 * l - 1] + p[2 * k - 1] * cj1[pos] + q[2 * k - 1] * dj1[pos]
      cj[2 * l] <- cj[2 * l] + p[2 * k] * cj1[pos] + q[2 * k] * dj1[pos]
    }
  }
  return(cj)
}
# コード5.6：blocks関数
blocks <- function(t) {
  n <- length(t)
  tj <- c(.1, .13, .15, .23, .25, .4, .44, .65, .76, .78, .81)
  hj <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
  f <- numeric(n)
  for (j in 1:length(tj)) {
    f <- f + hj[j] * (1 + sign(t-tj[j])) / 2
  }
  return(f)
}
# コード5.7：閾値関数
thresholding <- function(x, th, type) {
  idx <- which(abs(x) <= th)
  if (type == "s") {
    x[idx] <- 0
    x[-idx] <- x[-idx] - th * sign(x[-idx])
  }
  if (type == "h") {
    x[idx] <- 0
  }
  return(x)
}
# コード5.8
# サンプルの生成
set.seed(2017)
J <- 8; Jmin <- 4; NJ <- 2^J
sigma <- 0.4
t <- (0:(NJ - 1)) / (NJ - 1) 
f <- blocks(t) 
y <- f + rnorm(NJ, 0, sigma)
# Haar のウェーブレットのツースケール係数
fcoef <- c(1, 1) / sqrt(2) 
# 分解アルゴリズムの実行
cj <- y; w <- NULL
for (j in J:(Jmin + 1)) {
  cd <- wtf_1d(cj, fcoef)
  cj1 <- cd[1, ]
  dj1 <- cd[2, ]
  if (j == J) esigma <- median(abs(dj1)) / 0.6745 # 雑音分散の推定値
  w <- append(dj1, w)
  cj <- cj1
}
# 閾値処理の実行
tlev <- sqrt(2.0 * esigma^2 * log(NJ))
w <- thresholding(w, tlev, "h")
# 再構成アルゴリズムの実行
for (j in Jmin:(J - 1)) {
  Nj1 <- 2^j
  dj1 <- w[1:Nj1]
  w <- w[(Nj1 + 1):length(w)]
  cj <- iwtf_1d(cj1, dj1, fcoef)
  cj1 <- cj
}
# 結果の描画
par(omi = c(0, 0, 0, 0))
par(mai = c(0.8, 1.0, 0.2, 0.2))
par(ps = 14)
par(mfrow = c(2, 1))
xlim<-c(0, 1)
ylim<-c(-4, 6)
par(new = F)
plot(t, y, cex = 1, xlim = xlim, ylim = ylim, xlab = "t", ylab = "y and target")
par(new = T)
plot(t, f, type = "l", lwd = 0.5, xlim = xlim, ylim = ylim, ann = F)
par(new = F)
plot(t, f, type = "l", lwd = 0.5, xlim = xlim, ylim = ylim, xlab = "t", ylab = "target and estimate")
par(new = T)
plot(t, cj1, type = "l", lwd = 1.5, xlim = xlim, ylim = ylim, ann = F)
