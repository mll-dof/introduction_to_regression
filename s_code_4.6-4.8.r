# コード4.6
softmax <- function(X, beta) {
  ncls <- dim(beta)[2]
  comp <- exp(X %*% beta)
  eta <- comp / (apply(comp, 1, sum))
  return(eta)
}
# コード4.7
mlr_train <- function(X, y, beta, nu, itmax) {
  cat("Training...\n")
  ndim <- dim(beta)[1]
  ncls <- dim(beta)[2]
  n <- dim(X)[1]
  for(it in 1:itmax) { # 学習のループ
    cat("it=", it, "\n")
    eta <- softmax(X, beta)
    for(d in 1:ndim) {
      for(m in 1:ncls) {
        db <- sum(X[, d] * (y[, m] - eta[, m])) # 更新量の計算
        beta[d, m] <- beta[d, m] + nu * db # 更新式の計算
      }
    }
  }
  cat("Done!\n")
  return(beta)
}
# コード4.8
# サンプルの生成
set.seed(2017)
n1 <- 100; n <- 3 * n1
X1 <- cbind(rnorm(n1, 1, 1), rnorm(n1, 1, 1))
y1 <- cbind(rep(1, n1), rep(0, n1), rep(0, n1))
X2 <- cbind(rnorm(n1, -1, 1), rnorm(n1, -1, 1))
y2 <- cbind(rep(0, n1), rep(1, n1), rep(0, n1))
X3 <- cbind(rnorm(n1, 1, 0.5), rnorm(n1, -1, 0.5))
y3 <- cbind(rep(0, n1), rep(0, n1), rep(1, n1))
X <- rbind(X1, X2, X3); X <- cbind(rep(1, n), X)
y <- rbind(y1, y2, y3)
shuf <- sample(1:n, replace = F); X <- X[shuf, ]; y <- y[shuf, ]
# 学習
ndim <- 3; ncls <- 3; n <- 3 * n1
beta <- matrix(runif(ndim*ncls, -1, 1), ndim, ncls)
nu <- 0.01; itmax <- 100
beta <- mlr_train(X, y, beta, nu, itmax)
eta <- softmax(X, beta)# 出力（各クラスの生起確率）の計算
out <- apply(eta, 1, which.max) # 出力に基づきクラス番号を決定
tc <- apply(y, 1, which.max) # サンプルのクラス番号の特定
err <- sum(ifelse(out != tc, 1, 0)) / n # 誤判別率の計算
cat("err=", err, "\n")
# 分類結果の描画
par(new = F)
col <- c(grey(0), grey(0.5), grey(0.8))
mark <- numeric(n)
xlim<-c(-3,3); ylim<-c(-3,3)
for(k in 1:ncls) {
  idx <- which(y[, k] == 1)
  mark[idx] <- k
}
mark <- paste(mark)
for(k in 1:ncls) {
  idx <- which(out == k)
  plot(X[idx, 2], X[idx, 3], pch = mark[idx], col = col[k], xlim = xlim, ylim = ylim)
  par(new = T)
}


