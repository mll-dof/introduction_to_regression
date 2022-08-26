# コード4.2
lr_train <- function(X, y, beta, itmax) {
  for(it in 1:itmax) { # 学習のループ
    eta <- as.vector(1 / (1 + exp(-X %*% beta)))
    d <- (y - eta)
    score <- t(X) %*% d
    W <- diag(eta * (eta - 1))
    hessian <- t(X) %*% W %*%X
    dbeta <- solve(hessian) %*% score # 更新量の計算
    beta <- beta - dbeta # 更新式の計算
  }
  return(beta)
}
# コード4.3
# サンプルの生成
set.seed(2017)
n <- 100 
x <- (0:(n - 1)) / (n - 1) * 4 - 2 
tbeta <- c(0, 2) 
p <- 1/(1 + exp(-(tbeta[1] + tbeta[2] * x)))
y <- rbinom(n, 1, prob = p)
# lr_train による学習
beta <- runif(2, -1, 1) # 初期値
itmax <- 10 # 繰り返し回数
X <- cbind(rep(1, n), x) # デザイン行列
beta <- lr_train(X, y, beta, itmax) 
# サンプル，真の条件付き確率および推定した条件付き確率の描画
ep <- 1 / (1 + exp(-(beta[1] + beta[2] * x)))
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(-2, 2)
ylim <- c(-0.05, 1.05)
plot(x, y, xlim = xlim, ylim = ylim, xlab = "x", ylab = "Data and Prob.")
par(new = T)
plot(x, p, type = "l", lty = "dotted", lwd = 1.0, xlim = xlim, ylim = ylim, ann = F)
par(new = T)
plot(x, ep, type = "l", lwd = 1.5, xlim = xlim, ylim = ylim, ann = F)
