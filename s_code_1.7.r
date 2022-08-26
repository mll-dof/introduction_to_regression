# サンプルの生成
set.seed(2017)
n <- 50
x1 <- runif(n, -2, 2)
x2 <- runif(n, -1, 1)
X <- cbind(rep(1, n), x1, x2)
tb <- c(1, 0, -0.5)
tf <- X %*% tb
y <- tf + rnorm(n, 0, 0.4)
# 最小二乗推定量の計算
V <- solve(t(X) %*% X)
ecoef <- V %*% t(X) %*% y 
# 雑音分散の不偏推定量の計算
esigma2 <- sum((y - X %*% ecoef)^2) / (n - ncol(X))
# 最小二乗推定量の分散の計算
stderr <- sqrt(esigma2 * diag(V))
# t 値および p 値の計算
tval <- ecoef / stderr
pval <- 2 * pt(abs(tval), n - ncol(X), lower = F)
cat("pval : ", pval, "\n")
