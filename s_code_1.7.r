# �T���v���̐���
set.seed(2017)
n <- 50
x1 <- runif(n, -2, 2)
x2 <- runif(n, -1, 1)
X <- cbind(rep(1, n), x1, x2)
tb <- c(1, 0, -0.5)
tf <- X %*% tb
y <- tf + rnorm(n, 0, 0.4)
# �ŏ���搄��ʂ̌v�Z
V <- solve(t(X) %*% X)
ecoef <- V %*% t(X) %*% y 
# �G�����U�̕s�ΐ���ʂ̌v�Z
esigma2 <- sum((y - X %*% ecoef)^2) / (n - ncol(X))
# �ŏ���搄��ʂ̕��U�̌v�Z
stderr <- sqrt(esigma2 * diag(V))
# t �l����� p �l�̌v�Z
tval <- ecoef / stderr
pval <- 2 * pt(abs(tval), n - ncol(X), lower = F)
cat("pval : ", pval, "\n")
