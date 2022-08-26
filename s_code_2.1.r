# �T���v���̐���
set.seed(2017)
n <- 20
mmax <- 15
sigma <- 0.4
x <- 2 * (0:(n - 1)) / (n - 1) - 1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tb <- rep(0, mmax + 1)
tb[1:3] <- c(2.5, -0.5, -2)
tf <- X %*% tb
y <- tf + rnorm(n, 0, sigma)
# ������ς��Ȃ���ŏ���搄��ʂƎc����敽�ς��v�Z
fhat <- NULL
err <- NULL
for (m in 1:(mmax+1)) {
  Xm <- X[, 1:m]
# �ŏ���搄��ʂ̌v�Z
  bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
  fm <- Xm %*% bm
  fhat <- cbind(fhat, fm)
# �c����敽�ς̌v�Z
  err <- c(err, mean((y - fm)^2))
}
# �T���v���C�^�̉�A���̏o�͂���ѐ��肵���o�͂̕`��
par(omi = c(0, 0, 0, 0))
par(mai = c(0.5, 0.8, 0.1, 0.1))
par(ps = 14)
xlim <- c(-1, 1)
ylim <- c(0, 4)
par(mfrow = c(4, 2))
mcand <- c(0, 1, 2, 3, 6, 9, 12, 15)
for(m in mcand){
    plot(x, tf, type = "l", lty = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
    par(new = T)
    plot(x, fhat[, m+1], type = "l", xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "")
    par(new=T)
    plot(x, y, type = "p", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
    mtext("x", line = 2.3, side = 1, cex = 0.7)
    mtext("y & Output", line = 2.2, side = 2, cex = 0.7)
}

