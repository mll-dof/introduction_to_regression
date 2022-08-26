# �T���v���̐���
set.seed(2017)
n <- 50
tcoef <- c(-1, 1, 2)
x <- runif(n, -1, 1)
X <- cbind(rep(1, n), x, x^2)
tf <- X %*% tcoef
y <- tf + rnorm(n, 0, 0.5)
# �W���̍ŏ���搄��ʂ̌v�Z
ecoef <- solve(t(X) %*% X) %*% t(X) %*% y
# ����ɗp���Ă��Ȃ��_�ł̉�A���̌v�Z
nn <- 100
u <- seq(-1 , 1 , 2 / (nn - 1))
ef<-ecoef[1] + ecoef[2] * u + ecoef[3] * u^2
# �T���v���Ɛ��肵����A���̕`��
plot(x, y, xlim = c(-1, 1), ylim = c(-4, 2), type = "p", pch = 1, lwd = 0.5, cex = 2, xlab = "x", ylab = "y", col = gray(0))
par(new = T)
plot(u, ef, xlim = c(-1, 1), ylim = c(-4, 2), type = "l", lwd = 0.5, xlab = "", ylab = "", col = gray(0), ann = F)
