# �T���v���̐���
set.seed(202103)
n <- 10
l <- 100
sigma <- 0.4
x <- runif(n, -1, 1)
tf <- sin(pi * x)
y <- tf + rnorm(n, 0, sigma)
u <- 2 * (0:(l - 1)) / (l - 1) - 1
# �n�C�p�[�p�����[�^�̌��W���̍\��
Cans <- as.vector(c(1, 2, 5) %o% c(0.01, 0.1, 1))
Cant <- as.vector(c(1, 2, 5) %o% c(0.01, 0.1, 1))
# ���ׂẴn�C�p�[�p�����[�^�̑g�ɑ΂���ΐ����Ӗޓx�̌v�Z
mll <- NULL
for (esigma in Cans) {
  mlltmp <- NULL
  for (etau in Cant) {
    V11 <- exp(-as.matrix(dist(x))^2 / etau)
    V11 <- V11 + esigma^2 * diag(rep(1, n))
    V11i <- solve(V11)
    L <- -y %*% V11i %*% y / 2 - log(det(V11)) / 2 - n * log(2 * pi)
    mlltmp <- c(mlltmp, L)
  }
  mll <- rbind(mll, mlltmp, deparse.level = 0)
}
# �ΐ����Ӗޓx���ő剻����n�C�p�[�p�����[�^�l�̓���
pos <- which(mll == max(mll), arr.ind = T)
# ���߂��n�C�p�[�p�����[�^�l�ł̕��ϊ֐��̌v�Z
esigma <- Cans[pos[1]]
etau <- Cant[pos[2]]
V <- exp(-as.matrix(dist(c(x, u)))^2 / etau)
V11 <- V[1:n, 1:n]
V12 <- V[1:n, (n+1):ncol(V)]
V22 <- V[(n+1):nrow(V), (n+1):ncol(V)]
Vi <- solve(V11 + esigma^2 * diag(rep(1, n)))
mu <- t(V12) %*% Vi %*% y
Vy <- V22 - t(V12) %*% Vi %*% V12
mup <- mu + sqrt(diag(Vy))
mum <- mu - sqrt(diag(Vy))
# ���ʂ̕`��
par(omi = c(0, 0, 0, 0))
par(mai = c(0.8, 1.0, 0.2, 0.2))
par(ps = 14)
xlim <- c(-1, 1)
ylim <- c(-3, 3)
plot(x, y, xlim = xlim, ylim = ylim, lwd = 0.2, cex = 1.5, xlab = "x,  u", ylab = "y,  output")
par(new = T)
plot(u, mu, type = "l", xlim = xlim, ylim = ylim, xlab = "", ylab = "")
par(new = T)
plot(u, mup, type = "l", col = grey(0.8), xlim = xlim, ylim = ylim, xlab = "", ylab = "")
par(new = T)
plot(u, mum, type = "l", col = grey(0.8), xlim = xlim, ylim = ylim, xlab = "", ylab = "")

