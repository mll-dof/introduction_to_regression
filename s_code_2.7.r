# ���l�����̏���
set.seed(2017)
n <- 50
mmax <- 8
S1 <- 500
S2 <- 500
sigma <- 0.4
x <- 2*(0:(n - 1)) / (n - 1) - 1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tb <- c(2.5, -0.5, -2, 0, 0, 0, 0, 0, 0)
tf <- X %*% tb
# �\���p�T���v���Z�b�g�̐���
Z <- NULL
for (s2 in 1:S2) Z <- cbind(Z, tf + rnorm(n, 0, sigma))
# �e����p�T���v���ɑ΂��čŏ���搄����s���C�K���Ȃǂ��v�Z
pse <- fpe <- pe <- NULL
for (s1 in 1:S1) {
# ����p�T���v���̐���
  y <- tf + rnorm(n, 0, sigma)
# �G�����U�̕s�ΐ���ʂ̌v�Z
  bf <- solve(t(X) %*% X) %*% t(X) %*% y
  ff <- X %*% bf
  Sf <- sum((y - ff)^2)
  esigma2 <- Sf / (n - (mmax + 1))
# ������ς��Ȃ���ŏ���搄��ʂ��v�Z���C�K���Ȃǂ��v�Z
  PSEtmp <- FPEtmp <- PEtmp <- NULL
  for (m in 1:(mmax + 1)) {
    Xm <- X[, 1:m]
# �ŏ���搄��ʂ̌v�Z
    bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
    fm <- Xm %*% bm
    Sm <- mean((y - fm)^2)
# PSE��FPE�̌v�Z
    PSEtmp <- c(PSEtmp, Sm + 2 * esigma2 * m / n)
    FPEtmp <- c(FPEtmp, Sm * ((n + m) / (n - m)))
# �\�����덷�̌v�Z
    ptmp <- 0
    for (s2 in 1:S2) ptmp <- ptmp + mean((Z[, s2] - fm)^2)
    PEtmp <- c(PEtmp, ptmp / S2)
  }
  pse <- rbind(pse, PSEtmp)
  fpe <- rbind(fpe, FPEtmp)
  pe <- rbind(pe, PEtmp)
}
PSE <- apply(pse, 2, mean)
FPE <- apply(fpe, 2, mean)
PE <- apply(pe, 2, mean)
# �\�����덷�CPSE�����FPE�̕`��
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(0, mmax)
ylim <- c(0.1, 1)
plot(0:mmax, PE, log = "y", type = "l", xlim = xlim, ylim = ylim, xlab = "Order of Polynomial", ylab = "Error & Criteria", main = "")
par(new = T)
plot(0:mmax, PSE, log = "y", type = "p", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
par(new = T)
plot(0:mmax, FPE, log = "y", type = "p", pch = 16, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
labels <- c("Prediction Error", "PSE", "FPE")
pchs <- c(-1, 1, 16)
ltys <- c(1, -1, -1)
legend("topright",  legend  =  labels,  pch  =  pchs,  lty  =  ltys, bty = "n", y.intersp = 2.5, pt.cex = 2)
