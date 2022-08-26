# ���l�����̏���
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
# �\���p�T���v���Z�b�g�̐���
Z <- NULL
for (s2 in 1:S2) Z <- cbind(Z, tf + rnorm(n, 0, sigma))
# �e����p�T���v���ɑ΂��čŏ���搄����s���C�\�����덷�Ȃǂ��v�Z����
pe1 <- pe2 <- ee <- NULL
for (s1 in 1:S1) {
# ����p�T���v���̐���
  y <- tf + rnorm(n, 0, sigma)
  PEtmp1 <- PEtmp2 <- EEtmp <- NULL
# ������ς��Ȃ���ŏ���搄��ʂ��v�Z
  for (m in 1:(mmax + 1)) {
    Xm <- X[, 1:m]
    bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
    fm <- Xm %*% bm
    ptmp <- 0
    for (s2 in 1:S2) {
      ptmp <- ptmp + mean((Z[, s2] - fm)^2)
    }
# �\�����덷�̌v�Z
    PEtmp1 <- append(PEtmp1, ptmp / S2)
# ���X�N�̌v�Z
    PEtmp2 <- append(PEtmp2, mean((tf - fm)^2))
# �c����敽�ς̌v�Z
    EEtmp <- append(EEtmp, mean((y - fm)^2))
  }
  pe1 <- rbind(pe1, PEtmp1)
  pe2 <- rbind(pe2, PEtmp2)
  ee <- rbind(ee, EEtmp)
}
PE1 <- apply(pe1, 2, mean)
PE2 <- apply(pe2, 2, mean)
EE <- apply(ee, 2, mean)
# �c����敽�ρC�\�����덷����у��X�N�̕`��
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

