# �R�[�h3.2
greedy <- function(X, y) {
  mmaxp1 <- ncol(X)
  sidx <- NULL
  cidx <- 1:mmaxp1
# �X�e�b�v�Ɋւ��郋�[�v
  for (m in 1:mmaxp1) {
    errtmp <- NULL
# �I�΂�Ă��Ȃ��ϐ��̃C���f�b�N�X�Ɋւ��郋�[�v
    for (k in cidx) {
# cidx �̗v�f�� sidx �Ɉꎞ�I�� 1 �ǉ����čŏ���搄��ʂ��v�Z
      idx <- c(sidx, k)
      Xm <- X[, idx]
      bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
      fm <- Xm %*% bm
# cidx �̊e�v�f�ɑ΂��Ďc����敽�ς��i�[
      errtmp <- c(errtmp, mean((y - fm)^2))
    }
# �c����敽�ς��ŏ������� cidx �̗v�f�̃C���f�b�N�X�̌���
    pos <- which.min(errtmp)   
# sidx �� cidx �̃A�b�v�f�[�g
    sidx <- c(sidx, cidx[pos])
    cidx <- cidx[-pos]
  }
  return(sidx)
}
# �R�[�h3.3
set.seed(2017)
n <- 50; mmax <- 8; mmaxp1 <- mmax+1; sigma <- 0.4; S1 <- 500; S2 <- 500
x <- 2*(0:(n-1))/(n-1)-1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tf <- sin(pi * x)
Z <- NULL
for (s2 in 1:S2) Z <- cbind(Z, tf + rnorm(n, 0, sigma))
perr <- NULL; pse <- NULL
for (s1 in 1:S1) {
  y <- tf + rnorm(n, 0, sigma)
  bf <- solve(t(X) %*% X) %*% t(X) %*% y
  ff <- X %*% bf
  esigma2 <- sum((y - ff)^2)/(n - mmaxp1)
# �×~�@�ɂ���ĕ��ג������C���f�b�N�X sidx �̎擾
  sidx <- greedy(X, y)
  psem <- NULL
  perrm <- NULL 
# sidx �̏��Ԃōŏ���搄��ʂ��v�Z���C PSE ����ї\�����덷���v�Z
  for (m in 1:mmaxp1) {
    idx <- sidx[1:m]
    Xm <- matrix(X[, idx], n, m)
    bm <- solve(t(Xm) %*% Xm) %*% t(Xm) %*% y
    fm <- Xm %*% bm
    psetmp <- mean((y - fm)^2)+2*m * esigma2 / n
    psem <- c(psem, psetmp)
    perrtmp <- NULL
    for (s2 in 1:S2) perrtmp <- c(perrtmp, mean((Z[, s2] - fm)^2))
    perrm <- c(perrm, mean(perrtmp))
  }
# PSE �Ɨ\�����덷�̊i�[
  pse <- rbind(pse, psem)
  perr <- rbind(perr, perrm)
}
# PSE�Ɨ\�����덷�̕��ς̕`��
PSE <- apply(pse, 2, mean)
PERR <- apply(perr, 2, mean)
par(omi = c(0, 0, 0, 0))
par(mai = c(1.0, 1.5, 0.5, 0.5))
par(ps = 18)
xlim <- c(1, mmaxp1)
ylim <- c(0.15, 0.4)
plot(1:mmaxp1, PERR, log = "y", type = "l", xlim = xlim, ylim = ylim, xlab = "Number of Terms", ylab = "Errors", main = "")
par(new = T)
plot(1:mmaxp1, PSE, log = "y", type = "p", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "", ylab = "", ann = F)
labels <- c("Prediction Error", "PSE")
pchs <- c(-1, 1, 16)
ltys <- c(1, -1, -1)
legend("topright",  legend  =  labels,  pch  =  pchs,  lty  =  ltys, bty = "n", y.intersp = 2.5
, pt.cex = 2)

