# �R�[�h3.9
cal_cverr_l2reg <- function(X, y, folds, rpset) {
  n <- length(y)
# �T���v���̏��Ԃ��V���b�t������
  shuf <- sample(1:n, n, replace = F)
  X <- X[shuf, ]
  y <- y[shuf]
# �u���b�N���Ƃ̃T���v�����ƃu���b�N�̎n�_�ƏI�_���v�Z����
  n1b <- floor(n / folds)
  nB <- rep(n1b, folds)
  if (n %% folds != 0)  nB[1:(n - n1b * folds)] <- nB[1:(n - n1b * folds)] + 1
  st <- c(1, cumsum(nB)[-folds] + 1)
  ed <- cumsum(nB)
# ���،덷���v�Z����
  cve <- NULL
  for (rp in rpset) {
    cvetmp <- NULL
    for (j in 1:folds) {
# ����p�Z�b�g�ƌ��ؗp�Z�b�g���\������
      idx <- (st[j]:ed[j])
      Xe <- X[-idx, ]
      ye <- y[-idx]
      Xv <- X[idx, ]
      yv <- y[idx]
# ����p�Z�b�g�ɑ΂���p�����[�^����
      coef <- solve(t(Xe) %*% Xe + rp * diag(ncol(Xe))) %*% t(Xe) %*% ye
# ���ؗp�Z�b�g�ɑ΂���o�͂̌v�Z
      fv <- Xv%*%coef
# ���ؗp�Z�b�g�ɑ΂���덷�̌v�Z
      cvetmp <- c(cvetmp, (yv - fv)^2)
    }
    cve <- c(cve, mean(cvetmp))
  }
# ���،덷���ŏ������鐳�����p�����[�^�l�����߂�
  minpos <- which.min(cve) 
  rlist <- list(cve, minpos)
  names(rlist) <- c("cve", "minpos")
  return(rlist)
}
# �R�[�h3.10
# �����ϐ��̃T���v���̐���
set.seed(2017)
n <- 50; mmax <- 25; sigma <- 0.4; S <- 1000
x <- 2*(0:(n - 1)) / (n - 1) - 1
X <- rep(1, n)
for (m in 1:mmax) X <- cbind(X, x^m)
tf <- sin(pi * x)
# �������p�����[�^�̌��W���̐���
rpset <- seq(-3, 1, 1)
rpset <- c(1*10^rpset, 5*10^rpset)
rpset <- sort(rpset)
freq <- numeric(length(rpset))
folds <- 10
cve <- NULL
# �������p�����[�^�l���N���X�o���f�[�V�����ɂ���ċ��߂鐔�l����
for (s in 1:S) {
  y <- tf + rnorm(n, 0, sigma) # �ړI�ϐ��̃T���v���̐���
  rlist <- cal_cverr_l2reg(X, y, folds, rpset)
  pos <- rlist$minpos
  cve <- rbind(cve, rlist$cve)
}
# �]���덷�̕��ς̕`��
CVE <- apply(cve, 2, mean);
par(omi = c(0, 0, 0, 0));
par(mai = c(1.0, 1.5, 0.5, 0.5));
par(ps = 18);
xlim <- c(rpset[1], rpset[length(rpset)]);
ylim <- c(0.1, 1.0);
plot(rpset, CVE, type = "p", log = "xy", pch = 1, lwd = 0.5, cex = 2, xlim = xlim, ylim = ylim, xlab = "reguralization param.", ylab = "validation error");

