# �R�[�h4.2
lr_train <- function(X, y, beta, itmax) {
  for(it in 1:itmax) { # �w�K�̃��[�v
    eta <- as.vector(1 / (1 + exp(-X %*% beta)))
    d <- (y - eta)
    score <- t(X) %*% d
    W <- diag(eta * (eta - 1))
    hessian <- t(X) %*% W %*%X
    dbeta <- solve(hessian) %*% score # �X�V�ʂ̌v�Z
    beta <- beta - dbeta # �X�V���̌v�Z
  }
  return(beta)
}
# �R�[�h4.3
# �T���v���̐���
set.seed(2017)
n <- 100 
x <- (0:(n - 1)) / (n - 1) * 4 - 2 
tbeta <- c(0, 2) 
p <- 1/(1 + exp(-(tbeta[1] + tbeta[2] * x)))
y <- rbinom(n, 1, prob = p)
# lr_train �ɂ��w�K
beta <- runif(2, -1, 1) # �����l
itmax <- 10 # �J��Ԃ���
X <- cbind(rep(1, n), x) # �f�U�C���s��
beta <- lr_train(X, y, beta, itmax) 
# �T���v���C�^�̏����t���m������ѐ��肵�������t���m���̕`��
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
