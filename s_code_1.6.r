# �T���v���̐���
set.seed(2017)
n <- 100
x1 <- runif(n, -2, 2) 
x2 <- runif(n, -1, 1)
X <- cbind(rep(1, n), x1, x2)
tb <- c(2, 0, -1)
f <- X %*% tb
y <- f + rnorm(n, 0, 0.4)
# �ŏ���搄��ʂ̌v�Z
ecoef <- solve(t(X) %*% X) %*% t(X) %*% y
# �ΕW����A�W���̌v�Z
ecoef1 <- rbind(0, ecoef[2] * sd(x1) / sd(y), ecoef[3] * sd(x2) / sd(y)) 
cat("ecoef1 : ", ecoef1, "\n")

# �W�������ꂽ�T���v���̉��ł̍ŏ���搄��ʂ̌v�Z
u1 <- (x1 - mean(x1)) / sd(x1)
u2 <- (x2 - mean(x2)) / sd(x2)
z <- (y - mean(y)) / sd(y)
U <- cbind(rep(1, n), u1, u2)
ecoef2 <- solve(t(U) %*% U) %*% t(U) %*% z
cat("ecoef2 : ", ecoef2, "\n")
