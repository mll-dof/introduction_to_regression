set.seed(2017);
n <- 50;
tcoef <- c(-1, 2);
x <- runif(n, -1, 1);
tf <- tcoef[1] + tcoef[2] * x;
y <- tf + rnorm(n, 0, 0.5);

X <- cbind(rep(1, n), x);
ecoef <- solve(t(X) %*% X) %*% t(X) %*% y;

nn <- 100;
u <- seq(-1, 1, 2 / (nn - 1));
ef <- ecoef[1] + ecoef[2] * u;

par(omi = c(0, 0, 0, 0));
par(mai = c(1.0, 1.5, 0.5, 0.5));
par(ps = 24);
par(family = "Palatino");
plot(x, y, xlim = c(-1, 1), ylim = c(-4, 4), type = "p", pch = 1, lwd = 0.5, cex = 2, xlab = "x", ylab = "y", col = gray(0));
par(new = T);
plot(u, ef, xlim = c(-1, 1), ylim = c(-4, 4), type = "l", lwd = 0.5, xlab = "", ylab = "", col = gray(0), ann = F);
