# �R�[�h6.4
sigmoid <- function(x) {
  return(1.0  / (1.0 + exp(-x)))
}
dsigmoid <- function(x) {
  return(x * (1.0 - x))
}
initialize <- function(netsize, irng) {
    L <- length(netsize) - 1
    W <- NULL
    for(l in 1:L){
        k <- netsize[l]
        kn <- netsize[l + 1]
        Wtmp <- matrix(runif((k + 1) * kn, -irng, irng), k + 1, kn)
        W <- c(W, list(Wtmp))
    }
    return(W)
}
forward <- function(X, W) {
    N <- nrow(X)
    L <- length(W)
    out <- list(cbind(rep(1, N), X))
    for(l in 1:(L-1)){
        hout <- sigmoid(out[[l]] %*% W[[l]])
        out <- c(out, list(cbind(rep(1, N), hout)))
    }
    hout <- out[[L]] %*% W[[L]]
    out <- c(out, list(hout))
    return(out)
}
# �R�[�h6.5
bp <- function(X, y, W, eta, T, dT) {
    N <- nrow(X)
    L <- length(W)
    ehist <- NULL
# �w�K�̃��[�v
    for(t in 1:T){
# �덷�̌v�Z
        out <- forward(X, W)
        err <- (y - out[[L+1]])
        ehist <- c(ehist, mean(err^2))
# dT �񂲂ƂɌ덷��\��
        if(t %% dT == 0) cat("[",t,"] : ", mean(err^2), "\n") 
# �o�͑w�̎�O�̌����d�݂̍X�V�ʂ̌v�Z
        delta <- -err * rep(1, N) 
        dw <- t(t(delta) %*% out[[L]]) 
        for(l in L:2){
            W[[l]] <- W[[l]] - eta * dw # �����d�݂̍X�V
            dtmp <- delta
            delta <- NULL
# �����d�݂̍X�V�ʂ̌v�Z
            for(j in 2:nrow(W[[l]])){
                delta <- cbind(delta, dtmp %*% W[[l]][j, ] * dsigmoid(out[[l]][, j]))
            }
            dw <- t(t(delta) %*% out[[l - 1]]) 
        }
        W[[1]] <- W[[1]] - eta * dw # ���͑w�ƌ������Ă��錋���d�݂̍X�V
    }
# �덷�̌v�Z
    out <- forward(X, W)
    err <- (y - out[[L+1]])
    ehist <- c(ehist, mean(err^2))
    cat("[",t,"] : ", mean(err^2), "\n")

    return(list("weights" = W, "err_history" = ehist))
}
# �R�[�h6.6�Freadline�R�}���h�͓��͏�ԂȂ̂Ŏ��s���ɂ̓��^�[���L�[�������Ă�������
# �T���v���̐���
set.seed(2017)
n <- 100 
x <- runif(n, -1, 1)
X <- matrix(x, n, 1)    
y <- matrix(sin(pi * x), n, 1)
# �w�K�̏���
netsize <- c(1, 4, 1); irng <- 1; 
W <- initialize(netsize, irng)
eta <- 0.0012; T <- 100000; dT <- T / 100
# �w�K
ret <- bp(X, y, W, eta, T, dT)
# �w�K��̌����d�݂Ǝc�����a�̗����𒊏o
W <- ret$weights
ehist <- ret$err_history
L <- length(W)
# �c�����a�̗����̕`��
readline(">>> show Error history <<<")
plot(0:T, ehist, type = "l", log = "y", xlab = "iteration", ylab = "Error")
# ���t�f�[�^�Ɗw�K��̏o�͂̕`��
readline(">>> show output <<<")
out <- forward(X, W)
plot(x, y, xlim = c(-1, 1), ylim = c(-1.5, 1.5), lwd = 0.5, cex = 1.5, xlab = "x", ylab = "y, Output")
par(new = T)
xx <- sort(x, index.return = T)
plot(xx$x, out[[L+1]][xx$ix], type = "l", xlim = c(-1, 1), ylim = c(-1.5, 1.5), xlab = "", ylab = "")
