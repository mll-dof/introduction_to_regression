# コード6.4
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
# コード6.5
bp <- function(X, y, W, eta, T, dT) {
    N <- nrow(X)
    L <- length(W)
    ehist <- NULL
# 学習のループ
    for(t in 1:T){
# 誤差の計算
        out <- forward(X, W)
        err <- (y - out[[L+1]])
        ehist <- c(ehist, mean(err^2))
# dT 回ごとに誤差を表示
        if(t %% dT == 0) cat("[",t,"] : ", mean(err^2), "\n") 
# 出力層の手前の結合重みの更新量の計算
        delta <- -err * rep(1, N) 
        dw <- t(t(delta) %*% out[[L]]) 
        for(l in L:2){
            W[[l]] <- W[[l]] - eta * dw # 結合重みの更新
            dtmp <- delta
            delta <- NULL
# 結合重みの更新量の計算
            for(j in 2:nrow(W[[l]])){
                delta <- cbind(delta, dtmp %*% W[[l]][j, ] * dsigmoid(out[[l]][, j]))
            }
            dw <- t(t(delta) %*% out[[l - 1]]) 
        }
        W[[1]] <- W[[1]] - eta * dw # 入力層と結合している結合重みの更新
    }
# 誤差の計算
    out <- forward(X, W)
    err <- (y - out[[L+1]])
    ehist <- c(ehist, mean(err^2))
    cat("[",t,"] : ", mean(err^2), "\n")

    return(list("weights" = W, "err_history" = ehist))
}
# コード6.6：readlineコマンドは入力状態なので実行時にはリターンキーを押してください
# サンプルの生成
set.seed(2017)
n <- 100 
x <- runif(n, -1, 1)
X <- matrix(x, n, 1)    
y <- matrix(sin(pi * x), n, 1)
# 学習の準備
netsize <- c(1, 4, 1); irng <- 1; 
W <- initialize(netsize, irng)
eta <- 0.0012; T <- 100000; dT <- T / 100
# 学習
ret <- bp(X, y, W, eta, T, dT)
# 学習後の結合重みと残差二乗和の履歴を抽出
W <- ret$weights
ehist <- ret$err_history
L <- length(W)
# 残差二乗和の履歴の描画
readline(">>> show Error history <<<")
plot(0:T, ehist, type = "l", log = "y", xlab = "iteration", ylab = "Error")
# 教師データと学習後の出力の描画
readline(">>> show output <<<")
out <- forward(X, W)
plot(x, y, xlim = c(-1, 1), ylim = c(-1.5, 1.5), lwd = 0.5, cex = 1.5, xlab = "x", ylab = "y, Output")
par(new = T)
xx <- sort(x, index.return = T)
plot(xx$x, out[[L+1]][xx$ix], type = "l", xlim = c(-1, 1), ylim = c(-1.5, 1.5), xlab = "", ylab = "")
