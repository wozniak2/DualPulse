## Cluster-robust inference for one OLS coefficient with few, unbalanced clusters (sourced by 56 and 63).
## crinf(): CR2 standard error with Satterthwaite (Bell-McCaffrey) degrees of freedom (Pustejovsky & Tipton
## 2018; Imbens & Kolesar 2016), CV1 and CV3 (jackknife) standard errors, partial leverage of the clusters
## and the effective number of clusters G* (Carter, Schnepel & Steigerwald 2017, rho = 0).
## wcr(): restricted wild cluster bootstrap p-value (WCR-C, CV1 t, Webb six-point weights; MacKinnon,
## Nielsen & Webb 2023).
msq <- function(ev, pw) ev$vectors %*% diag(ifelse(ev$values > 1e-8, pmax(ev$values, 1e-8)^pw, 0), length(ev$values)) %*% t(ev$vectors)
## one coefficient k of an OLS fit, clustered by cl: CR2 with Satterthwaite df, CV1, CV3, partial leverage, G*
crinf <- function(fit, cl, k) {
  X <- model.matrix(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]; e <- resid(fit); n <- nrow(X); K <- ncol(X)
  B <- solve(crossprod(X)); w <- as.vector(X %*% B[match(k, colnames(X)), ])   # x_i'(X'X)^-1 l: the coefficient is sum(w * y)
  Mres <- diag(n) - X %*% B %*% t(X)
  gs <- split(seq_len(n), cl); G <- length(gs); Gm <- matrix(0, n, G); s1 <- s2 <- s3 <- 0; gam <- numeric(G)
  for (s in seq_len(G)) { i <- gs[[s]]; ev <- eigen(Mres[i, i, drop = FALSE], symmetric = TRUE)
    A2 <- msq(ev, -1/2); A3 <- msq(ev, -1); wi <- w[i]
    s1 <- s1 + sum(wi * e[i])^2; s2 <- s2 + sum(wi * (A2 %*% e[i]))^2; s3 <- s3 + sum(wi * (A3 %*% e[i]))^2
    Gm[, s] <- Mres[, i, drop = FALSE] %*% (A2 %*% wi); gam[s] <- sum(wi^2) }
  GG <- crossprod(Gm); df <- sum(diag(GG))^2 / sum(GG^2)
  lev <- gam / sum(gam)
  list(b = sum(w * fit$model[[1]]), se2 = sqrt(s2), df = df, se1 = sqrt(G / (G - 1) * (n - 1) / (n - K) * s1), se3 = sqrt((G - 1) / G * s3),
       G = G, Gstar = G / (1 + mean((gam - mean(gam))^2) / mean(gam)^2), lev_max = max(lev), lev_os = names(gs)[which.max(lev)])
}
## restricted wild cluster bootstrap (WCR-C): null imposed, CV1 t statistic, Webb six-point weights
wcr <- function(fit, cl, k, B = 9999) {
  X <- model.matrix(fit); keep <- !is.na(coef(fit)); X <- X[, keep, drop = FALSE]; y <- fit$model[[1]]; n <- nrow(X); K <- ncol(X)
  j <- match(k, colnames(X)); X0 <- X[, -j, drop = FALSE]
  f0 <- lm.fit(X0, y); yt <- f0$fitted.values; ut <- f0$residuals
  Bm <- solve(crossprod(X)); w <- as.vector(X %*% Bm[j, ]); gi <- as.integer(factor(cl)); G <- max(gi); cc <- G / (G - 1) * (n - 1) / (n - K)
  tstat <- function(Y) { b <- colSums(w * Y); E <- Y - X %*% (Bm %*% crossprod(X, Y)); S <- rowsum(w * E, gi); b / sqrt(cc * colSums(S^2)) }
  t0 <- tstat(matrix(y, ncol = 1))
  V <- matrix(sample(c(-sqrt(3/2), -1, -sqrt(1/2), sqrt(1/2), 1, sqrt(3/2)), G * B, replace = TRUE), G, B)
  ts <- tstat(yt + ut * V[gi, , drop = FALSE])
  mean(abs(ts) >= abs(t0))
}
