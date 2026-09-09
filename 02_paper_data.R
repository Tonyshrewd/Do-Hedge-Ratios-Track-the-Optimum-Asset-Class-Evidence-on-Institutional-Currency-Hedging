# Objects required by Project_updated.qmd
DATA <- "data/Replication/Data"

h  <- read.csv(file.path(DATA, "hedging_by_security.csv"))
h  <- h[!is.na(h$hedge_ratio), ]
h$q <- (h$yyyymm %/% 100) * 10 + ((h$yyyymm %% 100 - 1) %/% 3 + 1)

panel       <- read.csv("output/panel_R.csv")
panel_clean <- read.csv("output/panel_clean_R.csv")
m           <- read.csv("output/reit_merged.csv")

nw <- function(fit, L = 4) {
  X <- model.matrix(fit); u <- resid(fit); n <- nrow(X); k <- ncol(X)
  XtXi <- solve(crossprod(X))
  S <- crossprod(X * u)
  if (L > 0) for (l in 1:L) {
    w <- 1 - l / (L + 1)
    G <- crossprod(X[1:(n-l), , drop = FALSE] * u[1:(n-l)],
                   X[(l+1):n, , drop = FALSE] * u[(l+1):n])
    S <- S + w * (G + t(G))
  }
  V  <- XtXi %*% S %*% XtXi * (n / (n - k))
  se <- sqrt(diag(V))
  data.frame(coef = coef(fit), se = se, t = coef(fit)/se,
             p = 2 * pnorm(-abs(coef(fit)/se)))
}

