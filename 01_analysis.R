
####Setup####
setwd("~/Documents/Finance research")
dir.create("hedge_project/data", recursive = TRUE)
dir.create("hedge_project/output", showWarnings = FALSE)
setwd("~/Documents/Finance research/hedge_project")
getwd()
list.files()


pkg <- "~/Documents/Finance research/hedge_capacity_replication 3"
list.files(pkg)
file.exists(file.path(pkg, "00_raw/du_huber/Replication.zip"))


unzip(file.path(pkg, "00_raw/du_huber/Replication.zip"), exdir = "data")
list.files("data/Replication/Data")

list.files()

file.copy("~/Downloads/01_analysis.R", "01_analysis.R")
list.files()



# =============================================================
# Do Hedge Ratios Track the Optimum?
# Asset-class evidence on institutional currency hedging
# 01_analysis.R
# =============================================================

DATA <- "data/Replication/Data"
OUT  <- "output"

# --- 1. hedge ratios, split by asset class -------------------
h <- read.csv(file.path(DATA, "hedging_by_security.csv"))
h <- h[!is.na(h$hedge_ratio), ]

# --- 2. monthly asset and currency returns -------------------
r <- read.csv(file.path(DATA, "monthly_return.csv"))

# --- 3. hedging cost: the cross-currency basis ---------------
hd <- read.csv(file.path(DATA, "hedging_data.csv"))

#Data checks 
str(h)
table(h$entity, h$security)
head(r[, 1:5])

# Build hedge gap
# G = bond hedge ratio - equity hedge ratio, within investor.
# The cost of hedging is common to both sleeves and cancels in the
# difference, so G should reflect only the covariance differential.

h$q <- (h$yyyymm %/% 100) * 10 + ((h$yyyymm %% 100 - 1) %/% 3 + 1)

agg  <- aggregate(hedge_ratio ~ entity + q + security, data = h, FUN = mean)
bond <- agg[agg$security == "bond",   c("entity","q","hedge_ratio")]
eqty <- agg[agg$security == "equity", c("entity","q","hedge_ratio")]
names(bond)[3] <- "h_bond"
names(eqty)[3] <- "h_eqty"

gap <- merge(bond, eqty, by = c("entity","q"))
gap$G <- gap$h_bond - gap$h_eqty

aggregate(cbind(h_bond, h_eqty, G) ~ entity, data = gap, FUN = mean)


# Rolling betas
# beta_a = Cov(USD asset return, USD currency return) / Var(currency return),
# estimated over rolling windows. dbeta = beta_bond - beta_equity is the
# optimal hedge gap the theory predicts G should equal.

r$q <- (r$yyyymm %/% 100) * 10 + ((r$yyyymm %% 100 - 1) %/% 3 + 1)

roll_beta <- function(y, x, win) {
  n <- length(y); out <- rep(NA_real_, n)
  for (i in win:n) {
    ys <- y[(i - win + 1):i]; xs <- x[(i - win + 1):i]
    ok <- !is.na(ys) & !is.na(xs)
    if (sum(ok) >= win * 0.8) out[i] <- cov(ys[ok], xs[ok]) / var(xs[ok])
  }
  out
}

make_dbeta <- function(fx, win) {
  bb <- roll_beta(r$USD_excess_10Y_holding_1M,   fx, win)
  be <- roll_beta(r$USD_excess_stock_holding_1M, fx, win)
  data.frame(q = r$q, dbeta = bb - be)
}

last_in_q <- function(d) {
  d <- d[order(d$q), ]
  do.call(rbind, lapply(split(d, d$q), function(g) g[nrow(g), ]))
}

#Checking 
test <- make_dbeta(r$AUD_fx_holding_effective, 36)
summary(test$dbeta)


# Assemble panel
# Merge hedge gaps with rolling betas and the cross-currency basis.
# Australian pensions -> AUD, Dutch pensions -> EUR.

hd$q <- (hd$yyyymm %/% 100) * 10 + ((hd$yyyymm %% 100 - 1) %/% 3 + 1)
cost <- aggregate(cbind(basis_3M_qavg, diff_USD_less_CCY_ibor_3M, vol_1M_qavg) ~ currency + q,
                  data = hd, FUN = mean, na.action = na.pass)
names(cost) <- c("currency","q","basis","ibor","fxvol")

INV <- list(auspension = "AUD", dutchpension = "EUR")

panel <- do.call(rbind, lapply(names(INV), function(ent) {
  ccy <- INV[[ent]]
  g   <- gap[gap$entity == ent, ]
  fx  <- r[[paste0(ccy, "_fx_holding_effective")]]
  for (w in c(36, 60)) {
    b <- last_in_q(make_dbeta(fx, w))
    names(b)[2] <- paste0("dbeta", w)
    g <- merge(g, b, by = "q", all.x = TRUE)
  }
  merge(g, cost[cost$currency == ccy, c("q","basis","ibor","fxvol")], by = "q", all.x = TRUE)
}))


#Checking 

aggregate(cbind(G, dbeta36, basis) ~ entity, data = panel, FUN = mean, na.action = na.omit)



# Add mutual funds
# Mutual funds are not tied to one currency, so their FX return and
# hedging cost are USD-holdings-weighted averages across currency areas.

CCYS <- c("AUD","CAD","CHF","DKK","EUR","GBP","JPY","NOK","SEK","CLP","ILS")

wt <- aggregate(bondequity_usdtrillion ~ currency + q, data = hd, FUN = sum)
wt <- wt[wt$currency %in% CCYS, ]
tot <- aggregate(bondequity_usdtrillion ~ q, data = wt, FUN = sum)
names(tot)[2] <- "tot"
wt <- merge(wt, tot, by = "q")
wt$w <- wt$bondequity_usdtrillion / wt$tot

fx_mf <- rep(NA_real_, nrow(r))
for (i in seq_len(nrow(r))) {
  wq <- wt[wt$q == r$q[i], c("currency","w")]
  if (nrow(wq) == 0) next
  vals <- sapply(wq$currency, function(cc) r[[paste0(cc, "_fx_holding_effective")]][i])
  ok <- !is.na(vals)
  if (any(ok)) fx_mf[i] <- sum(vals[ok] * wq$w[ok]) / sum(wq$w[ok])
}

cost_mf <- do.call(rbind, lapply(split(wt, wt$q), function(g) {
  cc <- merge(g, cost, by = c("currency","q"))
  if (nrow(cc) == 0) return(NULL)
  data.frame(q = g$q[1],
             basis = weighted.mean(cc$basis, cc$w, na.rm = TRUE),
             ibor  = weighted.mean(cc$ibor,  cc$w, na.rm = TRUE),
             fxvol = weighted.mean(cc$fxvol, cc$w, na.rm = TRUE))
}))

g <- gap[gap$entity == "mf", ]
for (w in c(36, 60)) {
  b <- last_in_q(make_dbeta(fx_mf, w))
  names(b)[2] <- paste0("dbeta", w)
  g <- merge(g, b, by = "q", all.x = TRUE)
}
g <- merge(g, cost_mf, by = "q", all.x = TRUE)
panel <- rbind(panel, g[, names(panel)])


#Checking 

aggregate(cbind(G, dbeta36, basis) ~ entity, data = panel, FUN = mean, na.action = na.omit)


# Newey-West standard errors
# Hedge ratios and the basis are persistent, so OLS standard errors
# understate uncertainty. This corrects for serial correlation up to L lags.

nw <- function(fit, L = 4) {
  X <- model.matrix(fit); u <- resid(fit); n <- nrow(X); k <- ncol(X)
  XtXi <- solve(crossprod(X))
  S <- crossprod(X * u)
  if (L > 0) for (l in 1:L) {
    w  <- 1 - l / (L + 1)
    G  <- crossprod(X[1:(n-l), , drop = FALSE] * u[1:(n-l)],
                    X[(l+1):n, , drop = FALSE] * u[(l+1):n])
    S  <- S + w * (G + t(G))
  }
  V  <- XtXi %*% S %*% XtXi * (n / (n - k))
  se <- sqrt(diag(V))
  data.frame(coef = coef(fit), se = se, t = coef(fit)/se,
             p = 2 * pnorm(-abs(coef(fit)/se)))
}


# Main test
# G = a_i + b*dbeta + c*basis + d*ibor + e*fxvol
# Sleeve-level optimisation predicts b = 1 and c = 0.

for (w in c(36, 60)) {
  d <- panel[complete.cases(panel[, c("G", paste0("dbeta", w), "basis","ibor","fxvol")]), ]
  f <- as.formula(paste("G ~", paste0("dbeta", w),
                        "+ basis + ibor + fxvol + factor(entity)"))
  res <- nw(lm(f, data = d))
  b   <- res[paste0("dbeta", w), ]
  t1  <- (b$coef - 1) / b$se
  cat(sprintf("win=%2dm  b=%6.3f (%.3f) p=%.3f | H0:b=1 t=%7.2f p=%.4f | basis=%6.3f (%.3f) p=%.3f | N=%d\n",
              w, b$coef, b$se, b$p, t1, 2*pnorm(-abs(t1)),
              res["basis","coef"], res["basis","se"], res["basis","p"], nrow(d)))
}


# Reliability and the corrected null
# dbeta is estimated, not observed. Classical measurement error attenuates
# b toward zero by the reliability factor lambda, so sleeve-level
# optimisation predicts b = lambda, not b = 1.

d   <- panel[complete.cases(panel[, c("dbeta36","dbeta60")]), ]
cat(sprintf("pooled lambda = %.3f\n\n", cor(d$dbeta36, d$dbeta60)))

for (e in unique(panel$entity)) {
  s  <- panel[panel$entity == e & !is.na(panel$dbeta36) & !is.na(panel$dbeta60), ]
  vr <- sd(s$G, na.rm = TRUE) / sd(s$dbeta36, na.rm = TRUE)
  cat(sprintf("%-14s lambda = %6.3f   var ratio = %.3f   N=%d\n",
              e, cor(s$dbeta36, s$dbeta60), vr, nrow(s)))
}

###Summary####
# SO FAR
#
# Loaded three files from the Du & Huber (2026) archive: hedge ratios split by
# bond vs equity, monthly asset and FX returns, and the 3M cross-currency basis.
#
# Built G = bond hedge ratio - equity hedge ratio, by investor, quarterly.
#   auspension 0.345 | dutchpension 0.219 | mf 0.314
#   Bonds hedged more than equities in all three, as theory predicts.
#
# Built dbeta = beta_bond - beta_equity from rolling 36m and 60m regressions of
# USD bond and USD equity excess returns on the investor's USD currency return.
# This is the optimal gap. G should equal it exactly.
#   Mean dbeta36: 0.413 | 0.395 | 0.541 - all above the observed G.
#
# Main test: G on dbeta, basis, ibor, fxvol, with investor FE and Newey-West SEs.
#   b = 0.097 (0.027) at 36m, 0.059 (0.036) at 60m. Theory says 1.
#   H0: b = 1 rejected at p < 0.0001 in both. N = 108.
#   Basis insignificant (p = 0.40, 0.52).
#   => b ~ 0 and c ~ 0: hedge gaps look like fixed policy, not optimisation.
#
# Reliability of dbeta (correlation between the 36m and 60m estimates):
#   auspension   -0.347   var ratio 0.094   unusable
#   dutchpension  0.176   var ratio 0.522   barely usable
#   mf            0.720   var ratio 0.119   usable, and rejects
#   Measurement error attenuates b, so the null is lambda, not 1. The rejection
#   rests on mutual funds, where dbeta is well measured and G is about a sixth
#   as variable as the optimum.
#
# NEXT: (1) check whether hedge ratios are smoothed in the source data;
#       (2) re-estimate betas from daily returns to raise lambda.






code <- readLines("data/Replication/Code/Hedge_1_Process.R")
grep("hedge_ratio", code, value = TRUE)
grep("approx|na.locf|fill|interp|rollmean|rollapply|smooth|spline", code, value = TRUE)

# Smoothing check
# hedging_by_security.csv is only read and plotted in the Du & Huber code, never
# constructed, so it was hand-assembled from filings and regulator statistics.
# There is no code to audit, so test the series directly for interpolation:
# repeated values suggest carry-forward, repeated consecutive differences
# suggest linear interpolation between annual anchor points. If the observed
# flatness is a data artefact rather than behaviour, the result does not stand.

for (e in c("auspension","dutchpension","mf")) {
  for (s in c("bond","equity")) {
    x <- h$hedge_ratio[h$entity == e & h$security == s]
    d <- diff(x)
    cat(sprintf("%-13s %-6s N=%3d  unique=%3d  zero-diffs=%4.1f%%  repeated-diffs=%4.1f%%\n",
                e, s, length(x), length(unique(round(x, 6))),
                100*mean(d == 0), 100*mean(duplicated(round(d, 8)) & d != 0)))
  }
}


# Anchor frequency
# Locate the months where the change in the series actually changes.
# The spacing between those points is the true reporting frequency.
for (e in c("auspension","dutchpension","mf")) {
  x <- h$hedge_ratio[h$entity == e & h$security == "bond"]
  ym <- h$yyyymm[h$entity == e & h$security == "bond"]
  d <- round(diff(x), 8)
  brk <- which(c(TRUE, diff(d) != 0))
  cat(sprintf("%-13s breaks at months: %s ...\n", e,
              paste(head(ym[brk], 8), collapse = " ")))
  cat(sprintf("%-13s median spacing = %.1f months\n\n", e, median(diff(brk))))
}

# Drop interpolated series
# Dutch pensions are annual (December anchors) interpolated to monthly, so the
# quarterly panel contains 24 observations built from ~7 real ones. Australian
# pensions and mutual funds are natively quarterly and unaffected.

panel_clean <- panel[panel$entity != "dutchpension", ]

for (w in c(36, 60)) {
  d <- panel_clean[complete.cases(panel_clean[, c("G", paste0("dbeta", w), "basis","ibor","fxvol")]), ]
  f <- as.formula(paste("G ~", paste0("dbeta", w), "+ basis + ibor + fxvol + factor(entity)"))
  res <- nw(lm(f, data = d))
  b <- res[paste0("dbeta", w), ]
  t1 <- (b$coef - 1) / b$se
  cat(sprintf("win=%2dm  b=%6.3f (%.3f) | H0:b=1 t=%7.2f p=%.4f | basis=%6.3f (%.3f) p=%.3f | N=%d\n",
              w, b$coef, b$se, t1, 2*pnorm(-abs(t1)),
              res["basis","coef"], res["basis","se"], res["basis","p"], nrow(d)))
}


write.csv(panel, "output/panel_R.csv", row.names = FALSE)
write.csv(panel_clean, "output/panel_clean_R.csv", row.names = FALSE)




####Figures####
####Figures####

qdate <- function(q) (q %/% 10) + ((q %% 10) - 1)/4
cols  <- c(auspension = "#1f4e79", dutchpension = "#c0504d", mf = "#4f7942")

# Figure 1: observed gap vs optimal gap
png("output/fig1_gap_vs_optimum.png", width = 1500, height = 480, res = 130)
par(mfrow = c(1,3), mar = c(4,4,3,1))
for (e in c("auspension","dutchpension","mf")) {
  d <- panel[panel$entity == e, ]; d <- d[order(d$q), ]; x <- qdate(d$q)
  plot(x, d$dbeta36, type = "l", col = "grey55", lwd = 1.6,
       ylim = range(c(d$G, d$dbeta36), na.rm = TRUE),
       xlab = "", ylab = if (e == "auspension") "Bond - equity hedge ratio" else "",
       main = e)
  lines(x, d$G, col = cols[e], lwd = 3)
  abline(h = 0, col = "grey80")
  if (e == "auspension")
    legend("topleft", c("Observed gap G", "Optimal gap"),
           col = c(cols[e], "grey55"), lwd = c(3,1.6), bty = "n", cex = .85)
}
dev.off()


# Figure 2: scatter of G on dbeta, with the 45-degree line
png("output/fig2_scatter.png", width = 800, height = 800, res = 130)
par(mar = c(4.5,4.5,3,1))
d <- panel[complete.cases(panel[, c("G","dbeta36")]), ]
lim <- range(c(d$G, d$dbeta36), na.rm = TRUE)
plot(d$dbeta36, d$G, col = cols[d$entity], pch = 19, cex = .8,
     xlim = lim, ylim = lim, asp = 1,
     xlab = expression("Optimal gap  " * Delta * beta),
     ylab = "Observed gap  G",
     main = "Theory predicts every point on the 45-degree line")
abline(0, 1, lty = 2, lwd = 2, col = "black")
abline(lm(G ~ dbeta36, data = d), lwd = 2, col = "grey30")
legend("topleft",
       legend = c("45-degree (theory)", "fitted"),
       lty = c(2, 1), lwd = 2, col = c("black", "grey30"), bty = "n", cex = .8)
legend("bottomright", names(cols), col = cols, pch = 19, bty = "n", cex = .8)
dev.off()

# Figure 3: how variable is each series
png("output/fig3_variance.png", width = 850, height = 620, res = 130)
par(mar = c(4,4.5,3,1))
sds <- sapply(c("auspension","dutchpension","mf"), function(e) {
  s <- panel[panel$entity == e, ]
  c(G = sd(s$G, na.rm = TRUE), optimum = sd(s$dbeta36, na.rm = TRUE))
})
bp <- barplot(sds, beside = TRUE, col = c("#1f4e79","grey65"),
              ylim = c(0, max(sds, na.rm = TRUE) * 1.25),
              ylab = "Standard deviation",
              main = "Observed gap is far less variable than the optimum")
text(colMeans(bp), sds[2,] * 1.10,
     sprintf("ratio %.2f", sds[1,]/sds[2,]), cex = .85)
legend("topleft", c("Observed gap G","Optimal gap"),
       fill = c("#1f4e79","grey65"), bty = "n")
dev.off()


# Figure 4: interpolation diagnostic
png("output/fig4_interpolation.png", width = 1500, height = 480, res = 130)
par(mfrow = c(1,3), mar = c(4,4,3,1))
for (e in c("auspension","dutchpension","mf")) {
  s  <- h[h$entity == e & h$security == "bond", ]
  s  <- s[order(s$yyyymm), ]
  tt <- (s$yyyymm %/% 100) + ((s$yyyymm %% 100) - 1)/12
  dd <- round(diff(s$hedge_ratio), 8)
  brk <- which(c(TRUE, diff(dd) != 0))
  plot(tt, s$hedge_ratio, type = "l", col = "grey60", lwd = 1.5,
       xlab = "", ylab = if (e == "auspension") "Bond hedge ratio" else "",
       main = sprintf("%s (anchors every %.0f months)", e, median(diff(brk))))
  points(tt[brk], s$hedge_ratio[brk], pch = 19, col = cols[e], cex = 1.1)
}
dev.off()



# Figure 5: reliability of the covariance estimate
png("output/fig5_reliability.png", width = 1500, height = 480, res = 130)
par(mfrow = c(1,3), mar = c(4,4,3,1))
for (e in c("auspension","dutchpension","mf")) {
  s <- panel[panel$entity == e & complete.cases(panel[, c("dbeta36","dbeta60")]), ]
  plot(s$dbeta36, s$dbeta60, pch = 19, col = cols[e], cex = .9,
       xlab = "36-month estimate",
       ylab = if (e == "auspension") "60-month estimate" else "",
       main = sprintf("%s  (lambda = %.2f)", e, cor(s$dbeta36, s$dbeta60)))
  abline(0, 1, lty = 2, col = "grey50")
}
dev.off()


####REAL ESTATE####
install.packages(c("RPostgres", "DBI", "dplyr"))
library(DBI); library(RPostgres)

wrds <- dbConnect(Postgres(),
                  host    = "wrds-pgdata.wharton.upenn.edu",
                  port    = 9737,
                  dbname  = "wrds",
                  sslmode = "require",
                  user    = rstudioapi::askForPassword("anthony123"),
                  password = rstudioapi::askForPassword("Sutescop123!"))
dbGetQuery(wrds, "
  SELECT table_name FROM information_schema.tables
  WHERE table_schema = 'crsp'
    AND (table_name LIKE 'msf%' OR table_name LIKE 'msenames%' OR table_name LIKE 'stksecurityinfo%')
  ORDER BY table_name")

dbGetQuery(wrds, "
  SELECT column_name FROM information_schema.columns
  WHERE table_schema='crsp' AND table_name='msf'
  ORDER BY column_name")

dbGetQuery(wrds, "
  SELECT column_name FROM information_schema.columns
  WHERE table_schema='crsp' AND table_name='msenames'
  ORDER BY column_name")



q <- "
  SELECT a.permno, a.date, a.ret, a.prc, a.shrout, b.siccd, b.shrcd
  FROM   crsp.msf AS a
  LEFT JOIN crsp.msenames AS b
         ON a.permno = b.permno
        AND b.namedt <= a.date
        AND a.date   <= b.nameendt
  WHERE  a.date BETWEEN '2002-01-01' AND '2021-12-31'
    AND  b.siccd = 6798
    AND  b.exchcd IN (1,2,3)
    AND  b.shrcd IN (11,18)
    AND  a.ret IS NOT NULL
"
raw <- dbGetQuery(wrds, q)

nrow(raw)
length(unique(raw$permno))
table(raw$shrcd)


raw$mcap <- abs(raw$prc) * raw$shrout
raw <- raw[order(raw$permno, raw$date), ]

# lag market cap within firm
raw$mcap_lag <- ave(raw$mcap, raw$permno,
                    FUN = function(z) c(NA, z[-length(z)]))

d <- raw[!is.na(raw$mcap_lag) & raw$mcap_lag > 0 & !is.na(raw$ret), ]

reit <- do.call(rbind, lapply(split(d, d$date), function(g) {
  data.frame(date     = g$date[1],
             reit_ret = weighted.mean(g$ret, g$mcap_lag),
             n_firms  = nrow(g))
}))
reit$yyyymm <- as.integer(format(reit$date, "%Y%m"))
reit <- reit[order(reit$date), ]

nrow(reit)
range(reit$n_firms)
summary(reit$reit_ret)
write.csv(reit, "output/reit_monthly.csv", row.names = FALSE)


####Matching####
code <- unlist(lapply(list.files("data/Replication/Code", pattern="\\.R$", full.names=TRUE),
                      function(f) paste0(basename(f), ": ", readLines(f))))
grep("excess_stock|excess_10Y|holding_1M", code, value = TRUE)

p <- readLines("data/Replication/Code/Hedge_1_Process.R")
i <- grep("excess_stock_holding_1M", p)[1]
cat(p[i:(i+3)], sep = "\n")


j <- grep("excess_10Y_holding_1M", p)[1]
cat(p[j:(j+3)], sep = "\n")


# 1-month IBOR from their data, lagged, in percent
r$ibor1m_lag <- c(NA, head(r$USD_ibor_1M, -1))

m <- merge(reit[, c("yyyymm","reit_ret")],
           r[, c("yyyymm","ibor1m_lag")], by = "yyyymm")

# their convention: 12 x log capital gain, less log(1 + lagged short rate)
m$reit_excess <- 12 * log(1 + m$reit_ret) - log(1 + m$ibor1m_lag/100)

summary(m$reit_excess)
sd(m$reit_excess, na.rm = TRUE)



grep("ibor", names(r), value = TRUE, ignore.case = TRUE)
r$ibor1m_lag <- NULL

grep("ibor", names(hd), value = TRUE, ignore.case = TRUE)

m <- merge(reit[, c("yyyymm","reit_ret")], r[, c("yyyymm","q")], by = "yyyymm")
m$reit_excess <- 12 * log(1 + m$reit_ret)

summary(m$reit_excess)
sd(m$reit_excess, na.rm = TRUE)
sd(r$USD_excess_stock_holding_1M, na.rm = TRUE)
sd(r$USD_excess_10Y_holding_1M, na.rm = TRUE)



m <- merge(m, r[, c("yyyymm","AUD_fx_holding_effective","EUR_fx_holding_effective",
                    "USD_excess_10Y_holding_1M","USD_excess_stock_holding_1M")],
           by = "yyyymm")

for (ccy in c("AUD","EUR")) {
  fx <- m[[paste0(ccy, "_fx_holding_effective")]]
  bb <- cov(m$USD_excess_10Y_holding_1M,   fx, use="complete.obs") / var(fx, na.rm=TRUE)
  be <- cov(m$USD_excess_stock_holding_1M, fx, use="complete.obs") / var(fx, na.rm=TRUE)
  br <- cov(m$reit_excess,                 fx, use="complete.obs") / var(fx, na.rm=TRUE)
  cat(sprintf("%s   beta_bond=%6.3f   beta_equity=%6.3f   beta_REIT=%6.3f   REIT-equity=%6.3f\n",
              ccy, bb, be, br, br - be))
}


m$year <- m$yyyymm %/% 100

for (ccy in c("AUD","EUR")) {
  for (lab in c("full","ex 2020")) {
    s  <- if (lab == "full") m else m[m$year != 2020, ]
    fx <- s[[paste0(ccy, "_fx_holding_effective")]]
    be <- cov(s$USD_excess_stock_holding_1M, fx, use="complete.obs") / var(fx, na.rm=TRUE)
    br <- cov(s$reit_excess,                 fx, use="complete.obs") / var(fx, na.rm=TRUE)
    cat(sprintf("%s  %-8s  beta_equity=%6.3f  beta_REIT=%6.3f  REIT-equity=%6.3f  N=%d\n",
                ccy, lab, be, br, br - be, sum(!is.na(fx))))
  }
}


write.csv(m, "output/reit_merged.csv", row.names = FALSE)

#### plot####
png("output/fig6_reit_ordering.png", width = 900, height = 500, res = 130)
par(mfrow = c(1,2), mar = c(4,6,3,1))
for (ccy in c("AUD","EUR")) {
  fx <- m[[paste0(ccy, "_fx_holding_effective")]]
  b  <- c(Bond   = cov(m$USD_excess_10Y_holding_1M,   fx, use="complete.obs")/var(fx, na.rm=TRUE),
          Equity = cov(m$USD_excess_stock_holding_1M, fx, use="complete.obs")/var(fx, na.rm=TRUE),
          REIT   = cov(m$reit_excess,                 fx, use="complete.obs")/var(fx, na.rm=TRUE))
  plot(b, seq_along(b), pch = 19, cex = 1.6, yaxt = "n",
       xlim = c(-1, 0.2), ylim = c(0.5, 3.5),
       col = c("#1f4e79","#4f7942","#c0504d"),
       xlab = expression(beta), ylab = "", main = ccy)
  axis(2, at = seq_along(b), labels = names(b), las = 1)
  abline(v = 0, col = "grey70", lty = 2)
  segments(0, seq_along(b), b, seq_along(b), col = "grey60")
}
dev.off()

