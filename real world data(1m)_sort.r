
# --- Step 0: Clean Slate ---
rm(list = ls())
graphics.off()
cat("\014") # Clear Console

suppressPackageStartupMessages({
  library(tidyquant); library(JGL); library(clime); library(Matrix)
  library(R.utils); library(dplyr); library(tidyr); library(zoo); library(lubridate)
})

Sys.setenv(TZ = "UTC")
RNGkind("Mersenne-Twister","Inversion","Rejection")
set.seed(123) # Seed for reproducibility

# =========================================================
# 0. Helper: Safe scaling
# =========================================================
safe_scale <- function(m) {
  s <- scale(m, center = TRUE, scale = TRUE)
  s[is.nan(s)] <- 0
  as.matrix(s)
}

# =========================================================
# 1. Configuration & Data
# =========================================================

TICKERS <- c(
  "AAPL","MSFT","AMZN","TSLA","GOOGL","BRK-B","NVDA","META","V",
  "JPM","JNJ","XOM","WMT","BAC","PFE","HD","DIS","MA","ABT",
  "T","CMCSA","VZ","UNH","PG","INTC","NFLX","ADBE","CSCO","CRM",
  "PEP","COST","AVGO","CVX","KO","ABBV","MRK","TMO","ACN","NKE",
  "QCOM","WFC","DHR","LIN","UPS","AMD","MS","BA","ORCL","HON",
  "LOW","CAT","IBM","MMM","RTX","GS","AXP","BK","CL","C"
)

START_DATE <- "2020-01-01"
END_DATE   <- "2025-12-31"

# --- Parameters ---
W <- 12  # Rolling Window Size (in MONTH blocks)
S <- 1   # STEP size (in MONTH blocks). Set S=1 for monthly update, S=2 for every-2-months update, etc.

# --- StARS Grid (Robust Range) ---
L1_GRID_INIT <- seq(0.10, 0.50, length.out = 10)
L2_GRID_INIT <- seq(0.01, 0.20, length.out = 10)

STARS_SUBSAMPLES  <- 10
STARS_THRESH      <- 0.05

# --- Strategy Params ---
LMAX_ROLL   <- 1.0
LMAX_STATIC <- 1.0
LMAX_CLIME  <- 1.0
EMA_ALPHA   <- 0.2

# --- Download & Process ---
cat(">> [1/7] Downloading data...\n")
prices <- tq_get(TICKERS, get = "stock.prices", from = START_DATE, to = END_DATE)

returns_wide <- prices %>%
  group_by(symbol) %>% arrange(date) %>%
  mutate(ret = adjusted/lag(adjusted) - 1) %>%
  ungroup() %>% select(date, symbol, ret) %>%
  tidyr::drop_na() %>%
  pivot_wider(names_from = symbol, values_from = ret) %>%
  arrange(date) %>% tidyr::drop_na()

dates <- returns_wide$date
Rmat  <- as.matrix(select(returns_wide, -date))
p     <- ncol(Rmat)

# Create Calendar Blocks (MONTH)
build_calendar_classes <- function(dates){
  ym <- zoo::as.yearmon(dates)
  split(seq_along(dates), ym)
}
class_indices <- build_calendar_classes(dates)
data_list <- lapply(class_indices, function(ix) Rmat[ix, , drop = FALSE])
K <- length(data_list)
cat(sprintf(">> Data Ready: %d Assets, %d Months.\n", p, K))

# =========================================================
# 2. Robust StARS Tuning (Split-Sample Fix)
# =========================================================

compute_stars_instability <- function(data_list, l1_grid, l2_grid, n_subs = 10, beta = 0.8){
  n_l1 <- length(l1_grid); n_l2 <- length(l2_grid)
  instability <- matrix(NA, n_l1, n_l2)
  
  X_full <- do.call(rbind, data_list)
  n_total <- nrow(X_full)
  b <- floor(n_total * beta)
  if(b %% 2 != 0) b <- b - 1
  
  cat(">> [2/7] Running StARS (Split-Sample Mode)...\n")
  
  for(i in 1:n_l1){
    cat(sprintf("   Checking L1=%.2f ... ", l1_grid[i]))
    
    for(j in 1:n_l2){
      edges_bucket <- list()
      
      for(k in 1:n_subs){
        idx_sub <- sort(sample(1:n_total, b, replace = FALSE))
        mat_sub <- X_full[idx_sub, , drop = FALSE]
        
        half <- b / 2
        mat1 <- mat_sub[1:half, , drop = FALSE]
        mat2 <- mat_sub[(half + 1):b, , drop = FALSE]
        
        X_split <- list(safe_scale(mat1), safe_scale(mat2))
        
        try_fit <- tryCatch(
          JGL(Y = X_split, penalty = "fused",
              lambda1 = l1_grid[i], lambda2 = l2_grid[j],
              tol = 1e-3, maxiter = 100, return.whole.theta = TRUE),
          error = function(e) NULL
        )
        
        if(!is.null(try_fit)){
          theta <- try_fit$theta[[1]]
          adj <- (abs(theta) > 1e-6) * 1
          diag(adj) <- 0
          edges_bucket[[length(edges_bucket) + 1]] <- adj
        }
      }
      
      if(length(edges_bucket) > 1){
        avg_theta <- Reduce("+", edges_bucket) / length(edges_bucket)
        inst <- 2 * avg_theta * (1 - avg_theta)
        instability[i,j] <- mean(inst[upper.tri(inst)])
      } else {
        instability[i,j] <- 1.0
      }
    }
    cat(sprintf("Min Inst: %.4f\n", min(instability[i,], na.rm = TRUE)))
  }
  instability
}

inst_mat <- compute_stars_instability(data_list, L1_GRID_INIT, L2_GRID_INIT, n_subs = STARS_SUBSAMPLES)

select_stars_params <- function(instability, l1_grid, l2_grid, threshold = 0.05){
  valid_indices <- which(instability <= threshold, arr.ind = TRUE)
  
  if(nrow(valid_indices) == 0){
    min_val <- min(instability, na.rm = TRUE)
    cat(sprintf("\n>> WARNING: Min instability is %.4f (Threshold %.2f). Picking min.\n", min_val, threshold))
    best_idx <- which.min(instability)
    r <- row(instability)[best_idx]; c <- col(instability)[best_idx]
    return(list(l1 = l1_grid[r], l2 = l2_grid[c]))
  }
  
  sorted_valid <- valid_indices[order(l1_grid[valid_indices[,1]]), , drop = FALSE]
  best_r <- sorted_valid[1, 1]
  best_c <- sorted_valid[1, 2]
  list(l1 = l1_grid[best_r], l2 = l2_grid[best_c])
}

best_params <- select_stars_params(inst_mat, L1_GRID_INIT, L2_GRID_INIT, threshold = STARS_THRESH)
cat(sprintf("\n>> [3/7] FINAL PARAMETERS: Lambda1=%.4f, Lambda2=%.4f\n", best_params$l1, best_params$l2))

# =========================================================
# 3. Rolling JGL Estimation (UPDATED: S truly works)
# =========================================================

rolling_jgl_estimation <- function(data_list, l1, l2, W, S){
  K <- length(data_list)
  p_local <- ncol(data_list[[1]])
  theta_seq <- vector("list", K)
  
  pb <- txtProgressBar(min = 0, max = K, style = 3)
  
  last_theta <- diag(p_local)
  
  for (t in 1:K) {
    
    # Update only every S months; otherwise carry forward
    if (((t - 1) %% S) != 0) {
      theta_seq[[t]] <- last_theta
      setTxtProgressBar(pb, t)
      next
    }
    
    idx_start <- max(1, t - W + 1)
    idx_range <- idx_start:t
    
    Y_sub <- lapply(data_list[idx_range], safe_scale)
    
    if(length(idx_range) < 3) {
      theta_seq[[t]] <- last_theta
      setTxtProgressBar(pb, t)
      next
    }
    
    try_fit <- tryCatch(
      JGL(Y = Y_sub, penalty = "fused",
          lambda1 = l1, lambda2 = l2,
          tol = 1e-3, maxiter = 200, return.whole.theta = TRUE),
      error = function(e) NULL
    )
    
    if(!is.null(try_fit)){
      mat <- try_fit$theta[[length(try_fit$theta)]]
      last_theta <- (mat + t(mat)) / 2
    }
    
    theta_seq[[t]] <- last_theta
    setTxtProgressBar(pb, t)
  }
  
  close(pb)
  theta_seq
}

cat(">> [4/7] Estimating Rolling JGL...\n")
theta_roll_FINAL <- rolling_jgl_estimation(data_list, best_params$l1, best_params$l2, W, S)

# =========================================================
# 4. Benchmark Estimations
# =========================================================

cat(">> [5/7] Running Benchmarks...\n")

# Static JGL (Scale=TRUE)
Y_full_list <- list(safe_scale(do.call(rbind, data_list)))
fit_static <- tryCatch(
  JGL(Y = Y_full_list, penalty = "fused", lambda1 = best_params$l1, lambda2 = best_params$l2),
  error = function(e) NULL
)
theta_static <- if(!is.null(fit_static)) fit_static$theta[[1]] else diag(p)

# CLIME (Robust)
run_clime_safe <- function(X) {
  if(nrow(X) < p/2) return(diag(p))
  fit <- tryCatch(clime(X, standardize = TRUE, nlambda = 5), error = function(e) NULL)
  if(is.null(fit)) return(diag(p))
  idx <- ceiling(length(fit$lambda)/2)
  mat <- as.matrix(fit$Omegalist[[idx]])
  mat + diag(0.05, p)
}

# =========================================================
# 5. Backtesting (GMV + Momentum Tilt) (UPDATED: S truly works)
# =========================================================

compute_weights_gmv <- function(Theta, mu, Lmax){
  ones <- rep(1, p)
  Theta_safe <- Theta + diag(1e-4, p)
  w_gmv <- tryCatch(solve(Theta_safe, ones), error = function(e) ones)
  w_gmv <- as.numeric(w_gmv)
  w_gmv <- w_gmv / sum(w_gmv)
  
  w_final <- w_gmv * (1 + 0.5 * sign(mu))
  w_final <- w_final / sum(abs(w_final)) * Lmax
  w_final
}

# Signal Generation (monthly)
mu_raw <- lapply(data_list, colMeans)
mu_sig <- vector("list", K)
last_mu <- rep(0, p)
for(i in 1:K){
  curr <- mu_raw[[i]]
  last_mu <- last_mu * (1 - EMA_ALPHA) + curr * EMA_ALPHA
  mu_sig[[i]] <- last_mu
}

rets_roll <- c(); rets_stat <- c(); rets_clim <- c(); rets_bnch <- c()

w_r_last <- rep(1/p, p)
w_s_last <- rep(1/p, p)
w_c_last <- rep(1/p, p)

cat(">> [6/7] Backtesting...\n")
for(t in 1:K){
  X_real <- data_list[[t]]
  mu_t   <- mu_sig[[t]]
  
  # Rebalance only every S months; otherwise carry forward weights
  if (((t - 1) %% S) == 0) {
    
    # Rolling uses theta_roll_FINAL[t] (which is also carried forward)
    w_r_last <- compute_weights_gmv(theta_roll_FINAL[[t]], mu_t, LMAX_ROLL)
    
    # Static: rebalance on same schedule for fairness
    w_s_last <- compute_weights_gmv(theta_static, mu_t, LMAX_STATIC)
    
    # CLIME: update theta and rebalance on same schedule
    idx_win <- max(1, t - W + 1):t
    X_win <- do.call(rbind, data_list[idx_win])
    theta_c <- run_clime_safe(X_win)
    w_c_last <- compute_weights_gmv(theta_c, mu_t, LMAX_CLIME)
  }
  
  rets_roll <- c(rets_roll, as.numeric(X_real %*% w_r_last))
  rets_stat <- c(rets_stat, as.numeric(X_real %*% w_s_last))
  rets_clim <- c(rets_clim, as.numeric(X_real %*% w_c_last))
  
  # Equal Weight benchmark always "current month average return"
  rets_bnch <- c(rets_bnch, rowMeans(X_real))
}

# =========================================================
# 6. Visualization & Metrics
# =========================================================

cat(">> [7/7] Plotting Results...\n")

cum_roll <- cumprod(1 + rets_roll)
cum_stat <- cumprod(1 + rets_stat)
cum_clim <- cumprod(1 + rets_clim)
cum_bnch <- cumprod(1 + rets_bnch)

min_len <- min(length(cum_roll), length(cum_stat), length(cum_clim), length(cum_bnch))
dates_plot <- tail(dates, min_len)

sr_r <- mean(rets_roll, na.rm = TRUE) / sd(rets_roll, na.rm = TRUE) * sqrt(252)
sr_s <- mean(rets_stat, na.rm = TRUE) / sd(rets_stat, na.rm = TRUE) * sqrt(252)
sr_c <- mean(rets_clim, na.rm = TRUE) / sd(rets_clim, na.rm = TRUE) * sqrt(252)

# --- Plot 1: Wealth Index ---
op <- par(no.readonly = TRUE); par(mar = c(5, 5, 4, 2))
plot(dates_plot, tail(cum_roll, min_len), type = "l", lwd = 2, col = "blue",
     main = sprintf("FINAL Performance (W=%d months, S=%d months)\nSharpe: RW-JGL=%.2f, Static JGL=%.2f, CLIME=%.2f", W, S, sr_r, sr_s, sr_c),
     ylab = "Wealth Index", xlab = "Date",
     ylim = range(c(tail(cum_roll, min_len), tail(cum_stat, min_len), tail(cum_clim, min_len), tail(cum_bnch, min_len)), na.rm = TRUE))
lines(dates_plot, tail(cum_stat, min_len), col = "red", lwd = 2)
lines(dates_plot, tail(cum_clim, min_len), col = "green", lwd = 1)
lines(dates_plot, tail(cum_bnch, min_len), col = "gray", lwd = 2, lty = 2)
legend("topleft", c("Rolling Window JGL", "Static JGL", "CLIME", "Benchmark"),
       col = c("blue", "red", "green", "gray"), lwd = 2, bty = "n")
par(op)

# --- Plot 2: Box Plot (monthly endpoints) ---
class_end_idx <- vapply(class_indices, function(ix) tail(ix, 1), integer(1))
valid_idx <- class_end_idx[class_end_idx <= length(cum_roll)]

df_boxplot <- data.frame(
  Rolling_Window_JGL = cum_roll[valid_idx],
  Static_JGL         = cum_stat[valid_idx],
  Window_CLIME       = cum_clim[valid_idx],
  Equal_Weight       = cum_bnch[valid_idx]
)

final_vals <- as.numeric(tail(df_boxplot, 1))

op <- par(no.readonly = TRUE)
ymax <- max(df_boxplot, final_vals, na.rm = TRUE)
par(mar = c(5, 5, 4, 2), xpd = NA)

boxplot(df_boxplot,
        col = c("blue", "red", "darkgreen", "gray"),
        main = sprintf("Wealth Index Distribution (W=%d, S=%d)", W, S),
        ylab = "Cumulative Wealth",
        names = c("Rolling Window JGL", "Static JGL", "CLIME", "Benchmark"),
        border = "black",
        medlwd = 1.5,
        ylim = c(0, ymax * 1.12))

points(1:4, final_vals, pch = 19, cex = 1.5, col = "gold")
text(1:4, final_vals, labels = round(final_vals, 1),
     pos = 3, offset = 0.6, cex = 0.8, col = "black")
grid(nx = NA, ny = NULL)
par(op)

# --- Final Printout ---
cat("\n=== FINAL RESULTS ===\n")
cat("Rolling Window JGL Final:", tail(cum_roll, 1), "\n")
cat("Static JGL Final :", tail(cum_stat, 1), "\n")
cat("CLIME Final      :", tail(cum_clim, 1), "\n")
cat("Benchmark Final  :", tail(cum_bnch, 1), "\n")
cat(sprintf("Sharpe: RW-JGL=%.3f | Static=%.3f | CLIME=%.3f\n", sr_r, sr_s, sr_c))