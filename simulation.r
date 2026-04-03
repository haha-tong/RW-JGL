# --- Libraries ---
library(JGL)
library(clime)
library(MASS)
library(Matrix)
library(doParallel)
library(foreach)
library(doRNG)
library(R.utils) 
# --- 1. Generate Sparse Precision Matrix with Drift ---
generate_sparse_precision <- function(p, sparsity = 0.2, min_eig = 0.1, noise = 0.24, prev = NULL) {
  if (is.null(prev)) {
    repeat {
      mat <- matrix(0, p, p)
      for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
          if (runif(1) < sparsity) {
            val <- runif(1, -0.5, 0.5)
            mat[i, j] <- val
            mat[j, i] <- val
          }
        }
      }
      diag(mat) <- 1
      eigvals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
      if (all(eigvals > min_eig)) return(mat)
    }
  } else {
    noisy <- prev + matrix(rnorm(p^2, 0, noise), p, p)
    noisy <- (noisy + t(noisy)) / 2
    diag(noisy) <- 1
    noisy <- as.matrix(nearPD(noisy)$mat)
    return(noisy)
  }
}

# --- 2. Simulate Data ---
simulate_multiclass_data <- function(K, p, n_per_class) {
  data_list <- vector("list", K)
  mu_list <- vector("list", K)
  prev_Omega <- NULL
  for (k in 1:K) {
    Omega <- generate_sparse_precision(p, prev = prev_Omega)
    Sigma <- solve(Omega)
    Sigma <- as.matrix(nearPD(Sigma)$mat)
    mu <- rep(0, p)
    mu[1:(p / 2)]     <- 0.36 * sin(pi * k / 10) + 0.08
    mu[(p / 2 + 1):p] <- 0.36 * cos(pi * k / 8)
    if (k > K / 2) mu <- mu + 0.20
    X <- MASS::mvrnorm(n = n_per_class, mu = mu, Sigma = Sigma)
    data_list[[k]] <- X
    mu_list[[k]] <- mu
    prev_Omega <- Omega
  }
  list(data = data_list, mu = mu_list)
}

# --- 3. Estimation ---
estimate_jgl <- function(data_list, lambda1, lambda2) {
  fit <- JGL(Y = data_list, penalty = "fused", lambda1 = lambda1, lambda2 = lambda2)
  fit$theta
}
estimate_clime <- function(X, ridge = 1e-2) {
  fit <- clime(X)
  Theta <- as.matrix(fit$Omegalist[[1]])
  Theta <- Theta + diag(ridge, nrow(Theta))
  Theta
}

# --- 4. Portfolio ---
compute_weights <- function(Theta, mu) {
  w <- tryCatch({ solve(Theta, mu) }, error = function(e) rep(1 / length(mu), length(mu)))
  w / sum(w)
}
simulate_class_portfolio <- function(X, Theta, mu) {
  w <- compute_weights(Theta, mu)
  as.vector(X %*% w)
}

# --- 5. Rolling-window JGL ---
rolling_jgl <- function(data_list, lambda1, lambda2, W = 8, S = 1, center = TRUE, scale. = TRUE) {
  K <- length(data_list)
  theta_seq <- vector("list", K)
  last_theta <- NULL
  for (t in seq_len(K)) {
    if ((t - 1) %% S != 0 && !is.null(last_theta)) {
      theta_seq[[t]] <- last_theta
      next
    }
    i_start <- max(1, t - W + 1)
    idx <- i_start:t
    Ywin <- lapply(data_list[idx], function(X) scale(X, center = center, scale = scale.))
    fit <- tryCatch(
      JGL(Y = Ywin, penalty = "fused", lambda1 = lambda1, lambda2 = lambda2),
      error = function(e) NULL
    )
    if (is.null(fit) || is.null(fit$theta) || length(fit$theta) != length(Ywin)) {
      fit_fallback <- tryCatch(
        JGL(Y = list(Ywin[[length(Ywin)]]), penalty = "fused", lambda1 = lambda1, lambda2 = 0),
        error = function(e) NULL
      )
      theta_t <- if (!is.null(fit_fallback) && length(fit_fallback$theta) == 1) fit_fallback$theta[[1]] else diag(ncol(data_list[[1]]))
    } else {
      theta_t <- fit$theta[[length(fit$theta)]]
    }
    theta_seq[[t]] <- theta_t
    last_theta <- theta_t
  }
  theta_seq
}

# --- 6. One Replication ---
run_experiment <- function(K, p, n, lambda1, lambda2, max_time, W = 8, S = 1) {
  sim_data <- simulate_multiclass_data(K, p, n)
  data_list <- sim_data$data
  mu_list <- sim_data$mu
  
  theta_roll <- rolling_jgl(data_list, lambda1, lambda2, W = W, S = S)
  theta_static <- estimate_jgl(list(data_list[[1]], data_list[[1]]), lambda1, lambda2)[[1]]
  theta_clime_list <- lapply(data_list, estimate_clime)
  
  rets_roll <- rets_static <- rets_clime <- c()
  for (i in 1:K) {
    rets_roll    <- c(rets_roll,    simulate_class_portfolio(data_list[[i]], theta_roll[[i]],        mu_list[[i]]))
    rets_static  <- c(rets_static,  simulate_class_portfolio(data_list[[i]], theta_static,          mu_list[[i]]))
    rets_clime   <- c(rets_clime,   simulate_class_portfolio(data_list[[i]], theta_clime_list[[i]], mu_list[[i]]))
  }
  
  list(
    jgl_roll = cumprod(1 + rets_roll  )[1:max_time],
    static   = cumprod(1 + rets_static)[1:max_time],
    clime    = cumprod(1 + rets_clime )[1:max_time]
  )
}

# --- 6.5 StARS ---
compute_stars_instability <- function(data_list, lambda1_grid, lambda2_grid,
                                      n_subsamples = 10,
                                      subsample_frac = 0.8,
                                      timeout_sec = 120,
                                      k_per_stars = 8,
                                      min_success = 2,
                                      tol = 1e-4, maxiter = 5000) {
  p <- ncol(data_list[[1]])
  stars_matrix <- matrix(NA_real_, length(lambda1_grid), length(lambda2_grid))
  
  for (i in seq_along(lambda1_grid)) {
    for (j in seq_along(lambda2_grid)) {
      cat(sprintf(">> [StARS] Evaluating lambda1 = %.3f , lambda2 = %.3f ...\n",
                  lambda1_grid[i], lambda2_grid[j]))
      edges_list <- list()
      
      for (b in 1:n_subsamples) {
        subsample <- lapply(data_list, function(X) {
          X[sample(seq_len(nrow(X)), floor(subsample_frac * nrow(X))), , drop = FALSE]
        })
        subsample <- lapply(subsample, function(X) scale(X, center = TRUE, scale = TRUE))
        k_take <- min(k_per_stars, length(subsample))
        k_idx  <- sort(sample.int(length(subsample), size = k_take, replace = FALSE))
        subsample_k <- subsample[k_idx]
        
        fit <- tryCatch({
          withTimeout({
            JGL(Y = subsample_k, penalty = "fused",
                lambda1 = lambda1_grid[i], lambda2 = lambda2_grid[j],
                tol = tol, maxiter = maxiter)
          }, timeout = timeout_sec, onTimeout = "error")
        }, error = function(e) NULL)
        
        if (!is.null(fit) && !is.null(fit$theta) && length(fit$theta) == length(subsample_k)) {
          edges <- lapply(fit$theta, function(mat) {
            m <- as.matrix(mat)
            if (!is.matrix(m) || !all(dim(m) == c(p, p))) return(NULL)
            (abs(m) > 1e-6) * 1L
          })
          if (all(vapply(edges, function(M) is.matrix(M) && all(dim(M) == c(p, p)), TRUE))) {
            edges_list[[length(edges_list) + 1L]] <- edges
          }
        }
      }
      
      if (length(edges_list) < min_success) {
        stars_matrix[i, j] <- Inf
        next
      }
      k_lengths <- vapply(edges_list, length, 1L)
      if (length(unique(k_lengths)) != 1L) {
        stars_matrix[i, j] <- Inf
        next
      }
      ok_shape <- all(vapply(unlist(edges_list, recursive = FALSE),
                             function(M) is.matrix(M) && all(dim(M) == c(p, p)), TRUE))
      if (!ok_shape) {
        stars_matrix[i, j] <- Inf
        next
      }
      
      avg_edges_per_b <- lapply(edges_list, function(edge_set) {
        Reduce(`+`, edge_set) / length(edge_set)
      })
      edge_counts <- Reduce(`+`, avg_edges_per_b) / length(avg_edges_per_b)
      
      P <- edge_counts
      instability <- 4 * P * (1 - P)
      stars_matrix[i, j] <- mean(instability[lower.tri(instability)], na.rm = TRUE)
    }
  }
  
  if (all(!is.finite(stars_matrix))) {
    warning("All instability evaluations failed. Consider increasing timeout or reducing k_per_stars.")
    return(list(lambda1 = lambda1_grid[1], lambda2 = lambda2_grid[1], stars = stars_matrix))
  }
  
  idx <- which(stars_matrix == min(stars_matrix, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(lambda1 = lambda1_grid[idx[1]], lambda2 = lambda2_grid[idx[2]], stars = stars_matrix)
}

# --- 7. Parameters and Run ---
replications <- 50
K <- 30
p <- 20
n <- 100
max_time <- 100
set.seed(9375)

lambda1_grid <- seq(0.2, 0.40, length.out = 4)
lambda2_grid <- seq(0.10, 0.16, length.out = 4)

init_data <- simulate_multiclass_data(K, p, n)$data

stars_result <- compute_stars_instability(
  data_list = init_data,
  lambda1_grid = lambda1_grid,
  lambda2_grid = lambda2_grid,
  n_subsamples = 10,
  subsample_frac = 0.8,
  timeout_sec = 120,
  k_per_stars = 8,
  min_success = 2
)

lambda1 <- stars_result$lambda1
lambda2 <- stars_result$lambda2
cat(">> [StARS] Selected Lambda1:", lambda1, " Lambda2:", lambda2, "\n")

W <- 8
S <- 1

# --- 8. Replications (parallel) ---
cl <- makeCluster(max(1, detectCores() - 1))
registerDoParallel(cl)
results <- foreach(r = 1:replications, .combine = rbind,
                   .packages = c("JGL", "MASS", "Matrix", "clime")) %dorng% {
                     cat(sprintf(">> Running replication %d...\n", r))
                     res <- run_experiment(K, p, n, lambda1, lambda2, max_time, W = W, S = S)
                     c(res$jgl_roll, res$static, res$clime)
                   }
stopCluster(cl)

# --- 9. Plot Time Series ---
results_jgl    <- results[, 1:max_time, drop = FALSE]
results_static <- results[, (max_time + 1):(2 * max_time), drop = FALSE]
results_clime  <- results[, (2 * max_time + 1):(3 * max_time), drop = FALSE]

mean_jgl    <- colMeans(results_jgl)
mean_static <- colMeans(results_static)
mean_clime  <- colMeans(results_clime)

time <- 1:max_time

col_jgl    <- "blue"
col_static <- "red"
col_clime  <- "green"

plot(time, mean_jgl, type = "l", lwd = 2, col = col_jgl,
     ylim = range(c(mean_jgl, mean_static, mean_clime)),
     ylab = "Cumulative Return", xlab = "Time",
     main = sprintf("Average Portfolio Return with Rolling Window JGL (W=%d, S=%d)", W, S))
lines(time, mean_static, lwd = 2, col = col_static)
lines(time, mean_clime,  lwd = 2, col = col_clime)

legend("topleft",
       legend = c("Rolling Window JGL", "Static JGL", "CLIME"),
       col = c(col_jgl, col_static, col_clime),
       lty = 1, lwd = 2, bty = "n")

# --- 10. Boxplot of Final Returns ---
final_returns <- data.frame(
  Rolling_Window_JGL = results[, max_time],
  Static_JGL  = results[, 2 * max_time],
  CLIME       = results[, 3 * max_time]
)
boxplot(final_returns,
        col = c(col_jgl, col_static, col_clime),
        main = "Final Cumulative Return",
        ylab = "Cumulative Return")
points(x = 1:3, y = colMeans(final_returns),
       col = c(col_jgl, col_static, col_clime), pch = 19)
