
group_sequential_test_bivariate_transformed <- function(
    lambda_A, theta_A, 
    lambda_B, theta_B,
    phi_A = 0.1, phi_B = 0.1,
    mA1 = 50, mB1 = 50, 
    mA2 = 50, mB2 = 50,
    mA3 = 50, mB3 = 50,
    alpha1 = 0.02, beta1 = 0.015,
    alpha2 = 0.02, beta2 = 0.015,
    alpha3 = 0.01,
    num_simulations = 1000,
    statistic = "T1"  # Can be "T1", "T2", or "T3"
) {
  if (!statistic %in% c("T1", "T2", "T3")) {
    stop("statistic must be one of: 'T1', 'T2', 'T3'")
  }
  nA1 <- mA1
  nB1 <- mB1
  nA2 <- mA1 + mA2
  nB2 <- mB1 + mB2
  nA3 <- nA2 + mA3
  nB3 <- nB2 + mB3
  compute_statistic <- function(X_A, Y_A, X_B, Y_B, n_A, n_B) {
    X_A_bar <- mean(X_A)
    Y_A_bar <- mean(Y_A)
    X_B_bar <- mean(X_B)
    Y_B_bar <- mean(Y_B)
    W1 <- X_A_bar^(2/3) - X_B_bar^(2/3)
    W2 <- Y_A_bar^(2/3) - Y_B_bar^(2/3)
    if (statistic == "T1") {
      return(min(W1, W2))
    } else if (statistic == "T2") {
      return(max(W1, W2))
    } else if (statistic == "T3") {
      return((min(W1, W2) + max(W1, W2))/2)
    }
  }
  set.seed(123)
  T_values_H0 <- matrix(0, nrow = num_simulations, ncol = 3)
  T_values_H1 <- matrix(0, nrow = num_simulations, ncol = 3)
  for (i in 1:num_simulations) {
    X_A1_H0 <- rpois(mA1, lambda_A + phi_A)
    Y_A1_H0 <- rpois(mA1, theta_A + phi_A)
    X_B1_H0 <- rpois(mB1, lambda_A + phi_B)
    Y_B1_H0 <- rpois(mB1, theta_A + phi_B)
    T_values_H0[i, 1] <- compute_statistic(X_A1_H0, Y_A1_H0,
X_B1_H0, Y_B1_H0, nA1, nB1)
    X_A2_H0 <- c(X_A1_H0, rpois(mA2, mean(X_A1_H0) + phi_A))
    Y_A2_H0 <- c(Y_A1_H0, rpois(mA2, mean(Y_A1_H0) + phi_A))
    X_B2_H0 <- c(X_B1_H0, rpois(mB2, mean(X_B1_H0) + phi_B))
    Y_B2_H0 <- c(Y_B1_H0, rpois(mB2, mean(Y_B1_H0) + phi_B))
    T_values_H0[i, 2] <- compute_statistic(X_A2_H0, Y_A2_H0, 
X_B2_H0, Y_B2_H0, nA2, nB2)
    X_A1_H1 <- rpois(mA1, lambda_A + phi_A)
    Y_A1_H1 <- rpois(mA1, theta_A + phi_A)
    X_B1_H1 <- rpois(mB1, lambda_B + phi_B)
    Y_B1_H1 <- rpois(mB1, theta_B + phi_B)
    T_values_H1[i, 1] <- compute_statistic(X_A1_H1, Y_A1_H1, 
X_B1_H1, Y_B1_H1, nA1, nB1)
    X_A2_H1 <- c(X_A1_H1, rpois(mA2, mean(X_A1_H1) + phi_A))
    Y_A2_H1 <- c(Y_A1_H1, rpois(mA2, mean(Y_A1_H1) + phi_A))
    X_B2_H1 <- c(X_B1_H1, rpois(mB2, mean(X_B1_H1) + phi_B))
    Y_B2_H1 <- c(Y_B1_H1, rpois(mB2, mean(Y_B1_H1) + phi_B))
    T_values_H1[i, 2] <- compute_statistic(X_A2_H1, Y_A2_H1, 
X_B2_H1, Y_B2_H1, nA2, nB2)
    X_A3_H1 <- c(X_A2_H1, rpois(mA3, mean(X_A2_H1) + phi_A))
    Y_A3_H1 <- c(Y_A2_H1, rpois(mA3, mean(Y_A2_H1) + phi_A))
    X_B3_H1 <- c(X_B2_H1, rpois(mB3, mean(X_B2_H1) + phi_B))
    Y_B3_H1 <- c(Y_B2_H1, rpois(mB3, mean(Y_B2_H1) + phi_B))
    T_values_H1[i, 3] <- compute_statistic(X_A3_H1, Y_A3_H1, 
X_B3_H1, Y_B3_H1, nA3, nB3)
  }
  prob_function_a1 <- function(a1) mean(
T_values_H0[, 1] < a1) - alpha1
  prob_function_b1 <- function(b1) mean(
T_values_H1[, 1] > b1) - beta1
  a1_optimal <- uniroot(prob_function_a1, 
interval = c(-1000, 1000))$root
  b1_optimal <- uniroot(prob_function_b1, 
interval = c(-1000, 1000))$root
  prob_function_a2 <- function(a2) mean(
T_values_H0[, 2] < a2) - alpha2
  prob_function_b2 <- function(b2) mean(
T_values_H1[, 2] > b2) - beta2
  a2_optimal <- uniroot(prob_function_a2, 
interval = c(-1000, 1000))$root
  b2_optimal <- uniroot(prob_function_b2, 
interval = c(-1000, 1000))$root
  prob_function_b3 <- function(b3) mean(
T_values_H1[, 3] > b3) - alpha3
  b3_optimal <- uniroot(prob_function_b3, 
interval = c(-1000, 1000))$root
  p1_H1 <- 1 - mean(T_values_H1[, 1] >= a1_optimal & 
T_values_H1[, 1] <= b1_optimal)
  p2_H1 <- mean(T_values_H1[, 1] >= a1_optimal & 
T_values_H1[, 1] <= b1_optimal &
                (T_values_H1[, 2] > b2_optimal | 
T_values_H1[, 2] < a2_optimal))
  p3_H1 <- 1 - p1_H1 - p2_H1
  p1_H0 <- mean(T_values_H0[, 1] > b1_optimal |
T_values_H0[, 1] < a1_optimal)
  p2_H0 <- mean(T_values_H0[, 1] >= a1_optimal & 
T_values_H0[, 1] <= b1_optimal &
                (T_values_H0[, 2] > b2_optimal | 
T_values_H0[, 2] < a2_optimal))
  p3_H0 <- 1 - p1_H0 - p2_H0
  n1 <- mA1 + mB1
  n2 <- nA2 + nB2
  n3 <- nA3 + nB3
  ASN_H1 <- p1_H1 * n1 + p2_H1 * n2 + p3_H1 * n3
  ASN_H0 <- p1_H0 * n1 + p2_H0 * n2 + p3_H0 * n3
  list(
    statistic_used = statistic,
    boundaries = c(a1 = a1_optimal, b1 = b1_optimal,
                   a2 = a2_optimal, b2 = b2_optimal,
                   b3 = b3_optimal),
    ASN_H0 = ASN_H0,
    ASN_H1 = ASN_H1,
    stopping_probabilities_H0 = c(p1_H0, p2_H0, p3_H0),
    stopping_probabilities_H1 = c(p1_H1, p2_H1, p3_H1),
    mean_stat_H0 = colMeans(T_values_H0),
    mean_stat_H1 = colMeans(T_values_H1)
  )
}
compute_ASN_table <- function(
    lambda_A = 2,
    theta_A = 2.5,
    lambda_B_start = 2,
    lambda_B_end = 4,
    theta_B_start = 2.5,
    theta_B_end = 4.5,
    step_size = 0.2,
    phi_A = 0.1,
    phi_B = 0.1,
    mA1 = 50, mB1 = 50, 
    mA2 = 50, mB2 = 50,
    mA3 = 50, mB3 = 50,
    alpha1 = 0.02, beta1 = 0.015,
    alpha2 = 0.02, beta2 = 0.015,
    alpha3 = 0.01,
    num_simulations = 1000,
    statistic = "T1"
) {
  lambda_B_seq <- seq(lambda_B_start, lambda_B_end, by = step_size)
  theta_B_seq <- seq(theta_B_start, theta_B_end, by = step_size)
  n_values <- min(length(lambda_B_seq), length(theta_B_seq))
  lambda_B_seq <- lambda_B_seq[1:n_values]
  theta_B_seq <- theta_B_seq[1:n_values]
  results <- data.frame(
    lambda_B = lambda_B_seq,
    theta_B = theta_B_seq,
    ASN_H0 = numeric(n_values),
    ASN_H1 = numeric(n_values),
    stringsAsFactors = FALSE
  )
  for(i in 1:n_values) {
    sink(tempfile())
    result <- group_sequential_test_bivariate_transformed(
      lambda_A = lambda_A,
      theta_A = theta_A,
      lambda_B = lambda_B_seq[i],
      theta_B = theta_B_seq[i],
      phi_A = phi_A,
      phi_B = phi_B,
      mA1 = mA1, mB1 = mB1,
      mA2 = mA2, mB2 = mB2,
      mA3 = mA3, mB3 = mB3,
      alpha1 = alpha1, beta1 = beta1,
      alpha2 = alpha2, beta2 = beta2,
      alpha3 = alpha3,
      num_simulations = num_simulations,
      statistic = statistic
    )
    sink()
    results$ASN_H0[i] <- result$ASN_H0
    results$ASN_H1[i] <- result$ASN_H1
  }
  return(results)
}
resultT1 <- compute_ASN_table(
  lambda_B_start = 2, lambda_B_end = 4,
  theta_B_start = 2.5, theta_B_end = 4.5,
  step_size = 0.2,
  statistic = "T1"  # Can be "T1", "T2", or "T3"
)
resultT1 #min
