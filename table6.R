
es_adaptive_groupseq <- function(lambda_A, lambda_B, mA1 = 50, mB1 = 50, mA2 = 50, mB2 = 50, mA3 = 50, mB3 = 50,
                                 alpha1 = 0.03, beta1 = 0.02, alpha2 = 0.015, beta2 = 0.025, alpha3 = 0.005,
                                 num_simulations = 1000, compute_Z = NULL) {
  nA2 <- mA1 + mA2  
  nB2 <- mB1 + mB2 
  nA3 <- nA2 + mA3  
  nB3 <- nB2 + mB3 
  if (is.null(compute_Z)) {
    compute_Z <- function(X_A, X_B, n_A, n_B) {
      X_A_bar <- mean(X_A)
      X_B_bar <- mean(X_B)
      lambda_hat <- (n_A * X_A_bar + n_B * X_B_bar) / (n_A + n_B)
      numerator <- X_A_bar^(2/3) - X_B_bar^(2/3)
      denominator <- sqrt((4/9) * lambda_hat^(1/3) * (1/n_A + 1/n_B))
      numerator / denominator
    }
  }
  set.seed(123)  
  Z1_values_H0 <- numeric(num_simulations)
  Z2_values_H0 <- numeric(num_simulations)
  Z1_values_H1 <- numeric(num_simulations)
  Z2_values_H1 <- numeric(num_simulations)
  Z3_values_H1 <- numeric(num_simulations)
  for (i in 1:num_simulations) {
    X_A1_H0 <- rpois(mA1, lambda_A)
    X_B1_H0 <- rpois(mB1, lambda_A)
    Z1_values_H0[i] <- compute_Z(X_A1_H0, X_B1_H0, mA1, mB1)
    X_A2_H0 <- c(X_A1_H0, rpois(mA2, mean(X_A1_H0)))
    X_B2_H0 <- c(X_B1_H0, rpois(mB2, mean(X_A1_H0)))
    Z2_values_H0[i] <- compute_Z(X_A2_H0, X_B2_H0, nA2, nB2)
    X_A1_H1 <- rpois(mA1, lambda_A)
    X_B1_H1 <- rpois(mB1, lambda_B)
    Z1_values_H1[i] <- compute_Z(X_A1_H1, X_B1_H1, mA1, mB1)
    X_A2_H1 <- c(X_A1_H1, rpois(mA2, mean(X_A1_H1)))
    X_B2_H1 <- c(X_B1_H1, rpois(mB2, mean(X_B1_H1)))
    Z2_values_H1[i] <- compute_Z(X_A2_H1, X_B2_H1, nA2, nB2)
    X_A3_H1 <- c(X_A2_H1, rpois(mA3, mean(X_A1_H1)))
    X_B3_H1 <- c(X_B2_H1, rpois(mB3, mean(X_B1_H1)))
    Z3_values_H1[i] <- compute_Z(X_A3_H1, X_B3_H1, nA3, nB3)
  }
  prob_function_a1 <- function(a1) mean(
Z1_values_H0 < a1) - alpha1
  prob_function_b1 <- function(b1) mean(
Z1_values_H0 > b1) - beta1
  a1_optimal <- uniroot(prob_function_a1, 
interval = c(-100, 100))$root
  b1_optimal <- uniroot(prob_function_b1, 
interval = c(-100, 100))$root 
  prob_function_a2 <- function(a2) mean(
Z2_values_H0 < a2) - alpha2
  prob_function_b2 <- function(b2) mean(
Z2_values_H1 > b2) - beta2
  a2_optimal <- uniroot(prob_function_a2,
interval = c(-100, 100))$root
  b2_optimal <- uniroot(prob_function_b2, 
interval = c(-100, 100))$root
  prob_function_b3 <- function(b3) mean(
Z3_values_H1 > b3) - alpha3
  b3_optimal <- uniroot(prob_function_b3, 
interval = c(-100, 100))$root
  p1_H1 <- 1 - mean((Z1_values_H1 >= a1_optimal) & (
Z1_values_H1 <= b1_optimal))
  p2_H1 <- mean((Z1_values_H1 >= a1_optimal) & (
Z1_values_H1 <= b1_optimal) & 
                  ((Z2_values_H1 < a2_optimal) | (
Z2_values_H1 > b2_optimal)))
  p3_H1 <- 1 - p1_H1 - p2_H1
  p1_H0 <- 1 - mean((Z1_values_H0 >= a1_optimal) & (
Z1_values_H0 <= b1_optimal))
  p2_H0 <- mean((Z1_values_H0 >= a1_optimal) & (
Z1_values_H0 <= b1_optimal) & 
                  ((Z2_values_H0 < a2_optimal) | (
Z2_values_H0 > b2_optimal)))
  p3_H0 <- 1 - p1_H0 - p2_H0
  n1 <- mA1 + mB1  
  n2 <- nA2 + nB2  
  n3 <- nA3 + nB3  
  ASN_H1 <- p1_H1 * n1 + p2_H1 * n2 + p3_H1 * n3
  ASN_H0 <- p1_H0 * n1 + p2_H0 * n2 + p3_H0 * n3
  ASN_table <- data.frame(
    Stage = c("Stage 1", "Stage 2", "Stage 3"),
    Sample_Size = c(n1, n2, n3),
    Probability_H1 = c(p1_H1, p2_H1, p3_H1),
    Probability_H0 = c(p1_H0, p2_H0, p3_H0)
  )  
  cat("Optimal a1:", a1_optimal, "\n")
  cat("Optimal b1:", b1_optimal, "\n")
  cat("Optimal a2:", a2_optimal, "\n")
  cat("Optimal b2:", b2_optimal, "\n")
  cat("Optimal b3:", b3_optimal, "\n")
  print(ASN_table)
  cat("Average Sample Number (ASN_H1):", ASN_H1, "\n")
  cat("Average Sample Number (ASN_H0):", ASN_H0, "\n")
  return(list(ASN_H0 = ASN_H0, ASN_H1 = ASN_H1))
}
compute_ASN_table_ad <- function(
lambda_A = 2, lambda_B_seq = seq(2, 4, by = 0.2), 
mA1 = 50, mB1 = 50, mA2 = 50, mB2 = 50, mA3 = 50, mB3 = 50,
alpha1 = 0.03, beta1 = 0.02, alpha2 = 0.015, beta2 = 0.025, alpha3 = 0.005,
num_simulations = 10000, compute_Z = NULL) {
  results <- data.frame(
    lambda_B = numeric(),
    ASN_H0 = numeric(),
    ASN_H1 = numeric()
  )
  for (lambda_B in lambda_B_seq) {
    sink(tempfile())  
    result <- es_adaptive_groupseq(
lambda_A, lambda_B, mA1, mB1, mA2, mB2, mA3, mB3,
alpha1, beta1, alpha2, beta2, alpha3, num_simulations, compute_Z)
    sink()
    results <- rbind(results, data.frame(
      lambda_B = lambda_B,
      ASN_H0 = result$ASN_H0,
      ASN_H1 = result$ASN_H1
    ))
  }
  return(results)
}
compute_ASN_table_ad(lambda_A = 2, lambda_B_seq = seq(2, 4, by = 0.2),
                     alpha1 = 0.03, beta1 = 0.02, alpha2 = 0.015, 
                     beta2 = 0.025, alpha3 = 0.005)
