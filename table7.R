
fully_adaptive_group_seq <- function(lambdaA, lambdaB, 
                               RA_fun = RA,
                               sample_size_interim = 100,
                               alpha1 = 0.03, beta1 = 0.02,
                               alpha2 = 0.015, beta2 = 0.025, 
                               alpha3 = 0.005,
                               num_simulations = 1000) {
  set.seed(123)
  compute_Z <- function(X_A, X_B) {
    n_A <- length(X_A)
    n_B <- length(X_B)
    X_A_bar <- mean(X_A)
    X_B_bar <- mean(X_B)
    lambda_hat <- (n_A * X_A_bar + n_B * X_B_bar) / (n_A + n_B)
    numerator <- X_A_bar^(2/3) - X_B_bar^(2/3)
    denominator <- sqrt((4/9) * lambda_hat^(1/3) * (1/n_A + 1/n_B))
    numerator / denominator
  }
  adap_alloc <- function(lambdaA, lambdaB, samp_size = 100) {
    X_A <- rpois(4, lambdaA)
    X_B <- rpois(4, lambdaB)
    n_remaining <- samp_size - 8
    new_A <- numeric(n_remaining)
    new_B <- numeric(n_remaining)
    for (i in 1:n_remaining) {
      lambdaA_hat <- mean(X_A)
      lambdaB_hat <- mean(X_B)
      RA_val <- RA_fun(lambdaA_hat, lambdaB_hat)
      if (runif(1) <= RA_val) {
        new_A[i] <- rpois(1, lambdaA_hat)
        X_A <- c(X_A, new_A[i])
      } else {
        new_B[i] <- rpois(1, lambdaB_hat)
        X_B <- c(X_B, new_B[i])
      }
    }
    list(SampleA = X_A, 
         SampleB = X_B, 
         CountA = length(X_A), 
         PropA = length(X_A)/samp_size)
  }
  results <- list(
    Z1_H0 = numeric(num_simulations),
    Z2_H0 = numeric(num_simulations),
    Z1_H1 = numeric(num_simulations),
    Z2_H1 = numeric(num_simulations),
    Z3_H1 = numeric(num_simulations),
    RA1 = numeric(num_simulations),
    RA2 = numeric(num_simulations),
    RA3 = numeric(num_simulations),
    PropA1 = numeric(num_simulations),
    PropA2 = numeric(num_simulations),
    PropA3 = numeric(num_simulations)
  )
  for (i in 1:num_simulations) {
    alloc1_H0 <- adap_alloc(lambdaA, lambdaA, 
sample_size_interim)
    results$Z1_H0[i] <- compute_Z(alloc1_H0$SampleA, 
alloc1_H0$SampleB)
    alloc2_H0 <- adap_alloc(
mean(alloc1_H0$SampleA), mean(alloc1_H0$SampleA), 
                                 sample_size_interim)
    results$Z2_H0[i] <- compute_Z(
c(alloc1_H0$SampleA, alloc2_H0$SampleA),
c(alloc1_H0$SampleB, alloc2_H0$SampleB))
    alloc1_H1 <- adap_alloc(lambdaA, lambdaB, sample_size_interim)
    results$Z1_H1[i] <- compute_Z(
alloc1_H1$SampleA, alloc1_H1$SampleB)
    results$RA1[i] <- RA_fun(
mean(alloc1_H1$SampleA), mean(alloc1_H1$SampleB))
    results$PropA1[i] <- alloc1_H1$PropA
    alloc2_H1 <- adap_alloc(
mean(alloc1_H1$SampleA), mean(alloc1_H1$SampleB), 
                                 sample_size_interim)
    combined_A2 <- c(alloc1_H1$SampleA, alloc2_H1$SampleA)
    combined_B2 <- c(alloc1_H1$SampleB, alloc2_H1$SampleB)
    results$Z2_H1[i] <- compute_Z(combined_A2, combined_B2)
    results$RA2[i] <- RA_fun(
mean(combined_A2), mean(combined_B2))
    results$PropA2[i] <- length(
combined_A2)/(2*sample_size_interim)
    alloc3_H1 <- adap_alloc(
mean(combined_A2), mean(combined_B2), 
                                 sample_size_interim)
    combined_A3 <- c(combined_A2, alloc3_H1$SampleA)
    combined_B3 <- c(combined_B2, alloc3_H1$SampleB)
    results$Z3_H1[i] <- compute_Z(
combined_A3, combined_B3)
    results$RA3[i] <- RA_fun(
mean(combined_A3), mean(combined_B3))
    results$PropA3[i] <- length(
combined_A3)/(3*sample_size_interim)
  }
  prob_function_a1 <- function(a1) mean(
results$Z1_H0 < a1) - alpha1
  prob_function_b1 <- function(b1) mean(
results$Z1_H0 > b1) - beta1
  a1_optimal <- uniroot(prob_function_a1, 
interval = c(-100, 100))$root
  b1_optimal <- uniroot(prob_function_b1, 
interval = c(-100, 100))$root
  prob_function_a2 <- function(a2) mean(
results$Z2_H0 < a2) - alpha2
  prob_function_b2 <- function(b2) mean(
results$Z2_H1 > b2) - beta2
  a2_optimal <- uniroot(prob_function_a2, 
interval = c(-100, 100))$root
  b2_optimal <- uniroot(prob_function_b2, 
interval = c(-100, 1000))$root
  prob_function_b3 <- function(b3) mean(
results$Z3_H1 > b3) - alpha3
  b3_optimal <- uniroot(prob_function_b3,
interval = c(-100, 100))$root
  p1_H1 <- 1 - mean((results$Z1_H1 >= a1_optimal) & (
results$Z1_H1 <= b1_optimal))
  p2_H1 <- mean((results$Z1_H1 >= a1_optimal) & (
results$Z1_H1 <= b1_optimal) & 
                  ((results$Z2_H1 < a2_optimal) | (
results$Z2_H1 > b2_optimal)))
  p3_H1 <- 1 - p1_H1 - p2_H1
  p1_H0 <- 1 - mean((results$Z1_H0 >= a1_optimal) & (
results$Z1_H0 <= b1_optimal))
  p2_H0 <- mean((results$Z1_H0 >= a1_optimal) & (
results$Z1_H0 <= b1_optimal) & 
                  ((results$Z2_H0 < a2_optimal) | (
results$Z2_H0 > b2_optimal)))
  p3_H0 <- 1 - p1_H0 - p2_H0
  n1 <- sample_size_interim
  n2 <- 2 * sample_size_interim
  n3 <- 3 * sample_size_interim
  ASN_H1 <- p1_H1 * n1 + p2_H1 * n2 + p3_H1 * n3
  ASN_H0 <- p1_H0 * n1 + p2_H0 * n2 + p3_H0 * n3
  ASN_table <- data.frame(
    Stage = c("Stage 1", "Stage 2", "Stage 3"),
    Sample_Size = c(n1, n2, n3),
    Probability_H1 = c(p1_H1, p2_H1, p3_H1),
    Probability_H0 = c(p1_H0, p2_H0, p3_H0),
    RA_values = c(mean(results$RA1), mean(
results$RA2), mean(results$RA3)),
    Proportions_A = c(mean(results$PropA1), mean(
results$PropA2), mean(results$PropA3)),
    PropA_se = c(sd(results$PropA1),
sd(results$PropA2), sd(results$PropA3))  
  )
  list(
    boundaries = c(a1 = a1_optimal, b1 = b1_optimal, 
                   a2 = a2_optimal, b2 = b2_optimal, 
                   b3 = b3_optimal),
    ASN = list(H0 = ASN_H0, H1 = ASN_H1),
    ASN_table = ASN_table
  )
}
RA <- function(lambdaA, lambdaB) {
  (lambdaB + 0.5) / (lambdaA + lambdaB + 1) 
}
library(dplyr)
lambda_A <- 2
lambda_B_values <- seq(2, 6, 0.5)
results_table <- lapply(lambda_B_values, function(lambda_B) {
  result <- fully_adaptive_group_seq(lambdaA = lambda_A, 
            lambdaB = lambda_B, num_simulations = 1000, RA_fun = RA)
  data.frame(
    lambda_B = lambda_B,
    ASN_H0 = result$ASN$H0,
    ASN_H1 = result$ASN$H1,
    RA = result$ASN_table$RA_values[3],
    Prop_A = result$ASN_table$Proportions_A[3],
    Prop_A_SE = result$ASN_table$PropA_se[3]  
  )
}) %>% bind_rows()
print(results_table)
