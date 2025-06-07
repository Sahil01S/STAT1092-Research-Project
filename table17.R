
compute_stats <- function(X_A, Y_A, X_B, Y_B, statistic = "T1") {
  X_A_bar <- mean(X_A)
  Y_A_bar <- mean(Y_A)
  X_B_bar <- mean(X_B)
  Y_B_bar <- mean(Y_B)
  
  W1 <- X_A_bar^(2/3) - X_B_bar^(2/3)
  W2 <- Y_A_bar^(2/3) - Y_B_bar^(2/3)
  
  if (statistic == "T1") return(min(W1, W2))
  if (statistic == "T2") return(max(W1, W2))
  if (statistic == "T3") return((min(W1, W2) + max(W1, W2))/2)
}





# ==============================================================
# Bivariate Adaptive Group Sequential Design
# ==============================================================
# This implementation features:
# 1. Response-adaptive randomization based on Poisson outcomes
# 2. Three-stage design with interim analyses
# 3. T1/T2/T3 test statistics for bivariate endpoints
# 4. Tracking of allocation proportions and RA values

# --------------------------
# Allocation Rules
# --------------------------

RA <- function(lambdaA, lambdaB, thetaA, thetaB) {
  (lambdaB + thetaB + 0.5) / 
    (lambdaA + thetaA + lambdaB + thetaB + 1)
}


# --------------------------
# Core Adaptive Design Function  
# --------------------------

fully_adaptive_bivariate <- function(
    lambdaA, thetaA, lambdaB, thetaB,
    RA_fun = RA,
    sample_size_interim = 100,
    alpha1 = 0.03, beta1 = 0.02,
    alpha2 = 0.015, beta2 = 0.025,
    alpha3 = 0.005,
    num_simulations = 1000,
    statistic = "T3") {
  
  set.seed(123) # For reproducibility
  
  # --------------------------
  # Helper Functions
  # --------------------------
  
  # Adaptive allocation procedure
  adap_alloc <- function(lambdaA, thetaA, lambdaB, thetaB, samp_size = 100) {
    # Initial samples (4 per group for stability)
    X_A <- rpois(4, lambdaA)
    Y_A <- rpois(4, thetaA)
    X_B <- rpois(4, lambdaB)
    Y_B <- rpois(4, thetaB)
    
    # Adaptive allocation of remaining samples
    n_remaining <- samp_size - 8
    for (i in 1:n_remaining) {
      # Get current estimates
      lambdaA_hat <- mean(X_A)
      thetaA_hat <- mean(Y_A)
      lambdaB_hat <- mean(X_B)
      thetaB_hat <- mean(Y_B)
      
      # Calculate allocation probability
      RA_val <- RA_fun(lambdaA_hat, lambdaB_hat, thetaA_hat, thetaB_hat)
      
      # Assign next subject
      if (runif(1) <= RA_val) {
        X_A <- c(X_A, rpois(1, lambdaA_hat))
        Y_A <- c(Y_A, rpois(1, thetaA_hat))
      } else {
        X_B <- c(X_B, rpois(1, lambdaB_hat))
        Y_B <- c(Y_B, rpois(1, thetaB_hat))
      }
    }
    
    # Return allocation results
    list(
      SampleA_X = X_A, SampleA_Y = Y_A,
      SampleB_X = X_B, SampleB_Y = Y_B,
      CountA = length(X_A),
      PropA = length(X_A)/samp_size,
      RA = RA_fun(mean(X_A), mean(X_B), mean(Y_A), mean(Y_B))
    )
  }
  
  # --------------------------
  # Simulation Setup
  # --------------------------
  
  # Initialize results storage
  results <- list(
    T1_H0 = numeric(num_simulations),
    T2_H0 = numeric(num_simulations),
    T1_H1 = numeric(num_simulations),
    T2_H1 = numeric(num_simulations),
    T3_H1 = numeric(num_simulations),
    RA1 = numeric(num_simulations),
    RA2 = numeric(num_simulations),
    RA3 = numeric(num_simulations),
    PropA1 = numeric(num_simulations),
    PropA2 = numeric(num_simulations),
    PropA3 = numeric(num_simulations)
  )
  
  # --------------------------
  # Main Simulation Loop
  # --------------------------
  
  for (i in 1:num_simulations) {
    # Null hypothesis simulations (lambdaA = lambdaB, thetaA = thetaB)
    alloc1_H0 <- adap_alloc(lambdaA, thetaA, lambdaA, thetaA, sample_size_interim)
    results$T1_H0[i] <- compute_stats(alloc1_H0$SampleA_X, alloc1_H0$SampleA_Y,
                                     alloc1_H0$SampleB_X, alloc1_H0$SampleB_Y, statistic)
    
    alloc2_H0 <- adap_alloc(mean(alloc1_H0$SampleA_X), mean(alloc1_H0$SampleA_Y),
                           mean(alloc1_H0$SampleA_X), mean(alloc1_H0$SampleA_Y), sample_size_interim)
    results$T2_H0[i] <- compute_stats(
      c(alloc1_H0$SampleA_X, alloc2_H0$SampleA_X),
      c(alloc1_H0$SampleA_Y, alloc2_H0$SampleA_Y),
      c(alloc1_H0$SampleB_X, alloc2_H0$SampleB_X),
      c(alloc1_H0$SampleB_Y, alloc2_H0$SampleB_Y), 
      statistic
    )
    
    # Alternative hypothesis simulations
    alloc1_H1 <- adap_alloc(lambdaA, thetaA, lambdaB, thetaB, sample_size_interim)
    results$T1_H1[i] <- compute_stats(
      alloc1_H1$SampleA_X, alloc1_H1$SampleA_Y,
      alloc1_H1$SampleB_X, alloc1_H1$SampleB_Y, 
      statistic
    )
    results$RA1[i] <- alloc1_H1$RA
    results$PropA1[i] <- alloc1_H1$PropA
    
    # Stage 2 allocation
    alloc2_H1 <- adap_alloc(
      mean(alloc1_H1$SampleA_X), mean(alloc1_H1$SampleA_Y),
      mean(alloc1_H1$SampleB_X), mean(alloc1_H1$SampleB_Y), 
      sample_size_interim
    )
    combined_A_X <- c(alloc1_H1$SampleA_X, alloc2_H1$SampleA_X)
    combined_A_Y <- c(alloc1_H1$SampleA_Y, alloc2_H1$SampleA_Y)
    combined_B_X <- c(alloc1_H1$SampleB_X, alloc2_H1$SampleB_X)
    combined_B_Y <- c(alloc1_H1$SampleB_Y, alloc2_H1$SampleB_Y)
    
    results$T2_H1[i] <- compute_stats(
      combined_A_X, combined_A_Y, 
      combined_B_X, combined_B_Y, 
      statistic
    )
    results$RA2[i] <- alloc2_H1$RA
    results$PropA2[i] <- length(combined_A_X)/(2*sample_size_interim)
    
    # Stage 3 allocation
    alloc3_H1 <- adap_alloc(
      mean(combined_A_X), mean(combined_A_Y),
      mean(combined_B_X), mean(combined_B_Y),
      sample_size_interim
    )
    combined_A_X3 <- c(combined_A_X, alloc3_H1$SampleA_X)
    combined_A_Y3 <- c(combined_A_Y, alloc3_H1$SampleA_Y)
    combined_B_X3 <- c(combined_B_X, alloc3_H1$SampleB_X)
    combined_B_Y3 <- c(combined_B_Y, alloc3_H1$SampleB_Y)
    
    results$T3_H1[i] <- compute_stats(
      combined_A_X3, combined_A_Y3,
      combined_B_X3, combined_B_Y3,
      statistic
    )
    results$RA3[i] <- alloc3_H1$RA
    results$PropA3[i] <- length(combined_A_X3)/(3*sample_size_interim)
  }
  
  # --------------------------
  # Boundary Determination
  # --------------------------
  
  # Stage 1 boundaries
  prob_function_a1 <- function(a1) mean(results$T1_H0 < a1) - alpha1
  prob_function_b1 <- function(b1) mean(results$T1_H0 > b1) - beta1
  a1_optimal <- uniroot(prob_function_a1, interval = c(-100, 100))$root
  b1_optimal <- uniroot(prob_function_b1, interval = c(-100, 100))$root
  
  # Stage 2 boundaries
  prob_function_a2 <- function(a2) mean(results$T2_H0 < a2) - alpha2
  prob_function_b2 <- function(b2) mean(results$T2_H1 > b2) - beta2
  a2_optimal <- uniroot(prob_function_a2, interval = c(-100, 100))$root
  b2_optimal <- uniroot(prob_function_b2, interval = c(-100, 100))$root
  
  # Final boundary
  prob_function_b3 <- function(b3) mean(results$T3_H1 > b3) - alpha3
  b3_optimal <- uniroot(prob_function_b3, interval = c(-100, 100))$root
  
  # --------------------------
  # Performance Metrics
  # --------------------------
  
  # Stopping probabilities under H1
  p1_H1 <- 1 - mean((results$T1_H1 >= a1_optimal) & (results$T1_H1 <= b1_optimal))
  p2_H1 <- mean(
    (results$T1_H1 >= a1_optimal) & (results$T1_H1 <= b1_optimal) &
    ((results$T2_H1 < a2_optimal) | (results$T2_H1 > b2_optimal))
  )
  p3_H1 <- 1 - p1_H1 - p2_H1
  
  # Stopping probabilities under H0
  p1_H0 <- 1 - mean((results$T1_H0 >= a1_optimal) & (results$T1_H0 <= b1_optimal))
  p2_H0 <- mean(
    (results$T1_H0 >= a1_optimal) & (results$T1_H0 <= b1_optimal) &
    ((results$T2_H0 < a2_optimal) | (results$T2_H0 > b2_optimal))
  )
  p3_H0 <- 1 - p1_H0 - p2_H0
  
  # Average sample numbers
  n1 <- sample_size_interim
  n2 <- 2 * sample_size_interim
  n3 <- 3 * sample_size_interim
  ASN_H1 <- p1_H1 * n1 + p2_H1 * n2 + p3_H1 * n3
  ASN_H0 <- p1_H0 * n1 + p2_H0 * n2 + p3_H0 * n3
  
  # --------------------------
  # Results Compilation
  # --------------------------
  
  ASN_table <- data.frame(
    Stage = c("Stage 1", "Stage 2", "Stage 3"),
    Sample_Size = c(n1, n2, n3),
    Probability_H1 = c(p1_H1, p2_H1, p3_H1),
    Probability_H0 = c(p1_H0, p2_H0, p3_H0),
    RA_mean = c(mean(results$RA1), mean(results$RA2), mean(results$RA3)),
    RA_se = c(sd(results$RA1), sd(results$RA2), sd(results$RA3)),
    PropA_mean = c(mean(results$PropA1), mean(results$PropA2), mean(results$PropA3)),
    PropA_se = c(sd(results$PropA1), sd(results$PropA2), sd(results$PropA3))
  )
  
  list(
    boundaries = c(
      a1 = a1_optimal, b1 = b1_optimal,
      a2 = a2_optimal, b2 = b2_optimal,
      b3 = b3_optimal
    ),
    ASN = list(H0 = ASN_H0, H1 = ASN_H1),
    ASN_table = ASN_table,
    statistic_used = statistic,
    allocation_rule = deparse(substitute(RA_fun)),
    raw_results = results
  )
}

# --------------------------
# Table Generation Function
# --------------------------

#' Generate Performance Table Across Parameter Space
#' @param lambda_A Baseline rate for X in group A
#' @param theta_A Baseline rate for Y in group A
#' @param lambda_B_start Starting value for lambda_B sequence
#' @param lambda_B_end Ending value for lambda_B sequence  
#' @param theta_B_start Starting value for theta_B sequence
#' @param theta_B_end Ending value for theta_B sequence
#' @param step_size Increment for parameter sequences
#' @param ... Other parameters passed to fully_adaptive_bivariate
#' @return Data frame with performance metrics across parameter combinations
compute_adaptive_table <- function(
    lambda_A = 2,
    theta_A = 2.5,
    lambda_B_start = 2,
    lambda_B_end = 4,
    theta_B_start = 2.5,
    theta_B_end = 4.5,
    step_size = 0.2,
    sample_size_interim = 100,
    alpha1 = 0.03, beta1 = 0.02,
    alpha2 = 0.015, beta2 = 0.025,
    alpha3 = 0.005,
    num_simulations = 1000,
    statistic = "T3",
    allocation_rule = "RA") {
  
  # Create parameter sequences
  lambda_B_seq <- seq(lambda_B_start, lambda_B_end, by = step_size)
  theta_B_seq <- seq(theta_B_start, theta_B_end, by = step_size)
  n_values <- min(length(lambda_B_seq), length(theta_B_seq))
  
  # Initialize results
  results <- data.frame(
    lambda_B = lambda_B_seq[1:n_values],
    theta_B = theta_B_seq[1:n_values],
    ASN_H0 = numeric(n_values),
    ASN_H1 = numeric(n_values),
    RA = numeric(n_values),
    Prop_A = numeric(n_values),
    Prop_A_SE = numeric(n_values),
    stringsAsFactors = FALSE
  )
  
  # Select allocation rule
  alloc_func <- if (allocation_rule == "RA") RA else RA_star
  
  # Main computation loop
  for (i in 1:n_values) {
    design_result <- suppressMessages(
      fully_adaptive_bivariate(
        lambdaA = lambda_A,
        thetaA = theta_A,
        lambdaB = lambda_B_seq[i],
        thetaB = theta_B_seq[i],
        RA_fun = alloc_func,
        sample_size_interim = sample_size_interim,
        alpha1 = alpha1, beta1 = beta1,
        alpha2 = alpha2, beta2 = beta2,
        alpha3 = alpha3,
        num_simulations = num_simulations,
        statistic = statistic
      )
    )
    
    # Store results
    results$ASN_H0[i] <- design_result$ASN$H0
    results$ASN_H1[i] <- design_result$ASN$H1
    results$RA[i] <- design_result$ASN_table$RA_mean[3]
    results$Prop_A[i] <- design_result$ASN_table$PropA_mean[3]
    results$Prop_A_SE[i] <- design_result$ASN_table$PropA_se[3]
  }
  
  return(results)
}

# --------------------------
# Example Usage
# --------------------------

# Run with reduced parameters for demonstration
adaptive_results <- compute_adaptive_table(
  lambda_B_start = 2,
  lambda_B_end = 4,
  theta_B_start = 2.5,
  theta_B_end = 4.5,
  step_size = 0.5,
  num_simulations = 1000,
  allocation_rule = "RA"
)

# View results
print(adaptive_results)
