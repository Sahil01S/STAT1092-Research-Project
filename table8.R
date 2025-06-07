
RA_star <- function(lambdaA, lambdaB, n_iter = 50) {
  X_A <- rpois(n_iter, lambda = lambdaA)
  X_B <- rpois(n_iter, lambda = lambdaB)
  compar_mat_1 <- outer(X_A, X_B, FUN = "<")
  compar_mat_2 <- outer(X_A, X_B, FUN = "==")
  prob_quan_1 <- sum(compar_mat_1) / (n_iter^2)
  prob_quan_2 <- sum(compar_mat_2) / (n_iter^2)
  prob_quan_1 + (0.5 * prob_quan_2)
}
lambda_A <- 2
lambda_B_values <- seq(2, 6, 0.5)
results_table2 <- lapply(lambda_B_values, function(lambda_B) {
  result <- fully_adaptive_group_seq(lambdaA = lambda_A, 
            lambdaB = lambda_B, num_simulations = 1000, RA_fun = RA_star)
  data.frame(
    lambda_B = lambda_B,
    ASN_H0 = result$ASN$H0,
    ASN_H1 = result$ASN$H1,
    RA_Star = result$ASN_table$RA_values[3],
    Prop_A = result$ASN_table$Proportions_A[3],
    Prop_A_SE = result$ASN_table$PropA_se[3]  
  )
}) %>% bind_rows()
print(results_table2)
