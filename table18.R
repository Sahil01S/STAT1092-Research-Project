


RA_star <- function(lambdaA, lambdaB, thetaA, thetaB, n_iter = 50) {
  X_A <- rpois(n_iter, lambda = lambdaA)
  X_B <- rpois(n_iter, lambda = lambdaB)
  Y_A <- rpois(n_iter, lambda = thetaA)
  Y_B <- rpois(n_iter, lambda = thetaB)
  value <- 0
  for (i in 1:n_iter) {
    if(X_A[i] < X_B[i] && Y_A[i] < Y_B[i]) value <- value + 1
    else if(X_A[i] < X_B[i] && Y_A[i] > Y_B[i]) value <- value + 0.5
    else if(X_A[i] > X_B[i] && Y_A[i] < Y_B[i]) value <- value + 0.5
    else if(X_A[i] == X_B[i] && Y_A[i] < Y_B[i]) value <- value + 0.75
    else if(X_A[i] == X_B[i] && Y_A[i] > Y_B[i]) value <- value + 0.25
    else if(X_A[i] < X_B[i] && Y_A[i] == Y_B[i]) value <- value + 0.75
    else if(X_A[i] > X_B[i] && Y_A[i] == Y_B[i]) value <- value + 0.25
    else if(X_A[i] == X_B[i] && Y_A[i] == Y_B[i]) value <- value + 0.5

  }
  value/n_iter
}

adaptive_results2 <- compute_adaptive_table(
  lambda_B_start = 2,
  lambda_B_end = 4,
  theta_B_start = 2.5,
  theta_B_end = 4.5,
  step_size = 0.2,
  num_simulations = 1000,
  allocation_rule = "RA_star"
)

# View results
print(adaptive_results2)


