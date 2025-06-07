

phi_seq <- seq(0.2, 2, by = 0.2)
lambda_B_seq <- seq(2, 4, by = 0.2)
theta_B_seq <- seq(2.5, 4.5, by = 0.2)
alpha_choice <- c(0.03, 0.015, 0.005)
beta_choice  <- c(0.02, 0.025, 0.005)
phi_results <- list()

for (phi_val in phi_seq) {
  result <- compute_ASN_table(
    lambda_A = 2,
    theta_A = 2.5,
    lambda_B_start = 2,
    lambda_B_end = 4,
    theta_B_start = 2.5,
    theta_B_end = 4.5,
    step_size = 0.2,
    phi_A = phi_val,
    phi_B = phi_val,
    alpha1 = alpha_choice[1],
    beta1 = beta_choice[1],
    alpha2 = alpha_choice[2],
    beta2 = beta_choice[2],
    alpha3 = alpha_choice[3],
    num_simulations = 1000,
    statistic = "T3"
  )
  
  phi_results[[paste0("phi_", phi_val)]] <- result
}
colors <- rainbow(length(phi_seq))
# --- ASN under H0 ---
plot(NA,
     xlim = range(lambda_B_seq),
     ylim = c(100, max(sapply(phi_results, function(x) max(x$ASN_H0)))),
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "ASN under H0",
     main = "ASN under H0 for Different phi Values",
     type = "n")
grid()

for (i in seq_along(phi_seq)) {
  lines(phi_results[[i]]$lambda_B,
        phi_results[[i]]$ASN_H0,
        col = colors[i], lwd = 2)
}

legend("topright",
       legend = parse(text = paste("phi == ", phi_seq)),
       col = colors, lwd = 2, cex = 0.7, ncol = 2)

