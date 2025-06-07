
choices <- list(
  CHOICE1 = list(alpha = c(0.02, 0.02, 0.01),
                beta = c(0.015, 0.015, 0.02)),
  CHOICE2 = list(alpha = c(0.03, 0.01, 0.01),
                beta = c(0.02, 0.02, 0.01)),
  CHOICE3 = list(alpha = c(0.015, 0.025, 0.01),
                beta = c(0.03, 0.01, 0.01)),
  CHOICE4 = list(alpha = c(0.015, 0.015, 0.02),
                beta = c(0.025, 0.01, 0.015)),
  CHOICE5 = list(alpha = c(0.03, 0.015, 0.005),
                beta = c(0.02, 0.025, 0.005))
)

# Parameter sequences
lambda_B_seq <- seq(2, 4, by = 0.2)
theta_B_seq <- seq(2.5, 4.5, by = 0.2)  # Matching length with lambda_B_seq

# Initialize results storage
all_results <- list()

# Run simulations for each choice
for (choice_name in names(choices)) {
  choice <- choices[[choice_name]]
  
  result <- compute_ASN_table(
    lambda_A = 2,
    theta_A = 2.5,
    lambda_B_start = 2,
    lambda_B_end = 4,
    theta_B_start = 2.5,
    theta_B_end = 4.5,
    step_size = 0.2,
    alpha1 = choice$alpha[1],
    beta1 = choice$beta[1],
    alpha2 = choice$alpha[2],
    beta2 = choice$beta[2],
    alpha3 = choice$alpha[3],
    num_simulations = 1000,
    statistic = "T1"  # Can change to "T2" or "T3" as needed
  )
  
  all_results[[choice_name]] <- result
}

colors <- c("black", "red", "blue", "green", "purple", "orange")

plot(NA,
     xlim = range(lambda_B_seq),
     ylim = c(100, max(sapply(all_results, function(x) max(x$ASN_H0)))),
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "ASN under H0",
     main = "ASN under H0 across Different Alpha/Beta Choices",
     type = "n")
grid()

# Loop to add each line
for (i in seq_along(all_results)) {
  lines(all_results[[i]]$lambda_B,
        all_results[[i]]$ASN_H0,
        col = colors[i], lwd = 2)
}

legend("topright", legend = names(choices), col = colors, lwd = 2, cex = 0.8)


