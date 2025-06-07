
choices <- list(
  CHOICE1 = list(alpha = c(0.02, 0.02, 0.01), beta = c(0.015, 0.015, 0.02)),
  CHOICE2 = list(alpha = c(0.03, 0.01, 0.01), beta = c(0.02, 0.02, 0.01)),
  CHOICE3 = list(alpha = c(0.015, 0.025, 0.01), beta = c(0.03, 0.01, 0.01)),
  CHOICE4 = list(alpha = c(0.015, 0.015, 0.02), beta = c(0.025, 0.01, 0.015)),
  CHOICE5 = list(alpha = c(0.03, 0.015, 0.005), beta = c(0.02, 0.025, 0.005))
)
lambda_B_seq <- seq(2, 4, by = 0.2)
all_results <- list()
for (choice_name in names(choices)) {
  choice <- choices[[choice_name]]
  result <- compute_ASN_table(
    lambda_A = 2,
    lambda_B_seq = lambda_B_seq,
    alpha1 = choice$alpha[1],
    beta1 = choice$beta[1],
    alpha2 = choice$alpha[2],
    beta2 = choice$beta[2],
    alpha3 = choice$alpha[3]
  )
  all_results[[choice_name]] <- result
}
colors <- c("black", "red", "blue", "green", "purple", "orange")
plot(NA, 
     xlim = range(lambda_B_seq), 
     ylim = c(100, max(sapply(all_results, function(x) max(x$ASN_H1)))),
     xlab = expression(lambda[B]), 
     ylab = "ASN under H1", 
     main = "ASN Comparison Across Different Choices",
     type = "n")  # 'n' means no plotting
grid()
for (i in seq_along(all_results)) {
  lines(all_results[[i]]$lambda_B, all_results[[i]]$ASN_H1, 
        col = colors[i], lwd = 2)
}
legend("topright", 
       legend = names(choices), 
       col = colors,  
       lwd = 2, 
       cex = 0.8)
