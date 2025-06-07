
compile_ASN_tables <- function(     
                      lambda_A = 2, 
                      lambda_B_seq = seq(2, 4, by = 0.2),
                      mA1 = 50, mB1 = 50, mA2 = 50, 
                      mB2 = 50, mA3 = 50, mB3 = 50,
                      alpha1 = 0.02, beta1 = 0.015, 
                      alpha2 = 0.02, beta2 = 0.015, 
                      alpha3 = 0.01,
                      num_simulations = 10000) {
  non_adaptive_results <- compute_ASN_table(
    lambda_A = lambda_A, 
    lambda_B_seq = lambda_B_seq,
    mA1 = mA1, mB1 = mB1, mA2 = mA2,
    mB2 = mB2, mA3 = mA3, mB3 = mB3,
    alpha1 = alpha1, beta1 = beta1, 
    alpha2 = alpha2, beta2 = beta2, 
    alpha3 = alpha3,
    num_simulations = num_simulations
  )
  adaptive_results <- compute_ASN_table_ad(
    lambda_A = lambda_A, 
    lambda_B_seq = lambda_B_seq,
    mA1 = mA1, mB1 = mB1, mA2 = mA2, 
    mB2 = mB2, mA3 = mA3, mB3 = mB3,
    alpha1 = alpha1, beta1 = beta1, 
    alpha2 = alpha2, beta2 = beta2,
    alpha3 = alpha3,
    num_simulations = num_simulations
  )
  combined_results <- data.frame(
    lambda_B = lambda_B_seq,
    NonAdaptive_ASN_H0 = non_adaptive_results$ASN_H0,
    NonAdaptive_ASN_H1 = non_adaptive_results$ASN_H1,
    Adaptive_ASN_H0 = adaptive_results$ASN_H0,
    Adaptive_ASN_H1 = adaptive_results$ASN_H1
  )
  return(combined_results)
}
combined_table <- compile_ASN_tables()
plot(combined_table$lambda_B, 
     combined_table$NonAdaptive_ASN_H1, 
     type = 'b', col = 'blue', pch = 16,
     xlab = expression(lambda[B]), ylab = "ASN_H1",
     ylim = c(100,300),
     main = "Comparison between GSD and Adaptive GSD")
lines(combined_table$lambda_B, combined_table$Adaptive_ASN_H1, 
      type = 'b', col = 'red', pch = 16)
legend("topright",legend = c("Non Adaptive", "Adaptive"),
       fill = c("blue","red"),
       cex = 0.7)
