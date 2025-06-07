

plot(adaptive_results2$lambda_B,
     adaptive_results2$Prop_A,
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "Allocation for A",
     main = "Allocation Comparison for RA and RA_star",
     type = 'l', col = "orange", lwd = 1.5)
lines(adaptive_results$lambda_B,
      adaptive_results$Prop_A,
      col = 'purple',lty = 1, lwd = 1.5)
legend('bottomright', col = c('orange', 'purple'),
       lwd = c(1.5,1.5), cex = 0.8,
       legend = c("RA_Star","RA"))

