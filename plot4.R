plot(results_table2$lambda_B,
     results_table2$Prop_A,
     xlab = expression(lambda[B]), ylab = "Prop_A",
     main = "Allocation comparison for RA and RA_star",
     type = 'l', lwd = 2, col ='orange')
lines(results_table$lambda_B,
      results_table$Prop_A,
      col = "purple", lwd = 2)
legend('topleft', col = c('orange','purple'), cex = 0.8,
       lwd = c(2,2), legend = c("RA_star","RA"))