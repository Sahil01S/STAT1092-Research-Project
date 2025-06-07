
plot(result1$lambda_B, result1$ASN_H0,
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "ASN Under H0",
     ylim = c(100, 300),
     main = "ASN Comparison between non adaptive and adaptive",
     type = 'l', col = "red", lwd = 1.5)
lines(result1$lambda_B, result_ada1$ASN_H0,
      col = 'blue', lwd = 1.5)
lines(result1$lambda_B, result1$ASN_H1,
      col = 'red', lwd = 1.5, lty=2)
lines(result1$lambda_B, result_ada1$ASN_H1,
      col = 'blue', lwd = 1.5, lty=2)

legend('topright', col = c('red', 'blue','red', 'blue'),
       lwd = c(1.5,1.5,1.5,1.5), cex = 1,
       lty = c(1,1,2,2),
       legend = c("Non adaptive (H0)","Entry adaptive (H0)",
                "Non adaptive (H1)","Entry adaptive (H1)"))