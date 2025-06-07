<<echo=F,fig.height=5,fig.show='hold', comment=NA>>=
H0table <- data.frame(
  lambda_B = resultT1$lambda_B,
  theta_B  = resultT1$theta_B,
  ASN_T1   = resultT1$ASN_H0,
  ASN_T2   = resultT2$ASN_H0,
  ASN_T3   = resultT3$ASN_H0
)


plot(H0table$lambda_B, H0table$ASN_T1,
     type = "l",
     col = "orange", lwd = 2,
     ylim = c(100, 300),
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "Average Sample Number under H0",
     main = "Comparison of ASN under H0 for Different Statistics")

# Add lines for T2 and T3 statistics
lines(H0table$lambda_B, H0table$ASN_T2, col = "purple", lwd = 2)
lines(H0table$lambda_B, H0table$ASN_T3, col = "red", lwd = 2)

# Add a legend
legend("topright",
       legend = c(
         expression(T[1] * ":" ~ "min(" * Z[1] * "," ~ Z[2] * ")"),
         expression(T[2] * ":" ~ "max(" * Z[1] * "," ~ Z[2] * ")"),
         expression(T[3] * ":" ~ "avg(" * T[1] * "," ~ T[2] * ")")
       ),
       col = c("orange", "purple", "red"),
       lwd = 2,
       cex = 0.85)


