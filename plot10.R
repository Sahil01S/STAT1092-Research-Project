
plot(NA,
     xlim = range(lambda_B_seq),
     ylim = c(100, max(sapply(phi_results, function(x) max(x$ASN_H1)))),
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "ASN under H1",
     main = "ASN under H1 for Different phi Values",
     type = "n")
grid()

for (i in seq_along(phi_seq)) {
  lines(phi_results[[i]]$lambda_B,
        phi_results[[i]]$ASN_H1,
        col = colors[i], lwd = 2)
}

legend("topright",
       legend = parse(text = paste("phi == ", phi_seq)),
       col = colors, lwd = 2, cex = 0.7, ncol = 2)

