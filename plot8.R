
plot(NA,
     xlim = range(lambda_B_seq),
     ylim = c(100, max(sapply(all_results, function(x) max(x$ASN_H1)))),
     xlab = expression(lambda[B] ~ "(with increasing " ~ theta[B] ~ ")"),
     ylab = "ASN under H1",
     main = "ASN under H1 across Different Alpha/Beta Choices",
     type = "n")
grid()

# Loop to add each line
for (i in seq_along(all_results)) {
  lines(all_results[[i]]$lambda_B,
        all_results[[i]]$ASN_H1,
        col = colors[i], lwd = 2)
}

legend("topright", legend = names(choices), col = colors, lwd = 2, cex = 0.8)

