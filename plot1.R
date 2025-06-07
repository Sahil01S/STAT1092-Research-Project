plot_titles <- list(
  expression("Histogram of" ~ X),
  expression("Histogram of" ~ sqrt(X)),
  expression("Histogram of" ~ X^{2/3}),
  expression("Histogram of" ~ log(X)),
  expression("Histogram of" ~ log(sqrt(X)))
)
par(mfrow = c(3,2))
for (i in seq_along(all_transformations)) {
  x <- all_transformations[[i]]
  hist(x, prob = TRUE,
       main = plot_titles[[i]],
       xlab = "x values",
       col = "skyblue", border = "black")
  curve(dnorm(x, mean = mean(x), sd = sd(x)),
        col = "red", lwd = 2, add = TRUE)   
}
par(mfrow = c(1,1))