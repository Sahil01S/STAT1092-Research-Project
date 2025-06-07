resultT2 <- compute_ASN_table(
  lambda_B_start = 2, lambda_B_end = 4,
  theta_B_start = 2.5, theta_B_end = 4.5,
  step_size = 0.2,
  statistic = "T2"  # Can be "T1", "T2", or "T3"
)

resultT2 #max