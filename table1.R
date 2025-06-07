set.seed(100)
OrgSample <- rpois(n = 500, lambda = 15)
all_transformations <- data.frame(
  "OrgSample" = OrgSample,
  "sqrt" = sqrt(OrgSample),
  "two_thirds" = OrgSample^(2/3),
  "log" = log(OrgSample),
  "log_of_sqrt" = log(sqrt(OrgSample))
)
KSTest_p.value <- suppressWarnings(
  apply(all_transformations, 
        2, 
        function(x) ks.test(x, 'pnorm', 
        mean = mean(x), sd = sd(x))$p.value)
)
KSTest_table <- data.frame(
  Transformation = names(KSTest_p.value),
  P_Value = KSTest_p.value
)
print(format(KSTest_table, scientific=FALSE), 
row.names = FALSE)