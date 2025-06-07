comp <- cbind("Lambda_B" = result1[,1],
              "Theta_B" = result1[,2], 
              "Non Ada" = result1[,4],
              "Entry Ada" = result_ada1[,4], 
              "RA" = adaptive_results[,4],
              "RA_star" = adaptive_results2[,4])
cat("Comparison of ASN (Under H1) for different methods :\n")
comp