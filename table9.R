cat("Comparison of ASN (Under H1) for different methods\n")
comp <- cbind("Lambda_B" = result_nonada,1],
              "Non Adaptive" = result_nonada[,3],
              "Entry Adaptive" = result_es[,3], 
              "Adaptive (RA)" = results_table[,3],
              "Adaptive (RA_star)" = results_table2[,3])
comp