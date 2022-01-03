deterministic_twolocus_selection <- function(r, sA, sB, tA, tB, eAB, eAb, eaB, eab, h1, h2, h3, h4, gens)
{
  fit <- matrix(nrow=3, ncol=3, c(1-sA-sB-eAB, 1-sA, 1-sA-tB-eaB,  # col 1
                                  1-sB, 1, 1-tB,					 # col 2
                                  1-tA-sB-eAb, 1-tA, 1-tA-tB-eab)) # col 3
  print(fit);
  
  # initial haplotype frequencies
  h <- c(h1, h2, h3, h4) # f_A-B, f_A-b, f_a-B, f_a-b
  
  w <- vector(); #vector for marginal fitnesses
  anchor <- list(a1 = c(1,1), a2 = c(2,1), a3 = c(1,2), a4 = c(2,2)) # upper-left cell of each submatrix
  d <- setNames(data.frame(matrix(ncol = 9)), c("gen", "h1", "h2", "h3", "h4", "pA", "pB", "wbar", "D")) # initialize data matrix
  
  for (gen in 1:gens) {
    # calc marginal fitness of each haplotype
    for (i in 1:4) { 
      subfit <- fit[ anchor[[i]][1]:(anchor[[i]][1]+1), anchor[[i]][2]:(anchor[[i]][2]+1)]
      subfit <- c(subfit[,1], subfit[,2])
      print(subfit)
      w[i] <- sum(h * subfit)
    }
    
    # calc mean population fitness and D
    wbar <- sum(h*w)
    D <- h[1] - (h[1]+h[2]) * (h[1]+h[3])
    
    # calc new frequencies of haplotypes
    for (i in 1:4) {
      if (i ==2 | i == 3) {
        h[i] <- ( (h[i]*w[i]) + (r*D*fit[2,2]) ) / wbar
      } else {
        h[i] <- ( (h[i]*w[i]) - (r*D*fit[2,2]) ) / wbar
      }
    }
    d[gen,] = c(gen, h[1], h[2], h[3], h[4], h[1]+h[2], h[1]+h[3], wbar, D)
  }
  return(list("data" = d, "fitness_matrix" = fit))
}