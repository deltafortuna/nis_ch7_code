equilibrium_gridsearch <- function(r, sA, sB, tA, tB, eAB, eAb, eaB, eab)
{
	fit <- matrix(nrow=3, ncol=3, c(1-sA-sB-eAB, 1-sA, 1-sA-tB-eaB,  # col 1
									1-sB, 1, 1-tB,					 # col 2
									1-tA-sB-eAb, 1-tA, 1-tA-tB-eab));# col 3
	print(fit);

	anchor <- list(a1 = c(1,1), a2 = c(2,1), a3 = c(1,2), a4 = c(2,2)); # upper-left cell of each submatrix
	w <- vector();
	start <- c(0.01, 0.01, 0.01, 0.01);
	end <- c(0.99, 0.99, 0.99, 0.99);

	for (it in c(0.01)) {    #, 0.001, 0.0001, 0.00001)) {
		d <- list();
		combos <- data.frame(cbind(seq(start[1], end[1], it), seq(start[2], end[2], it), seq(start[3], end[3], it), seq(start[4], end[4], it)) );
		combos <- expand.grid(combos);
		combos$sum <- rowSums(combos);
		hgrid <- combos[combos$sum ==1,1:4];

		gridlist <- split(hgrid, seq(nrow(hgrid)));
		for (it2 in 1:nrow(hgrid)) {

			h <- unname(unlist(gridlist[it2])); # convert gridlist row to numeric vector w/o names

			# calc marginal fit of each haplotype
			for (i in 1:4) {
				subfit <- fit[ anchor[[i]][1]:(anchor[[i]][1]+1), anchor[[i]][2]:(anchor[[i]][2]+1)];
				subfit <- c(subfit[,1], subfit[,2]);  # don't need to recalc subfit every time # can calc before entering nested loops
				w[i] <- sum(h * subfit);
			}

			# calc strength of epistasis
			E <- w[1] + w[4] - w[2] - w[3];

			# calc mean population fitness and D
			wbar <- sum(h*w);
			D <- h[1] - (h[1]+h[2]) * (h[1]+h[3]);

			# calc new frequencies of haplotypes
			nexth <- vector();
			for (i in 1:4) {
				if (i ==2 | i == 3) {
					nexth[i] <- ( (h[i]*w[i]) + (r*D*fit[2,2]) ) / wbar;
				} else {
					nexth[i] <- ( (h[i]*w[i]) - (r*D*fit[2,2]) ) / wbar;
				}
			}
			d[[it2]] <- c(h[1], h[2], h[3], h[4], nexth[1]-h[1], nexth[2]-h[2], nexth[3]-h[3], nexth[4]-h[4], h[1]+h[2], h[1]+h[3], wbar, D, E);
		}
		df <- data.frame(matrix(unlist(d), nrow=length(d), byrow=T),stringsAsFactors=FALSE);
		names(df) <- c("h1", "h2", "h3", "h4", "dh1", "dh2", "dh3", "dh4", "p1", "p2", "wbar", "D", "E");
		df$totalchange <- abs(df$dh1) + abs(df$dh2) + abs(df$dh3) + abs(df$dh4);
		df <- df[order(df$totalchange),];
		print (df[1:10,]);

		# reset start and end vectors for higher resolution search grid
		for (i in 1:4) {
			start[i] <- mean(df[1:10,i]) - (5*it);
			end[i] <- mean(df[1:10,i]) + (5*it);
		}
	}
}
