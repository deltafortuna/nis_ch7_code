analyze_selective_regime <- function(s,h,relfitAA, relfitAa, relfitaa)
{
	p <- seq(0,1,0.001)
	q <- 1-p
	relfit <- list()

	relfit[[1]] <- eval(substitute(relfitAA))
	relfit[[2]] <- eval(substitute(relfitAa))
	relfit[[3]] <- eval(substitute(relfitaa))

	wA <- p*relfit[[1]] + q*relfit[[2]]
	wa <- p*relfit[[2]] + q*relfit[[3]]
	wbar <- p*wA + q*wa

	pt <- p * wA / wbar
	change <- pt-p
	vf <- data.frame("p" = p, "pdot" = change)
	popfit <- data.frame("p" = p, "wbar" = wbar)
	q <- list()
	q[[1]] <- vf
	q[[2]] <- popfit
	return(q)
}

plot_vector_field <- function(q)
{
	plot(q[[1]]$p, q[[1]]$pdot, type = "l", xlab = "p", ylab = "dp/dt")
	abline(h = 0., lty = 2)
}

plot_wbar <- function(q)
{
	plot(q[[2]]$p, q[[2]]$wbar, type = "l", xlab = "p", ylab = "wbar")
	m <- q[[2]][q[[2]]$wbar == max(q[[2]]$wbar),][1]
	print(m)
	abline(v = m, lty = 2)
}
