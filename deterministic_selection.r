deterministic_selection <- function(N, s, h, epsilon)
{
	timecount <- 0
	p <- 1 / (2*N)
	fitness <- c(1-s, 1-h*s, 1)
	marg_fitness <- calc_marginal_fitness(fitness, p)

	gen <- c(timecount)
	frequency <- c(p)
	popfit <- c(calc_population_fitness(p, marg_fitness))

	while (p <= 1-epsilon) {
		gen <- c(gen, timecount <- timecount + 1)
		marg_fitness <- calc_marginal_fitness(fitness, p)
		wbar <- calc_population_fitness(p, marg_fitness)
		frequency <- c(frequency, p <- p * ( marg_fitness[2] / wbar ))
		popfit <- c(popfit, wbar)
	}

	d <- data.frame("gen" = gen, "frequency" =  frequency, "popfit" = popfit)
	return(d)
}

calc_marginal_fitness <- function(fitness, p)
{
	q <- 1-p
	mfit <-  c(q*fitness[1] + p*fitness[2], p*fitness[3] + q*fitness[2])
	return(mfit)
}

calc_population_fitness <- function(p, mfit)
{
	q <- 1-p
	wbar <- q*mfit[1] + p*mfit[2]
}
