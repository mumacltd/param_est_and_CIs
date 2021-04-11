# prerequisites
library(Dark)
# a <- p0 <- dark$opt[1:5]
# ct <- p0[1]
# cc <- p0[2]
# tau <- p0[3]
# S2 <- p0[4]
# alpha <- p0[5]
# idx <- dark$time < dark$opt[7]
# x <- dark$time[idx]
# y <- dark$thrs[idx]
# plot(x,y)
fn <- P5c

# optim(p0, fn)
### step function 
H <- function(x, t) {
	round(1/(1 + exp(-20 * (x - t))), 0)
}


### exponential model
P3 <- function(a, x) {
	a[1] + a[2] * exp(-x/a[3])
}

### Abbreviated MLP
P5c <- function(a, x) {
	a[1] + a[2] * exp(-x/a[3]) + a[4] * (x - a[5]) * H(x, a[5])
}

### Full MLP
P7c <- function() {
	a[1] + a[2] * exp(-x/a[3]) + a[4] * (x - a[5]) * H(x, a[5]) + +a[6] * (x - a[7]) * H(x, a[7])
}
# Residuals
Resid <- function(a, fn, x, y) {
	y - fn(a, x)
}

# LeastSqs
LeastSq <- function(a, fn = P5c, x, y) {
	Resid(a, fn, x, y) %*% Resid(a, fn, x, y)
}

# Partial Differential Equations
pdf1 <- function(a, fn, x, y) {
	root <- x/x
	-2 * (root %*% Resid(a, fn, x, y))
}
# pdf1(a,P5c,x,y)

pdf2 <- function(a, fn, x, y) {
	root <- exp(-x/a[3])
	-2 * (root %*% Resid(a, fn, x, y))
}
# pdf2(a,P5c,x,y)

pdf3 <- function(a, fn, x, y) {
	root <- exp(-x/a[3]) * a[2] * x/(a[3] * a[3])

	-2 * (root %*% Resid(a, fn, x, y))
}
# pdf3(a,P5c,x,y)

pdf4 <- function(a, fn, x, y) {
	root <- (x - a[5])/(1 + exp(-200 * (x - a[5])))
	-2 * (root %*% Resid(a, fn, x, y))
}
# pdf4(a,P5c,x,y)

pdf5 <- function(a, fn, x, y) {

	k <- 20
	root = (a[4]/(1 + exp(-2 * k * (x - a[5]))) - a[4] * (x - a[5]) * (exp(-2 * k * (x - a[5])) * (-2 * k))/(1 + exp(-2 * k * (x - a[5])))^2)

	2 * (root %*% Resid(a, fn, x, y))
}

# pdf5(a,P5c,x,y)


# Gradient Vector
Grad <- function(a, fn, x, y) {
	ct <- a[1]
	cc <- a[2]
	tau <- a[3]
	S2 <- a[4]
	alpha <- a[5]
	k = 20
	L <- expression((y - (ct + cc * exp(-x/tau) + S2 * (x - alpha)/(1 + exp(-2 * k * (x - alpha)))))^2)
	Params <- c("ct", "cc", "tau", "S2", "alpha")
	pGrad <- c()
	for (ii in 1:5) {
		pGrad <- cbind(pGrad, sum(eval(D(L, Params[ii])), na.rm = TRUE))

	}
	#  c(pdf1(a, fn, x, y), pdf2(a, fn, x, y), pdf3(a, fn, x, y), pdf4(a, fn, x, y), pdf5(a, fn, x, y))
	pGrad
}
# Grad(a,P5c,x,y)

progGrad <- function(p0, fn = P5c, x, y) {
	ct <- p0[1]
	cc <- p0[2]
	tau <- p0[3]
	S2 <- p0[4]
	alpha <- p0[5]
	k = 10

	L <- expression((y - (ct + cc * exp(-x/tau) + S2 * (x - alpha)/(1 + exp(-2 * k * (x - alpha)))))^2)


	Params <- c("ct", "cc", "tau", "S2", "alpha")

	pGrad <- c()
	for (ii in 1:5) {
		pGrad <- cbind(pGrad, (eval(D(L, Params[ii]))))

	}
	idx <- is.na(pGrad)
	pGrad[idx] <- 0

	GradMat <- t(pGrad) %*% pGrad
	GradMat
}
# progGrad(p0, P5c, x, y)


# Batch Gradient Descent
findParametersGDPlain <- function(a, fn, x, y) {
	val = list()
	tol <- 1e-06
	delta <- NULL
	eta <- 1/c(1000, 1000, 1000, 10000, 100)
	LS0 <- 1e+07
	for (ii in 1:10000) {
		gradient <- Grad(a, fn, x, y)
		a <- a - eta * gradient
		LS1 <- LeastSq(a, fn, x, y)
		del = LS0 - LS1
		LS0 = LS1
		delta <- c(delta, del)
		if (abs(del) < tol) 
			break
	}
	resV <- Resid(a, fn, x, y)
	val$params = a
	val$runs = ii
	val$noise = sd(resV)
	val$MSE = resV %*% resV/(length(x) - 5)
	return(val)
}
# findParametersGDPlain(a, fn, x, y)


## GD with Momentum

findParametersGDMomentum <- function(a, fn, x, y) {
	val = list()
	tol <- 1e-06
	# delta <- NULL
	eta <- 1/c(1000, 1000, 1000, 10000, 100)
	LS0 <- 1e+07
	velocity = 0
	gamma = 0.9
	for (ii in 1:10000) {
		gradient <- Grad(a, fn, x, y)
		velocity = gamma * velocity + eta * gradient
		a <- a - velocity
		LS1 <- LeastSq(a, fn, x, y)
		del = LS0 - LS1
		LS0 = LS1
		# delta <- c(delta, del)
		if (abs(del) < tol) 
			break
	}
	resV <- Resid(a, fn, x, y)
	val$params = a
	val$runs = ii
	val$noise = sd(resV)
	val$MSE = resV %*% resV/(length(x) - 5)
	return(val)
}
# findParametersGDMomentum(a, fn, x, y)

## GD with Momentum Bounded

findParametersGDMomentumBounded <- function(a, fn, x, y) {
	val = list()
	tol <- 1e-06
	# delta <- NULL
	eta <- 1/c(1000, 1000, 1000, 10000, 100)
	LS0 <- 1e+07
	velocity = 0
	gamma = 0.9
	for (ii in 1:10000) {
		gradient <- Grad(a, fn, x, y)
		velocity = gamma * velocity + eta * gradient
		a <- a - velocity
		if(a[5]<0 | a[5]>max(x)){
			a[5]<-0
			a[4]<-0
		}
		LS1 <- LeastSq(a, fn, x, y)
		del = LS0 - LS1
		LS0 = LS1
		# delta <- c(delta, del)
		if (abs(del) < tol) 
			break
	}
	resV <- Resid(a, fn, x, y)
	val$params = a
	val$runs = ii
	val$noise = sd(resV)
	val$MSE = resV %*% resV/(length(x) - 5)
	return(val)
}
# findParametersGDMomentumBounded(a, fn, x, y)

## Nesterov Accelerated GD
findParametersGDNAG <- function(a, fn, x, y) {
	val = list()
	tol <- 1e-06
	# delta <- NULL
	eta <- 1/c(1000, 1000, 1000, 10000, 10)
	LS0 <- 1e+07
	velocity = 0
	gamma = 0.9
	for (ii in 1:10000) {
		gradient <- Grad(a - gamma * velocity, fn, x, y)
		velocity = gamma * velocity + eta * gradient
		a <- a - velocity
		LS1 <- LeastSq(a, fn, x, y)
		del = LS0 - LS1
		LS0 = LS1
		# delta <- c(delta, del)
		if (abs(del) < tol) 
			break
	}
	resV <- Resid(a, fn, x, y)
	val$params = a
	val$runs = ii
	val$noise = sd(resV)
	val$MSE = resV %*% resV/(length(x) - 5)
	return(val)
}
# findParametersGDNAG(a, fn, x, y) 

Hessian <- function(a,fn,x,y) {
	ct <- a[1]
	cc <- a[2]
	tau <- a[3]
	S2 <- a[4]
	alpha <- a[5]
	k = 20
	time <- x
	threshold <- y
	L <- expression((threshold - (ct + cc * exp(-time/tau) + S2 * (time - alpha)/(1 + exp(-20 * (time - alpha)))))^2)
	Params <- c("ct", "cc", "tau", "S2", "alpha")

	Hess = matrix(0, 5, 5)
	Hess[1, 1] <- 2*length(x)
	for (ii in 2:5) {
		for (jj in 2:5) {
			Hess[ii, jj] = sum(eval(D(D(L, Params[ii]), Params[jj])))
		}
	}
	Hess
}

# Hessian(a,fn,x,y)
