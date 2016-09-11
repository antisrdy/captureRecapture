# Data generation
f.generate <- function(nobs, alpha = 1, beta = 3) {
  ## Generation of X_1,...,X_n following a Beta(alpha, beta)
  data <- rbeta(nobs, alpha, beta)
  return(data)
}

# Method of moments
f.methodOfMoments <- function(data) {
  ## Empirical mean and variance
  xBar <- mean(data)
  v <- var(data)
  ## alpha_hat and beta_hat estimates
  alpha_hat <- xBar*(xBar*(1-xBar)/v-1)
  beta_hat <- (1-xBar)*(xBar*(1-xBar)/v-1)
  ## Result
  return(c(alpha_hat, beta_hat))
}

# Newton-Raphson algorithm
f.newtonRaphson <- function(data, precisionThreshold = 0.0001) {
  ## Initialization: theta_hat from method of moments
  ## Vector of precision to compare with, at each iteration
  ## (termination if the difference between 2 iterations is under precisionThreshold)
  theta_hat = f.methodOfMoments(data)
  precision <- c(10000, 10000) # Initialization of alpha and beta precision to huge values
  ## Iterations
  while (abs(precision[1] > precisionThreshold) | precision[2] > precisionThreshold) {
    alpha_hat <- theta_hat[1]
    beta_hat <- theta_hat[2]
    ### Let's compute the gradient of the likelihood function and its hessian matrix
    ### using digamma and trigamma functions for simpler calculus
    ### We get rid of n (number of observations) as they go when we take the product of gradient^(-1) with hessian
    gradient <- matrix(c(digamma(alpha_hat+beta_hat) - digamma(alpha_hat) + mean(log(data)),
                         digamma(alpha_hat+beta_hat) - digamma(beta_hat) + mean(log(1-data))),
                       2,1)
    hessian <- solve(matrix(c(trigamma(alpha_hat+beta_hat) - trigamma(alpha_hat), # solve = inverse matrix
                        trigamma(alpha_hat+beta_hat),
                        trigamma(alpha_hat+beta_hat),
                        trigamma(alpha_hat+beta_hat) - trigamma(beta_hat)),
                      2,2))
    ### Iteration updates
    precision <- hessian %*% gradient
    theta_hat <- theta_hat - precision
  }
  return(theta_hat)
}

# Comparison of the maximum likelihood and the method of moments
# Alpha and Beta properties are the same. Let's focus on Alpha
## How does the sample size play a role?
## Chart 1
### Retrieve estimates for different sample sizes
nobs <- seq(100, 100000, by=50)
alpha_methodOfMoments <- vector("numeric", length(nobs))
alpha_newtonRaphson <- vector("numeric", length(nobs))
for (i in seq(1:length(nobs))) {
  data <- f.generate(nobs = nobs[i])
  alpha_methodOfMoments[i] <- f.methodOfMoments(data)[1]
  alpha_newtonRaphson[i] <- f.newtonRaphson(data)[1]
}
### Plot
plot.new()
par(mfrow=c(1,1))
plot(nobs,
     alpha_newtonRaphson,
     type="p",
     xlim=range(nobs),
     ylim=range(c(alpha_newtonRaphson,alpha_methodOfMoments)),
     pch=19, cex = 0.1, axes = TRUE, col="red ", xlab = "sample size", ylab = "alpha estimation")
par(new = TRUE)
abline(h = 1, col = "black")
par(new = TRUE)
plot(nobs,
     alpha_methodOfMoments,
     type="p",
     xlim=range(nobs),
     ylim=range(c(alpha_newtonRaphson,alpha_methodOfMoments)),
     pch=19, cex = 0.1, axes = TRUE, col="green ", xlab = "", ylab = "")
legend("topleft", legend = c("Maximum likelihood", "Method of moments"), col=c("red","green"),
        pch=c(19,19) , bty="p", pt.cex = 1, cex = 0.8)
dev.off()

## What does the alpha distribution look like?
## Charts 2&3
f.drawHistogram <- function(nobs, nSamples = 3000) {
  ### Retrieve estimates
  alpha_methodOfMoments <- vector("numeric", length(nSamples))
  alpha_newtonRaphson <- vector("numeric", length(nSamples))
  for (i in seq(1:nSamples)) {
    data <- f.generate(nobs)
    alpha_methodOfMoments[i] <- f.methodOfMoments(data)[1]
    alpha_newtonRaphson[i] <- f.newtonRaphson(data)[1]
  }
  ### Plot
  histogramMOM <- hist(alpha_methodOfMoments)
  histogramNR <- hist(alpha_newtonRaphson)
  histogramRange <- range(c(alpha_methodOfMoments, alpha_newtonRaphson))
  plot.new()
  par(mfrow = c(1,1))
  plot(histogramMOM,
       col= "blue",
       xlim= histogramRange,
       xlab = "alpha estimates",
       ylab = "frequency")
  plot(histogramNR,
       col="lightblue",
       xlim=histogramRange,
       add=T)
  legend("topright", c("Method of moments","Maximum likelihood"), col=c("blue","lightblue") , pch=15:15)
}
f.drawHistogram(nobs = 1000)
f.drawHistogram(nobs = 100000)

# Adapated Newton Raphson and Wald test statistics
f.newtonRaphsonAdapted <- function(data, precisionThreshold = 0.0001) {
  ## Newton Raphson subsection
  theta_hat <- f.newtonRaphson(data = data, precisionThreshold = precisionThreshold)
  alpha_hat <- theta_hat[1]
  beta_hat <- theta_hat[2]
  ## Test subsection
  ### Test Statistics
  fisherInformation <- matrix(c(trigamma(alpha_hat + beta_hat) - trigamma(alpha_hat),
                                trigamma(alpha_hat + beta_hat),
                                trigamma(alpha_hat + beta_hat),
                                trigamma(alpha_hat + beta_hat) - trigamma(beta_hat)),
                              2,2)
  varCov <- 1/length(data) * solve(fisherInformation)
  statChiSquare <- (alpha_hat - beta_hat)^2 / (varCov[1,1] + varCov[2,2] - 2*varCov[1,2])
  return(abs(statChiSquare))
}

# Wald Test
f.waldTest <- function(alpha, nobs) {
  data <- f.generate(nobs = nobs, alpha = alpha, beta = alpha)
  return(f.newtonRaphsonAdapted(data = data))
}
qchisq(0.95, df=1, ncp=0, lower.tail=TRUE, log.p=FALSE)
f.waldTest(alpha = 3, 1000000)

# Asymptotic threshold value
## Generate statistics distribution under H0: alpha = beta
### genNumber: how many times do we compute the statistics for each number of observations?
f.statisticsDistribution <- function(nobs, genNumber = 100000) {
  statisticsValues <- vector("numeric", genNumber)
  for (i in seq(1:genNumber)) {
    data <- f.generate(nobs = nobs, alpha = 2, beta = 2)
    statisticsValues[i] <- f.newtonRaphsonAdapted(data = data)
  }
  result <- list("sampleSize" = nobs, "statisticsValues" = statisticsValues)
  return(result)
}

distribution1000 <- f.statisticsDistribution(1000) # Asymptotic distribution of the statistics for 100000 obs
chiSquareDistrib <- rchisq(1000, 1) # True chi square distribution under one df and for 100000 obs

## Plot results for small sizes
## Plot results for a large sample size (1000)
## Plot a true chi square distribution with 1 df
nobs = seq(10, 60, by = 10)
plot.new()
par(mfrow = c(4,2))
for (i in seq(1, length(nobs))) {
  statistics <- f.statisticsDistribution(nobs[i])
  print(paste(c("Quantiles for", nobs[i], "observations"), collapse = " "))
  print(quantile(statistics$statisticsValues))
  hist(statistics$statisticsValues, main = paste(c(nobs[i], "observations"), collapse = " "), xlab = "statistics values distribution")
}
hist(distribution1000$statisticsValues, main = paste(c(1000, "observations"), collapse = " "), xlab = "statistics values")
hist(chiSquareDistrib, main = paste(c(1000, "observations"), collapse = " "), xlab = "chi square distribution")
dev.off()
quantile(distribution1000$statisticsValues)
quantile(chiSquareDistrib)

finalTest <- f.generate(nobs = 20, alpha = 2, beta = 3)
f.newtonRaphsonAdapted(data = finalTest)
qchisq(0.95, df=1, ncp=0, lower.tail=TRUE, log.p=FALSE)