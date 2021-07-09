######### Home assignment 3 ######### 
# Philip Dahlqvist-Sjöberg
# 23-04-2020

library(dplyr)
# Write Latex matrices
write_matex <- function(x) {
  begin <- "\\begin{bmatrix}"
  end <- "\\end{bmatrix}"
  X <-
    apply(x, 1, function(x) {
      paste(
        paste(x, collapse = "&"),
        "\\\\"
      )
    })
  writeLines(c(begin, X, end))
}

wd <- "C:/Users/Eva/Downloads/Statistical Computation/HW3/"
######### 1 #########

# a)
f_prior <- function(mu){
  if(mu==0){1/(2*sqrt(2*pi))}
  else{(1-exp(-mu^(2)/2))/(sqrt(2*pi)*mu^(2))}
} 

g <- function(mu, xbar = 0.93, tao2 = 1){f_prior(mu) * exp(-(mu-xbar)^(2)/(2*tao2))}

trapezoidal <- function(f, interval, tol = 10^-6){
  nodes <- 4 # we start with 4 nodes
  h <- (interval[2]-interval[1])/nodes # the length of each subinterval
  x.i <- interval[1] + (1:(nodes-1))*h # the values of every x_i
  
  old <- h/2*f(interval[1])+sum(h*sapply(x.i, f))+h/2*f(interval[2])
  result.vector <- c(old) #save the first value of the integral
  i <- 3 # index
  conv <- 1
  nodes.vector <- c(nodes) # create a vector of difefrent nodes
  while (conv > tol )
  {
    nodes <- 2^i
    h <- (interval[2]-interval[1])/nodes
    x.i <- interval[1] + (1:(nodes-1))*h 
    result <- h/2*f(interval[1])+sum(h*sapply(x.i, f))+h/2*f(interval[2])
    result.vector <- c(result.vector,result)
    nodes.vector <- c(nodes.vector,nodes)
    conv <- abs(result/old-1) # stopping criterion
    old <- result
    i <- i+1
  }
  return (list(result,result.vector,nodes.vector))
}
interval <- c(-1000, 1000)
int <- trapezoidal(g, interval)
a <- int[[1]]
# build-in function
int2 <- integrate(g, -Inf, Inf, mu = 1)$value
cat("Comparison of the integrals:", a, int2)

# b)
grid <- seq(-10, 10, length = 1000)
plot(grid,c(sapply(grid[-1], g),0), type='l', ylab = "g(x)",
     xlab = expression(mu))
abline(v=c(-4,5), col = "red")

interval <- c(-4, 5)
int <- trapezoidal(g, interval)
a <- int[[1]]
# build-in function
int2 <- integrate(g, -Inf, Inf, mu = 1)$value
cat("Comparison of the integrals:", a, int2)

# c)
c <- a
round(c,4)

# d)
f_post <- function(mu){1/c*g(mu)}

f_exp <- function(mu){f_post(mu)*mu}
f_var <- function(mu){f_post(mu)*mu^2}

E_exp <- trapezoidal(f_exp, interval)
E_var <- trapezoidal(f_var, interval)
cat("Mean:", round(E_exp[[1]],4), "Variance:", round(E_var[[1]] - E_exp[[1]]^2,4))

######### 2 #########
n <- 10000
# a)
# f(x)
f <- function(x){
  if(x>1)         {return(0)}
  if(x<(-1))      {return(0)}
  if(-1<=x & x<=0){return(x + 1)}
  if(0<x & x<=1)  {return(1 - x)}
}
grid <- seq(-2, 2, length= 100)
plot(grid, sapply(grid, f), type = "l", ylab = "f(x)")

# F(x)
F_x <- function(x){
  if(x>1)         {return(1)}
  if(x<(-1))      {return(0)}
  if(-1<=x & x<=0){return((1/2)*x^2 + x + .5)}
  if(0<x & x<=1)  {return(.5 + x - (1/2)*x^2)}
}
plot(grid, sapply(grid, F_x), type = "l", ylab = "F(x)")
# F-1(y)
F1_y <- function(u){ 
  if(u>=1)        {return(1)}
  if(u<=(0))      {return(-1)}
  if(0<u & u<=0.5){return(-1 + sqrt(2)*sqrt(u))}
  if(0.5<u & u<1) {return(-((-10 + sqrt(-200*u + 200))/ 10))}
}

u <- runif(n, 0, 1)
x <- sapply(u, F1_y)
hist(x, main = "Histogram, sampled Inverse Transformation Method")

# b) 
plot(grid, sapply(grid, f), type = "l", ylab = "f(x)")

# Envelope
a <- 1/0.83 # Scaling parameter a < 1
e <- function(x){a*dnorm(x, 0, 0.45)} # normal distribution 
lines(grid, e(grid), col = 2) # a * e(x) > f(x) , all x
legend("topleft", legend = c("e(x)","f(x)"), col = c(2,1), lty = 19)

rejection <- function(n, envelope, funct, a = 1/0.83){
  # n - the size of the sample we want to create
  result <- c()
  for (i in 1:n){
    x.init <- rnorm(1,0, 0.45)
    u.init <- runif(1)
  if(u.init < f(x.init)/e(x.init))
    result <- c(result,x.init)
  }
  return(result)
}
set.seed(101);output <-rejection(n,e,f)
length(output)/n #percent of draws excepted

# c)
# Triangle 1
f_1y <- function(y){
  if(0<=y & y <=1){return(2-2*y)}
  else{return(0)}
  }

# Triangle 2
f_2y <- function(y){
  if(y<=0 & y >=-1){return((2+2*y))}
  else{return(0)}
}
grid <- seq(-1,1,length = 1000)
plot(grid, sapply(grid, f_1y), type = "l", ylab = "f(Y)")
lines(grid, sapply(grid, f_2y), type = "l", col = 2) # Looks like f(x) for sure
legend("topleft", legend = c("Y","-Y"), col = c(1,2), lty = 19)

library(triangle) # has a random generator for triangle distributions
composition <- function(n){
  x_c_vector <- c()
  for (i in 1:n) {
    p <- rbinom(1, size = 1, prob = 0.5)
    temp <- ifelse(p < 1, -rtriangle(1, 0, 1, 0), rtriangle(1, 0, 1, 0))
    x_c_vector[i] <- temp
  }
  return(x_c_vector)
}
x_c <- composition(n)
hist(x_c, breaks = 50, main = "Histogram", xlab = "Composition sampling")


# d)
Uni_diff <- function(n){
  u1 <- runif(n)
  u2 <- runif(n)
  x_lc <- u1-u2
  return(x_lc)
}
x_u <- Uni_diff(n)
hist(x_u, main = "Histogram", xlab = "Uniform diff sampling")

# e)
plot(grid, sapply(grid, f), type = "l", ylab = "f(x)", main = "X")
par(mfrow=c(2,2))
hist(x, breaks = 50, main = "Histogram a)", xlab = "Inverse Transformation Method")
hist(output, breaks = 50, main = "Histogram b)", xlab = "Envelope sampling")
hist(x_c, breaks = 50, main = "Histogram c)", xlab = "Composition sampling")
hist(x_u, breaks = 50, main = "Histogram d)", xlab = "Uniform diff sampling")
par(mfrow=c(1,1))

######### 3 #########

# a)
f <-function(x){
  exp(-x^2)*(2+cos(x*(64/pi)))/3.544909
}

sim <-1000000
set.seed(101);y <-rnorm(sim,sd=2)
w <- f(y)/dnorm(y,sd=2)

# Observations, 1 below times the weight or, 0, if condition is not met.
z <- (abs(y)>1)*w

# Mean of z.
p <- mean(z)/2
p

# Standard deviation of the importance sampling.
sdIS <-sqrt(var((y>1)*w)/sim)
sdIS

# Correlation for the antithetic sampling.
rho <- -p/(1-p)

# Standard deviation for the antithetic importance sampling.
sd  <- sdIS*(1+rho)/2
sd

# Confidence interval.
ci95_p <- round(c(p-1.96*sd, p, p+1.96*sd),4)
ci95_p

# Length of the CI
ci95_p_length <- ci95_p[3] - ci95_p[1]
ci95_p_length

# b)
B <- 500
p_bvector <- c()
time.start <- Sys.time() # 4 min run time
for(i in 1:B){
  y_b <- sample(y, sim, replace = T)
  w_b <- f(y_b)/dnorm(y_b,sd=2)
  z_b <- (abs(y_b)>1)*w_b
  p_b <- mean(z_b)/2
  p_bvector <- c(p_bvector, p_b)
}
Sys.time() - time.start
hist(p_bvector, breaks = 50, main = "Bootstrap", xlab = "P-sample distribution")
p_bvector2 <- sort(p_bvector)
p95c <- round(c(p_bvector2[round(B*0.025)], mean(p_bvector), p_bvector2[round(B*0.975)]),4)
p95c

p95c_length <- p95c[3] - p95c[1]
p95c_length

# c)
# The confidence in a) is shorter than in b). However, the bootstrap method requires less assumptions, 
# i.e., is a non-parametric method (see pdf for more analysis)

######### 4 #########

# a)  
library(BSDA)
mu <- seq(0,1, by = 0.1)
a <- 0.1
n <- 21
s <- 1000

t_power <- function(mu, a, n, s){
  
  talph<-qt(1-a, df= n-1)
  output <- c()
  
  for(i in 1:length(mu)){
    mu_temp <- mu[i]
    
    x <- matrix(rnorm(n*s, mu_temp, 1), ncol = n)
    meanv <- rowMeans(x)
    sdv <- sqrt(rowSums((x-meanv)^2)/(n-1))
    reject <-(sqrt(n)*meanv/sdv > talph)
    rre <- mean(reject) #Rejection rate estimate
    output <- c(output, rre)
    
  }
  output1 <- rbind(mu, output)
  return(output1)
}
set.seed(101);t_test <- t_power(mu, a = 0.1, n = 21, s = 10000)


sign_power <- function(mu, a, n, s){
  
  output <- c()
  
  for(i in 1:length(mu)){
    mu_temp <- mu[i]
    tmp_rre <- c()
    for (j in 1:s) {
      x <- rnorm(n, mu_temp, 1)
      tmp <- SIGN.test(x, md=0, conf.level = 1-a, alternative = "greater")
      tmp_rre <- c(tmp_rre, sum(tmp$p.value<=a))
    }
    
    rre <- mean(tmp_rre)
    output <- c(output, rre)
    
  }
  output1 <- rbind(mu, output)
  return(output1)
}
set.seed(101);sign_test <- sign_power(mu, a = 0.1, n = 21, s = 1000)

# b)
plot(t_test[1,], t_test[2,], type = "l", col = 3, ylab = "Power", xlab = expression(mu))
lines(sign_test[1,], sign_test[2,], type = "l", col = 2)
legend("topleft", legend = c("T-test","SIGN.test"), lty = 19 , col = c(3,2))

write_matex(sign_test)
#####################