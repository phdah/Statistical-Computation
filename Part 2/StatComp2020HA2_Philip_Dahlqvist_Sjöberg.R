######### Home assignment 2 ######### 
# Philip Dahlqvist-Sjöberg
# 09-04-2020

library(dplyr)
library(mvtnorm)
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

wd <- "C:/Users/"

######### 1 ######### 

eq <- function(x,y){sin(x + y) + (x - y)^2 - 1.5*x + 2.5*y + 1}
eq(2,4)

# a) 
x_grid <- seq(-1.5, 4,  length = 100)
y_grid <- seq(-3, 4, length = 100)
grid <- expand.grid(x_grid, y_grid)

eq_value <- as.vector(eq(grid[1],grid[2]))
eq_plot <- matrix(eq_value$Var1, nrow = 100, ncol = 100)

contour(x_grid, y_grid, eq_plot)

# b) 

# Gradient
gradient <- function(x, y) {
  gradient_x <- -1.5 + 2*x - 2*y + cos(x + y)
  gradient_y <- 2.5 - 2*x + 2*y + cos(x + y)
  matrix(c(gradient_x, gradient_y),2,1)
}
gradient(1,2)

# Hessian 
hessian <- function(x, y) {
  hessian_xx <- 2-sin(x+y)
  hessian_yy <- 2-sin(x+y)
  hessian_xy <- -2-sin(x+y)
  matrix(c(hessian_xx,hessian_xy,hessian_xy,hessian_yy),2,2,byrow=TRUE)
}
hessian(1,2)
newton.raphson <- function(start,gradient,hessian,col,type, tol = 0.000001) {
  conv <- 100
  iter <- 0
  while (conv > tol)
  {
    x <- start[1]
    y <- start[2]
    start.new <- c(x, y) - solve(hessian(x,y)) %*% gradient(x, y)
    conv = sum(abs(start.new-start))/sum(abs(start-0))
    iter <- iter +1
    start <- start.new
    a=lines(c(start[1],start.new[1]),c(start[2],start.new[2]),type=type,col=col)
    print(paste0("Function value: ", eq(start[1], start[2])))
  }
  return(list(start,iter))
}
contour(x_grid, y_grid, eq_plot, nlevels = 20)
# Starting point
lines(0,-1, col = 4, type = "b") 
# Best starting point c(0,-1)
newton.raphson(c(0,-1), gradient, hessian, col = 3, type ="p")
# Optimal point
lines(-0.5471976, -1.5471976, col = 2, type = "b")
legend(-0.5, -2, legend = "Local minimum:  -1.9132", fill = "red")

# c)
# Gradient
gradient <- function(x, y) {
  gradient_x <- -1.5 + 2*x - 2*y + cos(x + y)
  gradient_y <- 2.5 - 2*x + 2*y + cos(x + y)
  matrix(c(gradient_x, gradient_y), 2, 1)
}

# Hessian 
hessian <- function(x, y) {
  hessian_xx <- 2-sin(x+y)
  hessian_yy <- 2-sin(x+y)
  matrix(c(hessian_xx, hessian_yy), 2, 1)
}

newton.raphson <- function(start,gradient,hessian,col,type, tol = 0.000001) {
  conv_left <- 100
  iter_left <- 0
  # Left boundary (y)
  while (conv_left > tol){
    x <- -1.5
    y <- start[1]
    start.new_Left <- y - gradient(x, y)[2] / hessian(x,y)[2]
    conv_left = sum(abs(start.new_Left-start[1]))/sum(abs(start[1]-0))
    iter_left <- iter_left +1
    
    if(start.new_Left > 4){
      print("***Left search converged to large***")
      break
    }
    if(start.new_Left < -3){
      print("***Left search converged to small***")
      break
    }
    
    start[1] <- start.new_Left
    a=lines(x,c(start[1]),type=type,col=col[1])
    print(paste0("Function value left: ", eq(x, start[1])))
  }
  conv_right <- 100
  iter_right <- 0
  # Right boundary (y)
  while (conv_right > tol){
    x <- 4
    y <- start[2]
    start.new_right <- y - gradient(x, y)[2] / hessian(x,y)[2]
    conv_right = sum(abs(start.new_right-start[2]))/sum(abs(start[2]-0))
    iter_right <- iter_right +1
    
    if(start.new_right > 4){
      print("***Right search converged to large***")
      break
    }
    if(start.new_right < -3){
      print("***Right search converged to small***")
      break
    }
    
    start[2] <- start.new_right
    a=lines(x,c(start[2]),type=type,col=col[2])
    print(paste0("Function value right: ", eq(x, start[2])))
  }
  conv_top <- 100
  iter_top <- 0
  # Right boundary (x)
  while (conv_top > tol){
    x <- start[3]
    y <- 4
    start.new_top <- x - gradient(x, y)[1] / hessian(x,y)[1]
    conv_top = sum(abs(start.new_top-start[3]))/sum(abs(start[3]-0))
    iter_top <- iter_top +1
    
    if(start.new_top > 4){
      print("***Top search converged to large***")
      break
    }
    if(start.new_top < -1.5){
      print("***Top search converged to small***")
      break
    } 
    
    start[3] <- start.new_top
    a=lines(c(start[3]), y,type=type,col=col[3])
    print(paste0("Function value top: ", eq(start[3],y)))
  }
  conv_bot <- 100
  iter_bot <- 0
  # Right boundary (x)
  while (conv_bot > tol){
    x <- start[4]
    y <- -3
    start.new_bot <- x - gradient(x, y)[1] / hessian(x,y)[1]
    conv_bot = sum(abs(start.new_bot-start[4]))/sum(abs(start[4]-0))
    iter_bot <- iter_bot +1
    
    if(start.new_bot > 4){
      print("***Bottom search converged to large***")
      break
    }
    if(start.new_bot < -1.5){
      print("***Bottom search converged to small***")
      break
    } 
    
    start[4] <- start.new_bot
    a=lines(c(start[4]), y,type=type,col=col[4])
    print(paste0("Function value bottom: ", eq(start[4],y)))
  }
  return(matrix(c("Y left", "Y right", "X top", "X bottom",start[1], start[2], start[3], start[4], iter_left, iter_right, iter_top, iter_bot)
                , ncol = 3))
}

contour(x_grid, y_grid, eq_plot, nlevels = 20)
# start = (left, right, top, bottom)
newton.raphson(c(-2, -2, 1, 1), gradient, hessian, col = c(2,3,4,5), type ="p")
# No minimum on the border 

# # Plot 
# par(mfrow=c(1,2))
# contour(x_grid, y_grid, eq_plot, nlevels = 20)
# # Starting point
# lines(0,-1, col = 4, type = "b") 
# # Best starting point c(0,-1)
# newton.raphson(c(0,-1), gradient, hessian, col = 3, type ="p")
# # Optimal point
# lines(-0.5471976, -1.5471976, col = 2, type = "b")
# legend(-0.5, -2, legend = "Local minimum:  -1.9132", fill = "red")
# 
# contour(x_grid, y_grid, eq_plot, nlevels = 20)
# # start = (left, right, top, bottom)
# newton.raphson(c(-2, -2, 1, 1), gradient, hessian, col = c(2,3,4,5), type ="p")
# par(mfrow=c(1,1))

######### 2 ######### 
# a) 
data <- read.csv(paste0(wd,"threepops.csv"), header = F) 
data <- data$V1
emalg <- function(dat,eps=0.00001, type){
  n      <- length(dat)
  pi     <- matrix(c(rep(NA,n*4)), ncol = 4)    #initialize vector for prob. to belong to group 1  
  p1     <- 1/3          #Starting value for mixing parameter 
  p2     <- 1/3 
  p3     <- 1/3 
  sigma1 <- sd(dat)*1/3  #Starting value for variances
  sigma2 <- sigma1
  sigma3 <- sigma2
  mu1    <- mean(dat)    #Starting values for means
  mu2    <- mean(dat)-sigma1/2
  mu3    <- mean(dat)+sigma1/2
  pv     <- c(p1, p2, p3 , mu1, mu2, mu3, sigma1, sigma2, sigma3)  #parameter vector
  pv_plot <- round(pv, 2)
  pi_plot <- rep(2, n)
  cc     <- eps + 100    #initialize conv. crit. not to stop directly 
  while(cc > eps){
    pv1  <- pv           #Save previous parameter vector
    ### E step ###
    for (j in 1:n){
      pi1   <- p1*dnorm(dat[j],mean=mu1,sd=sigma1)
      pi2   <- p2*dnorm(dat[j],mean=mu2,sd=sigma2)
      pi3   <- p3*dnorm(dat[j],mean=mu3,sd=sigma3)
      pall  <- c(pi1/(pi1+pi2+pi3), pi2/(pi1+pi2+pi3), pi3/(pi1+pi2+pi3))
      pi[j, 1] <- which.max(pall)
      pi[j, 2:4] <- pall
    }
    ### M step ###
    p1      <- mean(pi[,2])
    p2      <- mean(pi[,3])
    p3      <- mean(pi[,4])
    mu1    <- sum(pi[,2]*dat)/(p1*n)
    mu2    <- sum(pi[,3]*dat)/(p2*n)
    mu3    <- sum(pi[,4]*dat)/(p3*n)
    sigma1 <- sqrt(sum(pi[,2]*(dat-mu1)*(dat-mu1)/(p1*n)))
    sigma2 <- sqrt(sum(pi[,3]*(dat-mu2)*(dat-mu2)/(p2*n)))
    sigma3 <- sqrt(sum(pi[,4]*(dat-mu3)*(dat-mu3)/(p3*n)))
    ######
    pv     <- c(p1, p2, p3 , mu1, mu2, mu3, sigma1, sigma2, sigma3)
    pv_ratio <- c(1, 1, p3 , mu1, mu2, mu3, sigma1, sigma2, sigma3)
    pv_plot <- rbind(pv_plot, pv)
    pi_plot <- rbind(pi_plot, pi[,2:4])
    if(type == "ratio"){
      cc <- sum(abs((pv-pv1)/(pv_ratio)))
    }else{
    cc     <- sum(abs(pv-pv1))
    }
  }
  pv <- round(pv, 2)
  names(pv) <- c("p_1", "p_2", "p_3" , "mu1", "mu2", "mu3", "sigma1", "sigma2", "sigma3")
  print(paste0("X1~N(",pv[4],",",pv[7],") , X2~N(",pv[5],",",pv[8],") , X3~N(",pv[6],",",pv[9],")", ": After ", nrow(pv_plot), " iterations"))
  # if(plot = "yes"){
  par(mfrow=c(1,2))
  plot(1, type="l", xlab="", ylab="", xlim=c(0, nrow(pv_plot)), ylim=c(min(data), max(data)), main = "Mean parameter")
  lines(pv_plot[,4], col = 2)
  lines(pv_plot[,5], col = 3)
  lines(pv_plot[,6], col = 4)
  plot(1, type="l", xlab="", ylab="", xlim=c(0, nrow(pv_plot)), ylim=c(min(data), max(data)), main = "Sigma parameter")
  lines(pv_plot[,7], col = 2)
  lines(pv_plot[,8], col = 3)
  lines(pv_plot[,9], col = 4)
  
  plot(1, type="l", xlab="", ylab="", xlim=c(1, nrow(pv_plot)), ylim=c(0,1), main = "P1 parameter")
  lines(pi_plot[,1], col = 2)
  plot(1, type="l", xlab="", ylab="", xlim=c(1, nrow(pv_plot)), ylim=c(0,1), main = "P2 parameter")
  lines(pi_plot[,2], col = 3)
  plot(1, type="l", xlab="", ylab="", xlim=c(1, nrow(pv_plot)), ylim=c(0,1), main = "P3 parameter")
  lines(pi_plot[,3], col = 4)
  
  par(mfrow=c(1,1))
  # }
  return(pv)
}
emalg(data, type = "normal")

# b)
emalg(data, type = "ratio") # Added ratio criteria

# c)
data <- read.csv(paste0(wd,"threepops.csv"), header = F) 
data <- data$V1
hist(data, nclass = 50) # "X1~N(3.83,0.72) , X2~N(1.07,1.2) , X3~N(6.58,0.6): After 200 iterations"
points(c(3.83, 1.07, 6.58), c(0,0,0), col = c(2,3,4), pch = 19)
legend(-2, 15, "Marked means from EM-alg.")
emalg(data, type = "normal")

# d) 
# Check plots
emalg(data, type = "normal")
par(mfrow=c(1,1))
######### 3 ######### 
# a)
emalg_multi <- function(dat, p_start = 0.5, sigma1_start = matrix(c(1,0,0,1),nrow = 2), sigma2_start = matrix(c(1,0,0,1),nrow = 2), 
                        mu1_start = c(0,0), mu2_start = c(0,0), eps=0.000001){
  iter   <- 0 
  pv_q   <- NULL
  n      <- nrow(dat)
  pi     <- rep(NA,n)     #initialize vector for prob. to belong to group 1 
  Qlog   <- rep(NA,n)
  p      <- p_start       #Starting value for mixing parameter 
  sigma1 <- sigma1_start  #Starting value for variances
  sigma2 <- sigma2_start
  mu1    <- mu1_start     #Starting values for means
  mu2    <- mu2_start
  pv     <- c(p, mu1, mu2, sigma1, sigma2)  #parameter vector
  cc     <- eps + 100     #initialize conv crit
  while(cc>eps){
    pv1  <- pv            #Save previous parameter vector
    pv_q <- pv_q
    ### E step ###
    for (j in 1:n){
      pi1   <- p*dmvnorm(dat[j,], mean = mu1, sigma = sigma1)
      pi2   <- (1 - p)*dmvnorm(dat[j,], mean = mu2, sigma = sigma2)
      pi[j] <- pi1/(pi1 + pi2)
    }
    ### M step ###
    p      <- mean(pi)
    mu1    <- c(colSums(pi*dat)/(p*n))
    mu2    <- c(colSums((1 - pi)*dat)/((1 - p)*n))
    sigma1 <- matrix(
      c(sum(pi*(dat[,1] - mu1[1])*(dat[,1] - mu1[1]))/(p*n), 
        sum(pi*(dat[,1] - mu1[1])*(dat[,2] - mu1[2]))/(p*n),
        sum(pi*(dat[,1] - mu1[1])*(dat[,2] - mu1[2]))/(p*n), 
        sum(pi*(dat[,2] - mu1[2])*(dat[,2] - mu1[2]))/(p*n)), nrow=2)
    sigma2 <- matrix(
      c(sum((1-pi)*(dat[,1]-mu2[1])*(dat[,1]-mu2[1]))/((1-p)*n), 
        sum((1-pi)*(dat[,1]-mu2[1])*(dat[,2]-mu2[2]))/((1-p)*n),
        sum((1-pi)*(dat[,1]-mu2[1])*(dat[,2]-mu2[2]))/((1-p)*n), 
        sum((1-pi)*(dat[,2]-mu2[2])*(dat[,2]-mu2[2]))/((1-p)*n)), nrow=2)
    for (i in 1:n) {
      Qlog[i] <- pi[i]*(log(p) + dmvnorm(dat[i,], mean = mu1, sigma = sigma1,log = TRUE)) + 
        (1-pi[i])*(log(1 - p)+dmvnorm(dat[i,], mean = mu2, sigma = sigma2, log = TRUE))
    }
    ######
    iter <- iter + 1
    Q_expected <- mean(Qlog)
    pv     <- c(p, (1-p), mu1, mu2, sigma1, sigma2)
    pv_tmp <- c(Q_expected, pv)
    pv_tmp <- t(pv_tmp)
    pv_q <- rbind(pv_q, pv_tmp)
    cc     <- t(pv-pv1)%*%(pv-pv1)
  }
  
  pv_q <- as.data.frame(pv_q)
  names(pv_q) <- c("Expected Q", "p1", "p2", "mu11", "mu12", "mu21", "mu22", "sigma1_11", "sigma1_12", "sigma1_21", "sigma1_22", 
                                    "sigma2_11", "sigma2_12", "sigma2_21", "sigma2_22")
  output <- list(iter, round(pv_q, 2))
  names(output) <- c("Iteration", "Parameters")
  return(output)
}

data <- read.csv(paste0(wd, "bivardat.csv"), header = F)
output <- emalg_multi(data) # Default: Bivariate standard normal
output 

# b) 
data <- read.csv(paste0(wd, "bivardat.csv"), header = F)
plot(data$V1, data$V2, xlab = "X-axis", ylab = "Y-axis", main = "Two-dimensional scatter plot") # Looks to bee two peaks at c(2.2, 3.5) and c(3,3) and quite small variance. 
hist(data$V1) # Mean looks to be at 2.2
hist(data$V2) # Mean look to be at 3.5

# c) 
set.time <- Sys.time()
output <- emalg_multi(data, p_start = 0.7, sigma1_start = matrix(c(1.5,1/20,1/20,1.5),nrow = 2), 
                      sigma2_start = matrix(c(1,1/20,1/20,1),nrow = 2), 
                      mu1_start = c(2.2,3.5), mu2_start = c(3,3))
Sys.time()-set.time
output
sigma1_start = matrix(c(1.5,1/20,1/20,1.5),nrow = 2)
sigma2_start = matrix(c(1,1/20,1/20,1),nrow = 2)
mu1_start = matrix(c(2.2,3.5),nrow = 2)
mu2_start = matrix(c(3,3),nrow = 2)
p <- matrix(c(0.7,1-0.7),nrow = 2)
write_matex(sigma1_start)
write_matex(sigma2_start)
write_matex(mu1_start)
write_matex(mu2_start)
write_matex(p)

######### 4 ######### 
# a) 
node <- matrix(c(1,-3.66847 ,0.00000143956,
                  2,-2.78329 ,0.000346819,
                  3,-2.02595, 0.0119114,
                  4,-1.32656 ,0.117228,
                  5,-0.656810, 0.429360,
                  6,0, 0.654759,
                  7,0.656810, 0.429360,
                  8,1.32656 ,0.117228,
                  9,2.02595, 0.0119114,
                  10,2.78329 ,0.000346819,
                  11,3.66847 ,0.00000143956), ncol = 3, byrow = T)
colnames(node) <- c("Node i", "Node xi", "Weights Ai")
plot(node[,2], node[,3], xlab = "Node i", ylab = "Weight ai", type = "b")
write_matex(node)
# b)
wq <- function(x){exp(-x^2)}
# i)
eq <- function(x){dt(x, 3)}
eq_star <- function(x)(dt(x, 3))/(exp(-x^2))
int1 <- round(sum(node[,3]*eq_star(node[,2])),4)

# ii)
eq <- function(x){0.5*exp(-abs(x))}
eq_star <- function(x){((1/2)*exp(-abs(x)))/(exp(-x^2))}
int2 <- round(sum(node[,3]*eq_star(node[,2])),4)

# iii)
eq <- function(x){dnorm(x, mean = c(0,1,2,3), sd = 1)}
eq_star <- function(x, i)(dnorm(x, mean = i, sd = 1))/(exp(-x^2))
int3 <- NULL
for (i in 0:3) {
  int3[i+1] <- round(sum(node[,3]*eq_star(node[,2], i)),4)
}

integration_table <- matrix(c("i)", "ii)", "iii) mean 0", "iii) mean 1", "iii) mean 2",
                              "iii) mean 3", int1, int2, int3[1],int3[2], int3[3], int3[3]), ncol = 2)
integration_table
write_matex(integration_table)
