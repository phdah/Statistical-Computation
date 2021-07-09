######### Home assignment 1 ######### 
# Philip Dahlqvist-Sjöberg
# 02-04-2020

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

wd <- "C:/Users/"
######### 1 #########
data <- read.csv2(paste0(wd,"villorSW.csv")) %>%
  mutate(x0 = rep(1, 29))

# a)
x <- matrix(c(data$x0, data$x1, data$x2, data$x3, data$x4), ncol = 5)
colnames(x) <- c("x0","x1","x2","x3","x4")

xx <- t(x)%*%x
dimnames(xx) <- list(c("x0","x1","x2","x3","x4"),c("x0","x1","x2","x3","x4"))

non_singular <- round(solve(xx)%*%xx,2) # = the identity matrix, i.e., the X*X matrix is nonsingualar. 

xy <- t(x)%*%data$y
dimnames(xy) <- list(c("x0","x1","x2","x3","x4"),c("y"))

write_matex(xx)
write_matex(xy)
write_matex(non_singular)

# b) 
beta <- round(solve(xx)%*%xy,4) # beta coefficients, rounded at 4 decimal (same as lm)
dimnames(beta) <- list(c("x0","x1","x2","x3","x4"),c("Coefficients calc"))
write_matex(beta)

# c)
fit <- lm(y~x1+x2+x3+x4, data = data)

beta_lm <- round(matrix(fit$coefficients, ncol = 1),4)
dimnames(beta_lm) <- list(c("x0","x1","x2","x3","x4"),c("Coefficients lm"))
write_matex(beta_lm)
# d)
beta_diff <- beta-beta_lm
colnames(beta_diff) <- "Coefficients diff (4 dec)"

set.time <- Sys.time()
set.seed(100):for(i in 1:10000){
  solve(t(x)%*%x)%*%t(x)%*%data$y
}
time1 <- Sys.time() - set.time # Calculation: Time difference of 0.6950631 secs

set.time <- Sys.time()
set.seed(100):for(i in 1:10000){
  lm(y~x1+x2+x3+x4, data = data)
}
time2 <- Sys.time() - set.time # lm: Time difference of 13.30621 secs

round(as.vector(time1)/as.vector(time2),3)

######### 2 #########
# a)
par(mfrow = c(1,2))

grid <- seq(0, 4, length = 100)
eq = function(x){(log(x+1))/(x^(3/2)+1)}
plot(grid,eq(grid), type='l', ylab = "g(x)")
# The maximum looks to be at x = 1

# b)
u <- log(x+1)
du <- 1/(x+1)

v <- x^(3/2)+1
dv <- (3/2)*x^(3/2-1)

# deq_R <- function(x){D(expression((log(x+1))/(x^(3/2)+1)), "x")}
deq <- function(x){(1/(x+1)*(x^(3/2)+1)-(3/2)*x^(3/2-1)*log(x+1))/(x^(3/2)+1)^2}
plot(grid, deq(grid), type='l', ylab = "g'(x)")
abline(h = 0, col = "red")

par(mfrow = c(1,1))

# c) 
# Package
local_max <- optimise(eq,c(0,2),maximum =TRUE, tol = 0.000001)#we use the default accuracy tol = 0.0001220703
x_value <- local_max[[1]] # x-value
y_value <- eq(local_max[[1]]) # y-value
local_max_matrix <- round(matrix(c(x_value, y_value), ncol = 2),4)
dimnames(local_max_matrix) <- list("Local maximum", c("X value", "Y value"))

# Own algorithm
bisection_opt <- function(fun = deq, a = 0, b = 2, eps = 0.000001){
  if(fun(a)*fun(b)>0){break()}
  cc <- eps + 100
  while(cc > eps){
    x <- (a+b)/2
    
    if(fun(a)*fun(x)<=0)
      {b <- x} else {a <- x}
    
    cc <- abs(b-a)
   
  }
  return(mean(a,b))
}
time.start <- Sys.time()
x_value <- bisection_opt() # x-value
Sys.time()-time.start # Time difference of 0.06684709 secs
y_value <- eq(x_value) # y-value
local_max_section <- round(matrix(c(x_value, y_value), ncol = 2),4) # Same as package
dimnames(local_max_section) <- list("Local maximum", c("X value", "Y value"))

# d)

secant_opt <- function(fun = deq, x1 = 0.8, x2 = 0.7, eps = 0.000001){
  x_1 <- x1
  x_2 <- x2
  cc <- eps + 100
  while(cc > eps){
    x_new <- x_1-fun(x_1)*((x_1-x_2)/(fun(x_1)-fun(x_2)))
    
    cc <- abs(x_new-x_1)
    x_2 <- x_1
    x_1 <- x_new
  }
  return(x_new)
}
time.start <- Sys.time()
x_value <- secant_opt() # x-value
Sys.time()-time.start # Time difference of 0.01559687 secs
y_value <- eq(x_value) # y-value
local_max_secant <- round(matrix(c(x_value, y_value), ncol = 2),4) # Same as package
dimnames(local_max_secant) <- list("Local maximum", c("X value", "Y value"))

all.equal(local_max_section[1:2], local_max_secant[1:2])
# They reach the same local maximum, however secant method is slighly faster

######### 3 #########
# a)
# See pdf for calculations. 

# b)
# Manually
eq <- function(a){x <- matrix(c(1,-1,1,-1,1,-a,a^2,-a^3,1,a,a^2,a^3,1,1,1,1),nrow = 4, ncol = 4, byrow = T)
  xx <- t(x)%*%x
  col2 <- xx[, 3]
  col3 <- xx[, 2]
  xx[ , 2] <- col2
  xx[ , 3] <- col3
  row2 <- xx[3, ]
  row3 <- xx[2, ]
  xx[2, ] <- row2
  xx[3, ] <- row3
  
  det_matrix <- (-1)*(-1)*(xx[1,1]*xx[2,2]-xx[1,2]*xx[2,1])*
    (xx[3,3]*xx[4,4]-xx[3,4]*xx[4,3])
  return(det_matrix)
  }

# With det()

Req <- function(a){det(t(matrix(c(1,-1,1,-1,1,-a,a^2,-a^3,1,a,a^2,a^3,1,1,1,1),nrow = 4, ncol = 4, byrow = T))%*%
    matrix(c(1,-1,1,-1,1,-a,a^2,-a^3,1,a,a^2,a^3,1,1,1,1),nrow = 4, ncol = 4, byrow = T))}
# Req <- function(a){det(matrix(c(4,0,2+2*a^2,0,
#                             0,2+2*a^2,0,2+2*a^4,
#                             2+2*a^2,0,2+2*a^4,0,
#                             0,2+2*a^4,0,2+2*a^6), 
#                             nrow = 4, ncol = 4, byrow = T))} # This is the X'X matrix pre-calculated

# c)
grid <- seq(0, 1, length = 1000)
par(mfrow=c(1,2))

grid1 <- sapply(grid, eq)
plot(grid, grid1, type = "l", ylab = "Funktion value", xlab = "Interval", main = "Without det()")

grid2 <- sapply(grid, Req)
plot(grid, grid2, type = "l", col = "red", ylab = "Funktion value", xlab = "Interval", main = "With det()")
par(mfrow=c(1,1))
D <- grid1-grid2

plot(grid, D, type = "l", ylab = "Difference", xlab = "Interval", main = "Difference between functions")
abline(grid, 0, col = "red")

# d)
det_max <- optimize(eq, c(0,1), maximum = T)
plot(grid, grid1, type = "l", ylab = "Funktion value", xlab = "Interval", main = "Maximum")
abline(v=det_max[1], col = "red")
points(det_max[1], det_max[2], type = "p")

######### 4 #########

# Step 1
ar <- 1
ae <- 2
ac <- 0.5
as <- 0.5

# Plot 1
# First iteration
# Step 2
x_best  <- c(0.72, 0.31)
x_bad   <- c(0.71, 0.26)
x_worst <- c(0.63, 0.19)

g_x_best <- 150
g_x_bad <- 140
g_x_worst <- 130

# Step 3
c <- (x_best + x_bad)/2

# Step 4
x_r <- c + ar*(c-x_worst)
g_x_r <- 130
# g(x_worst) >= g(x_r): True, perform inner contraction

# Step 6 b
x_c <- c + ac*(x_worst - c)
g_x_c <- 140
# g(x_c) > g(x_worst): True, discard x_worst and accept x_c as new vertices

# Second iteration
# Step 2
x_worst <- x_bad
g_x_worst <- g_x_bad
x_bad <- x_c
g_x_bad <- g_x_c

# Step 3
c <- (x_best + x_bad)/2

# Step 4
x_r <- c + ar*(c-x_worst)
g_x_r <- 150
# g(x_best) >= g(x_r) > g(x_bad): True, discard x_worst and accept x_r (stopping step)
x_worst <- x_bad
g_x_worst <- g_x_bad
x_bad <- x_r
g_x_bad <- g_x_r
# Final vertices:
final1_vertices <- matrix(c(x_best, g_x_best, x_bad, g_x_bad, x_worst, g_x_worst), ncol = 3, byrow = T)
dimnames(final1_vertices) <- list(c("x_best", "x_bad", "x_worst"),c("X-axis","Y-axis", "g(X,Y)"))

# Both x_best and x_bad have the local maximum after two iterations. 

# Plot 2
# First iteration
# Step 2
x_best  <- c(0.38, 0.79)
x_bad   <- c(0.45, 0.66)
x_worst <- c(0.52, 0.72)

g_x_best <- 160
g_x_bad <- 150
g_x_worst <- 120

# Step 3
c <- (x_best + x_bad)/2

# Step 4
x_r <- c + ar*(c-x_worst)
g_x_r <- 180
# g(x_r) > g(x_best): True, perform expansion

# Step 5
x_e <- c + ae*(x_r - c)
g_x_e <- 160
# g(x_e) > g(x_r): Not true, discard x_worst and accept x_r

# Second iteration
# Step 2
x_worst <- x_bad
x_bad <- x_best
x_best <- x_r

g_x_worst <- g_x_bad
g_x_bad <- g_x_best
g_x_best <- g_x_r

# # Step 3
# c <- (x_best + x_bad)/2
# 
# # Step 4
# x_r <- c + ar*(c-x_worst)
# g_x_r <- 130
# # g(x_r) < g(x_bad) < g(x_best): True, perform inner contraction
# 
# # Step 6 b
# x_c <- c + ac*(x_worst - c)
# g_x_c <- 160
# # g(x_c) > g(x_worst): True, discard x_worst and accept x_c as new vertices
# 
# x_worst <- x_c
# g_x_worst <- g_x_c

# Final vertices:
final2_vertices <- matrix(c(x_best, g_x_best, x_bad, g_x_bad, x_worst, g_x_worst), ncol = 3, byrow = T)
dimnames(final2_vertices) <- list(c("x_best", "x_bad", "x_worst"),c("X-axis","Y-axis", "g(X,Y)"))
