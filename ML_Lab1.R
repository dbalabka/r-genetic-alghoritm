source("./my_ga.R");

f <- function(x)  (x^2 + x) * cos(x) #+ rnorm(length(x),1, 10)
min <- -10
max <- +10
curve(f, min, max)
fitness <- function(x) -f(x)
monitor <- function(obj) {
  curve(f, min, max, main = paste("iteration =", obj@iter), font.main = 1)
  points(obj@population, -obj@fitnessValues, pch = 20, col = 2)
  rug(obj@population, col = 2)
  #Sys.sleep(0.1)
}
GA <- my_ga(fitness, min = min, max = max, monitor = monitor);
plot(GA)
summary(GA)
opt.sol <- optimize(f, lower = min, upper = max, maximum = FALSE)
nlm.sol <- nlm(function(...) -f(...), 0, typsize = 0.1)
curve(f, min, max)
points(GA@solution, GA@fitnessValue*-1, col = 2, pch = 20)
points(opt.sol$minimum, opt.sol$objective, col = 3, pch = 8)
points(nlm.sol$estimate, -nlm.sol$minimum, col = 4, pch = 17)
legend(x = -5, y = -30, legend = c("ga", "optimize", "nlm"), title = "Solutions", pch = c(20,8,17), col = 2:4)

f <- function(x)  (x^2 + x) * cos(x) #+ rnorm(length(x),1, 10)
min <- -10
max <- +10
curve(f, min, max)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, min = min, max = max)
plot(GA)
summary(GA)
monitor <- function(obj) {
  curve(f, min, max, main = paste("iteration =", obj@iteration), font.main = 1)
  points(obj@population, -obj@fitnessValues, pch = 20, col = 2)
  rug(obj@population, col = 2)
  Sys.sleep(0.1)
}
GA <- ga(type = "real-valued", fitness, min = min, max = max, monitor = monitor)

opt.sol <- optimize(f, lower = min, upper = max, maximum = FALSE)
nlm.sol <- nlm(function(...) -f(...), 0, typsize = 0.1)
curve(f, min, max)
points(GA@solution, GA@fitnessValue*-1, col = 2, pch = 20)
points(opt.sol$minimum, opt.sol$objective, col = 3, pch = 8)
points(nlm.sol$estimate, -nlm.sol$minimum, col = 4, pch = 17)
legend(x = -5, y = -30, legend = c("ga", "optimize", "nlm"), title = "Solutions", pch = c(20,8,17), col = 2:4)

Rastrigin <- function(x1, x2) {
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}
x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20)
filled.contour(x1, x2, f, color.palette = jet.colors)
monitor <- function(obj) {
  contour(x1, x2, f, drawlabels = FALSE, col = gray(0.5))
  title(paste("iteration =", obj@iter), font.main = 1)
  points(obj@population, pch = 20, col = 2)
  Sys.sleep(0.2)
}

GA <- ga(type = "real-valued",
  fitness = function(x) -Rastrigin(x[1], x[2]),
  min = c(-5.12, -5.12), max = c(5.12, 5.12), popSize = 50,
  maxiter = 100, monitor = monitor)
summary(GA)
NLM <- nlm(function(x) Rastrigin(x[1], x[2]), GA@solution)
NLM[c("minimum", "estimate")]


Noisy <- function (x1, x2) {
  x1^4+2*x2^4 +  rnorm(length(x1),1, 10)
}
x1 <- x2 <- seq(-2, 2, by = 0.1)
f <- outer(x1, x2, Noisy)
persp3D(x1, x2, f, theta = 50, phi = 20)
filled.contour(x1, x2, f, color.palette = jet.colors)
monitor <- function(obj) {
  contour(x1, x2, f, drawlabels = FALSE, col = gray(0.5))
  title(paste("iteration =", obj@iter), font.main = 1)
  points(obj@population, pch = 20, col = 2)
  Sys.sleep(0.1)
}

GA <- ga(type = "real-valued",
         fitness = function(x) -Noisy(x[1], x[2]),
         min = c(-5.12, -5.12), max = c(5.12, 5.12), popSize = 50,
         maxiter = 100, monitor = monitor)

