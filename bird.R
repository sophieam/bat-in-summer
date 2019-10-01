## From table 5.1 (p. 110) in Dynamic state variable models in ecology.

# Maximum energy reserves.
x_max = 2.4              # g

# Discretizised values of x.
x_d <- seq(from = 0, to = x_max, length.out=100)

# Basic daily predation risk patch 1/2.
mu <- c(0, 0.001, 0.005) # /day

# Increase in predation risk with body mass.
lambda = 0.46            # /g

# Net daily forage inntake for patch 1/2.
e      = c(0, 0.6, 2.0)  # g/day

# Mass of bird with zero fat reserves.
m_0    = 10              # g

# Metabolic cost in a good and bad night respectively.
c_g    = 0.48            # g
c_b    = 1.20            # g

# Probability of a night being bad.
p_b    = 0.167

#Metabolic rate.
gamma  = 0.04             # /g day

# Time periods per day 
Time = 50


# F(x, t, d) - Fitness (probability of survival)
Fitness <- array(data = NA, dim = c(length(x_d), Time+1), dimnames = list(
  x_d,
  1:(Time+1)) )

# Optimal decision given time and state: 
Decision <- array(data = NA, dim = c(length(x_d), Time), dimnames = list(
  x_d,
  1:Time))


# Ensure that a value is between a minimum and a maximum.
cap <- function(maximum, avalue) {
  if (avalue > maximum) {
    return(maximum)
  } else
    return(avalue)
}

# Returns the index of the closest discrete x value 
closest_discrete_x <- function(x) {
  return(which(abs(x_d - x) == min(abs(x_d - x)))[1])
}


# Interpolate between two values a and b. dx is a value between 0 and 1.
linear_interpolation <- function (a, b, dx) {
  return((1-dx)*a+ b*dx)
}

# Interpolate a fitness value for a given combination of x, t and d.
interpolate <- function (x, t) {
  closest <- closest_discrete_x(x)
  
  if (x < 0) { return(0) }
  
  if (x < x_d[closest]) {
    j1 <- closest -1
    j2 <- closest
  } else if (x > x_d[closest]){
    j1 <- closest
    j2 <- closest +1
  } else {
    return(Fitness[closest, t])
  }
  
  delta_x <- (x-x_d[j1])/(x_d[j2]-x_d[j1])
  
  return(linear_interpolation(Fitness[j1, t], Fitness[j2, t], delta_x))
}
  
# Terminal reward - Equation 5.1
for (j in 1:length(x_d)) {
  x <- x_d[j]
  
  if (x < c_g) {
    Fitness[j, Time +1] <- 0
  } else if (x < c_b) {
    Fitness[j, Time +1] <- (1-p_b)
  } else {
    Fitness[j, Time +1] <- 1
  }
}

for (t in Time:1) {
  for (j in 1:length(x_d)) {
    x <- x_d[j]
    # Vector for storing fitness values for each of the patches.
    F_i <- vector(mode = 'numeric', length=3)
    
    #Calculate resulting fitness of choosing each patch
    for (h in 1:3) {
      # Equation 5.4
      x_mark <- cap(x_max, (x + (1/Time)*(e[h] - gamma*(m_0+x)) ))
      
      
      # Equation 5.3
      F_i[h] <- ( 1 - (1/Time)*(mu[h]*exp(lambda*x)) ) *interpolate(x_mark, t+1)
    }
    
    # Fitness is the fitness of the patch that maximizes fitness.
    Fitness[j, t] <- max(F_i)[1]
    
    # Optimal patch choice is the one that maximizes fitness.
    # In cases where more than one patch shares the same fitness, 
    # only the first one is considered.
    Decision[j, t] <- which(F_i == max(F_i))[1]
  }
}



View(Decision)
View(Fitness)
library('plot.matrix')
plot(Decision, breaks=c(0.5, 1.5, 2.5, 3.5), col=c("black", "yellow", "red"))
plot(Fitness)

library(plot3D)

persp3D(z = Fitness, theta = 225, phi = 45,
        xlab = "State (x)", 
        ylab = "Time (t)",
        zlab = "Fitness (F)")


