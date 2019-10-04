# This model is the model described in chapter 5.1 of Clark and Mangel 1999 
# (Dynamic State Variable Models in Ecology - Methods and Applications).

# Note that the patch indecies used in this script are 1-3 instead of 0-2.

library('plot.matrix')
library('ggplot2')
library('plot3D')

## Parameters of the model - From table 5.1 (p. 110) in Clark and Mangel (1999):

# Maximum fat reserves.
x_max = 2.4              # g

# Discretized values of x.
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

# Number of days
Days = 120

# F(x, t, d) - Fitness (probability of winter survival)
Fit <- array(data = NA, dim = c(length(x_d), Time+1, Days), dimnames = list(
  x_d,
  1:(Time+1),
  1:(Days)) )

# Optimal decision given time and state: 
H <- array(data = NA, dim = c(length(x_d), Time, Days), dimnames = list(
  x_d,
  1:Time,
  1:Days))

# Interpolate a fitness value for a given combination of x, t and d.
interpolate <- function (x, t, d) {
  
  # Returns the index of the closest discrete x value 
  closest_discrete_x <- function(x) {
    # Note: only the first value is returned if there are multiple 
    # equidistant values.
    return(which(abs(x_d - x) == min(abs(x_d - x)))[1])
  }
  
  # Interpolate between two values a and b. dx is a value between 0 and 1.
  linear_interpolation <- function (a, b, dx) {
    return((1-dx)*a+ b*dx)
  }
  
  # No point in doing interpolation if energy reserves are negative.
  # Bird is dead.
  if (x < 0) { return(0) }
  
  # Figure out between which two discretized values of x our x lies.
  closest <- closest_discrete_x(x)
  if (x < x_d[closest]) {
    j1 <- closest -1
    j2 <- closest
  } else if (x > x_d[closest]){
    j1 <- closest
    j2 <- closest +1
  } 
  # Fitness value for x is already present in the matrix. No need to interpolate.
  else { return( Fit[closest, t, d]) }
  
  # Calculate how x is positioned in relation to x_d[j1] and x_d[j2].
  # 0: Closer to x_d[j1]
  # 1: Closer to x_d[j2]
  delta_x <- (x-x_d[j1])/(x_d[j2]-x_d[j1])
  
  # Interpolate.
  return(linear_interpolation(Fit[j1, t, d], Fit[j2, t, d], delta_x))
}
  
# Terminal reward - Equation 5.1
# Note: I use j as the index of an x-value in x_d and x 
# as the value of x_d[j].
for (j in 1:length(x_d)) {
  x <- x_d[j]
  
  if (x < c_g) {
    # Survival is impossible.
    Fit[j, Time +1, Days] <- 0
  } else if (x < c_b) {
    # Survival is possible with probability p_b.
    Fit[j, Time +1, Days] <- (1-p_b)
  } else {
    # Survival is always possible.
    Fit[j, Time +1, Days] <- 1
  }
}

for (d in Days:1) {
  # Terminal fitness has already been calculated, 
  # so don't calculate fitness for T+1 on the last day.
  # Equation 5.2
  if (d != Days) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]

      # Equation 5.3
      # Fitness is (1-p_survive)*fitness the morning after.
      if (x < c_g) {
        # Bird is dead, it just doesn't know it yet.
        Fit[j, Time+1, d] <- 0
      } else if (x < c_b) {
        Fit[j, Time+1, d] <- (1-p_b)*interpolate(x-c_g, 1, d+1)
      } else {
        Fit[j, Time+1, d] <- (1-p_b)*interpolate(x-c_g, 1, d+1) +
                                 p_b*interpolate(x-c_b, 1, d+1)
      }
        
    }
  }
  
  for (t in Time:1) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]
      # Vector for storing fitness values for each of the patches.
      F_i <- vector(mode = 'numeric', length=3)
      
      #Calculate resulting fitness of choosing each patch
      for (i in 1:3) {
        # Equation 5.4
        # Calculating the expected state in the future.
        x_mark <- ( x + (1/Time)*(e[i] - gamma*(m_0+x)) )
        
        # Ensure that x_mark does not exceed x_max.
        x_mark <- min(c(x_max, x_mark))
        
        # Equation 5.3
        # Calculate the expected fitness given patch choice h.
        # i.e: Chance of surviving time expected future fitness in this patch.
        F_i[i] <- (1 - (1/Time)*(mu[i]*exp(lambda*x)) ) *interpolate(x_mark, t+1, d)
      }
      
      # Fitness is the fitness of the patch that maximizes fitness.
      Fit[j, t, d] <- max(F_i)[1]
      
      # Optimal patch choice is the one that maximizes fitness.
      # In cases where more than one patch shares the same fitness, 
      # the first one (i.e. lower risk) is chosen.
      H[j, t, d] <- which(F_i == max(F_i))[1]
    }
  }
}


# Reverse the x-dimension of the array it is ordered from high to low. 
# Nicer when plotting.
Fit <- Fit[length(x_d):1,,]
H <-   H[length(x_d):1,,]

# Plotting the optimal decision at any given time.
# black:  Patch 1 (0)
# yellow: Patch 2 (1)
# red:    Patch 3 (2)
# Should be equivalent to Figure 5.2 in Clark and Mangel (1999).
plot(H[,,Days-20], breaks=c(0.5, 1.5, 2.5, 3.5), col=c("black", "yellow", "red"))

# Fitness 
plot(Fit[,,Days-20])

# Fitness landscape plot.
persp3D(z = Fit[,,Days-20], theta = 135, phi = 45,
        xlab = "State (x)", 
        ylab = "Time (t)",
        zlab = "Fitness (F)")

