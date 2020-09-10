# STOCHASTIC DYNAMIC MODEL: A BAT IN SUMMER

# Building up on the code by Myranda and Thorbjørn

#install.packages('plot.matrix')
library(plot.matrix)
#install.packages('plot3D')
library(plot3D)

#---- State ---- 
# Defining state variables. Fat reserves are important as they determine hunting efficiency.
#There is a negative feedback loop here, as fatter bats will be able to hunt less.

m_0    <- 10    # Mass of bat with zero fat reserves (g)
x_max  <- 6   # Maximum fat reserves (g)
                
x_d    <- seq(from = 0, to = x_max, 
              length.out=100)   # Discretized values of x


#---- Actions ---- 
# Defining decision-related variables.
# Bats can choose between a maximum of three actions: forage, roost (no torpor), roost (torpor).
# These are coded as i = [1, 2, 3] 
# Each patch (H[i]) has a different energy cost denoted by mu (due to metabolism).
# Each patch also has a different energetic gain denoted by e (in practice only the foraging patch leads to an energetic gain).
# The bird code included predation, but in bats energy cost is likely to be a stronger influence on decisions.

mu <- c(0.0005, 0.0002, 0.0001)  # Basic daily metabolic cost per patch (/day) 
e  <- c(0, 0.6, 0.6)      # Net daily forage intake for patch 1/2 (g/day)

lambda <- 0.06     # Decrease in hunting efficiency with body mass (/g fat reserves)

# The metabolic costs per day are also mass-dependent, equal to gamma * (m0 * x), 
# where m0 is the mass of the bat with zero fat reserves (defined above)

gamma  <- 0.04     # Metabolic rate (g/day)



#---- Fitness ---- 
# Summer consists of 120 days (Days) and each day is divided into 50 time steps (Time). 
# We can account for this with a nested loop (see backwards iteration below). 
# Environmental stochasticity arises from daily metabolic costs. 
# A bad day has a higher metabolic costs (c_b) and occurs with probability p_b
# A good day has a lower metabolic cost (c_g) and occurs with probability 1-p_b

Time <- 50    # Number of time periods per day (where time+1 is the beginning of daytime)
Days <- 120   # Number of days in summer (where days+1 is the start of winter)
c_g  <- 0.48  # Metabolic cost of a good day (g)
c_b  <- 1.20  # Metabolic cost of a bad day (g)
p_b  <- 0.167 # Probability of a bad day

# We generate a calendar array that will code for days when torpor is or is not available.
# Make vector of length Days, fill it with 3s
calendar <- c(rep(3, times=Days))
# Then overwrite these for days 70 to 90
# We create the following variables to make it easier to change these later
no_torpor_start <- 70
no_torpor_end <- 80
calendar[no_torpor_start:no_torpor_end] <- c(2)

# Recall that terminal fitness is the probability of surviving the summer (i.e. d+1)
# We can calculate it using equation 5.1 in C&M

# First, an empty array to hold fitness values from the loop
Fit<- array(data = NA, dim = c(length(x_d), Time+1, Days), 
             dimnames = list(x_d, 1:(Time+1), 1:(Days)) )

# Calculate terminal fitness F(x, t+1, d) where x=state 
for (j in 1:length(x_d)) {
  # Note: I use j as the index of an x-value in x_d and x 
  # as the value of x_d[j].
  x <- x_d[j]
  # If your state does not allow you to survive a good day,
  # survival is impossible and your fitness is zero. 
  if (x < c_g) {
    Fit[j, Time+1, Days] <- 0
  # If your state allows you to survive a good day, but not a bad one,
  # you survive with probability p_g=1-p_b (the probability of a good night). 
  } else if (x < c_b) {
    Fit[j, Time+1, Days] <- (1-p_b)
  # If your state allows you to survive a bad day, you will always survive
  # i.e. your probability of survival is 1. 
  } else {
    Fit[j, Time+1, Days] <- 1
  } # end if-else loop
} # end for loop

as.matrix(Fit[, Time+1, Days])
# The left column shows state, the right column shows the terminal fitness we just calculated. 
# To survive a good day requires 0.48 g of fat reserves.
# When x < 0.48 the bat dies (it cannot survive the best-case scenario). 
# To survive a bad day requires 1.2 g of fat reserves. 
# When 0.48 < x < 1.2 the probability of survival =  1-p_b = 0.833
# When x > 1.2 (the cost of a bad day), the bird always survives (probability =1)

#---- Discrete state variable ----
# The computer discretizes the state variable 
# but in reality energetic reserves is a continuous variable. 
# We overcome this using interpolation (see C&M 2.1)

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

# ---- Backwards Iteration ---- 

# We already have a blank array for the fitness loop
# Create one for the decision loop (patch choice, H)
H <- array(data = NA, dim = c(length(x_d), Time, Days), dimnames = list(
  x_d, 1:Time, 1:Days))

# Iterate backwards across days at Time+1
for (d in Days:1) {
  # Terminal fitness has already been calculated, 
  # so don't calculate fitness for T+1 on the last day.
  if (d != Days) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]
      
      # Equation 5.2
      # Notice we use the interpolation function here
      # If a bird can't survive a good night... 
      if (x < c_g) {
        Fit[j, Time+1, d] <- 0      # Bird is dead, it just doesn't know it yet.
        
        # If a bird can survive a good night but not a bad night
        # its probability of survival = probability of a good night 
        # times its state (reward) after you subtract the overnight energetic cost
      } else if (x < c_b) {
        Fit[j, Time+1, d] <- (1-p_b)*interpolate(x-c_g, 1, d+1)  
        
        # Otherwise is can survive a bad night, in which case
        # fitness = prob. of a good night*state after a good night + 
        # prob. of a bad night * state after a bad night
        # this follows the same logic as the HK model
      } else {
        Fit[j, Time+1, d] <- (1-p_b)*interpolate(x-c_g, 1, d+1) +
          p_b*interpolate(x-c_b, 1, d+1)  #prob of good night*state after good night plus same for bad night
      }
    } # end if-else loop
  } # end j loop
  
  # Iterate backwards from max Time for each day
  for (t in Time:1) {
    for (j in 1:length(x_d)) {
      x <- x_d[j]
      # Vector for storing fitness values for each of the three patches.
      F_i <- vector(mode = 'numeric', length=3)
      
      # Calculate resulting fitness of choosing each patch
      # Equation 5.4
      
#Here, add a check that looks up the day of the year in the calendar dataframe or array. If... else...
#Based on this we choose between two loops: with (1:3) or without (1:2) torpor. 
    
#first we need to grab the number of the loop, should be =d
# then we return the corresponding calendar value
currentday<-calendar[d]
#then we can do our if... else... with the torpor loop under if
# if calendar[d]=3, run prewritten loop
#elseif calendar[d]=2, run shortened options
#else return error
      
      
      for (i in 1:currentday) {
        # Equation 5.4
        # Calculating the expected state in the future.
        x_mark <- ( x + (1/Time)*(e[i] - gamma*(m_0+x)) )
        
        # Ensure that x_mark does not exceed x_max.
        x_mark <- min(c(x_max, x_mark))
        
        # Plug x_mark (5.4) into equation 5.3
        # Calculate the expected fitness given patch choice h.
        # i.e: Chance of surviving time expected future fitness in this patch.
        F_i[i] <- (1 - (1/Time)*(mu[i]*exp(lambda*x)) ) *interpolate(x_mark, t+1, d)
      } # end i loop 
      
      # Which patch choice maximizes fitness?
      Fit[j, t, d] <- max(F_i)[1]
      
      # Optimal patch choice is the one that maximizes fitness.
      # In cases where more than one patch shares the same fitness, 
      # the first one (i.e. lower risk) is chosen.
      H[j, t, d] <- which(F_i == max(F_i))[1]
      
    } # end j loop 
  } # end t loop
} # end d loop

# We can explore the decision matrix at different days, similar to Fig 5.2
# Keeping in mind that our output shows a reverse order on the y-axis
H[,,100] 


#--- Plots! ----

# Would be nice to add a color code for days with and without torpor

# Fitness function at two different times on day Days-20, t=25 and t=49 respectively.
# This figure should be equivalent to 5.1 in Clark and Mangel (1999).

# Create an empty plot for plotting lines.
plot(NA, type="n", 
     xlab="Fat reserves (g)",
     ylab="Fitness, F[x, t]", 
     xlim=c(0, x_max), ylim=c(0, 1))

# Plotting the lines.
lines(x_d, as.vector(Fit[,25,Days-20]), col = "black", lty = 3)
lines(x_d, as.vector(Fit[,49,Days-20]), col = "black", lty = 1)


# Reverse the x-dimension of the array so it is ordered from high to low. 
# Nicer when plotting.
Fit.rev <- Fit[length(x_d):1,,]
H.rev <-   H  [length(x_d):1,,]

# Plotting the optimal decision at any given time.
# black:  Patch 1 (0)
# yellow: Patch 2 (1)
# red:    Patch 3 (2)
# Should be equivalent to Figure 5.2 in Clark and Mangel (1999).
plot(H.rev[,,Days-20], breaks=c(0.5, 1.5, 2.5, 3.5), col=c("black", "yellow", "red"),
     xlab = "Time of day",
     ylab = "Fat reserves")

# Fitness 
plot(Fit.rev[,,Days-20],
     xlab = "Time of day",
     ylab = "Fat reserves")

# Fitness landscape plot.
persp3D(z = Fit.rev[,,Days-20], theta = 135, phi = 45,
        xlab = "State (x)", 
        ylab = "Time (t)",
        zlab = "Fitness (F)")



###### FORWARD ITERATION ######
# The below section implements forward iteration of the above model (i.e. 
# looking at the fate of individuals in the model).

# Number of days to forward iterate for.
Days  <- 120

# Time/day is given from the model above.

# Initial states for individual (A vector of indecies in x_d).
j_0 <- 1 # seq(from=1, to = length(x_d), by = 20)

# Number of individuals to iterate for each x_0. Total number of individuals
# is length(x_0)*N.ind:
N.ind <- 500

# X[j_0, N, D, T]:
# State of individual N with initial state j_0, at day D and time T.

X <- array(data = NA, dim = c(length(j_0), N.ind, Days*Time +1), dimnames = list(
  x_d[j_0],
  1:N.ind,
  1:(Days*Time+1)
))

for (j in 1:length(j_0)) {
  x_0 <- x_d[j_0[j]]
  for (n in 1:N.ind) {
    # Set the initial state.
    X[j, n, 1] <- x_0
    
    # Iterate through the time.
    for (z in 1:(Days*Time) ) {
      t = (z-1) %%  Time  +1 # Time of day
      d = (z-1) %/% Time  +1 # Day
      # The current state.
      x <- X[j, n, z]
      if (x < 0) {
        # Bird is dead. It will remain dead.
        x_new <- x
      } else {
        # The best decision given the current state (Found by rounding x to the closest item
        # in x_d).
        # TODO: interpolate the decision between the two closest optimal decisions.
        h <- H[which(abs(x_d - x) == min(abs(x_d - x)))[1], t, d]
        
        if (runif(1) <= (1/Time)*mu[h]*exp(lambda*x)) { # Predation risk.
          # Negative x = dead.
          x_new <- -1
        } else {
          metabolism <- (1/Time)*gamma*(m_0+x)
          foraging   <- (1/Time)*e[h]
          
          x_new <- x + foraging - metabolism
        }
      }
      
      if (t == Time) {
        # If it is the end of the day, nightly costs need to be applied as well.
        if (runif(1) < p_b) {
          # Bad night.
          x_new <- x_new - c_b
        } else {
          # Good night.
          x_new <- x_new - c_g
        }
        
        # Set the state for the start of the following day.
      } 
      
      # Set new state.
      X[j, n, z+1] <- x_new
    }
  }
}

## The below plots multiple individuals on the same plot.
# Red vertical lines:     Separates days.
# Light grey solid line:  Fat reserves necessary for surviving a bad night.
# Light grey dotted line: Fat reserves necessary for surviving a good night.

# Create an empty plot, with the necessary xlim and ylim.
plot(NA, type="n", 
     xlab="Time",
     ylab="Fat reserves (g)", xlim=c(1, (Days*Time +1)), ylim=c(0, x_max))


# Horizontal lines to indicate cost of good and bad night respectively.
abline(h=c_g, col = "blue", lty = 3)
abline(h=c_b, col = "blue", lty = 1)

# Plot each individual for the chosen x_0.
for (n in 1:N.ind) {
  lines(1:(Days*Time +1), X[1, n, ], col = "black")
}

# Plot vertical lines to indicate position of the night.
for (i in 1:Days) {
  abline(v=i*Time, col = "red", lty = 3)
}

# Calculate the proportion of individuals still alive at the start 
# of each day.

alive = data.frame(day = NULL, n_alive = NULL, p_alive = NULL)
for (d in 1:(Days+1)) {
  
  tmp <- data.frame(day = d,
                    n_alive = length(which(as.vector(X[1, , (d-1)*Time +1]) >= 0))
                    
                    
  )
  tmp$p_alive <- tmp$n_alive/N.ind
  
  alive <- rbind(alive, tmp)
}

# Plotting the proportion of individuals that are still alive 
# at the beginning of each day.
plot(NA, type="n", 
     xlab="Day",
     ylab="Proportion alive", xlim=c(1, Days+1), ylim=c(0, 1))

lines(alive$day, alive$p_alive)

####END####





