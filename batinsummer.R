#### STOCHASTIC DYNAMIC MODEL: A BAT IN SUMMER ####
# Building on the code by Myranda and Thorbjørn

#---- Variables (can be changed) ----
#### Required packages and data ####
#install.packages('plot.matrix')
library(plot.matrix)
#install.packages('plot3D')
library(plot3D)
library(ggplot2)
#install.packages('tidyr')
library(tidyr)
#install.packages('viridis')
library(viridis) # colourblind palette

predationdata<-read.csv("Data/PredationAndForagingEquation.csv")

#### State (fat reserves) ####
#the original bird values are commented out
mass_zero_fat    <- 7 #10    # Mass of bat with zero fat reserves (g)
fat_max  <- 7 #2.4   # Maximum fat reserves (g) (NOT mass of bat with max fat reserves!)
fat_discretized    <- seq(from = 0, to = fat_max, length.out=100)   # Discretized values of fat_state
predation_risk_increase <- 0.46     # Increase in predation risk due to body mass (/g fat reserves)

#### Fitness (time, probability of good or bad night) ####

# A season is divided in days, themselves divided in timesteps.
nb_timesteps <- 70    # Number of time periods per day (where nb_timesteps+1 is the beginning of daytime)
nb_days <- 153   # Number of days in summer (where nb_days+1 is the start of winter)
nb_hours <- 24 #Number of hours in a day (e.g. 12 or 24). Used later to calculate mebabolic cost per timestep instead of per hour.

#Metabolic costs
cost_flight_hourly <- 0.615  #Hourly cost of flight in g

# Environmental stochasticity (bad or good day) affects daily metabolic costs.
cost_good_day  <- 0.48  # Metabolic cost of a good day (g)
cost_bad_day  <- 1.20  # Metabolic cost of a bad day (g)
probability_bad_day  <- 0.167 # Probability of a bad day

#### Temperature, predation risk, prey availability and metabolism functions ####

#External temperature in degrees celsius, affects prey availability
get_temperature <- function(time_current){
  temperature_current <- 17.45 + 1.101*time_current - 0.05206*time_current^2 + 0.0007522*time_current^3 - 0.000003164*time_current^4
  return(temperature_current)
}
curve(expr = get_temperature, from = 1, to = nb_timesteps)

#Roost temperature in degrees celsius, affects metabolism for patches 1 and 2
get_temperature_roost <- function(time_current){
  temperature_roost_current <- 22.92 + 2.180*time_current - 0.1382*time_current^2 + 0.002574*time_current^3 - 0.00001465*time_current^4
  return(temperature_roost_current)
}
curve(expr = get_temperature_roost, from = 1, to = nb_timesteps)

#Prey availability
get_prey <- function(temperature_current){
    reward_prey_current <- 1 / (1 + exp( -0.524* (temperature_current - 11)))
    return(reward_prey_current)
 }

curve(expr = get_prey, from = 0, to = 20)

#Predation
#we are using existing data for this, rather than calculating it

risk_predation <- predationdata[,5]
plot(risk_predation)

get_predation <- function(time_current){
  predation_current <- risk_predation[time_current]
  return(predation_current)
}

#...but we could make a normal distribution instead? here is an attempt although too leptokurtic
# get_predation <- function(time_current){
#   mu <- nb_timesteps/2
#   sigm <- 5
#   e <- exp(1)
#   x <- time_current
#   predation_current <- 1-((1/(sigm*(sqrt(2*pi))))*e^(-((x-mu)^2)/(2*sigm^2))*12.5)
#   return(predation_current)
# }

test <- curve(expr = get_predation, from = 1, to = nb_timesteps)

#Metabolism, transform the hourly cost into a per timestep cost
standardize_metabo_cost <- function(cost_hourly){
  fraction <- nb_hours/nb_timesteps
  cost <- cost_hourly*fraction
  return(cost)
}

#Calculate metabolism per timestep, torpor (patch 1)

#Calculate metabolism per timestep, resting (patch 2)

#Calculate metabolism per timestep, foraging (patch 3)
cost_flight <- standardize_metabo_cost(cost_flight_hourly)

#make dataframe that stores timestep and corresponding prey availability
patch1 <- rep(0, times = nb_timesteps)
patch2 <- rep(0, times = nb_timesteps)
patch3 <- seq(1, nb_timesteps, by=1)
patch3 <- get_prey(get_temperature(patch3))
foraging_benefit <- data.frame(patch1, patch2, patch3)
rm(patch1)
rm(patch2)
rm(patch3)

#make dataframe that stores timestep and corresponding predation risk
patch1 <- rep(0, times = nb_timesteps)
patch2 <- rep(0, times = nb_timesteps)
patch3 <- seq(1, nb_timesteps, by=1)
patch3 <- get_predation(patch3)
predation_risk <- data.frame(patch1, patch2, patch3)
rm(patch1)
rm(patch2)
rm(patch3)

#make one for metabolic cost per patch; remember that patch 3 is foraging outdoors and uses external temperature
patch1 <- rep(0, times = nb_timesteps)
patch2 <- rep(0, times = nb_timesteps)
patch3 <- rep(cost_flight, times = nb_timesteps)
metabolic_cost_all <- data.frame(patch1, patch2, patch3)
rm(patch1)
rm(patch2)
rm(patch3)

#---- Empty arrays (do not change) ----
# We make an empty array to hold the fitness values  
fitness <- array(data = NA, dim = c(length(fat_discretized), nb_timesteps+1, nb_days), 
                 dimnames = list(fat_discretized, 1:(nb_timesteps+1), 1:(nb_days)) )

# And for the decision loop (patch choice)
patch_choice <- array(data = NA, dim = c(length(fat_discretized), nb_timesteps, nb_days), dimnames = list(
  fat_discretized, 1:nb_timesteps, 1:nb_days))


#---- Functions (do not change) ----
#### Calculating terminal fitness at the end of each day####
# Terminal fitness is the probability of surviving the last day (i.e. day_current+1)

#This function calculates the survival probability at a given state, when at the end of the night
#Is the individual fat enough to survive a bad day, a good day only, or not even a good day?
  calculate_survival_proba <- function(fat_state) {
    if (fat_state < cost_good_day) {
      survival_proba <- 0
    } else if (fat_state < cost_bad_day) {
      survival_proba <- (1-probability_bad_day)
    } else {
      survival_proba <- 1
    }
    return(survival_proba)
  }
  


# Calculate terminal fitness F(fat_state, timestep_current+1, day_current), for every day, where fat_state=state
# This for loop can be replaced by apply for cleaner code
for (j in 1:length(fat_discretized)) {
  fat_state <- fat_discretized[j]
    fitness[j, nb_timesteps+1, nb_days] <- calculate_survival_proba(fat_state)
}
rm(j)

  
# Visualise the state and corresponding terminal fitness
#as.matrix(fitness[, nb_timesteps+1, nb_days])

#---- Interpolating the discrete state variable ----
# The computer discretizes the state variable 
# but in reality energetic reserves is a continuous variable. 
# We overcome this using interpolation (see C&M 2.1)

interpolate <- function (fat_state, timestep_current, day_current) {
  
  # Returns the index of the closest discrete fat_state value 
  closest_discrete_x <- function(fat_state) {
    # Note: only the first value is returned if there are multiple 
    # equidistant values.
    return(which(abs(fat_discretized - fat_state) == min(abs(fat_discretized - fat_state)))[1])
  }
  
  # Interpolate between two values a and b. dx is a value between 0 and 1.
  linear_interpolation <- function (a, b, dx) {
    return((1-dx)*a+ b*dx)
  }
  
  # No point in doing interpolation if energy reserves are negative.
  # Bird is dead.
  if (fat_state < 0) { return(0) }
  
  # Figure out between which two discretized values of fat_state our fat_state lies.
  closest <- closest_discrete_x(fat_state)
  if (fat_state < fat_discretized[closest]) {
    j1 <- closest -1
    j2 <- closest
  } else if (fat_state > fat_discretized[closest]){
    j1 <- closest
    j2 <- closest +1
  } 
  # Fitness value for fat_state is already present in the matrix. No need to interpolate.
  else { return( fitness[closest, timestep_current, day_current]) }
  
  # Calculate how fat_state is positioned in relation to fat_discretized[j1] and fat_discretized[j2].
  # 0: Closer to fat_discretized[j1]
  # 1: Closer to fat_discretized[j2]
  delta_x <- (fat_state-fat_discretized[j1])/(fat_discretized[j2]-fat_discretized[j1])
  
  # Interpolate.
  return(linear_interpolation(fitness[j1, timestep_current, day_current], fitness[j2, timestep_current, day_current], delta_x))
}

# ---- Backwards Iteration ---- 



# Iterate backwards across days at nb_timesteps+1
for (day_current in nb_days:1) {
  # Terminal fitness has already been calculated, 
  # so don't calculate fitness for T+1 on the last day.
  if (day_current != nb_days) {
    for (j in 1:length(fat_discretized)) {
      fat_state <- fat_discretized[j]
      
      # Equation 5.2
      # Notice we use the interpolation function here
      # If a bird can't survive a good night... 
      if (fat_state < cost_good_day) {
        fitness[j, nb_timesteps+1, day_current] <- 0      # Bird is dead, it just doesn't know it yet.
        
        # If a bird can survive a good night but not a bad night
        # its probability of survival = probability of a good night 
        # times its state (reward) after you subtract the overnight energetic cost
      } else if (fat_state < cost_bad_day) {
        fitness[j, nb_timesteps+1, day_current] <- (1-probability_bad_day)*interpolate(fat_state-cost_good_day, 1, day_current+1)  
        
        # Otherwise it can survive a bad night, in which case
        # fitness = prob. of a good night*state after a good night + 
        # prob. of a bad night * state after a bad night
        # this follows the same logic as the HK model
      } else {
        fitness[j, nb_timesteps+1, day_current] <- (1-probability_bad_day)*interpolate(fat_state-cost_good_day, 1, day_current+1) +
          probability_bad_day*interpolate(fat_state-cost_bad_day, 1, day_current+1)  #prob of good night*state after good night plus same for bad night
      }
    } # end if-else loop
  } # end j loop
  
  # Iterate backwards from max nb_timesteps for each day
  for (timestep_current in nb_timesteps:1) {
    for (j in 1:length(fat_discretized)) {
      fat_state <- fat_discretized[j]
      # Vector for storing fitness values for each of the three patches.
      fitness_all_choices <- vector(mode = 'numeric', length=3)
      
      # Calculate resulting fitness of choosing each patch
      for (i in 1:3) {
        # Equation 5.4
        # Calculating the expected state in the future.
        # New calculation, not divided per timesteps as the values already are per timestep:
        fat_expected <- ( fat_state + foraging_benefit[timestep_current,i] - metabolic_cost_all[timestep_current,i]*(mass_zero_fat+fat_state))
        # instead of the old one:
        #fat_expected <- ( fat_state + (1/nb_timesteps)*(foraging_benefit[timestep_current,i] - metabolic_rate*(mass_zero_fat+fat_state)) )
        
        # Ensure that fat_expected does not exceed fat_max.
        fat_expected <- min(c(fat_max, fat_expected))
        
        # Plug fat_expected (5.4) into equation 5.3
        # Calculate the expected fitness given patch choice h.
        # And this is where we insert the predation value instead of predation_risk.
        # i.e: Chance of surviving time expected future fitness in this patch.
        fitness_all_choices[i] <- (1 - (predation_risk[timestep_current,i]*exp(predation_risk_increase*fat_state)) ) *interpolate(fat_expected, timestep_current+1, day_current)
      } # end i loop 
      
      # Which patch choice maximizes fitness?
      fitness[j, timestep_current, day_current] <- max(fitness_all_choices)[1]
      
      # Optimal patch choice is the one that maximizes fitness.
      # In cases where more than one patch shares the same fitness, 
      # the first one (i.e. lower risk) is chosen.
      patch_choice[j, timestep_current, day_current] <- which(fitness_all_choices == max(fitness_all_choices))[1]
      
    } # end j loop 
  } # end timestep_current loop
} # end day_current loop

# We can explore the decision matrix at different days, similar to Fig 5.2
# Keeping in mind that our output shows a reverse order on the y-axis
#patch_choice[,,100] 
#END OF BACKWARDS ITERATION#

#--- Plots ----

####Fitness function####
# Fitness function at two different times on day nb_days-20, timestep_current=25 and timestep_current=49 respectively.
# This figure should be equivalent to 5.1 in Clark and Mangel (1999).

# Create an empty plot for plotting lines.
plot(NA, type="n", 
     xlab="Fat reserves (g)",
     ylab="Fitness, F[fat_state, timestep_current]", 
     xlim=c(0, fat_max), ylim=c(0, 1))

# Plotting the lines.
lines(fat_discretized, as.vector(fitness[,25,nb_days-20]), col = "black", lty = 3)
lines(fat_discretized, as.vector(fitness[,49,nb_days-20]), col = "black", lty = 1)


# Reverse the x-dimension of the array so it is ordered from high to low. 
# Nicer when plotting.
fitness_reversed <- fitness[length(fat_discretized):1,,]
patch_choice_reversed <-   patch_choice  [length(fat_discretized):1,,]
####Optimal decision plot####
# Plotting the optimal decision at any given time.
# black:  Patch 1 (0)
# yellow: Patch 2 (1)
# red:    Patch 3 (2)
# Should be equivalent to Figure 5.2 in Clark and Mangel (1999).
plot(patch_choice_reversed[,,nb_days-20], breaks=c(0.5, 1.5, 2.5, 3.5), col=c("black", "yellow", "red"),
     xlab = "nb_timesteps of day",
     ylab = "Fat reserves",
)

#for more flexibility I will reshape the decision matrix into a long format dataframe. fat reserves, time of day and decision as columns.
#first, let's convert the array into a dataframe
optimaldecision <- as.data.frame(patch_choice_reversed)
optimaldecision <- cbind(rownames(optimaldecision), data.frame(optimaldecision, row.names=NULL, check.names=FALSE, stringsAsFactors=FALSE))
colnames(optimaldecision)[1] = 'fatreserves'
#now we do the reshaping using the tidy package (part of the tidyverse)
optimaldecision <- gather(optimaldecision,time,decision,-fatreserves)
#we can now plot our values



####Fitness plot####
plot(fitness_reversed[,,nb_days-20],
     xlab = "Time of day",
     ylab = "Fat reserves")

#### Fitness landscape plot####
persp3D(z = fitness_reversed[,,nb_days-20], theta = 135, phi = 45,
        xlab = "State (fat_state)", 
        ylab = "Timestep (timestep_current)",
        zlab = "Fitness (fitness)")




###### FORWARD ITERATION ######
# The below section implements forward iteration of the above model (i.e. 
# looking at the fate of individuals in the model).

# Number of days to forward iterate for.
nb_days_forward  <- 120

# nb_timesteps/day is given from the model above.

# Initial states for individual (A vector of indices in fat_discretized).
state_init_forward <- 1 # seq(from=1, to = length(fat_discretized), by = 20)

# Number of individuals to iterate for each state_init_forward_discretized. Total number of individuals
# is length(state_init_forward_discretized)*nb_indiv_forward:
nb_indiv_forward <- 500

# state_forward_all[state_init_forward, N, D, T]:
# State of individual N with initial state state_init_forward, at day D and time T.

state_forward_all <- array(data = NA, dim = c(length(state_init_forward), nb_indiv_forward, nb_days_forward*nb_timesteps +1), dimnames = list(
  fat_discretized[state_init_forward],
  1:nb_indiv_forward,
  1:(nb_days_forward*nb_timesteps+1)
))

for (j in 1:length(state_init_forward)) {
  state_init_forward_discretized <- fat_discretized[state_init_forward[j]]
  for (n in 1:nb_indiv_forward) {
    # Set the initial state.
    state_forward_all[j, n, 1] <- state_init_forward_discretized
    
    # Iterate through the time.
    for (z in 1:(nb_days_forward*nb_timesteps) ) {
      timestep_current = (z-1) %%  nb_timesteps  +1 # nb_timesteps of day
      day_current = (z-1) %/% nb_timesteps  +1 # Day
      # The current state.
      fat_state <- state_forward_all[j, n, z]
      if (fat_state < 0) {
        # Bird is dead. It will remain dead.
        state_new_forward <- fat_state
      } else {
        # The best decision given the current state (Found by rounding fat_state to the closest item
        # in fat_discretized).
        # TODO: interpolate the decision between the two closest optimal decisions.
        h <- patch_choice[which(abs(fat_discretized - fat_state) == min(abs(fat_discretized - fat_state)))[1], timestep_current, day_current]
        
        if (runif(1) <= predation_risk[timestep_current,h]*exp(predation_risk_increase*fat_state)) { # Predation risk.
          # Negative fat_state = dead.
          state_new_forward <- -1
        } else {
          metabolism <- metabolic_cost_all[timestep_current,h]*(mass_zero_fat+fat_state)
          foraging   <- foraging_benefit[timestep_current,h]
          
          state_new_forward <- fat_state + foraging - metabolism
        }
      }
      
      if (timestep_current == nb_timesteps) {
        # If it is the end of the day, nightly costs need to be applied as well.
        if (runif(1) < probability_bad_day) {
          # Bad night.
          state_new_forward <- state_new_forward - cost_bad_day
        } else {
          # Good night.
          state_new_forward <- state_new_forward - cost_good_day
        }
        
        # Set the state for the start of the following day.
      } 
      
      # Set new state.
      state_forward_all[j, n, z+1] <- state_new_forward
    }
  }
}

rm(h)
rm(i)
rm(j)
rm(z)

#END OF FORWARD ITERATION#

####Forward iteration survival plot####
## The below plots multiple individuals on the same plot.
# Red vertical lines:     Separates days.
# Light grey solid line:  Fat reserves necessary for surviving a bad night.
# Light grey dotted line: Fat reserves necessary for surviving a good night.

# Create an empty plot, with the necessary xlim and ylim.
plot(NA, type="n", 
     xlab="Timestep",
     ylab="Fat reserves (g)", xlim=c(1, (nb_days_forward*nb_timesteps +1)), ylim=c(0, fat_max))


# Horizontal lines to indicate cost of good and bad night respectively.
abline(h=cost_good_day, col = "blue", lty = 3)
abline(h=cost_bad_day, col = "blue", lty = 1)

# Plot each individual for the chosen state_init_forward_discretized.
for (n in 1:nb_indiv_forward) {
  lines(1:(nb_days_forward*nb_timesteps +1), state_forward_all[1, n, ], col = "black")
}

# Plot vertical lines to indicate position of the night.
for (i in 1:nb_days_forward) {
  abline(v=i*nb_timesteps, col = "red", lty = 3)
}

# Calculate the proportion of individuals still alive at the start 
# of each day.

alive = data.frame(day = NULL, n_alive = NULL, p_alive = NULL)
for (day_current in 1:(nb_days_forward+1)) {
  
  tmp <- data.frame(day = day_current,
                    n_alive = length(which(as.vector(state_forward_all[1, , (day_current-1)*nb_timesteps +1]) >= 0))
                    
                    
  )
  tmp$p_alive <- tmp$n_alive/nb_indiv_forward
  
  alive <- rbind(alive, tmp)
}

rm(n)
rm(i)

# Plotting the proportion of individuals that are still alive 
# at the beginning of each day.
plot(NA, type="n", 
     xlab="Day",
     ylab="Proportion alive", xlim=c(1, nb_days_forward+1), ylim=c(0, 1))

lines(alive$day, alive$p_alive)

###END OF SCRIPT###
