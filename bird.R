## From table 5.1 (p. 110) in Dynamic state variable models in ecology.

# Maximum energy reserves.
x_max = 2.4          # g

# Discretizised values of x.
x_discrete <- seq(from = 0, to = x_max, by = 0.05)
  
# Basic daily predation risk patch 1/2.
mu <- c(0, 0.001, 0.005) # /day

# Increase in predation rsik with body mass.
lambda = 0.46            # /g

# Net daily forage inntake for patch 1/2.
e      = c(0, 0.6, 2.0)  # g/day

# Nighttime metabolic cost, good/bad conditions respectively.
c_g    = 0.48            # g
c_b    = 1.20            # g

p_b    = 0.167

# Mass of bird with zero fat reserves.
m_0    = 10          # g

#Metabolic rate.
gamma  = 0.04        # /day

# Time periods per day 
Time = 50

# Days
Days = 120

# Possible actions at each time unit (i.e. patch choice):
H <- 1:3

# X(t, d) - Energy reserves at beginning of period t on day d (g)
X <- matrix(NA, nrow=Time, ncol=Days)

# F(x, t, d) - Fitness (probability of survival)
Fitness <- array(data = NA, dim = c(length(x_discrete), Time+1, Days), dimnames = NULL)

# Optimal decision given time and state: 
Decision <- array(data = 0L, dim = c(length(x_discrete), Time, Days), dimnames = NULL)




cap <- function(minimum, maximum, avalue) {
  if (avalue < minimum) {
    return(minimum)
  } else if (avalue > maximum) {
    return(maximum)
  } else
    return(avalue)
}

closest_discrete_x <- function(anX) {
  return(which(abs(x_discrete - anX) == min(abs(x_discrete - anX))))
}

# Terminal reward.
for (j in 1:length(x_discrete)) {
  if (x_discrete[j] <= c_g) {
    Fitness[j, Time +1, Days] <- 0
  } else if ( (x_discrete[j] > c_g) & (x_discrete[j] < c_b) ) {
    Fitness[j, Time +1, Days] <- 1 - p_b
  } else {
    Fitness[j, Time +1, Days] <- 1
  }
}


for (d in Days:1) {
  for (t in (Time+1):1) {
    for (j in (1:length(x_discrete))) {
      if ( (t == Time +1) & (d < Days) ) {
        if (x_discrete[j] < c_g) { # Not enough to survive good night.
          Fitness[j, t, d] <- 0
        } else if (x_discrete[j] < c_b) { # Not enough to survive bad night.
          Fitness[j, t, d] <- (1-p_b)* Fitness[closest_discrete_x(x_discrete[j]-c_g), 1, d +1]
        } else { 
          Fitness[j, t, d] <- (1-p_b)* Fitness[closest_discrete_x(x_discrete[j]-c_g), 1, d +1] + p_b*Fitness[closest_discrete_x(x_discrete[j] - c_b), 1, d +1]
        }
      } else if (t <= Time) {
        F_i <- vector(mode = 'numeric', length=3)
        
        for (h in H) {
          x_mark <- cap( 0, x_max, x_discrete[j] + (1/Time)*e[h] - gamma * (1/Time)*(m_0 + x_discrete[h]) )
          F_i[h] <- ( 1 - (1/Time)*(mu[h] + lambda*x_discrete[j]))*(Fitness[closest_discrete_x(x_mark), t +1, d] ) 
        }
        
        Fitness[j, t, d] <- max(F_i)
        Decision[j, t, d] <- which(F_i == max(F_i))[1]
      }
      
    }
  }
}
