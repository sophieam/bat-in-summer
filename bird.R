## From table 5.1 (p. 110) in Dynamic state variable models in ecology.

# Maximum energy reserves.
x_max = 2.4          # g

# Discretizised values of x.
x_discreet <- seq(from = 0, to = x_max, by = 0.05)
  
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
Fitness <- array(data = NA, dim = c(length(x_discreet), Time, Days), dimnames = NULL)


# Survival last night.
for (t in Time:1) {
  for (j in 1:length(x_discreet)) {
    if (x_discreet[j] < c_g) {
      Fitness[j, t, Days] <- 0
    } else if ( (x_discreet[j] > c_g) & (x_discreet[j] < c_b) ) {
      Fitness[j, t, Days] <- 1 - p_b
    } else {
      Fitness[j, t, Days] <- 1
    }
  }
}
  
# Other times.
for (d in (Days-1):1) {
  for (t in (Time:1)) {
    for (j in (1:length(x_discreet))) {
      F_i <- vector(mode = 'numeric', length=3)
      
      for (h in H) {
        x_mark <- x_discreet[j] + (1/Time)*e[h] - gamma * (1/Time)*(m_0 + x_discreet[h])
        F_i[h] <- ( 1 - (1/Time)*(mu[h] + lambda*x_discreet[j]))*(Fitness[which(abs(x_discreet-x_mark)==min(abs(x_discreet-x_mark))), t +1, d] ) 
      }
      
      Fitness[j, t, d] <- max(F_i)
    }
  }
}
