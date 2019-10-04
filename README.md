# Small bird in winter

This repository consists of r-implementations for different stochastic dynamic models described in chapter 5 (Winter survival strategies) of Clark and Mangel (2000)<sup>1</sup>

# 5.1: Foraging and fattening strategies for willow tits
This model considers what a bird should do during the winter to survive the night. The model considers multiple (120) days, each divided into 50 units of time. Each unit of time can be spent  in one of three patches:
* Patch 1: Rest - No predation risk and no foraging.
* Patch 2: Low predation risk (0.001/day) and low foraging rate (0.6g/day).
* Patch 3: High predation risk (0.005/day) and high foraging rate (2.0g/day).

The predation risk additionally increases with fat reserves, with a rate of increase of 0.46/g.

The energy reserves required at the end of the day to survive is stochastic. A bad night occurs with probability 0.167 and requires 1.2g of fat reserves. Otherwise it is a good night which requires 0.48g of fat reserves.

This model can be found in **5_1.R**.

# 5.2: Social interaction

# 5.4: Hoarding behaviour

# 5.5: Long-term hoarding



<sup>1</sup>: Clark, C. W. & Mangel, M. *Dynamic state variable models in ecology: methods and applications*. (Oxford University Press, 2000).

