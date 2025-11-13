rm(list = ls())

library(dplyr)
library(broom)

source("~/functions.R")

# Example --------------------------------------------------------------
Bc  <- 1.08e-07   
Ec  <- 3.7e-7         
m_c <- 1e-12         
phi <- 1.5
gamma <- 1.5
eta <- 1
dt <- 0.1

cv <- 0.2
n <- 100  # number of individuals
Nmax  <- 5e12
Nmax_mean <- rep(c(Nmax), each = n)
yrs <- rep(1, each = n)
n_individuals <- length(yrs)
 
res <- do.call(rbind, lapply(1:n_individuals, function(i){
    SOGM( Nmax_sd=cv, Nmax_mean=Nmax_mean[i] , 
                     years=yrs[i] ,  gamma=gamma, phi=phi,
                     dt = dt,
                     id = i, plot=F)
  }))

plot(log(res$B_mu)~log(res$mass))
abline(lm(log(res$B_mu)~log(res$mass) ))
summary(lm(log(res$B_mu)~log(res$mass)))

# comparison ------------------------------------------------------------------------
# Adult metabolic rate
mean(res$B_mu)
B_adult(Nmax, Bc, Ec, gamma, phi, eta, dt)

# Adult metabolic fluctuations
mean(res$B_lap_sd)
sigma_r_mc(Nmax, Bc, Ec, gamma, phi, eta, dt)

# Adult mass (kg)
mean(res$mass)
m_adult(Nmax, Bc, Ec, phi, m_c, eta, dt)                         

# Species-specific beta
coef(lm(log(res$B_mu)~log(res$mass) ))[2]
beta_species(Nmax, Bc, Ec, gamma, phi, eta, dt)                  

## SWEEP across phi and gamma -------------------------------------------------------
gamma_values <- seq(0.00001, 2, length.out = 500)  
phi_values <- seq(0.00001, 2, length.out = 500)  
Nmax_mean <- c(1e10, 5e12, 7e13, 5e14, 5e15, 5e16)
id <- expand.grid(Nmax_mean, gamma_values, phi_values )
colnames(id) <- c("Nmax", "gamma", "phi")
Bc  <- 1.08e-07   
Ec  <- 3.7e-7     
m_c <- 1e-12          
eta <- 0.7
dt <- 0.1

B <- M <- beta <- sd_lap <- list()
for(i in 1:nrow(id)){
       print(i)
       vars <- as.numeric(id[i,])
       B[[i]] <- B_adult(Nmax=vars[1], Bc, Ec, gamma=vars[2], phi=vars[3], eta, dt)  
       M[[i]] <- m_adult(Nmax=vars[1], Bc, Ec, phi=vars[3], m_c, eta, dt)                         
       beta[[i]] <- beta_species(Nmax=vars[1], Bc, Ec, gamma=vars[2], phi=vars[3], eta, dt)  
       sd_lap[[i]] <- sigma_r_mc(Nmax=vars[1], Bc, Ec, gamma=vars[2], phi=vars[3], eta, dt)
}

d <- data.frame( id, B=unlist(B), mass=unlist(M), beta=unlist(beta), 
                 sd_lap=unlist(sd_lap))

d <- d %>%
  mutate(
    animal = case_when(
      Nmax == 1e10  ~ "mouse",
      Nmax == 5e12  ~ "cat",
      Nmax == 7e13  ~ "human",
      Nmax == 5e14  ~ "cow",
      Nmax == 5e15  ~ "elephant",
      Nmax == 5e16  ~ "whale",
      TRUE ~ NA_character_
    ),
    animal = factor(animal, levels = c("mouse", "cat", "human", "cow", "elephant", "whale"))
  )

save(d, file ="~/analytical.RData")

