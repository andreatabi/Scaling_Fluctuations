library(dplyr)
library(extraDistr)
library(moments)
library(VGAM)

logLik_gaussian <- function(x) {
  mu <- mean(x)
  sigma <- sd(x)
  sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
}

logLik_laplace <- function(x) {
  mu <- median(x)
  b <- mean(abs(x - mu))
  sum(VGAM::dlaplace(x, location = mu, scale = b, log = TRUE))
}

# Stochastic ontogentic growth model 
SOGM <- function(Nmax_sd, Nmax_mean, id, years, gamma, phi, eta=0.7, dt=0.1, plot = FALSE) {
  message("Simulating ID: ", id)
  
  # Constants 
  Bc  <- 3e-11 * 3600  
  Ec  <- 3.7e-7         
  mc  <- 1e-12         
  Td <- years * 365 * 24
  times   <- seq(0, Td * 1.5, by = dt)
  a <- log(Nmax_mean) / Td

  b_Nmax <- sqrt(log(1 + Nmax_sd^2))
  Nmax   <- rlnorm(1, log(Nmax_mean) - 0.5 * b_Nmax^2, b_Nmax)
  
  # Initialise
  Nc <- Nmax * 0.01
  m  <- mc * Nc
  m_vec <- numeric(length(times))
  B_vec <- numeric(length(times))
  jumps <- numeric(length(times))
  
  for (i in seq_along(times)) {
    
    # Logistic growth
    det_term <- a * Nc * (1 - Nc / Nmax) * dt
    
    # Stochastic jump
    theta        <- 2
    dev_progress  <- Nc / Nmax
    jump_sd       <- Nc^0.5 * (phi + (1 - phi) * (1 - dev_progress)^theta)  * sqrt(dt/2)
    jump_term <- VGAM::rlaplace(1, 0, jump_sd) 
    
    dNc_dt <- det_term + jump_term
    
    # Energetics
    BM <- Nc * Bc
    if(dNc_dt >= 0){
      BG <- Ec * dNc_dt / eta          
      H_ineff <- (1 - eta) * Ec * dNc_dt  
    } else {
      BG <- 0                          
      H_ineff <- 0.3 * abs(dNc_dt) * Ec   
    }
    
    H <- H_ineff + gamma * abs(jump_term)
    
    B  <- max(BM + BG + H, 1e-20)  
    
    # Mass
    dm_dt <- (B - H) * mc / Ec - m * Bc / Ec
    m     <- max(m + dm_dt * dt, 1e-20)
    Nc    <- max(Nc + dNc_dt, 1e-20)
    
    m_vec[i] <- m
    B_vec[i] <- B
    jumps[i] <- jump_term
  }
  
  # Ontogenetic exponent beta
  growth_idx <- which(m_vec < 0.9 * max(m_vec))
  growth_idx <- unique(c(growth_idx, growth_idx + 1))
  beta <- coef(lm(log(B_vec[growth_idx]) ~ log(m_vec[growth_idx])))[2]
  
  #hist(unlist(jumps[growth_idx]), breaks = 50)
  
  # Fluctuation analysis
  tau <- 1
  
  # Adult phase
  BMR <- B_vec[floor(length(B_vec) * 0.9):length(B_vec)]
  r   <- log(BMR[(1 + tau):length(BMR)] / BMR[1:(length(BMR) - tau)])
  b_lap <- mean(abs(r - median(r)))
  b_sd <- sd(r)
  kurt  <- kurtosis(r)
  
  # Development phase
  BMR_dev <- B_vec[1:max(growth_idx)]
  r_dev   <- log(BMR_dev[(1 + tau):length(BMR_dev)] / BMR_dev[1:(length(BMR_dev) - tau)])
  b_lap_dev <- mean(abs(r_dev - median(r_dev)))
  kurt_dev  <- kurtosis(r_dev)
  
  # Model fit
  aic_gaussian <- -2 * logLik_gaussian(r) + 2 * 2  
  aic_laplace  <- -2 * logLik_laplace(r)  + 2 * 2  
  best_fit     <- ifelse(aic_gaussian < aic_laplace, "Gaussian", "Laplace")
  
  if (plot) {
    par(mfrow = c(3, 1))
    plot(m_vec, type = "l", main = "Body mass", ylab = "Mass", xlab = "Time step")
    plot(r, type = "l", main = "Metabolic fluctuations", ylab = "dB/dt", xlab = "Time step")
    hist(r, breaks = 100, main = "Histogram of metabolic fluctuations", xlab = "dB/dt")
  }
  
  data.frame(
    id          = id,
    mass        = m,
    B_mu        = mean(BMR) / 3600,
    mean_r      = mean(r),
    r_lap_sd    = b_lap,
    r_sd        = b_sd,
    r_lap_sd_dev= b_lap_dev,
    beta        = beta,
    kurt        = kurt,
    kurt_dev    = kurt_dev,
    Nmax        = formatC(Nmax, format = "e", digits = 2),
    Nmax_mean   = Nmax_mean,
    Nmax_sd     = Nmax_sd,
    gamma       = gamma,
    phi         = phi,
    jump        = mean(abs(jumps)),
    aic_gaussian= aic_gaussian,
    aic_laplace = aic_laplace,
    best_fit    = best_fit
  )
}

# Expected adult metabolic rate 
B_adult <- function(Nmax, Bc, Ec, gamma, phi, eta, dt = 0.1, kneg = 0.3) {
  s <- phi * sqrt(Nmax) * sqrt(dt / 2)
  add <- s * ( Ec * ( 1/(2*eta) + (1 - eta + kneg)/2 ) + gamma )
  (Nmax * Bc + add) / 3600
}

# Expected adult fluctuation MAD of r 
sigma_r_mc <- function(Nmax, Bc, Ec, gamma, phi, eta, dt = 0.1, kneg = 0.3, nsim = 2e5) {
  B0 <- Nmax * Bc
  s  <- phi * sqrt(Nmax) * sqrt(dt / 2)
  
  x1 <- VGAM::rlaplace(nsim, 0, s)
  x2 <- VGAM::rlaplace(nsim, 0, s)
  
  P1 <- pmax(x1, 0)            
  N1 <- pmax(-x1, 0)           
  P2 <- pmax(x2, 0)
  N2 <- pmax(-x2, 0)
  
  c_pos <- Ec/eta + (1 - eta) * Ec + gamma   
  c_neg <- kneg * Ec + gamma                 
  
  B1 <- B0 + c_pos * P1 + c_neg * N1
  B2 <- B0 + c_pos * P2 + c_neg * N2
  
  r <- log(B2 / B1)
  b_r <- mean(abs(r - stats::median(r)))
  
  aic_gaussian <- -2 * logLik_gaussian(r) + 2 * 2  
  aic_laplace  <- -2 * logLik_laplace(r)  + 2 * 2  
  best_fit     <- ifelse(aic_gaussian < aic_laplace, "Gaussian", "Laplace")
  
  return(list(sigma_r=b_r,  best_fit=best_fit ))
}

# Expected adult mass
m_adult <- function(Nmax, Bc, Ec, phi, m_c, eta, dt = 0.1) {
  s <- phi * sqrt(Nmax) * sqrt(dt / 2)
  m_c * ( Nmax + (Ec / (Bc * eta)) * (s / 2) )
}

# Expected species scaling exponent beta
beta_species <- function(Nmax, Bc, Ec, gamma, phi, eta, dt = 0.1, kneg = 0.3) {
  A  <- Bc
  C0 <- phi * sqrt(dt / 2) * ( Ec * ( 1/(2*eta) + (1 - eta + kneg)/2 ) + gamma )
  E0 <- phi * sqrt(dt / 2) * ( (Ec / Bc) * (1/(2*eta)) )
  
  n <- Nmax
  num <- (A + C0 / (2 * sqrt(n))) * (n + E0 * sqrt(n))
  den <- (A * n + C0 * sqrt(n)) * (1 + E0 / (2 * sqrt(n)))
  num / den
}



