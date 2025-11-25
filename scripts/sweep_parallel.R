library(foreach)
library(doSNOW)
source("~/functions.R")

n_cores <- parallel::detectCores() - 1  
cl <- makeCluster(n_cores, type = "SOCK")
registerDoSNOW(cl)

gamma_values <- seq(0.00001, 4, length.out = 400)  
phi_values   <- seq(0.00001, 1, length.out = 100)  
Nmax_mean    <- c(1e10, 5e12, 7e13, 5e14, 5e15, 5e16)

id <- expand.grid(Nmax_mean, gamma_values, phi_values)
colnames(id) <- c("Nmax", "gamma", "phi")

Bc  <- 3e-11 * 3600   
Ec  <- 3.7e-7         
m_c <- 1e-12          
eta <- 0.7
dt  <- 0.1

funs_to_export <- c("B_adult", "m_adult", "beta_species", "sigma_r_mc")

pb <- txtProgressBar(min = 0, max = nrow(id), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

res_sweep <- foreach(i = 1:nrow(id),
                     .combine = rbind,
                     .export  = funs_to_export,
                     .options.snow = opts) %dopar% {
                       
                       vars <- as.numeric(id[i, ])   # Nmax, gamma, phi
                       
                       B_i    <- B_adult(Nmax  = vars[1], Bc = Bc, Ec = Ec,
                                         gamma = vars[2], phi = vars[3],
                                         eta   = eta, dt  = dt)
                       
                       M_i    <- m_adult(Nmax  = vars[1], Bc = Bc, Ec = Ec,
                                         phi   = vars[3], m_c = m_c,
                                         eta   = eta, dt  = dt)
                       
                       beta_i <- beta_species(Nmax  = vars[1], Bc = Bc, Ec = Ec,
                                              gamma = vars[2], phi = vars[3],
                                              eta   = eta, dt  = dt)
                       
                       sim_sigma <- sigma_r_mc(Nmax  = vars[1], Bc = Bc, Ec = Ec,
                                                gamma = vars[2], phi = vars[3], eta   = eta, dt  = dt)
                       
                       sd_i   <- sim_sigma$sigma_r
                       
                       bf   <- sim_sigma$best_fit
                       
                       res <- data.frame(B = B_i, mass = M_i, beta = beta_i, sd_lap = sd_i, best_fit = bf)
                     }

close(pb)
stopCluster(cl)
registerDoSEQ()   

res_sweep <- as.data.frame(res_sweep)

d <- data.frame(id, res_sweep)

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

save(d, file = "/sims.RData")


