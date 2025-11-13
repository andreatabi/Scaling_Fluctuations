rm(list = ls())
library(dplyr)
library(broom)

source("~/functions.R")

gamma_values <- seq(0.00001, 2, length.out = 100)
phi_values    <- seq(0.00001, 2, length.out = 100)
Nmax_mean    <- c(1e10, 5e12, 7e13, 5e14, 5e15, 5e16)

param_grid <- expand.grid(Nmax = Nmax_mean,
                          gamma = gamma_values,
                          phi = phi_values)   
yrs <- rep(c(0.1, 1, 18, 4, 25, 30),   100)
param_grid$yrs <- yrs

Bc  <- 1.08e-07   
Ec  <- 3.7e-7         
m_c <- 1e-12          
eta <- 0.7
cv  <- 0.2           

n_individuals <- 10

all_results <- list()

for (row in seq_len(nrow(param_grid))) {
  this <- param_grid[row, ]
  
  res <- do.call(rbind, lapply(seq_len(n_individuals), function(i) {
    SOGM(
      Nmax_sd   = cv,
      Nmax_mean = this$Nmax,
      id        = paste0("id_", row, "_", i),
      years     = this$yrs,
      gamma     = this$gamma,
      phi       = this$phi,
      eta       = eta,
      plot      = FALSE
    )
  }))
  
  res$gamma <- this$gamma
  res$jf    <- this$jf
  
  all_results[[row]] <- res
}

df_results <- do.call(rbind, all_results)
head(df_results)

df_results <- df_results %>%
  mutate(
    animal = case_when(
      Nmax_mean == 1e10  ~ "mouse",
      Nmax_mean == 5e12  ~ "cat",
      Nmax_mean == 7e13  ~ "human",
      Nmax_mean == 5e14  ~ "cow",
      Nmax_mean == 5e15  ~ "elephant",
      Nmax_mean == 5e16  ~ "whale",
      TRUE ~ NA_character_
    ),
    animal = factor(animal, levels = c("mouse", "cat", "human", "cow", "elephant", "whale"))
  )

save(df_results, file = "~/sims.RData")




