Metabolic Scaling Fluctuations
==============================

Analysis scripts and data for the manuscript "Metabolic fluctuations explain allometric scaling diversity".

---------------------------------------
Overview
---------------------------------------
This repository contains code and data used to study scaling fluctuations in ecological systems.

The scripts implement:
• Data cleaning and preparation  
• Empirical analysis of trait variation across groups and species  
• Parameter sweeps over simulation models  
• Figure generation for manuscript or report development  


---------------------------------------
Data Description
---------------------------------------
data_ind.csv  
• Individual-level empirical dataset  
• Includes species, group labels, and trait measurements  
• Scripts compute number of individuals per species (n_ind)

sims.RData  
• Simulation results used in generating figures and comparisons.


---------------------------------------
Scripts Overview
---------------------------------------
functions.R  
• Utility and helper functions for simulations and analysis.

Figures.R  
• Generates all manuscript or report figures.

Fig5.R  
• Produces Figure 5 only.  
• Performs species filtering, group recoding, and scaling analysis.

sweep_parallel.R  
• Runs parallelized simulation sweeps across parameter sets.


---------------------------------------
Notes
---------------------------------------
• Some scripts may assume local file paths (e.g., '~/data_ind.csv'); update paths as needed.  
• Simulation scripts may require multi-core processing capabilities.  


---------------------------------------
License
---------------------------------------
MIT License  


