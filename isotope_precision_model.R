library(IsoplotR)
library(MonteCarlo)

# Moving mean function (credit: Matti Pastell, https://stackoverflow.com/questions/743812/calculating-moving-average)
ma <- function(x,n){filter(x,rep(1/n,n), sides=1)}


## Monte Carlo parameters ##

# Commented variables run model with more or fewer parameters. 

# Instrument precision
d_sd_grid <- c(3, 1.5, 0.5)  # Delta 13C standard deviations of of G2201-i and hypothetical instruments 
# d_sd_grid <- 3  # Delta 13C standard deviation of G2201-i
c_sd_grid <- 0.00039  # Concentration standard deviation of G2201-i
# c_sd_grid <- c(0.0020, 0.00039, 0.0001)  # Concentration standard deviations of UGGA, G2201-i, and hypothetical instrument

# Background air
d_bg_grid <- -47.56  # Based on 2015 mean Mace Head CH4 delta 13C concentration
c_bg_grid <- 1.91  # Based on 2015 mean Mace Head CH4 concentration

# Source
# d_sc_grid <- c(-62.2, -44.0, -22.2)  # Signatures of microbial, fossil fuel, and biomass burning from Schwietzke et al. 2016
d_sc_grid <- -44.0  # Estimate of fossil fuel from Schwietzke et al. 2016
# c_sc_grid <- c(0.5, 1, 2, 3.5, 6)  # Maximum peak height above background in ppm
c_sc_grid <- c(0.5, 1, 2.5, 5, 7.5, 10, 15, 20)  # Maximum peak height above background in ppm
sd_sc_grid <- 1 # Standard deviation/peak shape
# sd_sc_grid <- c(0.5 ,1, 2) # Standard deviation/peak shape

# Number of measurements per peak
n_grid <- c(100, 250, 500, 1000)

# Exchange rate, i.e. number of measurement cycles to replace gas in the instrument cavity.
exch_rt_grid <- c(20, 40, 60)
# exch_rt_grid <- 90

# Parameter grid for MC simulation
param_grid <- list(d_sd = d_sd_grid, c_sd = c_sd_grid, 
                   d_bg = d_bg_grid, c_bg = c_bg_grid, 
                   d_sc = d_sc_grid, c_sc = c_sc_grid, sd_sc = sd_sc_grid,
                   exch_rt = exch_rt_grid,
                   n = n_grid)


## Model ##

cavity_sim <- function(d_sd, c_sd, d_bg, c_bg, d_sc, c_sc, sd_sc, exch_rt, n) {
  
  # Construct peak
  peak_length <- seq(0 - 4 * sd_sc, 0 + 4 * sd_sc, length = n) # Include values up to 4 standard deviations from mean
  peak_conc <-dnorm(peak_length, mean = 0, sd = sd_sc)
  peak_conc <- peak_conc*c_sc/max(peak_conc)

  # Construct chamber
  cavity <- data.frame(timepoint = 1:(2 * length(peak_length) + exch_rt)) 
  cavity$peak_conc <- 0
  cavity$peak_conc[(exch_rt + 1):(length(peak_length) + exch_rt)] <- peak_conc
  cavity$CH4_true <- ma(cavity$peak_conc, exch_rt) + c_bg  # Concentration above background at timepoint (measurement cycle)
  cavity[is.na(cavity)] <- c_bg

  # Source signature calculated with two-pool mixing model
  cavity$d13C_true <- (d_sc * (cavity$CH4_true - c_bg) + d_bg * c_bg) / cavity$CH4_true 
  
  # Add random noise to concentration and isotope values
  cavity$c_sd <- c_sd
  cavity$c_noise <- rnorm(nrow(cavity), 0, c_sd)
  cavity$CH4_sim <- cavity$CH4_true + cavity$c_noise
  cavity$d_sd <- d_sd
  cavity$d_noise <- rnorm(nrow(cavity), 0, d_sd)
  cavity$d13C_sim <- cavity$d13C_true + cavity$d_noise
  
  # Set up data for York regression
  cavity <- cavity[cavity$CH4_true > c_bg, ]
  cavity$CH4_x_d13C <- cavity$CH4_sim * cavity$d13C_sim
  cavity$CH4_x_d13C_sd <- abs(
    cavity$CH4_x_d13C * sqrt(
    (cavity$c_sd / cavity$CH4_sim)^2 + (cavity$d_sd / cavity$d13C_sim)^2
  ))  # Gaussian error propagation for CH4 * d13C
  
  cavity$cor <- cor(cavity$CH4_x_d13C, cavity$CH4_sim)
  cavity <- cavity[c("CH4_sim", "c_sd", "CH4_x_d13C", "CH4_x_d13C_sd", "cor")]
  
  # Run York regression
  result <- york(cavity)
  result <- as.numeric(result$b[2])
  return(list(se = result))
  
}


## Monte Carlo simulation ##

# Run simulation
mc_results <- MonteCarlo(func = cavity_sim, 
                         param_list = param_grid, 
                         time_n_test = T, 
                         nrep = 1000,
                         save_res = T,
                         ncpus = 2)

# Compile results table in LaTeX format. Note that LaTeX packages "multirow" and "graphicx" are required.
rows <- c("sd_sc", "n", "exch_rt", "d_sd")
cols <- c("c_sc")

MakeTable(output=mc_results, rows=rows, cols=cols, digits=2, include_meta=FALSE)
