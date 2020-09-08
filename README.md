# Spatio-temporal-decoupling-water-splitting-
Data and code used for technoeconomic analysis in "Spatio-temporal decoupling of water electrolysis for dual-use grid energy storage and hydrogen generation"

## Data
Electricity price data from Seel et al. [1] were used for this analysis. The cleaned data for use in our analysis are in the csv files, caiso_data.csv and nyiso_data.csv.

## Code
Julia was used in this analysis due to its efficiency for large optimization problems. Two files are provided in the repository. The first, "Before_After_Data.jl", takes the electricity price data and clusters the data into a certain number of "representative time periods" (RTPs). The second, "Stochastic_optimizer.jl", performs the optimization of the operation and sizing of the device. 
