This contains the scripts used for the Master's Thesis "The relative Site Frequency Spectrum (rSFS) of subsamples versus a population". Their usage can be simplified in 3 steps: run the SLiM simulations to obtain your population data, obtain various statistics from this data using mstatspop, and finally graph the data of interest in R. There are 2 folders, the old scripts were the ones used originally at the start, while the updated ones had various later parameter changes and optimizations for better results.

# SLiM Simulations

Install SLiM from their website (https://messerlab.org/slim/), then download "sweeps_ind.slim" and then run the code from "0.run_selsweeps.sh" in the console.

# mstatspop calculations

Copy and compile mstatspop from their repository (https://github.com/CRAGENOMICA/mstatspop), download "columns_to_choose.txt" inside the bin folder in the main mstatspop folder, then run the code from "1.run_mstatspop.sh" in the console.

# Graphs in R

Download and install R (https://www.r-project.org/), download the "Sweep_SFSline_diversity_sims.R" R script, then run the code from "2.run_Rscript.sh" in the console.
