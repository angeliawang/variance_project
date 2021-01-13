# variance_project
 granule cells what do they do

Centralized repository for recurrent inhibition project, mostly so I can track changes.

Code parknotes:
spiking_network_rp.m generates and saves spiking data according to the model and experimental parameters specified in the script. Biological values are taken from visual system experiments (see Clopath et al. 2010). 

streamlined_sd.m is my (cleaned up) version of an old spike_discrimination.m code. Generates stimulus discrimination results based off the specified time windows. Need to include the same inhibitory parameters that are specified in spiking_network_rp.m

streamlined_sd_plotting.m makes plots based off the results from streamlined_sd.m. This particular script is definitely a mess (or shall we say... thoroughly customized to whatever plots I needed at the time).

power_spectrum.m is a function that generates the power spectral density of the specified spiking data

plot_things.m and plot_things2.m are scripts I wrote while investigating the analytical variance and PSD calculations in the 1- and 2-cell case, respectively. 