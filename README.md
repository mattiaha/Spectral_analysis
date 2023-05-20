Github repository containing the code and figures used in the master's thesis "Searching for Dark Matter with CTA". A Python environment containing the Gammapy v1.1 is needed in order to run the code. 

The simulated data from the DC-1 can be found at https://redmine.cta-observatory.org/ , but it requires permission from CTA to access it. The files are too large to share here.
The galactic center-models are stored in the gc_models.yaml file. The code is written such that this file should be in the same folder as the script.
The irf's used should also be accessible from the cta website, and the environmental variables CALDB and CTADATA should be defined as they suggest here.

read_fits.py is a script that selects observations for the galactic center and gives a list of the observational ID's fitting to our criteria. In this case our criteria is any pointing
position for an observation located between 0.5 and 1.0 degrees from the galactic center, but this can easily be changed. It saves these observational ID's as a list in the file obs_ids_dc.npy .

This list can then be used in the file on_off.py which is the code which performs the on-off analysis. It is based on the tutorial https://docs.gammapy.org/dev/tutorials/analysis-1d/spectral_analysis.html . Running this script performs a fit of a DM model defined in the code to the simulated observations which ID's are given by obs_ids_dc.npy. The coordinates of the on-region can be changed through the variable target_position, and the radius of the region can be changed in the variable on_region_radius. To change the energies we include in our fit we can change the values in the energy_axis variable. on line 100. The first value determines the minimum value and the second the maximum value. 
It then plots the predicted excess counts from our model after the fit against the counts in the dc-1 simulated observations.

find_sigma.py performs a similar task but the aim here is to find the scale values that corresponds to a delta test statistic of 3.8 for each mass value. This is then saved along with errors and mass values and used in plot_sigma.py to show the exclusion plot for mass. 

make_onoff_obs.py is the code used to create simulations on our own. These events are saved as FITS-files and can be loaded into the on_off.py script in the same way as the DC-1 simulations as long as the observational ID's are known.
