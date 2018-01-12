# Mass_model_grid_codes

gps_gen.py is a code to generate the gps mac file to be used to run the Geant4 simulation.
It creates the mac files with the name Tddd.dd_Pddd.dd_Edddd.dd_gps.mac in the directories 
Tddd.dd_Pddd.dd where T, P and E refer to theta, phi and energy for the particular run.
Copy over these directories to where you want the outputs to be written. The final fits files,
Tddd.dd_Pddd.dd_Edddd.dd_SingeEventFile.fits will be written in the respective Tddd.dd_Pddd.dd directory.

MAKE SURE THESE DIRECTORIES ARE COPIED TO THE RIGHT PLACE. OTHERWISE, NO PRODUCT WILL BE WRITTEN EVEN AFTER THE RUN!!

plot_basemap.py is still being developed so that we can visualise HTM pixels. Will be updated soon ;D
