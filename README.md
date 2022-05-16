# Snowfall_Scavenging_MCS
Repository of the cloud settling model, data and scripts used in analyzing water transport to the polar regions of Mars via scavenging by CO2 snowfall using MCS data

The original MCS data is obtained from NASA’s Planetary Data System (https://atmos.nmsu.edu/data_and_services/atmospheres_data/MARS/mcs.html![image](https://user-images.githubusercontent.com/40036308/168688055-24306332-0a99-44ce-9785-ad6ef7df0a3d.png)
The data obtained is the observation date, ls, local time, altitude, latitude, longitude, temperature, pressure , dust opacity, and h2oice opacity for both poles above 50 degrees latitude, for Mars years 29 to 35.

The MCS_modules.py file contains the cloud settling code, as well as all the functions used to bin and manipulate the data.

The North_Pole_Water_Scavenging_by_CO2_Snowfall_MY29-35.ipynb and South_Pole_Water_Scavenging_by_CO2_Snowfall_MY29-35.ipynb files are jupyter notebooks with step by step process of binning and processing the data, as well as running it through the cloud settling model. 

The results of the data processing and cloud settling results are stored in the two pickle files.

Snowfall_Scavenging_Plot_Maker.ipynb is a jupter notebook with additional plots and results generated from the results in the pickle files.
