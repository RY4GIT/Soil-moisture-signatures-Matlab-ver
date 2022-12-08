# Soil-moisture-signatures-Matlab-ver
This repository contains the Matlab soil moisture signature codes used for the analysis in the paper Araki et al., (2022). This  repository is in development: e.g., still missing complete instructions and the comments in the code. If you have any questions or feedback, or if you spotted an error or bug, please create an issue on Github. Alternatively, email Ryoko Araki (raraki8159 (at) sdsu.edu)

# Citations 
If you wish to analyze & publish using the codes, please cite the following paper. 
### A paper that developed original soil moisture signature sets
Branger F, McMillan HK. (2020). Deriving hydrological signatures from soil moisture data. Hydrological processes 34 (6): 1410–1427 DOI: https://doi.org/10.1002/hyp.13645
### A paper that developed Matlab codes based on Branger & McMillan (2020)
Araki, R., Branger, F., Wiekenkamp, I., & McMillan, H. (2022). A signature-based approach to quantify soil moisture dynamics under contrasting land-uses. Hydrological Processes 36 (4) e14553. DOI: https://doi.org/10.1002/hyp.14553
### A paper that developed streamflow Matlab signature codes, which are applied to the event-based soil moisture signature codes
Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. DOI: https://doi.org/10.1016/j.envsoft.2021.104983
### About the response order signature
For using the response type signature, please contact to Inge Wiekenkamp (inge.wiekenkamp@gfz-potsdam.de). 

# Instruction

## Codes:
All the codes are in [5_code_sig](https://github.com/RY4GIT/Soil-moisture-signatures-Matlab-ver/tree/main/5_code_sig)
- **sig_xxx.m: The primary signature codes**. 
- **util_xxx.m: The utility functions that are called by sig_xxx functions**
- main.m: The main script. It basically reads the input data, runs the signature codes, and saves the results with iterating networks, stations, and sensor depths. You might need to create your original main code suite for your own purpose. This is written messy (from me to me before the 3 years of graduate training), which I want to rewrite it in more modular/object-oriented style in future :() 
- io_siteinfo.m: Contains information specific to the soil moisture network. The same disclaimer to the main function applies. 

## Input data:
Example data is in [4_data_after_qc](https://github.com/RY4GIT/Soil-moisture-signatures-Matlab-ver/tree/main/4_data_after_qc/HB_F_cleaned_csv)
- **Soil moisture timeseries data** (The main script accepts files in CSV format, two columns of time and volumetric soil water content, in the hourly interval, per each sensor. You can tweak the code)
- **Rainfall timeseries data** (The same specification as the soil moisture timeseries. The unit is in mm/hr)
- [Optional] Station quality flag (The flags helps skip some sensors that you would want to skip. However, you don't need this if you are only interested in the signature functions, not the main script. Basically, the flag 100 means good data, and the flag 101 means bad data to skip. Organized by the sensor stations(row) x the sensor depths(column) per a soil moisture sensor network.)