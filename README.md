# Soil-moisture-signatures-Matlab-ver
Repository for the Matlab soil moisture signature scripts used for the analysis in the paper [Araki et al., (2022)](https://doi.org/10.1002/hyp.14553) at Hydrological Processes. This repository is in development: e.g., still missing complete instructions and the comments in the code.

If you have any questions or feedback, or if you spotted an error or bug, please create an issue on Github. Alternatively, email Ryoko Araki (raraki8159 (at) sdsu.edu)

# Citations 
If you wish to analyze & publish using the scripts, please cite the following paper. 
- **A paper that developed Matlab scripts based on Branger & McMillan (2020)**
    - Araki, R., Branger, F., Wiekenkamp, I., & McMillan, H. (2022). A signature-based approach to quantify soil moisture dynamics under contrasting land-uses. Hydrological Processes 36 (4) e14553. DOI: https://doi.org/10.1002/hyp.14553
- **A paper that developed original soil moisture signature sets**
    - Branger F, McMillan HK. (2020). Deriving hydrological signatures from soil moisture data. Hydrological processes 34 (6): 1410–1427 DOI: https://doi.org/10.1002/hyp.13645
- **A paper that developed streamflow Matlab signature scripts, which are applied to the event-based soil moisture signature scripts**
    - Gnann, S.J., Coxon, G., Woods, R.A., Howden, N.J.K., McMillan H.K., 2021. TOSSH: A Toolbox for Streamflow Signatures in Hydrology. Environmental Modelling & Software. DOI: https://doi.org/10.1016/j.envsoft.2021.104983
- **About the response order signature**
    - For using the response type signature, please contact to Inge Wiekenkamp (inge.wiekenkamp@gfz-potsdam.de). 

# How to use the scripts
## Scripts
All the scripts are in [5_code_sig](https://github.com/RY4GIT/Soil-moisture-signatures-Matlab-ver/tree/main/5_code_sig).
- **sig_xxx.m**
    -  **The primary signature scripts**. 
- **util_xxx.m**
    - **The utility functions that are called by sig_xxx functions**
- main.m
    - The main script. It basically reads the input data, runs the signature scripts, and saves the results with iterating networks, stations, and sensor depths.
    - You might need to create your original main code suite for your own purpose. This is written messy (from me to me before the 3 years of graduate training), which I want to rewrite it in more modular/object-oriented style in future :() 
- io_siteinfo.m
    - Contains information specific to the soil moisture network. 
    - The same disclaimer to the main function applies. 

## Input data
An example dataset is in [4_data_after_qc](https://github.com/RY4GIT/Soil-moisture-signatures-Matlab-ver/tree/main/4_data_after_qc/HB_F_cleaned_csv)
- **Soil moisture timeseries data**
    - The main script accepts files in CSV format, two columns of time and volumetric soil water content, in the hourly interval, per each sensor. You can tweak the code if you wish to use data in different specifications. 
- **Rainfall timeseries data**
    - The same specification as the soil moisture timeseries. The unit is in mm/hr. 
- [Optional] Station quality flag
    - You don't need this if you are only interested in the signature functions, not the main script. 
    - The flags helps skip some sensors that you would want to skip. Basically, the flag 100 means good data, and the flag 101 means bad data to skip. Organized by the sensor stations(row) x the sensor depths(column) per a soil moisture sensor network.
    

## List of signatures
| Signature                 | Unit                      | Corresponding function    | Output abbreviation   | Description |
| -------------             | -------------         | -------------             | -------------             | ------------- |
| Normalized amplitude      | [-]                       | sig_event.m               | amplitude             | Event amplitude was calculated as the difference between the soil moisture values at their maximum and at the start of the event, normalized using estimated field capacity and wilting point at the station  |
| Rising time               | [hour]                    | sig_event.m               | risingtime            | For each event, event rising time was calculated as the time-lag from the start of an event to the soil moisture peak  |
| No-response rate          | [-]                       | sig_event.m               | noresrate             | No-response rate was calculated as the number of events with no response divided by the number of all events  |
| Rising limb density       | [1/timestep]              | sig_RLD.m                 | RLD                   | Rising limb density was calculated as the inverse of the average rising time of all events |
| Response type             | [sequential/non-sequential/no response]               | sig_restype.m             | restype               | Response type was classified as ‘sequential’ when the response order was sequential from the shallow to the deeper sensor; as ‘non-sequential’ when the order of response times is non-sequential for at least one sensor; ‘No-response’ was assigned when none of the sensors responded |
| Duration of seasonal transition (dry to wet season)       | [days]            | sig_seasontrans_ForestFlow.m  | duration_dry2wet_p    | Seasonal transition signatures were calculated by fitting a piecewise linear model to the soil moisture timeseries for each wet-to-dry and dry-to-wet transition period. Transition duration was defined as the length of time between the start and the end day |
| Start date of seasonal transition (dry to wet season)     | [day of the year] | sig_seasontrans_ForestFlow.m  | sdate_dry2wet_p       | The start and end days of the transition were defined as the inflection points of the piecewise linear model, expressed in the day of the year.  |
| End date of seasonal transition (dry to wet season)       | [day of the year] | sig_seasontrans_ForestFlow.m  | edate_dry2wet_p       | Same as above  |
| Estimated field capacity  | [m3/m3]               | sig_fcwp.m                | fc                    | We calculated the estimated field capacity and wilting point as the peaks of the soil moisture PDF. First, peaks of the soil moisture PDF were detected. The peak with the largest volumetric soil moisture content was defined as the estimated field capacity |
| Estiamted wilting point   | [m3/m3]               | sig_fcwp.m                | wp                    |  The peak with the smallest volumetric soil moisture content was defined as the estimated wilting point |
| Distribution type         | [uni/bi/multi-modal]  | sig_pdf.m                 | disttype              | Soil moisture PDFs were classified according to the number of peaks into ‘unimodal’ (one peak), ‘bimodal’ (two peaks), or ‘multimodal’ (three or more peaks) |

![alt text](./readme/signature_schematics.png)
