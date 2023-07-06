# swu_dhe - Sectoral water use responses to droughts and heatwaves
The scripts here allow to reproduce the calculations and analysis of sectoral water use responses during droughts, heatwaves and compound events.

## Requirements
    - python 3.9.12
    - xarray 2022.3.0
    - pandas 1.4.2
    - geopandas 0.10.2
    - numpy 1.21.6
    - multiprocess 0.70.12.2
    - scipy 1.8.0
    - python-dateutil 2.8.2

## Content
    - extremes.py
    - heatwaves_W5E5.py
    - droughts_GSIM-GRDC.py
    - droughts_CCI-GLEAM.py
    - compound_events.py
    - analysis_global-scale.py
    - analysis_global-country-scale.py
    - analysis_country-scale.py
    - analysis_local-scale.py

## Configuration
### Input files
Input files can be obtain from Zenodo (**10.5281/zenodo.7657304**), in the folder '**/inputs**':
-    **gsim-grdc_discharges.csv**: GSIM and GRDC combined dataset
-    **datasets_global-scale_aquastat.csv**: AQUASTAT water withdrawal data at country level
-    **datasets_country-scale_usgs.csv**: USGS water withdrawal data at county level
-    **datasets_country-scale_useia.xlsx**: EIA data of power plants (cooling water withdrawals and net electricity production)
-    **datasets_local-scale.xlsx**: compilation of city-scale water use records
    
Other input files required:
-    Global reconstructed water use data for 1971-2010 v2.0 (Huang et al. (2018), https://doi.org/10.5281/zenodo.1209296)
-    WFDE5 over land merged with ERA5 over the ocean v2.0: Daily Maximum Near-Surface Air Temperature (Lange et al. (2021),
    https://data.isimip.org/datasets/38d4a8f4-12e8-44ff-afe3-0c7ce0e0dad6/)
-    European Space Agency - Climate Change Initiative: Soil Moisture v6.1 (Dorigo et al. (2017),
    https://www.esa-soilmoisture-cci.org/index.php?q=dataregistration)
-    Global Land Evaporation Amsterdam Model v3.5a: Root-zone soil moisture(Martens et al. (2017), https://www.gleam.eu/)- The Global Streamflow Indices and Metadata Archive - Part 1: Station catalogand Catchment boundary (Do et al. (2018),
    https://doi.pangaea.de/10.1594/PANGAEA.887477)
    
Download these files and place them in folder '**/inputs**'. Datasets downloaded sliced by time (**W5E5 tasmax**, **ESA CCI SM** and **GLEAM SMroot**) must be merged into single files that cover the entire time of analysis (e.g., 1990-2019).

### Drought-Heatwave events (DHE) identification
Scripts: **heatwaves_W5E5.py**, **droughts_GSIM-GRDC.py**, **droughts_CCI-GLEAM.py**, **compound_events.py**.
Modify the input and output directories accordingly to user's requirements.
Parameters used to identify extreme events could be modified in the section 'Input' in each script. Details on the parameters used can be found in the manuscript.
Other important information about the scripts:

-    first, run scripts that identify independent extreme events (**heatwaves_W5E5.py**, **droughts_GSIM-GRDC.py**, **droughts_CCI-GLEAM.py**).
-    '**compound_events.py**' must be run after the other independent extreme event files were generated.
-    '**extremes.py**' contain tools required for the analysis; it is not needed to be run.

Resulting datasets of identified extreme events contain tags in the file name indicating the parameters considered during calculation.
By default, generated files are stored in folder **'/dhe'**.

### Analysis of sectoral water use responses
Scripts: **analysis_global-scale.py**, **analysis_global-country-scale.py**, **analysis_country-scale.py**, **analysis_local-scale.py**.
Modify the input and output directories accordingly to user's requirements.
Only '**analysis_global-scale.py**' allows to modify the period of analysis.
Script '**analysis_global-country-scale.py**' allows analysis for AQUASTAT and USGS datasets.
By default, generated files are stored in folder '**/analysis**'.
File names include the water dimension, sector and period of interest evaluated, except for local-scale results, which varies by city and data availability.

## Contact
For further questions, please contact Gabriel Cardenas (g.a.cardenasbelleza@uu.nl).

swu_dhe v1.0.0: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8122382.svg)](https://doi.org/10.5281/zenodo.8122382)

## References
Huang, Z., Hejazi, M., Li, X., Tang, Q., Vernon, C., Leng, G., Liu, Y., Döll, P., Eisner, S., Gerten, D., Hanasaki, N., & Wada, Y. (2018). Reconstruction of global gridded monthly sectoral water withdrawals for 1971-2010 and analysis of their spatiotemporal patterns. Hydrology and Earth System Sciences, 22(4), 2117–2133. https://doi.org/10.5194/hess-22-2117-2018

Lange, S., Menz, C., Gleixner, S., Cucchi, M., Weedon, G. P., Amici, A., Bellouin, N., Müller Schmied, H., Hersbach, H., Buontempo, C., & Cagnazzo, C. (2021). WFDE5 over land merged with ERA5 over the ocean (W5E5 v2.0). ISIMIP Repository. https://doi.org/10.48364/ISIMIP.342217

Dorigo, W., Wagner, W., Albergel, C., Albrecht, F., Balsamo, G., Brocca, L., Chung, D., Ertl, M., Forkel, M., Gruber, A., Haas, E., Hamer, P. D., Hirschi, M., Ikonen, J., de Jeu, R., Kidd, R., Lahoz, W., Liu, Y. Y., Miralles, D., … Lecomte, P. (2017). Remote Sensing of Environment ESA CCI Soil Moisture for improved Earth system understanding : State-of-the art and future directions. Remote Sensing of Environment, 203, 185–215. https://doi.org/10.1016/j.rse.2017.07.001

Martens, B., Miralles, D. G., Lievens, H., Van Der Schalie, R., De Jeu, R. A. M., Fernández-Prieto, D., Beck, H. E., Dorigo, W. A., & Verhoest, N. E. C. (2017). GLEAM v3: Satellite-based land evaporation and root-zone soil moisture. Geoscientific Model Development, 10(5), 1903–1925. https://doi.org/10.5194/gmd-10-1903-2017

Do, H. X., Gudmundsson, L., Leonard, M., Westra, S. (2018): The Global Streamflow Indices and Metadata Archive - Part 1: Station catalog and Catchment boundary. PANGAEA, https://doi.org/10.1594/PANGAEA.887477 
