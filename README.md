# Summer Intership 2023

This repository contains all the codes that I created during my 2023 summer intership at Université du Québec à Montéal under the direction of professor Alejandro Di Luca. The main objective for this intership was to recreate figures 1 and 2 in [Chen & al (2022)](https://doi.org/10.1029/2022GL098776) about the seasonality of extratropical cyclone (ETC) trajectories and characteristics. In the paper, the region of interest is located over Northeastern North America. However, during my internship, the region of interest is the domain of the 6th version of the Canadian Regional Climate Model (CRCM6/GEM5).

The data source comes from the first version of the [North American Extratropical Cyclone (NAEC) catalogue](https://doi.org/10.5683/SP3/LH8OBV), which provides information about extratropical cyclone trajectories over North America from January 1979 to December 2020. The source data for this catalogue comes from hourly calculations of the fifth version of the European Center for Medium-Range Weather Forecasts [(ECMWF ERA5)](https://www.ecmwf.int/en/forecasts/dataset/ecmwf-reanalysis-v5) reanalysis with a resolution of 0.25x0.25°. Information about ETCs trajectories also includes impact variables such as near surface wind speed, wind gust and precipitation, averaged over various radii around the storm center. 

In addition to recreating the figures that represent the seasonnality of the storms trajectories as well as the relation between the near-surface wind speed and the 850hPa wind gust, I also created some code to represent other storm characteristics, such as the cyclone center density and the track density. Definitions for these two variables can be found in [Neu & al. (2013)](https://doi.org/10.1175/BAMS-D-11-00154.1). The relation between near-surface wind speed and precipitation was also explored. 


## How to use
Import python libraries before using `.ipynb`and `.py` files : 
```
$ module load python3
$ source activate base_plus
```

### Jupyter Folder
`jupyter` folder contains all `.ipynb` codes. I mainly use them to 

1. Run quick codes (less than 5 minutes)  
2. Create maps and figures  
3. Debug some bigger code  
   
### Python Folder
`python` folder contains all `.py` executable files. These files are mainly used to process heavier data, like the whole NAEC Catalogue. The codes are mainly for :  

1. Calculating track density and averaged intensity
2. Filter the catalogue to keep the storms that were active in the region of interest (CRCM6 domain)






