import pandas as pd
import numpy as np
import xarray as xr
import time

"""

Author : Maxine Cloutier-Gervais
Version : August 17th, 2023
Credits : Ting-Chen Chen /home/chen/codes/cir_disttr/Track_density_all_seasons.f90

This code takes the NAEC catalogue and calculates the track densities of all storms around 
a 250 km radius on ERA5 grid in a 0.5x0.5° resolution. 

Track density is the number of tracks passing a grid cell, with repeated entries of the same track being
counted as one (Neu & al, 2013)

To launch the code : 
python3 /home/cloutier/summer_2023/python/ERA5_den_track_1979_2020_250km.py

You need to import python libraries first with 
module load python3 
source activate base_plus

The code takes about two hours to run, so it's better to run it in the background
(with tmux, for example). 

"""

# start time
start = time.time()

# Read catalogue with season column
df = pd.read_csv('/pampa/cloutier/storm_tracks/NAEC/NAEC_1979-2020_max_season.csv', index_col=0)


def calc_dist(lat1, lon1, lat2, lon2) : 

"""

Calculate Haversine distance between two grid points in km 

Paramters : 
- lat1, lon1, lat2, lon2 : Coordinates (in °) of point 1 and point 2

Returns : 
dist12 : Distance (in km) between point 1 and point 2

"""
    # Earth's radius in meters
    r = 6371.22E3
    
    rlat1 = np.radians(lat1)
    rlon1 = np.radians(lon1)
    rlat2 = np.radians(lat2)
    rlon2 = np.radians(lon2)
    
    # Haversine formula
    dlat = abs(rlat2 - rlat1)
    dlon = abs(rlon2 - rlon1)
    a = np.sin(dlat / 2)**2 + np.cos(rlat1) * np.cos(rlat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    dist12 = r * c 
    
    return dist12


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### ####  PART ONE : CREATE 2D ARRAYS  #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####


# create 2d arrays so we can refer to season, latitude and longitude with 
# storm and lifetime 

# Step 1: Create a dictionary to store the storm information, including lifetimes, latitude, longitude, and season
storms_info = {}

for _, row in df.iterrows():
    storm = row['storm']
    lifetime = row['lifetime']
    season_value = row['season']
    lat_value = row['latitude']
    lon_value = row['longitude']
    
    if storm not in storms_info:
        storms_info[storm] = {'lifetimes': [], 'latitudes': [], 'longitudes': [], 'seasons': []}
    
    storms_info[storm]['lifetimes'].append(lifetime)
    storms_info[storm]['latitudes'].append(lat_value)
    storms_info[storm]['longitudes'].append(lon_value)
    storms_info[storm]['seasons'].append(season_value)

# Step 2: Get unique storm names
unique_storms = df['storm'].unique()

# Step 3: Find the maximum lifetime and the number of storm IDs 
# and use it to define the dimensions of the arrays
max_lifetimes = [max(storms_info[storm]['lifetimes']) for storm in unique_storms]
num_storms = len(unique_storms)
max_lifetime = max(max_lifetimes)

# Step 4: Create the arrays using the maximum lifetime
season = np.empty((num_storms, max_lifetime), dtype=object)
lat = np.empty((num_storms, max_lifetime), dtype=float)
lon = np.empty((num_storms, max_lifetime), dtype=float)
season.fill('')  # Initialize the array with empty strings
lat.fill(np.nan)
lon.fill(np.nan)

# Step 5: Fill the arrays with the corresponding values
for i, storm in enumerate(unique_storms):
    lifetimes = storms_info[storm]['lifetimes']
    latitudes = storms_info[storm]['latitudes']
    longitudes = storms_info[storm]['longitudes']
    seasons = storms_info[storm]['seasons']
    
    for j, lifetime in enumerate(lifetimes):
        lat[i, j] = latitudes[j]
        lon[i, j] = longitudes[j]
        season[i, j] = seasons[j]


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### PART TWO : CALCULATE THE STORM TRACK DENSITY #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

R = 2.5e5 # searching radius in meters

# Get the unique storm IDs and all lifetimes
storm = np.unique(df['storm'])
lifetime = df['lifetime'].to_numpy()

# Get the lifetime of each storm (pmax)
storm_max_lifetime = df.groupby('storm')['lifetime'].max()
pmax = storm_max_lifetime.to_numpy()

# virtual ERA5 grid (0.25° x 0.25°) of dim : 
# Lon : ((359.5 - 0) / 0.5) + 1(indices start at 0)
# Lat : ((90 - 0) / 0.5) + 1
vix = 720
vjx = 181

# Values of lat and lon in virtual grid
vlat = np.arange(0, 90.5, 0.5)
vlon = np.arange(0, 360 , 0.5)

# Track density at each virtual grid point for each season
# 1 = JJA and 4 = MAM
trackDen1 = np.zeros((vix, vjx))
trackDen2 = np.zeros((vix, vjx))
trackDen3 = np.zeros((vix, vjx))
trackDen4 = np.zeros((vix, vjx))

# Determining track density...

for num in storm :
    print('Processing storm #', num)
    # The virtual grid point is set to true if a storm track passes within a 250km radius around it.
    Trackpass1 = np.full((vix, vjx), False, dtype=bool)
    Trackpass2 = np.full((vix, vjx), False, dtype=bool)
    Trackpass3 = np.full((vix, vjx), False, dtype=bool)
    Trackpass4 = np.full((vix, vjx), False, dtype=bool)
    
    # Keep in mind that num starts at 1 whereas indices in arrays start with 0 
    if pmax[num - 1] > 1 : 
        
        # Iterate through all grid points of the storm
        for point in range(1, pmax[num - 1] + 1): 
            Tlat = lat[num - 1, point-1] # storm center lat
            Tlon = lon[num - 1, point-1] # storm center lon
            
            # Find the right II and JJ that points to the corresponding vlat and vlon
            # when lat = 0, vi = 0
            JJ = int(Tlat / 0.5) 
            II = int(Tlon / 0.5) 
            if Tlat%0.5 > 0.25 : 
                JJ += 1
            if Tlat%0.5 > 0.25 : 
                II += 1
            #print(Tlon, 'is II = ', II, 'and', Tlat, 'is JJ = ', JJ)
            
            # Determine searching distance in virtual grid points
            sgd = int(R/1.1e5/0.5) + 2 # 2 is a buffer since 1.1e5 is an approx. 
            maxi = max(0, II-sgd)
            mini = min(II+sgd, vix)
            maxj = max(0, JJ-sgd)
            minj = min(JJ+sgd, vjx)

            # Searching i and j that are within a distance R from II, JJ
            for i in np.arange(maxi, mini, 1) : 
                for j in np.arange(maxj, minj, 1) : 
                    dist = calc_dist(vlat[j], vlon[i], Tlat, Tlon)

                    if dist < R : 
                        # Affect True to the grid cell at II, JJ
                        if season[num-1, point-1] == 'JJA' : 
                            Trackpass1[i, j] = True
                        if season[num-1, point-1] == 'SON' : 
                            Trackpass2[i, j] = True
                        if season[num-1, point-1] == 'DJF' : 
                            Trackpass3[i, j] = True
                        if season[num-1, point-1] == 'MAM' : 
                            Trackpass4[i, j] = True
            
            if point > 1 : 
                # Get the distance with the previous grid point
                Dlat = lat[num-1, point-1] - lat[num-1, point-2]
                Dlon = lon[num-1, point-1] - lon[num-1, point-2]
                
                # Linear interpolation of the track path   
                if abs(Dlat) > abs(Dlon) : 
                    seg_point = int(abs(Dlat)/ 0.2)
                else : 
                    seg_point = int(abs(Dlon)/0.2)

                if seg_point > 1 : 
                    # Check where the trajectory line falls one the grid 
                    for tp in range (1, seg_point + 1) : 
                        Tlon = lon[num-1, point-2] + (tp * Dlon/seg_point)
                        Tlat = lat[num-1, point-2] + (tp / seg_point) * Dlat
                        JJ = int(Tlat / 0.5) 
                        II = int(Tlon / 0.5)
                        if Tlat%0.5 > 0.25 : 
                            JJ += 1
                        if Tlat%0.5 > 0.25 : 
                            II += 1

                        # Determine searching distance in virtual grid points
                        sgd = int(R/1.1e5/0.5) + 2 # 2 is a buffer 
                        maxi = max(0, II-sgd)
                        mini = min(II+sgd, vix)
                        maxj = max(0, JJ-sgd)
                        minj = min(JJ+sgd, vjx)
                        
                        # Searching i and j that are within a distance R from II, JJ
                        for i in np.arange(maxi, mini) : 
                            for j in np.arange(maxj, minj) : 
                                if calc_dist(vlat[j], vlon[i], Tlat, Tlon) < R : 
                                    # Affect True to the grid cell at JJ, II
                                    if season[num-1, point-1] == 'JJA' : 
                                        Trackpass1[i, j] = True
                                    if season[num-1, point-1] == 'SON' : 
                                        Trackpass2[i, j] = True
                                    if season[num-1, point-1] == 'DJF' : 
                                        Trackpass3[i, j] = True
                                    if season[num-1, point-1] == 'MAM' : 
                                        Trackpass4[i, j] = True

# Once all storms are looped through, we iterate through all the grid points in Trackpass 
# and determine how many time each grid point was set to true.    

    for vi in np.arange(0, vix-1, 1) : 
        for vj in range (0, vjx-1, 1) : 
            if Trackpass1[vi, vj] : 
                trackDen1[vi,vj] = trackDen1[vi,vj] + 1
            if Trackpass2[vi, vj] : 
                trackDen2[vi,vj] = trackDen2[vi,vj] + 1
            if Trackpass3[vi, vj] : 
                trackDen3[vi,vj] = trackDen3[vi,vj] + 1
            if Trackpass4[vi, vj] : 
                trackDen4[vi,vj] = trackDen4[vi,vj] + 1

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### ####  PART THREE : CREATE NETCDF FILE   #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

trackDen1_T = np.transpose(trackDen1)
trackDen2_T = np.transpose(trackDen2)
trackDen3_T = np.transpose(trackDen3)
trackDen4_T = np.transpose(trackDen4)

#Create xarray DataArrays for track densities for each season with the correct coordinates
trackDen_JJA_da = xr.DataArray(trackDen1_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_SON_da = xr.DataArray(trackDen2_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_DJF_da = xr.DataArray(trackDen3_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_MAM_da = xr.DataArray(trackDen4_T, coords=[('latitude', vlat), ('longitude', vlon)])

# Create a new xarray Dataset to store the variables
dataset = xr.Dataset({
    'longitude': ('longitude', vlon),
    'latitude': ('latitude', vlat),
    'trackDen_JJA': trackDen_JJA_da,
    'trackDen_SON': trackDen_SON_da,
    'trackDen_DJF': trackDen_DJF_da,
    'trackDen_MAM': trackDen_MAM_da,
})

# Save the dataset to a NetCDF file
dataset.to_netcdf('/pampa/cloutier/density/NAEC/Track_density_1979_2020_all_seasons_250km.nc')

# print time of execution
print('Time of execution : ', time.time() - start)