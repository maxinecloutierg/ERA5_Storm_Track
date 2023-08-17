import pandas as pd
import numpy as np
import xarray as xr
import time

"""

    Author : Maxine Cloutier-Gervais
    Version : August 17th, 2023
    Credits : Ting-Chen Chen (/home/chen/codes/cir_disttr/Cyclone_intensity_sesons.f90)

    This code calculates the average VORSmax around a 250km radius of a given grid point in ERA5 domain
    in a 0.5°x0.5° resolution

    To launch the code : 
    python3 /home/cloutier/summer_2023/python/ERA5_storm_inten_1979_2020_250.py

    You need to import python libraries first with 
    module load python3 
    source activate base_plus

    The code takes about 40 minutes to run, so it's better to run it in the background
    (with tmux, for example). 

"""

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

# create 2d arrays so we can refer to season, latitude, longitude and VORSmax  with 
# storm and lifetime 

# Step 1: Create a dictionary to store the storm information, including lifetimes, latitude, longitude, and season
storms_info = {}
print('creating dictionnary...')
for _, row in df.iterrows():
    storm = row['storm']
    lifetime = row['lifetime']
    season_value = row['season']
    lat_value = row['latitude']
    lon_value = row['longitude']
    vors_value = row['VORSmax']
    
    if storm not in storms_info:
        storms_info[storm] = {'lifetimes': [], 'latitudes': [], 'longitudes': [], 'seasons': [], 'vors' : []}
    
    storms_info[storm]['lifetimes'].append(lifetime)
    storms_info[storm]['latitudes'].append(lat_value)
    storms_info[storm]['longitudes'].append(lon_value)
    storms_info[storm]['seasons'].append(season_value)
    storms_info[storm]['vors'].append(vors_value)

# Step 2: Get unique storm names
unique_storms = df['storm'].unique()

# Step 3: Find the maximum lifetime for each storm and use it to define the dimensions of the arrays
max_lifetime = max(max(storms_info[storm]['lifetimes']) for storm in unique_storms)
num_storms = len(unique_storms)

# Step 4: Create the arrays using the maximum lifetime
season = np.empty((num_storms, max_lifetime), dtype=object)
lat = np.empty((num_storms, max_lifetime), dtype=float)
lon = np.empty((num_storms, max_lifetime), dtype=float)
vor = np.empty((num_storms, max_lifetime), dtype=float)
season.fill('')  # Initialize the array with empty strings
lat.fill(np.nan)
lon.fill(np.nan)
vor.fill(np.nan)

# Step 5: Fill the arrays with the corresponding values
for i, storm in enumerate(unique_storms):
    lifetimes = storms_info[storm]['lifetimes']
    latitudes = storms_info[storm]['latitudes']
    longitudes = storms_info[storm]['longitudes']
    seasons = storms_info[storm]['seasons']
    vors = storms_info[storm]['vors']
    
    for j, lifetime in enumerate(lifetimes):
        lat[i, j] = latitudes[j]
        lon[i, j] = longitudes[j]
        season[i, j] = seasons[j]
        vor[i,j] = vors[j]


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### PART TWO : CALCULATE AVERAGE VORTICITY  #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

R = 2.5e5 # searching radius in meters

# Get the unique storm IDs and all lifetimes. 
storm = np.unique(df['storm'])
lifetime = df['lifetime'].to_numpy()

# Get the max lifetime of each storm (pmax)
storm_max_lifetime = df.groupby('storm')['lifetime'].max()
pmax = storm_max_lifetime.to_numpy()

# virtual ERA5 grid (0.25° x 0.25°) of dim : 
# Lon : ((359.5 - 0) / 0.5) + 1 (indices start at 0)
# Lat : ((90 - 0) / 0.5) + 1
vix = 720
vjx = 181
# Values of lat and lon in virtual grid
vlat = np.arange(0, 90.5, 0.5)
vlon = np.arange(0, 360 , 0.5)

# Storm intensity at each virtual grid point for each season
Sinten1 = np.zeros((vix, vjx))
Sinten2 = np.zeros((vix, vjx))
Sinten3 = np.zeros((vix, vjx))
Sinten4 = np.zeros((vix, vjx))

# Storm count at each virtual grid point for each season 
Snum1 = np.zeros((vix, vjx))
Snum2 = np.zeros((vix, vjx))
Snum3 = np.zeros((vix, vjx))
Snum4 = np.zeros((vix, vjx))

# Determining track density...

for num in storm :
    print('Processing storm #', num)
    
    # Keep in mind that num starts at 1 whereas indices in python start at 0
        
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
        sgd = int(R/1.1e5/0.5) + 2 # 2 is a buffer 
        maxi = max(0, II-sgd)
        mini = min(II+sgd, vix)
        maxj = max(0, JJ-sgd)
        minj = min(JJ+sgd, vjx)

        # Searching i and j that are within a distance R from II, JJ
        for i in np.arange(maxi, mini, 1) : 
            for j in np.arange(maxj, minj, 1) : 
                dist = calc_dist(vlat[j], vlon[i], Tlat, Tlon)

                if dist < R : 
                    # Affect True to the grid cell at II, JJ and add 1 to the storm count
                    if season[num-1, point-1] == 'JJA' : 
                        Sinten1[i,j] = Sinten1[i,j]+vor[num-1,point-1]*1.e5
                        Snum1[i,j] = Snum1[i,j]+1

                    if season[num-1, point-1] == 'SON' : 
                        Sinten2[i,j] = Sinten2[i,j]+vor[num-1,point-1]*1.e5
                        Snum2[i,j] = Snum2[i,j]+1
                    
                    if season[num-1, point-1] == 'DJF' : 
                        Sinten3[i,j] = Sinten3[i,j]+vor[num-1,point-1]*1.e5
                        Snum3[i,j] = Snum3[i,j]+1
                    
                    if season[num-1, point-1] == 'MAM' : 
                        Sinten4[i,j] = Sinten4[i,j]+vor[num-1,point-1]*1.e5
                        Snum4[i,j] = Snum4[i,j]+1
            
# When all storms are looped through, we iterate through all the virtual grid points
# and calculate the average vorticity  
for vi in np.arange(0, vix-1, 1) : 
    for vj in range (0, vjx-1, 1) : 
        if Snum1[vi, vj] > 0 : 
            Sinten1[vi,vj] = Sinten1[vi,vj]/Snum1[vi,vj]
        if Snum2[vi, vj] > 0 : 
            Sinten2[vi,vj] = Sinten2[vi,vj]/Snum2[vi,vj]
        if Snum3[vi, vj] > 0 :
            Sinten3[vi,vj] = Sinten3[vi,vj]/Snum3[vi,vj]
        if Snum4[vi, vj] > 0 : 
            Sinten4[vi,vj] = Sinten4[vi,vj]/Snum4[vi,vj]


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### ####  PART THREE : CREATE NETCDF FILE   #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

Sinten1_T = np.transpose(Sinten1)
Sinten2_T = np.transpose(Sinten2)
Sinten3_T = np.transpose(Sinten3)
Sinten4_T = np.transpose(Sinten4)

#Create xarray DataArrays for track densities for each season with the correct coordinates
trackInt_JJA_da = xr.DataArray(Sinten1_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackInt_SON_da = xr.DataArray(Sinten2_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackInt_DJF_da = xr.DataArray(Sinten3_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackInt_MAM_da = xr.DataArray(Sinten4_T, coords=[('latitude', vlat), ('longitude', vlon)])


# Create a new xarray Dataset to store the variables
dataset = xr.Dataset({
    'longitude': ('longitude', vlon),
    'latitude': ('latitude', vlat),
    'trackInt_JJA': trackInt_JJA_da,
    'trackInt_SON': trackInt_SON_da,
    'trackInt_DJF': trackInt_DJF_da,
    'trackInt_MAM': trackInt_MAM_da,
})

# Save the dataset to a NetCDF file
dataset.to_netcdf('/pampa/cloutier/density/NAEC/Track_intensity_1979_2020_all_seasons_250km.nc')

print('Execution Time : ', time.time() - start)
