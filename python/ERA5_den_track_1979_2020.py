import pandas as pd
import numpy as np
import xarray as xr
import time

# start time (print time of execution at the end)
start = time.time()

# Read catalogue with season column
df = pd.read_csv('/pampa/cloutier/storm_tracks/NAEC/NAEC_1979-2020_max_season.csv', index_col=0)

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

# Step 3: Find the maximum lifetime for each storm and use it to define the dimensions of the arrays
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

R = 2.5e5 # searching radius in meters

storm = np.unique(df['storm'])
lifetime = df['lifetime'].to_numpy()

storm_max_lifetime = df.groupby('storm')['lifetime'].max()
pmax = storm_max_lifetime.to_numpy()

# virtual ERA5 grid (0.25° x 0.25°) of dim : 
# Lon : (375.75 - 0) / 0.25
# Lat : (0 - 90) / 0.25
vix = 1439
vjx = 360

# Values of lat and lon in virtual grid
vlat = np.arange(0, 90.25, 0.25)
vlon = np.arange(0, 360  , 0.25)

# Track density at each virtual grid point for each season
trackDen1 = np.zeros((vix+1, vjx+1))
trackDen2 = np.zeros((vix+1, vjx+1))
trackDen3 = np.zeros((vix+1, vjx+1))
trackDen4 = np.zeros((vix+1, vjx+1))

# Determining track density...

for num in storm :
    print('Processing Storm # ... ', num, '/ 24604')
    # The virtual grid point is set to true if a storm track 
    # passes within a 250km radius around it.
    #print('Initializing Trackpass...')
    Trackpass1 = np.full((vix+1, vjx+1), False, dtype=bool)
    Trackpass2 = np.full((vix+1, vjx+1), False, dtype=bool)
    Trackpass3 = np.full((vix+1, vjx+1), False, dtype=bool)
    Trackpass4 = np.full((vix+1, vjx+1), False, dtype=bool)
    
    # Keep in mind that num starts at 1 whereas 
    # indices in arrays start with 0 
    if pmax[num - 1] > 1 : 
        
        # Iterate through all grid points of the storm
        for point in range(1, pmax[num - 1] + 1): 
            Tlat = lat[num - 1, point-1] # storm center lat
            Tlon = lon[num - 1, point-1] # storm center lon
            
            # Find the right II and JJ that points to the corresponding vlat and vlon
            # when lat = 0, vi = 0
            JJ = int(Tlat / 0.25) 
            II = int(Tlon / 0.25)   
            #print(Tlon, 'is II = ', II, 'and', Tlat, 'is JJ = ', JJ)
            
            # Affect True to the grid cell at JJ, II
            if season[num-1, point-1] == 'JJA' : 
                Trackpass1[II, JJ] = True
            if season[num-1, point-1] == 'SON' : 
                Trackpass2[II, JJ] = True
            if season[num-1, point-1] == 'DJF' : 
                Trackpass3[II, JJ] = True
                #print ('Initializing Trackpass3[', II, ']', '[', JJ, '] to ', Trackpass3[II, JJ])
            if season[num-1, point-1] == 'MAM' : 
                Trackpass4[II, JJ] = True
            
            # Linear interpolation of the track path 
            if point > 1 : 
                Dlat = lat[num-1, point-1] - lat[num-1, point-2]
                Dlon = lon[num-1, point-1] - lon[num-1, point-2]
                
                if abs(Dlat) > abs(Dlon) : 
                    seg_point = int(abs(Dlat)/ 0.2)
                else : 
                    seg_point = int(abs(Dlon)/0.2)
                
                if seg_point > 1 : 
                    # Check where the trajectory line falls one the grid 
                    for tp in range (1, seg_point + 1) : 
                        Tlon = lon[num-1, point-2] + (tp * Dlon/seg_point)
                        Tlat = lat[num-1, point-2] + (tp / seg_point) * Dlat
                        JJ = int(Tlat / 0.25) 
                        II = int(Tlon / 0.25)
                        
                        if season[num-1, point-1] == 'JJA' : 
                            Trackpass1[II, JJ] = True
                        if season[num-1, point-1] == 'SON' : 
                            Trackpass2[II, JJ] = True
                        if season[num-1, point-1] == 'DJF' : 
                            Trackpass3[II, JJ] = True
                        if season[num-1, point-1] == 'MAM' : 
                            Trackpass4[II, JJ] = True
    
 
    for vi in range (0, vix + 1) : 
        for vj in range (0, vjx + 1) : 
            if Trackpass1[vi, vj] : 
                trackDen1[vi,vj] = trackDen1[vi,vj] + 1
            if Trackpass2[vi, vj] : 
                trackDen2[vi,vj] = trackDen2[vi,vj] + 1
            if Trackpass3[vi, vj] : 
                trackDen3[vi,vj] = trackDen3[vi,vj] + 1
            if Trackpass4[vi, vj] : 
                trackDen4[vi,vj] = trackDen4[vi,vj] + 1

trackDen1_T = np.transpose(trackDen1)
trackDen2_T = np.transpose(trackDen2)
trackDen3_T = np.transpose(trackDen3)
trackDen4_T = np.transpose(trackDen4)

# Create xarray DataArrays for track densities for each season with the correct coordinates
trackDen_JJA_da = xr.DataArray(trackDen1_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_SON_da = xr.DataArray(trackDen2_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_DJF_da = xr.DataArray(trackDen3_T, coords=[('latitude', vlat), ('longitude', vlon)])
trackDen_MAM_da = xr.DataArray(trackDen4_T, coords=[('latitude', vlat), ('longitude', vlon)])


# ... Your previous code to save the data ...

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
#dataset.to_netcdf('/pampa/cloutier/track_density/NAEC/Track_density_1979_2020_all_seasons.nc')

# print time of execution
print('Time of execution : ', time.time() - start)
