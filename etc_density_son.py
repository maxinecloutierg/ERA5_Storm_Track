import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import xarray as xr
import time

start_time = time.time()

# open catalogue dataset
df = pd.read_csv('/home/data/ReAnalysis/ERA5/Storm_analysis/NAECv1/NAEC_1979_2020_v1.csv')

# open catalogue that contains storm that were active for more than 24h
# more in CRCM6 domain
df24 = pd.read_csv('/pampa/cloutier/etc_24_nna.csv')

# open netcdf mask file
file = '/pampa/picart/Masks/mask_GEM5_ERA5grid'
data = xr.open_dataset(file)

# export netcdf to dataframe
mask = data.to_dataframe()

# drop index lat lon, but keep columns
mask = mask.reset_index()

# keep in mask values where HU = true
maskT = mask.loc[mask.HU != False]

# open catalogue that contains storm that were active for more than 24h in CRCM domain 
# for each season
djf = pd.read_csv('/pampa/cloutier/etc_24_nna_djf.csv')
mam = pd.read_csv('/pampa/cloutier/etc_24_nna_mam.csv')
jja = pd.read_csv('/pampa/cloutier/etc_24_nna_jja.csv')
son = pd.read_csv('/pampa/cloutier/etc_24_nna_son.csv')

# Define function to calculate distance between two points using Haversine formula
# Haversine package ?
def calculate_distance(lat1, lon1, lat2, lon2):
    # Earth radius in kilometers
    earth_radius = 6371

    # Convert latitude and longitude to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Calculate differences
    dlat = lat2_rad - lat1_rad
    dlon = lon2_rad - lon1_rad

    # Haversine formula
    a = math.sin(dlat / 2) ** 2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(dlon / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Calculate distance in kilometers
    distance = earth_radius * c

    return distance


# test de la densité avec seulement SON
r = 2.25
son_density = {} # will contain the grid point and density

# iterate through all grid point
for _, row in son.iterrows() : 
    track = {} #  will contain the year and storm id of the grid points found
    lat1 = row['latitude'] 
    lon1 = row['longitude'] 

    # restricting the research in a 2.5 x 2.5deg square
    neighbors = son[
        (son['latitude'] >= lat1 - r) & (son['latitude'] <= lat1 + r) &
        (son['longitude'] >= lon1 - r) & (son['longitude'] <= lon1 + r)
        ]

    # Find the grid points that are within a 250km radius 
    for _, row in neighbors.iterrows():
        lat2 = row['latitude'] 
        lon2 = row['longitude'] 
        year = str(row['datetime'])[:4]
        storm_id = row['storm']

        # calculate distance between grid points
        distance = calculate_distance(lat1, lon1, lat2, lon2)

        # add year and storm_id of grid points that are within the radius
        if distance <= 250 : 
            track['year'] = track.get('year', []) + [year]
            track['storm_id'] = track.get('storm_id', []) + [storm_id]

    # turn track dictionnary into dataframe and count the unique occurence of 
    # every storm for each year
    track_df = pd.DataFrame(data = track)
    track_df = track_df.groupby('year')['storm_id'].nunique()

    # number of storm per season in average
    average = track_df.mean()

    # add coord and average density
    son_density['lat'] = son_density.get('lat', []) + [lat1]
    son_density['lon'] = son_density.get('lon', []) + [lon1]
    son_density['avg_per_season'] = son_density.get('avg_per_season', []) + [average]

# transfort m_density into dataframe and save as csv
df_density = pd.DataFrame(data = son_density)

df_density.to_csv('/pampa/cloutier/test_son_density.csv')

end_time = time.time()

exec_time = end_time - start_time

print('temps pour rouler SON ', exec_time)
