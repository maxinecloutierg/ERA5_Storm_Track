import pandas as pd
import pdb

"""
    Maxine Cloutier-Gervais

Created : 

    June 9th, 2023

Info : 
    
    This script creates calculates the storm density of every grid point in etc24.csv. 
    The storm density represents the number of unique storm tracks that were active within 
    a given grid point. 

"""

def get_den(df_in) : 
    
    """
    Determines storm density of all the storms that passed within a grid point.

    Parameters : 
        df_in : Name of the dataframe variable 

    Returns : 
        df_dn : dataframe that contains all grid points and their track density

    """

    # Step 1 : Group by latitude, longitude and season to get the count of the unique 
    #          occurence of every storm that passed by the grid point. 

    df_dn = df_in.groupby(['latitude', 'longitude', 'season']).agg({'storm' : 'nunique'}).reset_index()
    df_dn = df_dn.rename(columns={'storm' : 'storm_count'})

    return df_dn

import pandas as pd
import xarray as xr

#df = pd.read_csv('/pampa/cloutier/etc24_consec_v4.csv')

# test avec les 6636 tempêtes
#df = pd.read_csv('/pampa/cloutier/etc24_consec.csv')
#dens = get_den(df)
#dens.to_csv('/pampa/cloutier/etc24_den_with_etc24_consec_v4.csv')


# TEST POUR CRÉER UN FICHIER NETCDF AU LIEU D'UN FICHIER CSV

import pandas as pd
import xarray as xr

# Read the CSV file
csv_file = '/pampa/cloutier/etc24_consec_v4.csv'
df = pd.read_csv(csv_file)

# Extract required columns
latitudes = df['latitude']
longitudes = df['longitude']
seasons = df['season']

# Calculate storm track density grouped by latitude, longitude, and season
storm_track_density = df.groupby(['latitude', 'longitude', 'season']).size().rename('storm_count')

# Determine unique latitude, longitude, and season values
unique_latitudes = latitudes.unique()
unique_longitudes = longitudes.unique()
unique_seasons = seasons.unique()

# Create a new netCDF file
output_file = '/pampa/cloutier/den_each_gp.nc'
nc = xr.Dataset()

# Define dimensions
nc['latitude'] = xr.DataArray(unique_latitudes, dims='latitude')
nc['longitude'] = xr.DataArray(unique_longitudes, dims='longitude')
nc['season'] = xr.DataArray(unique_seasons, dims='season')

# Create empty storm track density array
density_data = xr.DataArray(data=None, coords=[unique_latitudes, unique_longitudes, unique_seasons], dims=['latitude', 'longitude', 'season'])

# Assign storm track density data to the array
density_data.loc[:] = storm_track_density.values

# Define storm track density variable
nc['storm_track_density'] = density_data

# Write data to the netCDF file
nc.to_netcdf(output_file)

# Close the netCDF file
nc.close()

print(f"NetCDF file '{output_file}' created successfully.")




