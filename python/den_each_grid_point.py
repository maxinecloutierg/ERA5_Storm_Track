import pandas as pd
import xarray as xr
import pdb
import time

"""
    Maxine Cloutier-Gervais

Created : 

    June 9th, 2023

Info : 
    
    This script creates calculates the cyclone center density of every grid point. 
    The cyclone center density represents the number of centers that were active within 
    a given grid point. 

"""

def get_den(df_in) : 
    
    """
    Determines storm density of all the storms that passed within a grid point.

    Parameters : 
        df_in : Name of the dataframe variable 

    Returns : 
        df_dn : dataframe that contains all grid points and their center density

    """

    # Step 1 : Group by latitude, longitude and season to get the count of the unique 
    #          occurence of every storm that passed by the grid point. 

    #df_dn = df_in.groupby(['latitude', 'longitude', 'season']).agg({'storm' : 'nunique'}).reset_index()
    df_dn = df_in.groupby(['latitude', 'longitude', 'season'])['storm'].count().reset_index()
    df_dn = df_dn.rename(columns={'storm' : 'storm_count'})

    return df_dn

#df = pd.read_csv('/pampa/cloutier/etc24_consec_v4.csv')
start = time.time()
df = pd.read_csv('/pampa/cloutier/storm_tracks/NAEC/NAEC_1979-2020_month_to_season.csv')

# test for 2000 - 2020
#df_20 = df.loc[(df.lifetime // 1000000 >= 2000) | (df.lifetime // 1000000 <= 2020)]
dens = get_den(df)
dens.to_csv('/pampa/cloutier/density/NAEC/NAEC_month_to__season_center_den_each_gp.csv', index = False)

print('Time of execution : ' , time.time() - start)
