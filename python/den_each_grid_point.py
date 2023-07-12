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
df = pd.read_csv('/pampa/cloutier/storm_tracks/NAEC/NAEC_1979-2020_season.csv')

# test avec les 6636 tempêtes
#df = pd.read_csv('/pampa/cloutier/etc24_consec.csv')
dens = get_den(df)
dens.to_csv('/pampa/cloutier/track_density/NAEC/NAEC_track_den_each_gp.csv')
