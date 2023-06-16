import pandas as pd
import pdb

"""
    Maxine Cloutier-Gervais

Created : 

    June 9th, 2023

Info : 
    
    This script creates calculates the storm density and the average maximum vorticity of every grid 
    point in etc24.csv. The storm density represents the number of unique storm trackes that were active 
    within a 250km radius around a grid point and the average max vorticity represents the average value 
    (over their lifetime in the crcm6 domain and over a 800 km radius) of all the storms that passed 
    within a given grid point. 

"""

def get_vort_den(df_in) : 
    
    """
    Determine the average maximum vorticity and storm density of all the storms that passed within a grid point.
    The average is calculated only over the grid points within the crcm6 domain.
    The density is obtained by calculating the unique count of all storms that passed by a given grid point

    Parameters : 
        df_in : Name of the dataframe variable 

    Returns : 
        df_dv : dataframe that contains all grid points, their density and their average max vorticity

    """

    # Step 1 : Only keep grid points that are within crcm6 domain
    
    crcm6 = df_in.loc[df_in.HU == True]

    # Step 2 : Group by storm get VORSmax average. Drop the index to have a dataframe instead of 
    #          a pandas series as a result. Rename VORSmax as avgVORSmax

    storm_avgvort = crcm6.groupby(['storm'])['VORSmax'].mean().reset_index()
    storm_avgvort = storm_avgvort.rename(columns={'VORSmax' : 'avgVORSmax'})

    # Step 3 : Merge the average max values of every storm with the original dataframe

    merge = df_in.merge(storm_avgvort, how='left')

    # Step 4 : Group by latitude, longitude and season to get the average max vorticity of each storm as well as 
    #          the count of the unique occurence of every storm that passed by the grid point. 

    df_avgvort = merge.groupby(['latitude', 'longitude', 'season']).agg({'avgVORSmax' : 'mean', 'storm' : 'nunique'}).reset_index()
    df_avgvort = df_avgvort.rename(columns={'storm' : 'storm_count'})

    return df_avgvort

