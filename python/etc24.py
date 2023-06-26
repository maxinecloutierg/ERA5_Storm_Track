import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import xarray as xr

""" 

    Maxine Cloutier-Gervais

Created : 

    June 7th, 2023

Info : 
    
    This code creates and filters, a file that contains ETC that were active for 24 consecutives 
    hours or more in the crcm6 domain.

"""

def open_cat_mask(cat_in, mask_in) : 

    """
    Open NAEC catalogue and mask netCDF file and convert them into dataframes

    Parameters  : 
        cat_in  : Path of the NAEC catalogue
        mask_in : Path of the mask netCDF

    return   : 
        df   : Dataframe that contains NAEC catalogue
        mask : Dataframe that contains the mask

    """

    # Step 1 : Open catalogue

    df = pd.read_csv('/home/data/ReAnalysis/ERA5/Storm_analysis/NAECv1/NAEC_1979_2020_v1.csv')
    
    # Step 2 :  Open mask, convert into dataframe and rename columns

    mask = xr.open_dataset(mask_path)
    mask = mask.to_dataframe() # convert netCDF to dataframe
    mask = mask.reset_index() # drop lat and lon indexes
    # rename lat and lon columns
    mask = mask.rename(columns={'lat' : 'latitude', 'lon' : 'longitude'})

    return df, mask 




def create_df24_consec(df, mask, output_file): 
    
    """
    Extract ETC that were active for 24 consecutive hours or more in CRCM6 domain

    parameters: 
        df : name of the dataframe variable to extract the ETCs from
        mask : name of the mask variable
        output_file : path and name of the resulting csv file

    return : dataframe that contains the filtered ETCs (df24)

    """
    
    # Step 1 : Merge df and mask to add 'HU' mask column in df and initialise df24

    merge = df.merge(mask, how='left', on=['latitude', 'longitude'])
    merge = merge.fillna(value = False)
    df24_consec = pd.DataFrame(columns = df.columns)
    
    # Step 2 : Filter ETCs in merge and save result as csv

    #iterate through each storm 
    for storm_id in merge['storm'].unique() : 
        storm_data = merge[merge['storm'] == storm_id].copy()
        count = 0
        
        # count how many consecutive grid points are in crcm6 domain 
        # ('HU' == True) 
        for _, row in storm_data.iterrows() : 

            if row['HU'] == True : 
                count += 1

            # if we have count >= 24 and encounter False value, exit the 
            # current for loop
            if row['HU'] == False and count >= 24 : 
                break

            # If we don't have count >=24 yet but encounter False value, 
            # reset the count and continue to another grid point
            if row['HU'] == False and count < 24 : 
                count = 0
                continue

        if count >= 24 : 
            df24_consec = df24_consec.append(storm_data)
            # print the latest datetime to check where the code is at 
            print('Year in process : ', df24_consec['datetime'].iloc[-1])
        
    # move 'HU' column next to 'longitude'
    df24_consec.insert(5, 'HU', df24_consec.pop('HU'))
    
    # save df24 as csv
    df24_consec.to_csv(output_file, index = False)




def create_df24_noconsec(df, output_file) : 

    """
    Extract ETC that were active for 24 hours or more (total) in the CRCM6 domain

    parameters:
            df : name of the dataframe variable to extract the ETCs from
            mask : name of the msk variable
            output_file : path and name of the resulting csv file

    return : dataframe that contains the filtered ETCs (df24)

    """

    # Step 1 : Merge df and mask to add 'HU' mask column in df and initialise df24

    merge = df.merge(mask, how='left', on=['latitude', 'longitude'])
    merge = merge.fillna(value=False)
    df24_noconsec = pd.DataFrame(columns=df.columns)

    # Step 2 : Filter ETCs in merge and save result as csv

    # iterate through each storm
    for storm_id in merge['storm'].unique():
        storm_data = merge[merge['storm'] == storm_id].copy()
        count = 0
        
        # Get the number of grid points for the given storm 
        # that are in or out the crcm6 domain 
        group = storm_data.groupby(['storm','HU']).size()
        
        # Calculate the sum of 'HU' values that are True
        hu_sum = storm_data.loc[storm_data['HU'] == True, 'HU'].sum()  

        if hu_sum >= 24 : 
            df24_noconsec = df24_noconsec.append(storm_data)
            # check where the code is at 
            print('Year in process : ', df24_noconsec['datetime'].iloc[-1])
    
    # move 'HU' column next to 'longitude'
    df24_noconsec.insert(4, 'HU', df24_noconsec.pop('HU'))
    
    # save df24 as csv
    df24_noconsec.to_csv(output_file, index = False)




def add_season(df, output_file) : 

    """
    Add a column called 'season' in df24 that gives the season in which the ETC occured. 
    If the ETC occured in two or more season, the chosen season will be the one in which the ETC has the most grid point

    DJF : December, January & November
    MAM : March, April & May
    JJA : June, July & April
    SON : September, October and December
    
    Parameters : 
        df (dataframe) : dataframe variable name 
        output_file (string) : The path and name of csv file

    returns : overwrites dataframe with season column 
    
    """

    seasons = { 'SON': [9, 10, 11], 'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8] }

    # Step 1 : Add 'month' column in dataframe 

    df['month'] = (df.datetime // 10000) % 100

    # Step 2 : Group the storms by their ID and count the number of grid point 
    #          in each month

    storm_seasons = df.groupby(['storm', 'month']).size().unstack().fillna(0)

    # Step 3 : Determine the month with the maximum grid points for each storm

    storm_seasons['season'] = storm_seasons.idxmax(axis=1)
    
    # Step 4 : Transform month number into season
    
    storm_seasons['season'] = storm_seasons['season'].map(
    lambda month: next((season for season, months in seasons.items() if month in months), None)
    )
    
    # Step 5 : Merge the season column into original dataframe
    
    df = df.merge(storm_seasons['season'], on='storm', how='left')

    # Step 6 : Delete month and season_x columns

    df = df.drop(['month'], axis = 1)
    
    # Step 7 : move season column next to datetime (TODO)
    
    df.insert(3, 'season', df.pop('season'))

    # Step 8 : Overwrite dataframe with df
    df.to_csv(output_file, index = False)

#output_file = '/pampa/cloutier/etc24_consec_v2.csv'
#create_df24_consec(df, mask, output_file)
# test
