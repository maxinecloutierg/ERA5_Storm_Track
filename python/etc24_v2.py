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

def open_cat_mask(cat_in, bnd_in, mask_in) : 

    """
    Open NAEC catalogue, boundary file  and transform netCDF mask into dataframe 

    Parameters  : 
        cat_in  : path of the catalogue csv file
        mask_in : path of the mask netCDF file
        bnd_in  : path of the csv file that defines the boundary grid points
                  of CRCM6 domain
                    
    Returns : 
        cat : Dataframe containing NAEC catalogue data
        mk  : Dataframe containing mask data
        bnd : Dataframe containing boundary grid points 

    """
    
    # Step 1 : Open catalogue and boundary csv file

    cat = pd.read_csv(cat_in)
    bnd = pd.read_csv(bnd_in, index_col = 0)

    # Step 2 : Open mask netCDF file and convert into dataframe
    mk = xr.open_dataset(mask_in)
    mk = mk.to_dataframe()

    # Step 3 :  Drop index lat lon, but keep columns
    mk = mk.reset_index()

    # Step 4 : Rename lat & lon columns for latitude & longitude
    mk = mk.rename(columns={'lat' : 'latitude', 'lon' : 'longitude'})

    return cat, bnd, mk


def get_distance(latS, lonS, bnd) : 

    """
    Determine if a given grid point is at a minimal distance of 5deg from
    all CRCM6 boundary grid point domain

    Paramters : 
        latS  : Latitude of the catalogue grid point
        lonS  : Longitude of the catalogue grid point
        bnd   : Dataframe containing boundary grid points

    Returns       : 
        dist_cond : True if all grid points are within a minimal distance of 5deg
                    from all boundary layer grid points and False if not.
    """

    dist_cond = True

    for _, row1 in bnd.iterrows():
        latD = row1['lat']
        lonD = row1['lon']
        dist = ((latS-latD)**2 + (lonS - lonD)**2)**0.5
        
        if dist <= 5 : 
            dist_cond = False
            break
            
    return dist_cond




def add_season(df) : 

    """
    Add a column called 'season' in df24 that gives the season in which the ETC occured. 
    If the ETC occured in two or more season, the chosen season will be the one in which 
    the ETC has the most grid point

    DJF : December, January & November
    MAM : March, April & May
    JJA : June, July & April
    SON : September, October and December
    
    Parameters : 
        df (dataframe) : Dataframe to which we want to add the season column

    returns : 
        df_new : Dataframe with the season column
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
    
    df_new = df.merge(storm_seasons['season'], on='storm', how='left')

    # Step 6 : Delete month column

    df_new = df_new.drop(['month'], axis = 1)
    
    # Step 7 : move season column next to datetime (TODO)
    
    #df_new.insert(3, 'season', df_new.pop('season'))

    return df_new



""" MAIN PROGRAM """


# Step 1 : Open catalogue, boundary catalogue and mask
cat_in = ('/home/data/ReAnalysis/ERA5/Storm_analysis/NAECv1/NAEC_1979_2020_v1.csv')
bnd_in = ('/pampa/cloutier/outline_crcm6_domain.csv')
mask_in = ('/pampa/picart/Masks/mask_GEM5_ERA5grid')

cat, bnd, mk = open_cat_mask(cat_in, bnd_in, mask_in)

# Step 2 : Merge cat and mask to add HU column in cat

merge = cat.merge(mk, how='left', on=['latitude', 'longitude'])
merge = merge.fillna(value = False)

# Step 3 : Initialize empty dataframe that will contain the final result

df24 = pd.DataFrame(columns = cat.columns)

# Step 4 : Filter catalogue data

# Iterate through each storm 
for storm_id in merge['storm'].unique():
    storm_data = merge[merge['storm'] == storm_id].copy() # copy of merge for the given storm
    count = 0 # lifetime count

    # Iterate through each grid point of the storm
    for _, row in storm_data.iterrows() : 
        cond = False 
        hu = row['HU']
        latS = row['latitude']
        lonS= row['longitude']
        # check if storm center is within subdomain and at a 5° minimal 
        # distance from boundaries
        if hu : 
#            print('check distance ...')
            cond = get_distance(latS, lonS, bnd) 
        stInDom.append(cond)
    
    # add a new column that determines if each storm center agrees or not with the above condition
#    print('adding StInDom in storm_data ...')
    storm_data['StInDom'] = stInDom
#    print(storm_data['StInDom'])
    n = 23
    # exclude the last 23 lines in the search
    storm_data_rows = storm_data.head(len(storm_data) - n)
    
    count = 0
    
#    print('second for loop : ', row['storm'])
    for idx, row in storm_data_rows.iterrows() : 
        if row ['StInDom'] : 
#            print('Initializing count to one')
            count = 1
            
            # keep iterating for the next 23 grid points
            for _, row in islice(storm_data.iterrows(), idx, idx+23) : 
                if row['StInDom'] : 
                    count += 1
                
                else :  
#                    print(row['StInDom'], 'storm not in domain')
                    break
            
            if count >= 24 : 
                break
                
#    print(count)                           
    if count >= 24 :
        df24 = df24.append(storm_data)
        print('Year in process : ', df24['datetime'].iloc[-1])              

#        # For each grid point that is in the subdomain (HU = True), determine if the latter is at
#        # a minimal distance of 5deg from all grid points of the subdomain boundary 
#        if hu : 
#            dist_cond = get_distance(latS, lonS, bnd)

            # If the dist_cond is true, it also means that HU = True 
#            if dist_cond : 
#                count += 1 # Iterate the lifetime in subdomain count

        # If the grid point is not within the subdomain, but the count is already above 24, 
        # we can get out of the "for" loop
#        if hu == False and count >= 24 : 
#           break

        # If the grid point is not within the subdomain but the count is still under 24, 
        # we need to restart the count and keep iterating through the storms grid point
#        if hu == False and count < 24 : 
#            count = 0

    # Add to the new dataframe all etc that were alive for 24 consecutive hours or more 
    # in CRCM6 subdomain
#    if count >= 24 : 
#        df24 = df24.append(storm_data)
#        print('Year in process : ', df24['datetime'].iloc[-1])

# Step 5 : Add the season column in dataframe
df24_season = add_season(df24)

# Step 6 : Save dataframe as csv 
df24_season.to_csv('/pampa/cloutier/etc24_consec_v2.csv', index = False)
