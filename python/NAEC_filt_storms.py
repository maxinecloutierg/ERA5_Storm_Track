import pandas as pd
import math
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np
import xarray as xr
import time

""" 

    Author : Maxine Cloutier-Gervais
    Version : June 7th, 2023

    This code creates a file that contains ETC that has at least 24 grid points (lifetime of 24h) 
    in the CRCM6 domain at a minimal distance of 5° from the boundary. 

    To launch the code : 
    python3 /home/cloutier/summer_2023/python/ERA5_storm_inten_1979_2020_250.py

    You need to import python libraries first with 
    module load python3 
    source activate base_plus

    The code takes a while to run, so it's better to run it in the background
    (with tmux, for example). 

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

    # Step 3 :  Drop index lat lon
    mk = mk.reset_index()

    # Step 4 : Rename lat & lon columns for latitude & longitude
    mk = mk.rename(columns={'lat' : 'latitude', 'lon' : 'longitude'})

    return cat, bnd, mk

    
def eucl_dist(lat1, lon1, lat2, lon2) : 

    """
    Use eucledian formula to get the distance between two grid points. 

    Parameters : 
        lat1   : Latitude of the storm grid point
        lon1   : Longitude of the storm grid point
        lat2   : Latitude of the domain grid point
        lon2   : Longitude of the domain grid point

    Returns  : 
        dist : Distance (in degrees) between the two coordinates
    """ 

    dist=((lat2-lat1)**2+(lon2-lon1)**2)**0.5

    return dist

def get_cond(latS, lonS, bnd) : 

    """
    Determine if a given grid point is at a minimal distance of 5deg from
    CRCM6 boundary domain

    Parameters : 
        latS   : Latitude of the storm grid point
        lonS   : Longitude of the storm grid point
        bnd    : Dataframe containing boundary grid points

    Returns       : 
        dist_cond : True if all grid points are within a minimal distance of 5deg
                    from all boundary layer grid points and False if not.
    """

    dist_cond = True

    for _, row1 in bnd.iterrows():
        latD = row1['lat']
        lonD = row1['lon']
        dist = eucl_dist(latS, lonS, latD, lonD)
        
        if dist <= 5 : 
            dist_cond = False
            break
            
    return dist_cond


def add_season(df1) : 

    """
    Add a column called 'season' in df24 that gives the season in which the ETC occured. 
    If the ETC occured in two or more season, the chosen season will be the one in which 
    the ETC has the most grid point

    DJF : December, January & November
    MAM : March, April & May
    JJA : June, July & April
    SON : September, October and December
    
    Parameters : 
        df1 (dataframe) : Dataframe to which we want to add the season column

    returns : 
        df_new : Dataframe with the season column
    """

    seasons = { 'SON': [9, 10, 11], 'DJF': [12, 1, 2], 'MAM': [3, 4, 5], 'JJA': [6, 7, 8] }

    # Step 1 : Add 'month' column in dataframe 

    df1['month'] = (df1.datetime // 10000) % 100

    # Step 2 : Group the storms by their ID and count the number of grid point 
    #          in each month

    storm_seasons = df1.groupby(['storm', 'month']).size().unstack().fillna(0)

    # Step 3 : Determine the month with the maximum grid points for each storm

    storm_seasons['season'] = storm_seasons.idxmax(axis=1)
    
    # Step 4 : Transform month number into season

    """
    Steps for this line : 

        1. 'map' function is called on 'season' column to apply a function 
            on each element in the 'season' column. 
        2.  Inside 'map' function, there is a lambda function that takes 'month'
            as an input.
                a.  The function iterates over the 'seasons' dictionnary with 
                    season for season, months in seasons.items()
                b.  For each (season, months) pair in the dictionnary, the 
                    function checks if the given month is present in the list 
                    of months. 
                c.  If a match is found, it returns the season associated with the 
                    month list in the dictionnary. The (next) function is used to 
                    retreive the first match encountered. 
                d.  None is returned of no match is found. 
        3.  Because the lambda function is used on each value in the 'season' column with 
            'map', the resulting values are applied back to the 'season' column to change the 
            month number with the associated season. 
    """
    
    storm_seasons['season'] = storm_seasons['season'].map(
    lambda month: next((season for season, months in seasons.items() if month in months), None)
    )
    
    # Step 5 : Merge the season column into original dataframe
    
    df_new = df1.merge(storm_seasons['season'], on='storm', how='left')

    # Step 6 : Delete month column

    df_new = df_new.drop(['month'], axis = 1)
    
    # Step 7 : move season column next to datetime (TODO)
    
    #df_new.insert(3, 'season', df_new.pop('season'))

    return df_new



""" MAIN PROGRAM """

start = time.time()

# Step 1 : Open catalogue, boundary catalogue and mask

cat_in = ('/home/data/ReAnalysis/ERA5/Storm_analysis/NAECv1/NAEC_1979_2020_v1.csv')
bnd_in = ('/pampa/cloutier/outline_crcm6_domain.csv')
mask_in = ('/pampa/picart/Masks/mask_GEM5_ERA5grid')

print('reading files...')
cat, bnd, mk = open_cat_mask(cat_in, bnd_in, mask_in)
print('files opened')

# Step 2 : Merge cat and mask to add HU column in cat

merge = cat.merge(mk, how='left', on=['latitude', 'longitude'])
merge = merge.fillna(value = False)

# Step 3 : New dataframe that will contain filtered ETCs

df24 = pd.DataFrame(columns = cat.columns) # new dataframe with filtered etc

# create groups for each storm and iterate through them
for storm_id, group in merge.groupby('storm'):
   
    hu_count = group['HU'].sum()

    # We skip storms that have all HU == False values or less than 24 HU == True values in total
    if not group['HU'].any() or hu_count < 24:
        continue

    # within each group, iterate through each storm center (with apply function) 
    # to determine if the given grid point is within subdomain and 
    # has a minimal distance to the boundary > 5 degree
    
    stInDom = group['HU'] & group.apply(
                    lambda row: get_cond(row['latitude'], row['longitude'], bnd) 
                    if row['HU'] else False, axis=1
                        )

    count = 0
    
    # count the consecutive True values in stInDom. We want to keep ETC 
    # that were active for 24 CONSECUTIVE hours or more in CRCM6 subdomain
    for value in stInDom:
        if value:
            count += 1
            if count >= 24:
                df24 = df24.append(group)
                print('Year in process:', df24['datetime'].iloc[-1])
                break
        else:
            count = 0

# Step 5 : Add the season column in dataframe
df24_season = add_season(df24)

# Step 6 : Save dataframe as csv
df24.to_csv('/pampa/cloutier/ERA5_storm_tracks/filtered/etc24_consec_v3.csv', index = False)

print('Execution time = ', time.time() - start)
