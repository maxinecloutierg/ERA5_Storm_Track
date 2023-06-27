"""
Calculate density around a 250km radius

"""

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


df = pd.read_csv('/pampa/cloutier/density.csv')
mam = df.loc[df.season == 'MAM']

# test de la densité avec seulement MAM
# Avec la moyenne pas bonne de 250km
r = 2.25
mam_density = {} # will contain the grid point and density

# iterate through all grid point
for _, row in mam.iterrows() : 
    track = {} #  will contain the year and storm id of the grid points found
    lat1 = row['latitude'] 
    lon1 = row['longitude'] 

    # restricting the research in a 2.5 x 2.5deg square
    neighbors = mam[
        (mam['latitude'] >= lat1 - r) & (mam['latitude'] <= lat1 + r) &
        (mam['longitude'] >= lon1 - r) & (mam['longitude'] <= lon1 + r)
        ]

    # Find the grid points that are within a 250km radius 
    for _, row in neighbors.iterrows():
        lat2 = row['latitude'] 
        lon2 = row['longitude'] 
        storm_count = row['storm_count']
        avgvors = row['avgVORSmax']

        # calculate distance between grid points
        distance = calculate_distance(lat1, lon1, lat2, lon2)

        # add year and storm_id of grid points that are within the radius
        if distance <= 250 : 
            track['storm_count'] = track.get('storm_count', []) + [storm_count]

    # turn track dictionnary into dataframe and count the unique occurence of 
    # every storm for each year
    track_df = pd.DataFrame(data = track)
    #track_df = track_df.groupby('year')['storm_id'].nunique()

    # number of storm per season in average
    average = track_df['storm_count'].mean()

    # add coord and average density
    mam_density['lat'] = mam_density.get('lat', []) + [lat1]
    mam_density['lon'] = mam_density.get('lon', []) + [lon1]
    mam_density['storm_count'] = mam_density.get('storm_count', []) + [average]
    mam_density['avgvors'] = mam_density.get('avgvors', []) + [avgvors]

# transfort m_density into dataframe and save as csv
df_density = pd.DataFrame(data = mam_density)
df_density.to_csv('/pampa/cloutier/test_densite_mam.csv')
