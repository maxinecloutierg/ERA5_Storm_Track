# juste un test vite vite pour créer une map de densité test avec mam

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import Circle
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.feature import NaturalEarthFeature
from matplotlib.colors import ListedColormap

# Load the CSV file
data = pd.read_csv('/pampa/cloutier/test_densite_mam.csv')

# Extract necessary data
lats = data['lat']
lons = data['lon']
storms= data['storm_count']

grid_x, grid_y = np.mgrid[lons.min():lons.max():0.25, lats.min():lats.max():0.25]

# Perform nearest neighbor interpolation
storms_interp = griddata((lons, lats), storms, (grid_x, grid_y), method='nearest')

# Create a figure and axes for the map
fig, ax = plt.subplots(figsize=(12, 9), subplot_kw={'projection': ccrs.LambertConformal(
    central_longitude=-80, central_latitude=51)})

# Apply shading with a 250km radius
radius = 250 * 1000  # 250km radius in meters
for lat, lon in zip(lats, lons):
    circle = Circle((lon, lat), radius, fill=False, color='black', alpha=0.5)
    ax = plt.gca()
    ax.add_patch(circle)

ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.BORDERS, edgecolor='black')
ax.add_feature(NaturalEarthFeature(category='cultural', name='admin_0_boundary_lines_land', 
                                   scale='50m', edgecolor='black', facecolor='none'))
ax.add_feature(NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', 
                                   scale='50m', edgecolor='black', facecolor='none'))

# Add colorbar
#cbar = map.colorbar(contour, location='right')
#cbar.set_label('ETC Track Density')

# Show the map
plt.title('ETC Track Density Map')
plt.savefig('/pampa/cloutier/fig/test_densite_mam.png')


