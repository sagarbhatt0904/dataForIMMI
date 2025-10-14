# Code to find the most common grain Id in each block and assigning the grain that orientation.

#%% Load modules
import numpy as np
import pandas as pd
import netCDF4 as nc
from scipy.spatial.distance import cdist
from scipy.stats import mode
import time
import sys
#%% 
# Some functions
def normalize_coordinates(original_value, min_value, max_value, new_min, new_max):
    normalized_value = (original_value - min_value) / (max_value - min_value) * (new_max - new_min) + new_min
    return normalized_value

# Check if the correct number of command-line arguments is provided
if len(sys.argv) != 6:
    print("Usage: python orientationsForEachGrain.py <mesh-file.e> <exaca-input.csv> <ori-file-base> <starting-block-#> <ending-block-#>")
    sys.exit(1)

meshFile=sys.argv[1]
exacaFile=sys.argv[2]

#%% Read the mesh
mesh = nc.Dataset(meshFile)
for name, dimension in mesh.dimensions.items():
    if name in ['num_el_blk']:
        num_blocks = dimension.size
        break
for name, dimension in mesh.dimensions.items():
    if name in ['num_elem']:
        num_elem = dimension.size
        break

# First lets get the vertices of the element
coordinates = np.transpose(np.array([mesh.variables['coordx'][:].data, mesh.variables['coordy'][:].data, mesh.variables['coordz'][:].data]))
#%% 
# Read the point data
data = pd.read_csv(exacaFile)


#%%
# change the coordinates of the original data to lie in the range for the coordinates in the mesh

# Define min and max values for each coordinate
min_x, max_x = data['X (m)'].min(), data['X (m)'].max()
min_y, max_y = data['Y (m)'].min(), data['Y (m)'].max()
min_z, max_z = data['Z (m)'].min(), data['Z (m)'].max()

# Define new min and max values for normalization (0 and 100)
new_min_x, new_max_x = np.min(mesh.variables['coordx'][:].data), np.max(mesh.variables['coordx'][:].data)
new_min_y, new_max_y = np.min(mesh.variables['coordy'][:].data), np.max(mesh.variables['coordy'][:].data)
new_min_z, new_max_z = np.min(mesh.variables['coordz'][:].data), np.max(mesh.variables['coordz'][:].data)

# Normalize X, Y, Z coordinates
data['X (m)'] = data['X (m)'].apply(lambda x: normalize_coordinates(x, min_x, max_x, new_min_x, new_max_x))
data['Y (m)'] = data['Y (m)'].apply(lambda y: normalize_coordinates(y, min_y, max_y, new_min_y, new_max_y))
data['Z (m)'] = data['Z (m)'].apply(lambda z: normalize_coordinates(z, min_z, max_z, new_min_z, new_max_z))


#%% Get centroid of the elements by block
# Since the objective is to get one orientation for each block. A quick and dirty way, just compute the centroid of each element, find the closest point, assign that orientation. Get the mode for all elements in the block and assign the block that orientation. 

t0 = time.time()

# Specify the tolerance
tolerance = 1
euler = np.zeros((num_blocks,3))
el_index = 0
for blk in range(int(sys.argv[4]), int(sys.argv[5])):
    if blk==0:
        continue
    if blk>num_blocks:
        continue
    num_el_in_blk = 0
    connectivity = np.array(mesh.variables['connect'+str(blk)][:])
    num_elems_in_blk = connectivity.shape[0]
    gids_elem_in_blk = []
    for element_idx in range(0,num_elems_in_blk):
        element_vertices = coordinates[connectivity[element_idx] - 1, :]
        centroid = np.mean(element_vertices, axis=0)
         # Generate min and max values for filtering points
        min_values = centroid - tolerance
        max_values = centroid + tolerance

        # Filter points within the specified range
        filtered_points = data[
            (data['X (m)'] >= min_values[0]) & (data['X (m)'] <= max_values[0]) &
            (data['Y (m)'] >= min_values[1]) & (data['Y (m)'] <= max_values[1]) &
            (data['Z (m)'] >= min_values[2]) & (data['Z (m)'] <= max_values[2])
        ]
        if not filtered_points.empty:
            # Calculate distances only for the filtered points
            distances = cdist([centroid], filtered_points[['X (m)', 'Y (m)', 'Z (m)']])

            # Find the nearest point's index and gid among the filtered points
            nearest_point_index_within_filtered = np.argmin(distances)
            nearest_point_gid = np.abs(filtered_points.iloc[nearest_point_index_within_filtered]['gid'])%10000
            gids_elem_in_blk.append(nearest_point_gid)
    
    blkGid = mode(gids_elem_in_blk)[0]
    first_occurrence = data.loc[np.abs(data['gid'])%10000 == int(blkGid)].iloc[0]
    euler[blk-1,:] = np.rad2deg([first_occurrence['phi1'], first_occurrence['Phi'], first_occurrence['phi2']])

t1 = time.time()
print("Time taken: ", t1-t0)
print("Ori for ",int(sys.argv[4])," to ", int(sys.argv[5])," done!!!")
oriFileName=sys.argv[3]+"_"+sys.argv[4]+"_"+sys.argv[5]+".txt"
np.savetxt(oriFileName,euler, delimiter=' ',fmt='%.6f')

