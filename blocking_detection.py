
'''This file contains code to calculate blocking indices blocks based off geopotential height data'''

# Import necessary Python modules and libraries
import xarray as xr
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import calendar
import time
import psutil

plt.ion() # Turn on interactive mode

### This function should probably go after the blocking indicies have been calculated, although as they are calculated, I will want to place them into
### their respective grid cell.

def create_blocking_index_netcdf(filepath, lats, lons, year, month, day):
    '''Creates a new netcdf file containing the blocking indices for each day analyzed in the blocking detection function.
    These files will later be used to plot a blocking frequency.'''
    
    # Open the existing NetCDF file
    current_dataset = nc.Dataset(filepath, 'r')
    
    # Create a new NetCDF file
    blocking_index_directory = r"C:/Users/nweat/OneDrive/School/Blocking Research/Blocking Index Dataset"
    blocking_index_path = os.path.join(blocking_index_directory, f'blocking_indices_{year}_{month}_{day}.nc')
    blocking_index_dataset = nc.Dataset(blocking_index_path, 'w', format = 'NETCDF3_CLASSIC')
    
    # Define the dimensions in the new file
    lat_dim = blocking_index_dataset.createDimension('latitude', len(lats))
    lon_dim = blocking_index_dataset.createDimension('longitude', len(lons))

    #Fill value for coordinates outside of the range the blocking indices are being calculated over
    fill_value = np.nan
    
    # Define the coordinate variables in the new file
    lat_var = blocking_index_dataset.createVariable('latitude', np.float32, ('latitude',))
    lon_var = blocking_index_dataset.createVariable('longitude', np.float32, ('longitude',))
    B_raw = blocking_index_dataset.createVariable('B_raw', np.float32, ('latitude', 'longitude'), fill_value = fill_value)
    B_filt = blocking_index_dataset.createVariable('B_filt', np.float32, ('latitude', 'longitude'), fill_value = fill_value)

    # Copy the coordinate data to the new file
    lat_var[:] = lats
    lon_var[:] = lons
    
    return blocking_index_dataset


def blocking_indices(input_directory, num_consec_lons = 10, lon_range = [245, 315], lat_range = [60, 10]):
    ''' The only parameter for this function is the input directory which contains the files of days you want to detect blocks in. 
    Initially, I'll be looking at JJA only. This function will create/return multiple things, including new netcdf files for each day storing blocking indices.
    I also want to create a time series of blocking indices in the domain we are looking at for each period. Finally, there will be a folder dedicated to each 
    blocking event. 
    
    PARAMETERS
    ==========
    
    input_directory: string
        The path to the directory containing all of the geopotential height data
    
    num_consec_lons: integer
        The threshold of consecutive longitudes used to identify a block
    
    lat_range, lon_range: two-element lists
        Specify the range of latitude and longitudes to detect blocking across.
        Default is lon_range = [245, 315], lat_range = [60, 10] for current area of interest
    '''

    start = time.time() # start time

    # The change in latitude to calculate the north and south geopotential height gradients across
    delta_phi = 20  

    #The range of phi_0 I am interested in to search for blocking
    latitude_range = np.arange(lat_range[1] + delta_phi/2, lat_range[0] - delta_phi/2 + 0.25, 0.25)

    # Adjust longitude threshold for filtered based on data ERA5 data resolution
    num_points = num_consec_lons*4 + 1 # ensure that num_points is odd so that the search around it executes smoothly

    # floor division for number of cells to search on either side of the longitude
    filt_search = num_points // 2

    #Biggest loop will look through all the files in the directory.
    for file in os.listdir(input_directory):
        filepath = os.path.join(input_directory, file)
        data = xr.open_dataset(filepath) # open each file
        specified_data = data.sel(latitude = slice(lat_range[0], lat_range[1]), # slice each into specified domain
                                  longitude = slice(lon_range[0], lon_range[1])) 
        
        z500 = specified_data['z'][18]/9.81 #The 18th index will get data from 18z

        #Access the sliced longitudes and latitudes so that the new netcdf files have the correct domain
        lons = specified_data['longitude']
        lats = specified_data['latitude']

        #Get the different parts of the filename for naming the new netcdf files
        file_parts = filepath.split('_')
        year, month, day = file_parts[5], file_parts[6], file_parts[7]
    
        # Call the file creator function above
        blocking_index_dataset = create_blocking_index_netcdf(filepath, lats, lons, year, month, day)
        raw_indices = blocking_index_dataset['B_raw']  # Access raw blocking index variable
        filtered_indices = blocking_index_dataset['B_filt']  # Access filtered blocking index variable

        print("File has been created!")

        # The following set of loops is for calculating the raw blocking indices (prior to filtering)
    
        for phi0 in latitude_range: # Loop through the range of phi0 
            phi0_index = np.where(lats == phi0)[0][0] # Index where latitude is phi0
    
            #Nest a loop through the longitudes. For each line of latitude, every longitude is accounted for.
            for lon_index in range(lons.size):
                #Access all of the heights 10 degrees north and south of phi_0
                north_lat_indices = (lats >= phi0) & (lats <= phi0 + delta_phi/2) # the geopotential height indices north of phi_0
                south_lat_indices = (lats >= phi0 - delta_phi/2) & (lats <= phi0) # the geopotential height indices south of phi_0
                
                #Use the indices to index into z500
                north_z = z500[north_lat_indices, lon_index]
                south_z = z500[south_lat_indices, lon_index]
    
                #Take the sum of all the components north and south, then divide by N to get an average
                Zn = np.sum(north_z)/north_z.size
                Zs = np.sum(south_z)/south_z.size
    
                #Calculate the blocking index: the averaged meridional height gradient
                B = (Zn - Zs)/delta_phi

                raw_indices[phi0_index, lon_index] = B # Assign the calculated indices to their correct position 

        # The next set of loops is for filtering the blocking indices to narrow the blocked region

        # Setup a range of longitudes to search through, eliminates edgecases
        lon_indices = range(filt_search, lons.size - filt_search, 1) # range to loop through so that filt_search stays within the grid

        for phi0 in latitude_range: # Loop through range of phi0
            phi0_index = np.where(lats == phi0)[0][0] # Index where latitude is phi0
            
            ## OLD CODE TO PRESERVE ORIGINAL METHODS
            for i in lon_indices:  # For each line of latitude, every longitude is accounted for
                if np.all(raw_indices[phi0_index, i - filt_search: i + filt_search] > 0):
                    filtered_indices[phi0_index, i] = 1
                else:
                    filtered_indices[phi0_index, i] = 0
                
        blocking_index_dataset.close() # close the dataset
    print("All files complete!")   

    end = time.time()
    fxn_time = end - start
    print(f'Time it took the function to run: {round(fxn_time)} secs') 


def plot_indices(year, month, day):
    '''Plot both the raw and filtered binary blocking to verify that the functions above are working correctly
    
    PARAMETERS
    ==========

    year, month, day: integers
        used to access the proper file to access the indices from
    '''
    # local directory that contains all of the blocking indices
    directory = "C:/Users/nweat/OneDrive/School/Blocking Research/blocking indices"

    for file in os.listdir(directory):
        if f'{year}_{month:02d}_{day:02d}' in file:
            filepath = os.path.join(directory, file)
            data = xr.open_dataset(filepath)
            B_raw = data['B_raw']
            B_filt = data['B_filt']
            lats = data['latitude']
            lons = data['longitude']

    # Get the name of the month from the day specified for the figure
    month_name = calendar.month_name[month].lower()
    
    # Create a figure and projection for the plots
    projection = ccrs.PlateCarree(central_longitude=0)
    fig, axs = plt.subplots(1, 2, figsize = (16,7), subplot_kw={'projection': projection})

    # Add some cartopy details to each subplot
    for ax in axs:  # For each subplot
        ax.coastlines()  # Add Cartopy details
        ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='black')
        ax.add_feature(cfeature.LAKES, linewidth=0.5, edgecolor='black')
        ax.add_feature(cfeature.OCEAN, linewidth=0.5, edgecolor='black')
        ax.add_feature(cfeature.STATES, linewidth=0.5, edgecolor='black')

        # Add grid lines for latitude and longitude
        grid_lines = ax.gridlines(draw_labels=True, xlocs=[-105, -95, -85, -75, -65, -55], ylocs=[20, 30, 40, 50],
                                linestyle='--', color = 'white', zorder = 2)
        grid_lines.xlabel_style = {'size': 10, 'color': 'black'}
        grid_lines.ylabel_style = {'size': 10, 'color': 'black'}

        # Get rid of grid line labels on the top and right of the plot
        grid_lines.top_labels = False
        grid_lines.right_labels = False
    
    levels = np.arange(-16, 10, 2)

    # First subplot containing the raw indices
    cf0 = axs[0].contourf(lons, lats, B_raw, levels = levels, cmap='seismic', norm=colors.CenteredNorm())
    cbar0 = plt.colorbar(cf0, pad = 0.01, shrink = 0.7)
    cl0 = axs[0].contour(cf0, colors='k', linewidths = 0.5)
    cbar0.add_lines(cl0)
    cbar0.set_label("Blocking Strength ($\\frac{m}{^\circ\\text{lat}}$)")
    axs[0].set_title("Blocking Strength from 500mb Height Gradient")

    # Second subplot containing the binary indices
    cmap = ListedColormap(['blue', 'red'])
    cf1 = axs[1].contourf(lons, lats, B_filt, cmap=cmap)
    cl1 = axs[1].contour(cf1, colors='k', vmin=0, vmax=1, linewidths=0.5)
    axs[1].set_title("Binary Blocking Indices")

    axs[1].text(0.89, 0.93, "Blocked", transform = axs[1].transAxes, fontsize = 12, backgroundcolor = 'red',
                ha = 'center', va = 'center', zorder = 3)
    axs[1].text(0.89, 0.87, "Not Blocked", transform = axs[1].transAxes, fontsize = 12, backgroundcolor = 'blue',
                ha = 'center', va = 'center', zorder = 3)

    # Center the main title over the figure and make it larger
    fig.suptitle(f'{month_name.capitalize()} {day}, {year} at 18z', fontsize=16, ha ='center', y = 0.95)
    fig.tight_layout()

    return fig
    
#for date in np.arange(15,31):
#    fig = plot_indices(2023, 8, date)
#    fig.savefig(f"C:/Users/nweat/OneDrive/School/Blocking Research/Aug 19-25, 2023 Blocking Event/Aug{date}.png", dpi = 300)
#    plt.close()