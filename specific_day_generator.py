import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import cmasher
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import calendar

plt.ion()

def plot_specific_day(date, hour, great_lakes_region = False):
    '''Plot 500mb heights and winds, IVT, 500 mb height gradient strength (raw blocking indices), and filtered (binary) 
    blocking indices for any day at 18z from 1979-2023.

    Parameters:
    - date --> a string in the format 'yyyy_mm_dd' representing the day you want plots for.
    - hour --> an integer from 0 to 23 (UTC) that specifies the hour of the day at which you want plots.
    - great_lakes_region --> boolean, default is true. Specify whether you want to plot the full domain or just around the GL. '''

    #The directories to pull data from
    #gph_file = rf"X:/dewhirst/1979-2023_500mb_Heights_ERA5/era5_geopotential_{date}_500.nc"
    gph_anom_file = rf"C:/Users/nweat/OneDrive/School/Blocking Research/July 1995 gph anoms/gph_anoms_{date}.nc"
    u_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_u_component_of_wind_{date}_500.nc"
    v_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_v_component_of_wind_{date}_500.nc"
    east_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_eastward_water_vapour_flux_{date}_.nc"
    north_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_northward_water_vapour_flux_{date}_.nc"
    blocking_file = rf"X:/dewhirst/blocking_indices/blocking indices/blocking_indices_{date}.nc"

    # Convert meters per second to knots
    ms_to_knots = 1.94384449

    #gph = xr.open_dataset(gph_file)['z'][hour]/9.81 # open and convert gp to gph
    #gph = gph.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice

    gph_data = xr.open_dataset(gph_anom_file)
    gph_anom = gph_data['z_anom'].sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice 

    lats = gph_anom['latitude'].sel(latitude = slice(50,20)) # get latitudes and longitudes
    lons = gph_anom['longitude'].sel(longitude = slice(250,300))

    u = xr.open_dataset(u_file)['u'][hour]*ms_to_knots # get u component and convert to knots
    u = u.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice 
                                     
    v = xr.open_dataset(v_file)['v'][hour]*ms_to_knots # get v component and convert to knots
    v = v.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice
                                
    east_ivt = xr.open_dataset(east_ivt_file)['p71.162'][hour]  # get east/west component of ivt    
    east_ivt = east_ivt.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice                       

    north_ivt = xr.open_dataset(north_ivt_file)['p72.162'][hour] # get north/south component of ivt
    north_ivt = north_ivt.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice

    blocking_indices = xr.open_dataset(blocking_file)
    B_raw = blocking_indices['B_raw'].sel(longitude = slice(250, 300), latitude = slice(50, 20)) # get raw blocking indices
    B_filt = blocking_indices['B_filt'].sel(longitude = slice(250, 300), latitude = slice(50, 20)) # get filtered blocking_indices

    # Get the name of the month from the day specified to use to search for the correct geopotential height climatology.
    # The other parts of the date string will be used in the title of the figure.
    date_parts = date.split('_')
    year = date_parts[0]
    month = int(date_parts[1])
    day = date_parts[2]
    month_name = calendar.month_name[month].lower()

    # Calculate the magnitudes of 500mb wind and ivt
    winds = np.sqrt(u**2 + v**2)
    ivt = np.sqrt(east_ivt**2 + north_ivt**2)

    #Create the figure for the cartopy plot, specify the projection
    projection = ccrs.PlateCarree(central_longitude = 0)
    fig, axes = plt.subplots(2, 2, figsize = (25, 25), subplot_kw={'projection': ccrs.PlateCarree()})

    fig.tight_layout(pad = 5, h_pad = 8) 
    
    Lons, Lats = np.meshgrid(lons, lats)
    
    for ax in axes.flatten():
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth = 0.5, edgecolor = 'black', zorder = 1)
        ax.add_feature(cfeature.LAKES, linewidth = 0.5, edgecolor = 'black', zorder = 2)
        ax.add_feature(cfeature.OCEAN, linewidth = 0.5, edgecolor = 'black', zorder = 1)
        ax.add_feature(cfeature.STATES, linewidth = 0.5, edgecolor = 'black', zorder = 2)

    if great_lakes_region is True:
        #Domain around the Great Lakes to plot
        west = 265
        east = 285
        south = 37
        north = 53
        for ax in axes.flatten():
            #Add grid lines for latitude and longitude
            grid_lines = ax.gridlines(draw_labels = True, xlocs = [-90, -85, -80], ylocs = [40,45,50], linestyle = '--', color = 'white', zorder = 3)
            grid_lines.xlabel_style = {'size': 16, 'color': 'black'}
            grid_lines.ylabel_style = {'size': 16, 'color': 'black'}

            #Get rid of grid line labels on the top and right of the plot
            grid_lines.top_labels = False 
            grid_lines.right_labels = False 

            #Set the extent just around the GLR
            ax.set_extent([west, east, south, north], crs = projection)
            
        #Set the spacing between the quivers and barbs for winds and ivt
        spacing = 8

    if great_lakes_region is False:
        spacing = 14

        for ax in axes.flatten():
            #Add grid lines for latitude and longitude
            grid_lines = ax.gridlines(draw_labels = True, xlocs = [-105, -95, -85, -75, -65], ylocs = [25,35,45], 
                                      linestyle = '--', lw = 0.5, color = 'white', zorder = 3)
            grid_lines.xlabel_style = {'size': 16, 'color': 'black'}
            grid_lines.ylabel_style = {'size': 16, 'color': 'black'}

            #Get rid of grid line labels on the top and right of the plot
            grid_lines.top_labels = False 
            grid_lines.right_labels = False 

    #Subsample each field, besides heights, for barbs and quivers
    sub_lats = Lats[::spacing, ::spacing]
    sub_lons = Lons[::spacing, ::spacing]
    sub_u = u[::spacing, ::spacing]
    sub_v = v[::spacing, ::spacing]
    sub_east_ivt = east_ivt[::spacing, ::spacing]
    sub_north_ivt = north_ivt[::spacing, ::spacing]

    # Get the maximum blocking index for each day and plot as a star
    B_max = np.max(B_raw)

    # Get the coordinates of the maximum blocking index within the condition specified above
    B_max_lat_idx = np.where(B_raw == B_max)[0][0]
    B_max_lon_idx = np.where(B_raw == B_max)[1][0]

    #gph_levs = [540, 546, 552, 558, 564, 570, 576, 582, 588, 594, 600] # in dam
    gph_anom_levs = np.arange(-200, 250, 50) # in meters
    ivt_levs = [0, 75, 150, 225, 300, 375, 450, 525, 600, 675, 750] # in kg/ms
    block_levs = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6]
       
    #Setup the height contours, colors, and colorbar
    cf1 = axes[0,0].contourf(lons, lats, gph_anom, levels = gph_anom_levs, cmap = 'seismic', norm=colors.CenteredNorm(), 
                             zorder = 1, extend = 'both')
    cbar1 = plt.colorbar(cf1, pad = 0.01, extend = 'both')
    cl1 = axes[0,0].contour(cf1, colors = 'k', linewidths = 0.5, zorder = 1)
    cbar1.add_lines(cl1)
    cbar1.set_label("500 Anomaly (m)", fontsize = 16) 
    axes[0,0].set_title(f'Z500 Anomalies', fontsize = 22) # not plotting winds for now
    #axes[0,0].barbs(sub_lons, sub_lats, sub_u, sub_v, length = 6.5, color = 'white', zorder = 2)

    cf2 = axes[0,1].contourf(lons, lats, ivt, levels = ivt_levs, cmap = cmasher.torch, zorder = 1, extend = 'max')
    cbar2 = plt.colorbar(cf2, pad = 0.01, extend = 'max')
    #cl2 = ax.contour(cf2, colors = 'k', linewidths = 0.5, zorder = 1)
    #cbar2.add_lines(cl2)
    cbar2.set_label("IVT ($kgm^{-1}s^{-1}$)", fontsize = 16)  
    axes[0,1].set_title(f'Integrated Vapor Transport', fontsize = 22)
    axes[0,1].quiver(sub_lons, sub_lats, sub_east_ivt, sub_north_ivt, scale = 7500, color = 'white', zorder = 2)

    # Move the top right subplot to the left a little bit
    pos2 = axes[0, 1].get_position()  # Get the original position
    new_pos2 = [pos2.x0 - 0.045, pos2.y0, pos2.width, pos2.height]  # Shift left by 0.05
    axes[0, 1].set_position(new_pos2)  # Set the new position

    # Adjust the colorbar position of the top right plot also
    cbar2_pos = cbar2.ax.get_position()  # Get the original position of the colorbar
    new_cbar2_pos = [cbar2_pos.x0 - 0.045, cbar2_pos.y0, cbar2_pos.width, cbar2_pos.height]  # Shift left by 0.05
    cbar2.ax.set_position(new_cbar2_pos)  # Set the new position of the colorbar

    # First subplot containing the raw indices
    cf3 = axes[1, 0].contourf(lons, lats, B_raw, levels = block_levs, cmap='seismic', norm=colors.CenteredNorm(), extend = 'both')
    cbar3 = plt.colorbar(cf3, pad = 0.01, extend = 'both')
    cl3 = axes[1, 0].contour(cf3, colors='k', linewidths = 0.5)
    cbar3.add_lines(cl3)
    cbar3.set_label("Height Gradient ($\\frac{m}{^\circ\\text{lat}}$)", fontsize = 16)
    axes[1, 0].set_title("B", fontsize = 22)

    #Plot the maximum blocking index as a star
    axes[1,0].plot((360-lons[B_max_lon_idx])*-1, lats[B_max_lat_idx], 'w*', markersize = 10, zorder = 5)

    # Second subplot containing the binary indices
    cmap = ListedColormap(['blue', 'red'])
    cf4 = axes[1,1].contourf(lons, lats, B_filt, cmap=cmap)
    cl4 = axes[1,1].contour(cf4, colors='k', vmin=0, vmax=1, linewidths=0.5)
    axes[1,1].set_title(" Positive B of at least 10${^\circ}$Longitude", fontsize = 22)

    # Move the bottom right subplot to the left a little bit
    pos4 = axes[1, 1].get_position()  # Get the original position
    new_pos4 = [pos4.x0 - 0.045, pos4.y0, pos4.width, pos4.height]  # Shift left by 0.05
    axes[1, 1].set_position(new_pos4)  # Set the new position

    # Add a title
    fig.suptitle(f'{date}')

    fig.savefig(rf"C:/Users/nweat/OneDrive/School/Blocking Research/AMS_Poster_Figures/{date}.png", dpi = 500)




    return fig


for i in np.arange(1,32):
    fig = plot_specific_day(f'1995_07_{i:02d}', 18)
    fig.savefig(rf"C:/Users/nweat/OneDrive/School/Blocking Research/July 1995 Plots/Day {i}.png")

#fig = plot_specific_day('2023_08_21', 18)
