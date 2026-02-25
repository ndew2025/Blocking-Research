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

def ivt_block_comp():
    '''Plot 500mb heights, IVT, 500 mb height gradient strength (raw blocking indices), and filtered (binary) 
    blocking indices composite.'''

    #The directories to pull data from
    gph_file = r"X:/dewhirst/files_for_ivt_block_comp/gph_composite"
    east_ivt_file = r"X:/dewhirst/files_for_ivt_block_comp/east_ivt_composite"
    north_ivt_file = r"X:/dewhirst/files_for_ivt_block_comp/north_ivt_composite"
    raw_block_file = r"X:/dewhirst/files_for_ivt_block_comp/raw_blocking_composite"
    filt_block_file = r"X:/dewhirst/files_for_ivt_block_comp/filt_blocking_composite"

    # Convert meters per second to knots
    ms_to_knots = 1.94384449

    gph = xr.open_dataset(gph_file)['z'] # open 
    #gph = gph.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice
                                
    east_ivt = xr.open_dataset(east_ivt_file)['p71.162']  # get east/west component of ivt    
    #east_ivt = east_ivt.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice                       

    north_ivt = xr.open_dataset(north_ivt_file)['p72.162'] # get north/south component of ivt
    #north_ivt = north_ivt.sel(longitude = slice(250, 300), latitude = slice(50, 20)) # slice

    raw_blocking_indices = xr.open_dataset(raw_block_file)
    B_raw = raw_blocking_indices['B_raw'] # get raw blocking indices
 
    filt_blocking_indices = xr.open_dataset(filt_block_file)
    B_filt = filt_blocking_indices['B_filt']

    # Calculate the magnitudes of 500mb wind and ivt
    ivt = np.sqrt(east_ivt**2 + north_ivt**2)

    lons = np.arange(250, 300.25, 0.25)
    lats = np.arange(50, 19.75, -0.25)

    #Create the figure for the cartopy plot, specify the projection
    projection = ccrs.PlateCarree(central_longitude = 0)
    fig, axes = plt.subplots(2, 2, figsize = (25, 25), subplot_kw={'projection': ccrs.PlateCarree()})

    fig.tight_layout(pad = 5, h_pad = 8) 
    
    for ax in axes.flatten():
        ax.coastlines()
        ax.add_feature(cfeature.BORDERS, linewidth = 0.5, edgecolor = 'black', zorder = 1)
        ax.add_feature(cfeature.LAKES, linewidth = 0.5, edgecolor = 'black', zorder = 2)
        ax.add_feature(cfeature.OCEAN, linewidth = 0.5, edgecolor = 'black', zorder = 1)
        ax.add_feature(cfeature.STATES, linewidth = 0.5, edgecolor = 'black', zorder = 2)

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
    
    Lons, Lats = np.meshgrid(lons, lats)

    #Subsample each field, besides heights, for barbs and quivers
    sub_lats = Lats[::spacing, ::spacing]
    sub_lons = Lons[::spacing, ::spacing]
    sub_east_ivt = east_ivt[::spacing, ::spacing]
    sub_north_ivt = north_ivt[::spacing, ::spacing]

    gph_levs = [540, 546, 552, 558, 564, 570, 576, 582, 588, 594, 600] # in dam
    ivt_levs = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300] # in kg/ms
    raw_block_levs = [-12, -10, -8, -6, -4, -2, 0, 2, 4, 6]
    filt_block_levs = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
       
    #Setup the height contours, colors, and colorbar
    cf1 = axes[0,0].contourf(lons, lats, gph/10, levels = gph_levs, cmap = cmasher.freeze, zorder = 1, extend = 'both')
    cbar1 = plt.colorbar(cf1, pad = 0.01, extend = 'both')
    cl1 = axes[0,0].contour(cf1, colors = 'k', linewidths = 0.5, zorder = 1)
    cbar1.add_lines(cl1)
    cbar1.set_label("Geopotential Height (dam)", fontsize = 16) 
    axes[0,0].set_title(f'500 hPa Heights', fontsize = 22)

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
    cf3 = axes[1, 0].contourf(lons, lats, B_raw, levels = raw_block_levs, cmap='seismic', norm=colors.CenteredNorm(), extend = 'both')
    cbar3 = plt.colorbar(cf3, pad = 0.01, extend = 'both')
    cl3 = axes[1, 0].contour(cf3, colors='k', linewidths = 0.5)
    cbar3.add_lines(cl3)
    cbar3.set_label("Height Gradient ($\\frac{m}{^\circ\\text{lat}}$)", fontsize = 16)
    axes[1, 0].set_title("B", fontsize = 22)

    # Second subplot containing the binary indices
    cf4 = axes[1,1].contourf(lons, lats, B_filt*100, levels = filt_block_levs, cmap= 'Reds')
    cbar4 = plt.colorbar(cf4, pad = 0.01)
    cbar4.set_label("Frequency (%)", fontsize = 16)
    axes[1,1].set_title("Frequency of Positive B of at least 10${^\circ}$Longitude", fontsize = 21)

    # Move the bottom right subplot to the left a little bit
    pos4 = axes[1, 1].get_position()  # Get the original position
    new_pos4 = [pos4.x0 - 0.045, pos4.y0, pos4.width, pos4.height]  # Shift left by 0.05
    axes[1, 1].set_position(new_pos4)  # Set the new position

    # Adjust the colorbar position of the top right plot also
    cbar4_pos = cbar4.ax.get_position()  # Get the original position of the colorbar
    new_cbar4_pos = [cbar4_pos.x0 - 0.045, cbar4_pos.y0, cbar4_pos.width, cbar4_pos.height]  # Shift left by 0.05
    cbar4.ax.set_position(new_cbar4_pos)  # Set the new position of the colorbar


    # Mark the IVT search domain in  the GLR
    axes[0,1].plot([-92, -76], [49, 49], transform=ccrs.PlateCarree(), color = 'blue', linewidth = 5, label = 'IVT > 500 $kgm^{-1}s^{-1}$ Search Area')
    axes[0,1].plot([-92, -76], [41, 41], transform=ccrs.PlateCarree(), color = 'blue', linewidth = 5)
    axes[0,1].plot([-92, -92], [49, 41], transform=ccrs.PlateCarree(), color = 'blue', linewidth = 5)
    axes[0,1].plot([-76, -76], [49, 41], transform=ccrs.PlateCarree(), color = 'blue', linewidth = 5)

    # Mark the domain for the blocking being searched for
    axes[1,1].plot([-100, -80], [40, 40], transform=ccrs.PlateCarree(), color = 'black', linewidth = 5, label = 'Positive B Search Area')
    axes[1,1].plot([-100, -80], [33, 33], transform=ccrs.PlateCarree(), color = 'black', linewidth = 5)
    axes[1,1].plot([-100, -100], [40, 33], transform=ccrs.PlateCarree(), color = 'black', linewidth = 5)
    axes[1,1].plot([-80, -80], [40, 33], transform=ccrs.PlateCarree(), color = 'black', linewidth = 5)

    # Add legends to the subplots
    legend1 = axes[0,1].legend(loc = 'lower right', fontsize = 12)
    legend1.get_frame().set_alpha(1)

    legend2 = axes[1,1].legend(loc = 'lower right', fontsize = 12)
    legend2.get_frame().set_alpha(1)

    fig.savefig(rf"C:/Users/nweat/OneDrive/School/Blocking Research/AMS_Poster_Figures/ivt_block_climo.png", dpi = 500)