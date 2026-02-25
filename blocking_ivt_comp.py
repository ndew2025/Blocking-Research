import numpy as np
import xarray as xr
import os
import shutil
import datetime

def get_JJA_files():
    '''Get blocking and ivt files in JJA only and move to apropriate folders'''

    #blocking_dir = r"X:/dewhirst/blocking_indices/blocking indices"
    ivt_dir = r"X:/dewhirst/1979-2023 U, V, IVT Components ERA5"

    #blocking_JJA_dir = r"X:/dewhirst/JJA_blocking_indices"
    ivt_JJA_dir = r"X:/dewhirst/JJA_ivt"

    # List all NetCDF files in the folder
    #blocking_files = [os.path.join(blocking_dir, f) for f in os.listdir(blocking_dir)]
    ivt_files = [os.path.join(ivt_dir, f) for f in os.listdir(ivt_dir) if "water" in f]

    #for file in blocking_files:
    #    # Extract the month from the filename
    #    filename = os.path.basename(file)
    #    month = filename.split('_')[-2]
    #    if month in ['06', '07', '08']:
            # Copy the file to the destination directory
    #        shutil.copy(file, blocking_JJA_dir)
    
    for file in ivt_files:
        # Extract the month from the filename
        filename = os.path.basename(file)
        month = filename.split('_')[-3]
        if month in ['06', '07', '08']:
            # Copy the file to the destination directory
            shutil.copy(file, ivt_JJA_dir)


# look for blocking in the Southeast, southern US in JJA
# get instances where ivt is greater than 250 kg/ms in the Great Lakes region in JJA
#       Use a box that approximates the Great Lakes Basin 
# look for blocking two days before, one day before, or on the day of the ivt threshold instances at 18z
#       any of these days!

# first find all files with ivt > 250 kg/ms in the great lakes region and all files with blocking in the S/SE US
# by opening them, searching for the ivt and blocking (B_filt = 1) in a sliced dataset, then append all of them to separate lists

# Check the names of each file, split?, to see which ones match and how many. 

# Calculate the frequency of blocks with 250 ivt in GLR (divide by number of all blocking detection in S/SE), 
# and the frequency of 250 ivt in GLR with blocks (divide by number of ivt events above 250 in GLR)

# Climatology of this sort of event? Would be really cool, the 500mb heights, raw and filtered blocking indices, and ivt
#       Which day to use? 
#            Has to be the day of the ivt event, the blocking event also regardless if blocking is actually present 
#            on the same day as the ivt


def ivt_block_overlap():
    
    raw_blocking_datasets = []
    filt_blocking_datasets = []
    east_ivt_datasets = []
    north_ivt_datasets = []
    gph_datasets = []

    # Loop through a list of times, find the files that match the date and get the blocking files from the same date,
    # and one and two days before --> also will need to add blocking files for the last 2 days of may for it to work properly
    
    # Loop through each year
    for year in range(1979, 2024): # original loop started at 1979 for all data
        # Loop through June, July, and August
        for month in [6, 7, 8]:
            # Loop through each day of the current month
            for day in range(1, 32):
                try:
                    date = datetime.date(year, month, day)
                    one_day_bef = date - datetime.timedelta(days=1)
                    two_days_bef = date - datetime.timedelta(days=2)

                    # Format the date as yyyy_mm_dd
                    date = date.strftime('%Y_%m_%d')
                    one_day_bef = one_day_bef.strftime('%Y_%m_%d')
                    two_days_bef = two_days_bef.strftime('%Y_%m_%d')

                except ValueError:
                    # If day does not exist, skip to next iteration
                    continue

                # Open all of the necessary data: both ivt components and blocking data 0-2 days day before the current date
                # Get 18z data (needed for ivt) and slice into the subdomains for searching
                east_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_eastward_water_vapour_flux_{date}_.nc"
                east_ivt_data = xr.open_dataset(east_ivt_file)
                east_ivt = east_ivt_data['p71.162'][18].sel(latitude = slice(49, 41), longitude = slice(268, 284))

                north_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_northward_water_vapour_flux_{date}_.nc"
                north_ivt_data = xr.open_dataset(north_ivt_file)
                north_ivt = north_ivt_data['p72.162'][18].sel(latitude = slice(49, 41), longitude = slice(268, 284))

                blocking_file = rf"X:/dewhirst/blocking_indices/blocking indices/blocking_indices_{date}.nc"
                day_of_blocking_data = xr.open_dataset(blocking_file)
                day_of_blocking = day_of_blocking_data['B_filt'].sel(latitude = slice(40, 33), longitude = slice(260, 280))

                one_day_bef_blocking_data = xr.open_dataset(rf"X:/dewhirst/blocking_indices/blocking indices/blocking_indices_{one_day_bef}.nc")
                one_day_bef_blocking = one_day_bef_blocking_data['B_filt'].sel(latitude = slice(40, 33), longitude = slice(260, 280))

                two_days_bef_blocking_data = xr.open_dataset(rf"X:/dewhirst/blocking_indices/blocking indices/blocking_indices_{two_days_bef}.nc")
                two_days_bef_blocking = two_days_bef_blocking_data['B_filt'].sel(latitude = slice(40, 33), longitude = slice(260, 280))

                # Also get the gph file for the climatology later
                gph_file = rf"X:/dewhirst/1979-2023 500mb Heights ERA5/era5_geopotential_{date}_500.nc"
                gph = xr.open_dataset(gph_file)['z'][18].sel(longitude = slice(250, 300), latitude = slice(50, 20))/9.81 # divide by g to get gph
                
                # Calculate magnitude of ivt
                ivt_mag = np.sqrt(east_ivt**2 + north_ivt**2)

                # If IVT > 500 is present in the GLR and blocking is present 0-2 days (either of them) before
                if np.any(ivt_mag > 500) and (np.any(day_of_blocking == 1) or np.any(one_day_bef_blocking == 1) 
                                              or np.any(two_days_bef_blocking == 1)):
                    
                    # Grab a larger domain for the composite plots
                    raw_to_append = day_of_blocking_data['B_raw'].sel(latitude = slice(50, 20), longitude = slice(250, 300))
                    filt_to_append = day_of_blocking_data['B_filt'].sel(latitude = slice(50, 20), longitude = slice(250, 300))
                    east_ivt_append = east_ivt_data['p71.162'][18].sel(latitude = slice(50, 20), longitude = slice(250, 300))
                    north_ivt_append = north_ivt_data['p72.162'][18].sel(latitude = slice(50, 20), longitude = slice(250, 300))

                    # append the data to the lists
                    raw_blocking_datasets.append(raw_to_append)
                    filt_blocking_datasets.append(filt_to_append)
                    east_ivt_datasets.append(east_ivt_append)
                    north_ivt_datasets.append(north_ivt_append) 
                    gph_datasets.append(gph) # append the needed files into lists
                
                # Close all of the files at the end of an iteration
                east_ivt_data.close()
                north_ivt_data.close()
                day_of_blocking_data.close()
                one_day_bef_blocking_data.close()
                two_days_bef_blocking_data.close()
                gph.close()

                print(f"{date} is complete!")

    # Concatenate the raw blocking data along the desired dimension (e.g., 'time')
    combined_raw = xr.concat(raw_blocking_datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_raw.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/raw_blocking_concat")

    summed_raw = combined_raw.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_raw.sizes['time'] # get the number of timesteps
    raw_block_comp = summed_raw / num_time_steps # compute the composite
    raw_block_comp.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/raw_blocking_composite")



    # Concatenate the filtered blocking data along the desired dimension (e.g., 'time')
    combined_filt = xr.concat(filt_blocking_datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_filt.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/filt_blocking_concat")

    summed_filt = combined_filt.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_filt.sizes['time'] # get the number of timesteps
    filt_block_comp = summed_filt / num_time_steps # compute the composite
    filt_block_comp.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/filt_blocking_composite")



    # Concatenate the east ivt data along the desired dimension (e.g., 'time')
    combined_east_ivt = xr.concat(east_ivt_datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_east_ivt.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/east_ivt_concat")

    summed_east_ivt = combined_east_ivt.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_east_ivt.sizes['time'] # get the number of timesteps
    east_ivt_comp = summed_east_ivt / num_time_steps # compute the composite
    east_ivt_comp.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/east_ivt_composite")



    # Concatenate the north ivt data along the desired dimension (e.g., 'time')
    combined_north_ivt = xr.concat(north_ivt_datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_north_ivt.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/north_ivt_concat")

    summed_north_ivt = combined_north_ivt.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_north_ivt.sizes['time'] # get the number of timesteps
    north_ivt_comp = summed_north_ivt / num_time_steps # compute the composite
    north_ivt_comp.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/north_ivt_composite")



    # Concatenate the raw blocking data along the desired dimension (e.g., 'time')
    combined_gph = xr.concat(gph_datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_gph.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/gph_concat")

    summed_gph = combined_gph.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_gph.sizes['time'] # get the number of timesteps
    block_comp = summed_gph / num_time_steps # compute the composite
    block_comp.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/gph_composite")

    print("All files have been processed!")



def ivt_block_climo():
    '''Calculate a composite for blocking and ivt overlapping events in the eastern US'''

    block_dir = r"X:/dewhirst/blocking_with_ivt_JJA"
    ivt_dir = r"X:/dewhirst/ivt_with_blocking_JJA"
    gph_dir = r"X:/dewhirst/gph_for_ivt_JJA_climo"

    # List all NetCDF files in the folder
    blocking_files = [os.path.join(block_dir, f) for f in os.listdir(block_dir)]
    east_ivt_files = [os.path.join(ivt_dir, f) for f in os.listdir(ivt_dir) if "eastward" in f]
    north_ivt_files = [os.path.join(ivt_dir, f) for f in os.listdir(ivt_dir) if "northward" in f]
    gph_files = [os.path.join(gph_dir, f) for f in os.listdir(gph_dir)]

    print("Files have been grouped.")

    # # # Composite for raw blocking indices
    # datasets = []
    # for file in blocking_files:
    #     ds = xr.open_dataset(file)
    #     # Extract the variable of interest and assign it to a new DataArray
    #     variable = ds['B_raw'].sel(longitude = slice(250, 300), latitude = slice(50, 20))
    #     datasets.append(variable)

    #     ds.close() # close each file after opening
    
    # # Concatenate the DataArrays along the desired dimension (e.g., 'time')
    # combined_data = xr.concat(datasets, dim = 'time')

    # # Optionally, save the combined data to a new NetCDF file
    # combined_data.to_netcdf(r"X:/dewhirst/files_for_ivt_block_composite/raw_indices_concat")

    # summed_indices = combined_data.sum(dim = 'time') # sum up all of the indices
    # num_time_steps = combined_data.sizes['time'] # get the number of timesteps
    # freq = summed_indices / num_time_steps # compute normalized frequency (annual percentage)
    # freq.to_netcdf(r"X:/dewhirst/files_for_ivt_block_composite/raw_indices_composite")

    # print("Raw indices are complete!")


    # #Composite for filtered blocking indices
    # datasets = []
    # for file in blocking_files:
    #     print()
    #     ds = xr.open_dataset(file)
    #     # Extract the variable of interest and assign it to a new DataArray
    #     variable = ds['B_filt'].sel(longitude = slice(250, 300), latitude = slice(50, 20))
    #     datasets.append(variable)

    #     ds.close() # close each file after opening
    
    # # Concatenate the DataArrays along the desired dimension (e.g., 'time')
    # combined_data = xr.concat(datasets, dim = 'time')

    # # Optionally, save the combined data to a new NetCDF file
    # combined_data.to_netcdf(r"X:/dewhirst/files_for_ivt_block_composite/filtered_indices_concat")

    # summed_indices = combined_data.sum(dim = 'time') # sum up all of the indices
    # num_time_steps = combined_data.sizes['time'] # get the number of timesteps
    # freq = summed_indices / num_time_steps # compute normalized frequency (annual percentage)
    # freq.to_netcdf(r"X:/dewhirst/files_for_ivt_block_composite/filtered_indices_composite")

    # print("Filtered indices are complete!")


    # Composite for west/east ivt
    datasets = []
    for file in east_ivt_files:
        ds = xr.open_dataset(file)
        # Extract the variable of interest and assign it to a new DataArray
        variable = ds['p71.162'][18].sel(longitude = slice(250, 300), latitude = slice(50, 20))
        datasets.append(variable)

        ds.close() # close each file after opening

    # Concatenate the DataArrays along the desired dimension (e.g., 'time')
    combined_data = xr.concat(datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_data.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/east_ivt_concat")

    summed_indices = combined_data.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_data.sizes['time'] # get the number of timesteps
    freq = summed_indices / num_time_steps # compute normalized frequency (annual percentage)
    freq.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/east_ivt_composite")

    print("East/west ivt is complete!")


    # Composite for north/south ivt
    datasets = []
    for file in north_ivt_files:
        ds = xr.open_dataset(file)
        # Extract the variable of interest and assign it to a new DataArray
        variable = ds['p72.162'][18].sel(longitude = slice(250, 300), latitude = slice(50, 20))
        datasets.append(variable)

        ds.close() # close each file after opening
    
    # Concatenate the DataArrays along the desired dimension (e.g., 'time')
    combined_data = xr.concat(datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_data.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/north_ivt_concat")

    summed_indices = combined_data.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_data.sizes['time'] # get the number of timesteps
    freq = summed_indices / num_time_steps # compute normalized frequency (annual percentage)
    freq.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/north_ivt_composite")

    print("North/south ivt is complete!")


    # Composite for 500mb gph
    datasets = []
    for file in gph_files:
        ds = xr.open_dataset(file)
        # Extract the variable of interest and assign it to a new DataArray
        variable = ds['z'][18]
        gph = variable.sel(longitude = slice(250, 300), latitude = slice(50, 20))/9.81 # divide by g to get gph
        datasets.append(gph)

        ds.close() # close each file after opening
    
    # Concatenate the DataArrays along the desired dimension (e.g., 'time')
    combined_data = xr.concat(datasets, dim = 'time')

    # Optionally, save the combined data to a new NetCDF file
    combined_data.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/gph_concat")

    summed_indices = combined_data.sum(dim = 'time') # sum up all of the indices
    num_time_steps = combined_data.sizes['time'] # get the number of timesteps
    freq = summed_indices / num_time_steps # compute normalized frequency (annual percentage)
    freq.to_netcdf(r"X:/dewhirst/files_for_ivt_block_comp/gph_composite")

    print("GPHs are complete!")


def count_files(directory_path):
    try:
        entries = os.listdir(directory_path)
        files = [f for f in entries]
        return len(files)
    except Exception as e:
        print(f"Error counting files in directory {directory_path}: {e}")
        return 0


def freq_ivt_blocking():
    '''Calculate the percentage of ivt > 500 in the GLR that coincides with blocking in the S/SE US
    and the percentage of blocking in the S/SE US that coincides with IVT > 500'''

    blocking_with_ivt_JJA = r"X:/dewhirst/blocking_with_ivt_JJA"
    ivt_with_blocking_JJA = r"X:/dewhirst/ivt_with_blocking_JJA"

    all_blocking = r"X:/dewhirst/blocking_indices"
    all_ivt = r"X:/dewhirst/1979-2023 U, V, IVT Components ERA5"

    # Initialize variables to keep track of the blocking days in the S/SE and with ivt > 500 in the GLR
    blocking_count = 0
    ivt_count = 0

    # Loop through each year
    for year in range(1979, 2024):
        # Loop through June, July, and August
        for month in [6, 7, 8]:
            # Loop through each day of the current month
            for day in range(1, 32):
                try:
                    date = datetime.date(year, month, day)
                    # Format the date as yyyy_mm_dd
                    date = date.strftime('%Y_%m_%d')

                    # Print or process the formatted date
                    print(date)
                except ValueError:
                    # If day does not exist, skip to next iteration
                    continue

                # Open all of the necessary data: both ivt components and blocking data 0-2 days day before the current date
                # Get 18z data (needed for ivt) and slice into the subdomains for searching
                east_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_eastward_water_vapour_flux_{date}_.nc"
                east_ivt_data = xr.open_dataset(east_ivt_file)
                east_ivt = east_ivt_data['p71.162'][18].sel(latitude = slice(49, 41), longitude = slice(268, 284))

                north_ivt_file = rf"X:/dewhirst/1979-2023 U, V, IVT Components ERA5/era5_vertical_integral_of_northward_water_vapour_flux_{date}_.nc"
                north_ivt_data = xr.open_dataset(north_ivt_file)
                north_ivt = north_ivt_data['p72.162'][18].sel(latitude = slice(49, 41), longitude = slice(268, 284))

                blocking_file = rf"X:/dewhirst/blocking_indices/blocking indices/blocking_indices_{date}.nc"
                blocking_data = xr.open_dataset(blocking_file)
                blocking = blocking_data['B_filt'].sel(latitude = slice(38, 30), longitude = slice(268, 284))
                
                # Calculate magnitude of ivt
                ivt_mag = np.sqrt(east_ivt**2 + north_ivt**2)

                # If IVT > 500 is present in the GLR and blocking is present 0-2 days (either of them) before
                if np.any(ivt_mag > 500):
                    ivt_count += 1
                
                if np.any(blocking == 1):
                    blocking_count += 1
                    
                # Close all of the files at the end of an iteration
                east_ivt_data.close()
                north_ivt_data.close()
                blocking_data.close()
    
    num_ivt_files = count_files(ivt_with_blocking_JJA)
    num_blocking_files = count_files(blocking_with_ivt_JJA)

    block_freq_with_ivt = blocking_count/num_blocking_files*100
    ivt_freq_with_block = ivt_count/num_ivt_files/2*100 # divide by 2 because of the two ivt components

    print(f"The percentage of all IVT > 500 in the GLR assicated with blocking in the S/SE US: {ivt_freq_with_block}")
    print(f"The percentage of all blocking in the S/SE near the time of IVT > 500 in the GLR: {block_freq_with_ivt}")
    


    