import numpy as np
import xarray as xr
from collections import defaultdict
import os
import calendar
import subprocess

data_directory = r"X:/dewhirst/1979-2023 500mb Heights ERA5"
output_directory = r"X:/dewhirst/500_heights_monthly_climo"

def preprocess(ds):
    # Perform the desired selection on each individual dataset
    return ds.sel(longitude=slice(250, 310), latitude=slice(55, 15)).isel(time=slice(18, 19))

def save_climo(data_directory, output_directory):
    """
    Process geopotential height data files to compute monthly averages and save to a directory.
    
    Parameters:
    - data_directory: str, directory containing the input data files
    - output_directory: str, directory to save the output files
    
    The function will generate 24 files, two for each month. One contains all of the raw concatenated data, 
    and the second contains the averages.
    """
    # Step 1: Group files by month
    files_by_month = defaultdict(list)
    for filename in os.listdir(data_directory): # Loop through files in the data directory
        if 'geopotential' in filename: # Check if the file contains geopotential data
            parts = filename.split('_')
            if len(parts) > 3:  # Ensure filename contains date in the correct format
                month = int(parts[3])  # Extract the month from the filename, the fourth split part
                files_by_month[month].append(os.path.join(data_directory, filename))
    
    # Step 2: Calculate monthly averages and save to disk
    # Assume files_by_month is a dictionary mapping month numbers (1-12) to lists of file paths
    for month, files in files_by_month.items():
        # Open all files corresponding to this month as a single dataset
        monthly_data = xr.open_mfdataset(files, chunks={'time': 100}, preprocess = preprocess)

        # Calculate the monthly average
        monthly_average = monthly_data.mean(dim='time')

        # Convert month number to its corresponding name
        month_name = calendar.month_name[month].lower()

        # Save the concatenated dataset for the entire month
        concatenated_filename = os.path.join(output_directory, f'concatenated_geopotential_files_for_{month_name}.nc')
        monthly_data.to_netcdf(concatenated_filename)

        # Save the monthly averaged dataset
        output_filename = os.path.join(output_directory, f'geopotential_monthly_average_{month_name}.nc')
        monthly_average.to_netcdf(output_filename)

        # Close the datasets, although this is typically handled by xarray
        monthly_data.close()