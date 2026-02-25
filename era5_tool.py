import os
from datetime import datetime, timedelta
import gcsfs
import xarray as xr

def download_era5_data(first_year, end_year, first_month, end_month, first_day, end_day, 
                       variables, file_directory, lat_range = None, lon_range = None, pressure_levels = None):
    '''Downloads ERA5 data for any pressure level variables (or single level variables) and saves them into a single NetCDF file, 
    concatenated by time. Before concatentation, one can choose to specify both the time domain (data at a specified intervals) and the 
    spatial domain (a specified latitude and longitude range).
    Parameters: 
    - first_year: the first year you want data from.
    - end_year: the last year you want data from.
    - first_month: the first month you want data from.
    - end_month: the last month you want data from.
    - first_day: the first day you want data from.
    - end_day: the last day you want data from.
    - variables: a list of ERA5 variables you want to download, or a single variable.
    - pressure_levels: a list of pressure levels you want data taken from, or None if single level is desired.
    - file_directory: The directory where the files will be stored.
    - lat_range: tuple or list specifying latitude range (e.g., [min_lat, max_lat]). This must be defined from north to south.
    - lon_range: tuple or list specifying longitude range (e.g., [min_lon, max_lon]). This must have positive values starting from west to east.
    '''

    # Define the date range you want to download
    start_date = datetime(first_year, first_month, first_day)
    end_date = datetime(end_year, end_month, end_day)
    
    # Generate list of dates, by day, using the start and end dates
    date_range = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]

    # If pressure_levels is None, set it to an empty list
    if pressure_levels is None:
        pressure_levels = []

    # Initialize GCS filesystem
    fs = gcsfs.GCSFileSystem()

    # Loop through the date range in chunks
    for date in date_range:
        year = date.year
        month = date.month
        day = date.day

        #Loop through the parameters and pressure levels wanted
        for variable in variables:
            #If pressure levels are specified, the function will get the parameters from each.
            if pressure_levels:
                for level in pressure_levels:
                    #The constant url for raw era5 data. The source path is specified by the time, varible, and pressure level needed.
                    #The target path is wherver you want it stored on your computer.
                    base_url = "gs://gcp-public-data-arco-era5/raw/date-variable-pressure_level"
                    source_path = f'{base_url}/{year}/{month:02d}/{day:02d}/{variable}/{level}.nc'
                    target_path = f'{file_directory}/era5_{variable}_{year}_{month:02d}_{day:02d}_{level}.nc'

                    #Open the dataset to amend it. Here, I slice the dataset by taking out the lat/lon domain I need and by the hours I need.
                    #These are specified by the lat_range, lon_range, and hourly step, getting even time intervals throughout the day. 
                    with fs.open(source_path, 'rb') as fsrc:
                        ds = xr.open_dataset(fsrc)
                        ds = ds.sel(latitude = slice(lat_range[0], lat_range[1]), longitude = slice(lon_range[0], lon_range[1]))
                        ds.to_netcdf(target_path)
                            
            #If a pressure level is not specified, the function will assume you want a "single level" variable (such as IVT)               
            else:
                base_url = "gs://gcp-public-data-arco-era5/raw/date-variable-single_level"
                source_path = f'{base_url}/{year}/{month:02d}/{day:02d}/{variable}/surface.nc'
                target_path = f'{file_directory}/era5_{variable}_{year}_{month:02d}_{day:02d}_.nc'
                with fs.open(source_path, 'rb') as fsrc:
                    ds = xr.open_dataset(fsrc)
                    ds = ds.sel(latitude = slice(lat_range[0], lat_range[1]), longitude = slice(lon_range[0], lon_range[1]))
                    ds.to_netcdf(target_path)
    
    print("The unbearably long download process is complete!")