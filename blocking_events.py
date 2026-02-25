import xarray as xr
import netCDF4 as nc
import numpy as np
import os

# function to write the indices to a text file
def write_indices_to_file(indices, index_file):
    with open(index_file, 'w') as outfile:
        for index in indices:
            outfile.write(f"{index}\n")

# something else i need to add to this function is that an event also ends when the blocking maximum moves outside of the original box 
# of the onset that is 50% greater in both latitude and longitude than the daily box check

def blocking_events(consec_day_threshold):
    '''Using the previously calculated blocking indices, identify all of the events 
    and combine the days together for every event
    
    PARAMETERS
    ==========
    
    consec_day_threshold: integer, the number of consecutive days (time steps at 18z) 
    required for a series of blocked days to be considered an event. Probably start with 5 days.'''

    # blocking and gph anomaly directories. Change paths to Bucky
    blocking_dir = '/data/users/nickdew/test_blocking_indices'
    gph_anom_dir = '/data/users/nickdew/test_gph_anoms'

    # Get a list of all of the blocking files and sort them to ensure chronological order
    blocking_files = sorted(os.listdir(blocking_dir))

    # Get a list of all of the gph anomaly files
    gph_anom_files = sorted(os.listdir(gph_anom_dir))

    # text file to store the event day indices in
    index_file = '/data/users/nickdew/event_day_indices.txt'

    # Index to access the correct day. Begins at 0, will increase by one for each iteration in the loop
    current_index = 0

    # Index for consecutively blocked days to compare to consec_day_threshold.
    num_consec_days = 0

    # Index to compare to num_consec_days to see when an event has ended
    prev_num_consec_days = 0

    # boolean to serve as a check if the number of consecutive days is satisfied. Default is False
    is_greater_than_num_consec_days = False

    # Total number of days/iterations in the dataset
    num_days = len(blocking_files)

    # Initialize variables to store the indices of the coordinates of the maximum blocking index from the previous day
    prev_B_max_lat_idx = None
    prev_B_max_lon_idx = None

    # Size, in grid cells, for the latitude and longitude limits of the box to search inside
    # Multiply by 4 because the resolution of ERA5 data is 0.25 deg
    lat_box_cells = 27 * 4
    lon_box_cells = 36 * 4

    # Setup a while loop
    while current_index < num_days:
        # Get the current blocking file, open it, and access the variables. Slice accordingly to avoid issues with indices and nans
        current_blocking_file = os.path.join(blocking_dir, blocking_files[current_index])
        blocking_data = xr.open_dataset(current_blocking_file)
        B_raw = blocking_data['B_raw'].sel(latitude = slice(50,20), longitude = slice(250,310))
        B_filt = blocking_data['B_filt'].sel(latitude = slice(50,20), longitude = slice(250,310))

        # Get current GPH anomaly file and get the variable
        current_gph_anom_file = os.path.join(gph_anom_dir, gph_anom_files[current_index])
        gph_anom_data = xr.open_dataset(current_gph_anom_file)
        gph_anom = gph_anom_data['z_anom'].sel(latitude = slice(50,20), longitude = slice(250,310))

        # If all of the binary/filtered indices are 0 (not satisfying the longitude threshold), skip that day
        if np.all(B_filt == 0):
            # Reset the consecutive days count and the box limits
            num_consec_days = 0
            prev_B_max_lat_idx = None
            prev_B_max_lon_idx = None

            # If the previous number of consecutive days is greater than the current (since it has been reset),
            # this indicates that an event has ended and it needs to be stored.
            # Capture the correct indices of the file in the sorted directory and write them to the index file
            if (prev_num_consec_days > num_consec_days) and (is_greater_than_num_consec_days == True):
                indices = np.arange(current_index - prev_num_consec_days, current_index)
                write_indices_to_file(indices, index_file)
            
            # Reset after the above comparison
            prev_num_consec_days = 0
            is_greater_than_num_consec_days = False

        else:  # Otherwise, if there is a region of large blocked flow...
            # Find the area where the blocking indices are wide enough in longitude and where the gph anom is positive
            condition = (B_filt == 1) & (gph_anom > 0)

            # Get the area in B_raw where the above is true
            constrained_B = B_raw.where(condition)

            # Get the maximum inside that area
            B_max = np.max(constrained_B)

            # Get the coordinates of the maximum blocking index within the condition specified above
            B_max_lat_idx = np.where(B_raw == B_max)[0][0]
            B_max_lon_idx = np.where(B_raw == B_max)[1][0]

            # If the coordinates of B_max correspond with a neutral or negative height anomaly, move on
            if constrained_B.size == 0:
                # Reset the consecutive days count and the box limits
                num_consec_days = 0
                prev_B_max_lat_idx = None
                prev_B_max_lon_idx = None

                # If the previous number of consecutive days is greater than the current (since it has been reset),
                # This indicates that an event has ended and it needs to be stored.
                # Capture the correct indices of the file in the sorted directory and write them to the txt file
                if (prev_num_consec_days > num_consec_days) and (is_greater_than_num_consec_days == True):
                    indices = np.arange(current_index - prev_num_consec_days, current_index)
                    write_indices_to_file(indices, index_file)
                
                # Reset after the above comparison
                prev_num_consec_days = 0
                is_greater_than_num_consec_days = False

            else: # If constrained_B does have an actual size
                # Check if the previous blocking index maximum doesn't exist or if the maximum blocking index
                # is within the bounds of the box defined by the latitude and longitude limits
                if (prev_B_max_lat_idx is None) or (
                    (B_max_lat_idx >= prev_B_max_lat_idx - lat_box_cells / 2) and
                    (B_max_lat_idx <= prev_B_max_lat_idx + lat_box_cells / 2) and
                    (B_max_lon_idx >= prev_B_max_lon_idx - lon_box_cells / 2) and
                    (B_max_lon_idx <= prev_B_max_lon_idx + lon_box_cells / 2)):

                    # Update box position to the current day to be used for the next
                    prev_B_max_lat_idx = B_max_lat_idx
                    prev_B_max_lon_idx = B_max_lon_idx

                    # Update the previous number of consecutive before updating the current below
                    prev_num_consec_days = num_consec_days

                    # If the above criteria are satisfied, add 1 to num_consec_days
                    num_consec_days += 1

                    # Check if current block surpasses the threshold
                    if num_consec_days >= consec_day_threshold:
                        # Update to true if the above is satisfied
                        is_greater_than_num_consec_days = True

                else: # if the next day is not contained in the previous box, the event also ends here
                    # Reset the consecutive days count and the box limits
                    num_consec_days = 0
                    prev_B_max_lat_idx = None
                    prev_B_max_lon_idx = None

                    # If the previous number of consecutive days is greater than the current (since it has been reset),
                    # This indicates that an event has ended and it needs to be stored.
                    # Capture the correct indices of the file in the sorted directory and write them to the txt file
                    if (prev_num_consec_days > num_consec_days) and (is_greater_than_num_consec_days == True):
                        indices = np.arange(current_index - prev_num_consec_days, current_index)
                        write_indices_to_file(indices, index_file)
                    
                    # Reset after the above comparison
                    prev_num_consec_days = 0
                    is_greater_than_num_consec_days = False

        # print the index to see how the code is progressing
        print(current_index)

        # For every iteration, add 1 to the current_index so that it goes through every day
        current_index += 1


if __name__ == "__main__":
    blocking_events(5)