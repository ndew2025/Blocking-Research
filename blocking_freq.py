import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pandas as pd
import cmasher
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import os
import calendar

def blocking_freq(period = 'annual')
'''Calculate a normalized blocking frequency for the eastern United States from 1979-2023 based on the filtered blocking indices 
already calculated

PARAMETERS

period: string
    Specifies the type of frequency (annual or seasonal) to calculate. Choices are 'annual', 'DJF', and 'JJA'. Annual will use all
    files and find the frequency for an entire year. DJF is winter, JJA is summer, in general.
'''

data_dir = # Turbo drive directory here

if period == 'annual' # if you want an annual frequency

    

