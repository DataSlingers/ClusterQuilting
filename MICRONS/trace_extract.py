#!/usr/bin/env python3

import collections
import datajoint as dj
import numpy as np
import pandas as pd
from pathlib import Path 
import sys # if necessary for getting the right path

# Set credentials (these are defaults from a fresh download)
dj.config["database.user"] = "root"
dj.config["database.password"] = "microns123"

# Allows for accessing SQL
from microns_phase3 import nda, utils

# Creates matrix of activity traces that are aligned by timestamps.
## scan_key: dict containing desired session and scan id.
## save_name: save name of trace data frame
def fetch_function(scan_key, save_name):
    
    # Get experiment time stamps
    stim_times = nda.Trial & scan_key
    start_ind = stim_times.fetch('start_idx')
    end_ind = stim_times.fetch('end_idx')
    
    # Specific spontaneous time period at end of recording
    spont_start_time = start_ind[len(start_ind) - 1] + 1
    
    # Construct query to get table of neurons within scans
    query = nda.ScanUnit() & scan_key
    
    # Fetch table of neurons within scans
    neuron_table = query.fetch(as_dict=True)
    
    for ii in range(0, len(neuron_table)):
        unit_keys = {'session': neuron_table[ii]['session'], 
                    'scan_idx': neuron_table[ii]['scan_idx'], 
                    'unit_id': neuron_table[ii]['unit_id']} 
        trace_curr = (nda.Fluorescence & (nda.ScanUnit & unit_keys)).fetch1('trace')
        column_name = "u" + str(neuron_table[ii]['unit_id'])
        
        if ii == 0:
            df = pd.DataFrame({column_name: trace_curr})
        else:
            df_col = pd.DataFrame({column_name: trace_curr})
            # In case lengths differ, truncate to shortest time length
            len_new = min(len(df), len(df_col))
            df = pd.concat([df.reset_index(drop=True), 
                           df_col.reset_index(drop=True)],
                           axis=1)

    df.to_csv(save_name)


# Desired save location, change to desired path.
save_dir = Path.home()

# Scans included in data
scan_key = [{'session': 4, 'scan_idx': 7}] # specific to paper analysis
save_name = save_dir / "session4_scan7.csv"
fetch_function(scan_key, save_name)

scan_key = [{'session': 8, 'scan_idx': 5}] # specific to paper analysis
save_name = save_dir / "session8_scan5.csv"
fetch_function(scan_key, save_name)
