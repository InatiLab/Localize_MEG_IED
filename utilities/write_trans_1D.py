#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 2 14:02:40 2024

@author: Leela Srinivasan
"""

import os
import sys
import mne
import numpy as np

# Set variables
subj = sys.argv[1] # Put subject code here

# Set paths
dspm_run_dir =  sys.argv[2] # Put path to dspm meg run here
trans_file_dir = sys.argv[3] # Put path to .fif trans file created earlier here
            
if os.path.isfile(os.path.join(dspm_run_dir, 'trans_mat.1D')):
    print("Transformation 1D file found in dSPM run directory. Proceeding with clustering...")
else:
    written=False

    trans_arr = mne.read_trans(trans_file_dir)['trans']
                
    if np.array_equal(trans_arr[3,:], np.array([0,0,0,1])):
        trans_arr = np.delete(trans_arr, 3, 0)
        np.savetxt(os.path.join(dspm_run_dir, 'trans_mat.1D'), trans_arr)
        written=True
                    
    else:
        raise Exception('Manually verify fourth row of transformation matrix. Exiting...')
        
if written==False:
    raise Exception('Transformation matrix not found/written. Exiting...')
