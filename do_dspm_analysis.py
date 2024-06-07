import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import mne
mne.set_log_level(False)
import mne.channels
from mne.preprocessing import ICA
import sys

# User inputs
subj = '' # Put subject code here
meg_run = '' # Put path to desired MEG run here
spike_mark = '' # Put the mark designated by a clinician as a "spike peak" here
baseline_mark = '' # Put the mark designated by a clinician as a "quiet period" here

# User directories
trans_file = '' # Put path to MNE transformation file here 
out_dir = '' # Put path to desired dSPM output directory here
fs_dir = '' # Put path to folder with Freesurfer reconstructions for all subects

# Set output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Load specified run into mne
raw = mne.io.read_raw_ctf(meg_run, system_clock = 'ignore', preload = False) # reading in CTF file that has been marked for defined subject

if raw:
    print(f'{subj} exists with run {meg_run}') # change to get run number with split 

# Begin preprocess - pick meg channels, downsamples, filters, and plots for bad channels:
raw.load_data() 
# Picks only meg channels
raw.pick_types(meg = True, 
            eeg = False, 
            ref_meg = False
            ) 

# Downsample to 600 Hz, if needed
if not raw.info['sfreq'] == 600:
    raw.resample(600) 
Fs = raw.info['sfreq']

# Apply bandpass filter from 5-50 Hz
raw.filter(5,50) 

## In the limited cases that we required multiple runs to obtain at least 10 spikes, we ran the following line of code on the the second run to obtain the same head position as in the first run
## Spiking and baseline epochs from each run were then concatenated and used to compute the noise covariance matrix and averaged to create the source timecourse
# raw = mne.preprocessing.maxwell_filter(raw, origin='auto',coord_frame='head', destination=dev_head_t_ref)

# Read in marks
events, event_id = mne.events_from_annotations(raw) # adding in events from marker files, if no events exist will get error
raw.plot(clipping = None, n_channels = 50) # plot to remove channels - manually select them
out_text=input("Check raw for bad channels and press enter when finished.")
if out_text != '':
    sys.exit(1)

# Store bad meg channels that were manually chosen above into bad_meg structure for future reference
if len(raw.info['bads']) >= 1:
    bad_meg = raw.info['bads'] # store bad channels so they can be applied to patient's recording 
    raw.pick_types(meg = True, exclude = bad_meg) # excludes bad meg channels 

# Run ICA to remove artifact
ica = ICA(method = 'fastica',
    random_state = 97,
    n_components = 30,
    verbose=True
    )
ica.fit(raw,
    verbose = True,
    reject_by_annotation = True)

# Plot ICA components and manually select components containing phyiologic or other artifact 
ica.plot_sources(raw,title='ICA') # plots ICA component activitations as time series 
out_text=input('Select the components that you wish to exclude from analysis and type enter when finished.')
if out_text != '':
    sys.exit(1)
ica.apply(raw)

# Read in previously created BEM model; looks for folder named 'bem' within this folder that contains all bem model data
# Bem model is expected to be previously created by mne.make_bem_model
bem = mne.make_bem_model(
                        subject = subj, 
                        subjects_dir = fs_dir, 
                        conductivity = [0.3]
                        ) 
bem_sol = mne.make_bem_solution(bem) # Creates a solution that will return the input model

# Set up source space with all candidate dipole locations, created by dividing the Freesurfer model as an icosahedron 4 times (ico4)
src = mne.setup_source_space(
                            subject = subj, 
                            surface = 'white', 
                            spacing = 'ico4', # USER INPUT - change based on desired resolution of source space
                            subjects_dir = fs_dir, 
                            add_dist = True
                            ) 

# Using the file containing the transformation matrix, create the forward solution
n_jobs = 4
forward = mne.make_forward_solution(
                                    info = raw.info, 
                                    trans = trans_file, 
                                    src = src, 
                                    bem = bem_sol, 
                                    meg = True, 
                                    eeg = False, 
                                    # ignore_ref=True
                                    )
print(forward)
src = forward['src']

# Manually check alignment of model before computing outputs
mne.viz.plot_alignment(raw.info, trans=trans_file, subject=subj, subjects_dir=fs_dir, surfaces='outer_skin', show_axes=True, dig=True, eeg=[], meg='sensors',coord_frame='meg')
out_text = input('If bem looks good, hit return.')
if out_text != '':
    sys.exit(1)

# Create spike epochs from 1.5 seconds before your mark to 0.5 seconds after the marked peak
epochs = mne.Epochs(raw, 
                events, 
                event_id = event_id[spike_mark], 
                tmin = -1.5, 
                tmax = 0.5, 
                baseline = None,
                proj = False,
                reject_by_annotation = None,
                preload = True
                ) 

# Inspect all spiking epochs and click to remove if they contain significant artifact
epochs.plot()
out_text = input('Remove any bad epochs and hit return when finished.')
if out_text != '':
    sys.exit()

# Create baseline epochs from 1.5 seconds before your mark to 0.5 seconds after the marked peak
baseline_epochs = mne.Epochs(raw, 
                    events, 
                    event_id = event_id[baseline_mark],
                    tmin = -1.5, 
                    tmax = 0.5, 
                    baseline = None
                    ) 

# Compute the noise covariance matrix from baseline epochs
noise_cov=mne.compute_covariance(baseline_epochs)

# Create the dSPM inverse using the forward solution and noise covariance matrix selected above
method = "dSPM"
snr = 3.0 # Edit to fine-tune 
lambda2 = 1.0 / snr ** 2 
inv = mne.minimum_norm.make_inverse_operator(info=raw.info, 
                                                forward=forward, 
                                                noise_cov=noise_cov, 
                                                loose='auto', 
                                                depth=0.8, 
                                                fixed='auto',
                                                rank=None, 
                                                use_cps=True
                                                )

# Create an averaged epoch from all spiking epochs above
avg_epoch=epochs.average()

# Apply the dSPM inverse solution to the averaged epoch to obtain the source timecourse
stc=mne.minimum_norm.apply_inverse(avg_epoch,inv,lambda2,method,label=None)

# Create moving average of the averaged epoch
window=30
smoothed_stcs = np.zeros((1,1201)) 
data_arr=stc.data
row,column=data_arr.shape
for tp in np.arange(0,row):
    tp_data=pd.DataFrame(stc.data[tp])
    moving_avg_source=[]
    moving_avg_source = tp_data.rolling(window = window, min_periods = 1, center = True).mean()
    moving_avg_source = moving_avg_source.to_numpy().T
    if tp == 0:
        smoothed_stcs=moving_avg_source
    else:
        smoothed_stcs = np.vstack([smoothed_stcs, moving_avg_source])

# Save moving average to source timecourse
stc.data=smoothed_stcs

# Save source timecourse array for analysis
np.save(os.path.join(out_dir,'stc_array'),smoothed_stcs)
stc.save('full-stcs', overwrite=True)
src = inv['src']

# Save model for conversion to AFNI surface
mne.write_source_spaces(fname = 'mymodel.fif', 
                        src = src, 
                        overwrite = True
                    )
