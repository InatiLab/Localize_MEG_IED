import glob
import nibabel as nb
import numpy as np
import os
import pandas as pd
import subprocess
from read_dip_results import *

subj='' # Put your subject code here
run='' # Put the path to the subject's MEG run here

# Analysis window and accepted dipole fit error  
latency=25
error=70

# Set main directory here
ctf_run_dir = '' # Put path to CTF .ds file here
ctf_dir = ctf_run_dir.split('/'+(ctf_run_dir.split('/')[-1]))[0] # Get directory with CTF file in it
project_dir = '' # Folder with patient's MRI, and to write out project files

#====================================================================================================================
# DEFINE FUNCTIONS

def extract_dipole_per_metric(dipole_dframe,metric=None):
    if metric == 'first':
        first_idx=dipole_dframe['fit_seq_idx'].idxmin()
        dipole_dframe=dipole_dframe[first_idx:(first_idx+1)]
    elif metric == 'moving':
        dipole_dframe=dipole_dframe
    return dipole_dframe

def get_temp_dataframe(index):
    tmp_dip_disp=dip_disp.copy(deep=True)
    tmp_dip_disp
    for ind in tmp_dip_disp.index:
        seq_num=tmp_dip_disp['fit_seq_idx'][ind]
        if seq_num != index:
            tmp_dip_disp=tmp_dip_disp.drop(ind,axis=0)
    return tmp_dip_disp

def write_dipole_overlays(dipole_mat=None, mri_dir=None, open_afni=False, type=None, index=None, mark=None):
    # Set Filenames
    os.chdir(mri_dir)
    mri_input_filename=os.path.join(mri_dir,'ortho+orig.HEAD')
    dipole_mri_output_filename='DIPOLE_NII.nii'
    # Load MRI header
    mri=nb.load(mri_input_filename)
    # Dipole NIFTI image
    dipole_nii_output=nb.Nifti1Image(dipole_mat, mri.affine, header=mri.header)
    dipole_nii_output.to_filename(dipole_mri_output_filename)
    if type == 'moving':
        base='_dipole+orig'
        name=str(index)+base
        # Generate an afni readable file
        subprocess.run(('3dcopy {} {}'.format('',name)).split(' '))
        subprocess.run('3dcopy {} {}'.format(dipole_mri_output_filename, name).split(' '))
    else:
        # Generate an afni readable file
        name=type+'_'+mark+'_dipole+orig'
        subprocess.run('3dcopy {} dipole+orig'.split(' '))
        subprocess.run('3dcopy {} {}'.format(dipole_mri_output_filename, os.path.join(mri_dir,name)).split(' ')) #Change the data to afni format
        os.remove(dipole_mri_output_filename)

def transform_dipole_LPI_to_voxel(mri_dir,final_dipole_dframe):
    '''The MRI dataset is in LPI space Use the inverse affine matrix to map the LPI dipoles to the voxel coordinates'''
    ortho_mri=os.path.join(mri_dir,'ortho+orig.HEAD')
    mri=nb.load(ortho_mri)
    mri_matrix=mri.get_fdata()
    affine_inv=np.linalg.inv(mri.affine)
    # Create empty dipole array
    dipole_matrix=np.zeros(mri_matrix.shape[0:3]) #, dtype=np.int8
    # Loop over dipoles and write to the matrix.  Incrementing in overlapping areas
    for index,row in final_dipole_dframe.iterrows():
        x,y,z=row[['xLPI','yLPI','zLPI']]
        x,y,z=np.round(x),np.round(y),np.round(z)
        i,j,k=np.sum(affine_inv*(x,y,z,1),axis=1)[0:3]
        dipole_matrix[vox_fill(i,j,k, fill_vox=4)]+=1
    # May or may not be necessary to expand dimensions to 4
    dipole_matrix=np.expand_dims(dipole_matrix,axis=3)
    return dipole_matrix

def vox_fill(i,j,k, fill_vox=None): 
    '''Create padding around the dipoles for display'''
    i,j,k=np.round([i,j,k]).astype(int)
    if fill_vox==None:
        return i,j,k
    else:
        fill_vox=int(fill_vox)
        return slice(i-fill_vox, i+fill_vox,1),slice(j-fill_vox, j+fill_vox,1), slice(k-fill_vox, k+fill_vox,1) 

def initialize_and_fill_new_dframe(error,metric):
    tmp_dip_display=pd.DataFrame(columns=dipole_dframe.columns)
    tmp_dip_display=extract_dipole_per_metric(dipole_dframe,metric=metric)
    tmp_dip_display=tmp_dip_display[tmp_dip_display['Err(%)']<(100-int(error))]
    return tmp_dip_display

def get_dipole_idx(label):
    '''Return the dipole index from file fits
    This is important to track if multiple markers are present per run'''
    if '_' in label:
        return str(label.split('_')[1].split('.')[0])
    else:
        return str(label.split('.')[0])

#===================================================================================================
# START SCRIPT

# Rename t1
cmd="3dcopy {}/t1.nii {}/ortho+orig"
cmd=cmd.format(project_dir,project_dir)
subprocess.run(cmd,shell=True)

# Make averaged run
ctf_avg_dir = (ctf_run_dir.strip('.ds'))+'-avg.ds'
os.chdir(ctf_dir)

# Write average spike from all epochs
cmd="averageDs -marker S -includeBad -overlap 2.0 -time -1.5 0.5 {} {}"
cmd=cmd.format(ctf_run_dir,ctf_avg_dir)
subprocess.run(cmd,shell=True)

# Add a marker in the averaged epoch at the peak of the spike (0)
cmd="addMarker -n avgspike -t 0.0 {}"
cmd=cmd.format(ctf_avg_dir)
subprocess.run(cmd,shell=True)

# Check if file needs to be resampled to 600 Hz
os.chdir(ctf_avg_dir)
resampled_file=ctf_avg_dir.strip('-avg.ds')+'_resampled-avg.ds'
for file in os.listdir(os.getcwd()):
    hist_file=glob.glob('*hist')[0]
    d=open(hist_file)
    lines=d.readlines()
    for line in lines:
        line=line.strip()
        if 'Sample rate:' in line:
            sample_rate=line
    sample_rate=int((sample_rate.strip()).split(' ')[2])
    factor=int(sample_rate/600)
    # If the file needs to be resampled to 600, create a new resampled dataset
    if factor != 1 and not os.path.exists(resampled_file):
        cmd="newDs -resample {} {} {}"
        cmd=cmd.format(factor,ctf_avg_dir,resampled_file)
        subprocess.run(cmd,shell=True)

# If we created a resampled file, reset the averaged directory to the resampled file
if os.path.exists(resampled_file):
    ctf_avg_dir=(resampled_file)

## In the very limited cases that multiple runs were needed, each run was resampled and averaged, as above
## Then, the following lines of code were used to average the runs together, and the dipole was fit on the grand averaged run.
# cmd="grandAverageDs{} grandAv-c.ds"
# cmd=cmd.format(averaged_dsets)
# subprocess.run(cmd,shell=True)

# Write dipole over the course of the selected latency of the avgspike mark
dipole='avgspike'
os.chdir(ctf_avg_dir)
cmd="dfit -a -z -b -{} -e 0.001 -i 4 -h {}/default.hdm -f {}/default.dip -p /home/jstout/processing.cfg -m avgspike {} avgspike.dip"
cmd=cmd.format((latency/1000),ctf_avg_dir,ctf_avg_dir,ctf_avg_dir)
subprocess.run(cmd,shell=True)

# Assemble a dipole dataframe of position and error over the moving dipole fit
os.chdir(ctf_avg_dir)
dipole_dframe=read_dip_results.assemble_dipole_dframe('avgspike.dip')
dipole_dframe['Run']=ctf_avg_dir.split('/')[7]
dipole_dframe['Dipole_idx']=dipole_dframe.Label.apply(get_dipole_idx)

# Find dipole localizations with an error above 70%; get rid of sample after peak
filtered=(dipole_dframe[dipole_dframe['Err(%)']<(100-error)][0:-1]).reset_index(inplace=False)
if not filtered.empty == True:
    first_lat=float(filtered['Sample'][0])
    peak_lat=float(((filtered)[-1:])['Sample'])
else:
    first_lat=np.nan
    peak_lat=np.nan

# Write out a text file containing the selected dipole timing
out_df=pd.DataFrame({"dipole":'avg',"first":first_lat,"peak":peak_lat},index=[0])
if not os.path.exists(project_dir):
    os.mkdir(project_dir)
out_df.to_csv(os.path.join(project_dir,'dipole_timing.txt'),header=True)

# Write out dipoles for analysis
metrics="first","peak"
for metric in metrics:
    time=out_df[metric][0]
    for a in np.arange(0,len(dipole_dframe)):
        test=float(dipole_dframe['Sample'][a])
        if test == time:
            dip_disp=dipole_dframe[a:(a+1)]
            dipole_mat=transform_dipole_LPI_to_voxel(project_dir,final_dipole_dframe=dip_disp)
            write_dipole_overlays(dipole_mat=dipole_mat,mri_dir=project_dir,open_afni=False,index=[a],mark='avgspike',type=metric)
