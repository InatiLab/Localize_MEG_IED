import os
import subprocess
import pandas as pd
import shutil
import numpy as np
from sh import gunzip

# Set variables 
subj = '' # Set subject code here 

# Set paths
fs_dir = '' # Set path to subject's Freesurfer reconstruction here
trans_file = '' # Set path to transformation matrix, typically in the FS folder
dspm_dir = '' # Set path to subject's stc array 
out_dir = '' # Set path to desired output directory; ideally, this is the directory that you wrote the dipole timing file into
utilities_dir = os.path.join(os.getcwd(),'utilities')

# Gather all necessary surface reconstruction files from the patient's SUMA directory in a standard mesh
surfvol_file = '' # Set path to subject's SurfVol file, created by Freesurfer
rh_pial = '' # Set path to subject's right hemisphere pial surface reconstruction
rh_smoothwm = '' # Set path to subject's right hemisphere white matter surface reconstruction
lh_pial = '' # Set path to subject's left hemisphere pial surface reconstruction
lh_smoothwm = '' # Set path to subject's left hemisphere white matter surface reconstruction
rh_specfile = '' # Set path to subject's right hemisphere .spec file
lh_specfile = '' # Set path to subject's left hemisphere .spec file

# Read in selected timepoints, created by do_ECD_analysis.py
times_df=pd.read_csv(os.path.join(out_dir,'dipole_timing.txt'))

# Read in stc file, created by do_dspm_analysis.py
data_array=np.load(os.path.join(dspm_dir,'stc_array.npy'))

#====================================================================================================================
#====================================================================================================================
# START SCRIPT

# Send dSPM outputs to surface

os.chdir(utilities_dir)
subprocess.run(("./get_new_stc2gii.sh {} {} {} {} {}".format(subj,fs_dir,dspm_dir,trans_file,utilities_dir)),shell=True)

#====================================================================================================================
# Create clusters

os.chdir(dspm_dir)
for a in 'first','peak':
    for hemi in "rh","lh":

        # Read in suma2mne file to relate VS and SUMA nodes
        infile=subj+'_'+hemi+'_suma2mne.csv'
        master_df=pd.read_csv(os.path.join(dspm_dir,infile),delimiter=' ')

        # Add an empty column to the suma2mne dataframe to contain the dSPM values
        zeros_df=pd.DataFrame(np.zeros(len(master_df)))
        print('Now clustering at '+a)
        master_df=pd.concat([master_df,zeros_df],axis=1)
        # Rename the column to the timepointg and hemisphere
        in_col = a+'_'+hemi
        master_df=master_df.rename(columns={0:a+'_'+hemi})

        # Get the timepoint from the dipole writing; subtract 1, for the difference in sampling points between CTF and MNE
        select_tp=int(times_df[a][0])-1
        
        # Make a new array that contains the dSPM value for every virtual sensor at that timepoint
        vs_df=pd.DataFrame((np.arange(len(data_array[:,select_tp])))) # Get index for every virtual sensor
        tmp_df=pd.DataFrame(data_array[:,select_tp]) # Get the dSPM values for every virtual sensor at that timepoint
        tp_df=(pd.concat([vs_df,tmp_df],axis=1,ignore_index=True)).rename(columns={0:'vs',1:'val'})

        # For each virtual sensor, find the instances of it being related to a SUMA vertex and insert its dSPM value there
        for vs in tp_df['vs']:
            indices=master_df[master_df['virtual_sensor_node'] == vs].index.tolist()
            for index in indices:
                master_df[in_col][index] = tp_df['val'][vs]
        
        # Make a new empty dataframe with each of the dSPM values for the SUMA vertices
        tmp_df=(pd.DataFrame(np.arange(0,36002))).rename(columns={0:'SUMA_vertex'})
        new_df=(pd.concat([tmp_df,pd.DataFrame(np.zeros(36002))],axis=1)).rename(columns={0:'val'})
        for x in tmp_df['SUMA_vertex']:
            indices=master_df[master_df['SUMA_vertex'] == x].index.to_list() 
            if len(indices) != 0:
                new_df['val'][x] = master_df[in_col][indices[0]]

        # Write out the dataframe with the dSPM value for every SUMA vertex
        new_df['val'].to_csv(os.path.join(fs_dir,'SUMA','tp.1D.dset'),sep=' ',index=False,header=False)
        os.chdir(os.path.join(fs_dir,'SUMA'))

        # Smooth dSPM values on the white matter surface
        smooth_cmd="SurfSmooth -spec {} -surf {} -input tp.1D.dset -met HEAT_07 -target_fwhm 20 -Niter -1 -output {}/tp_smoothed_{}.niml.dset"
        if hemi == 'rh':
            smooth_cmd=smooth_cmd.format(rh_specfile,rh_smoothwm,out_dir,hemi)
        else:
            smooth_cmd=smooth_cmd.format(lh_specfile,lh_smoothwm,out_dir,hemi)
        subprocess.run(smooth_cmd,shell=True)

        # Convert the values to a 1D for each hemisphere
        subprocess.run(("ConvertDset -o_1D -input {}/tp_smoothed_{}.niml.dset -prefix {}_tp_smoothed".format(out_dir,hemi,hemi)),shell=True)
        
        # Remove the index column so AFNI can read the file
        nodes=pd.read_csv(os.path.join(fs_dir,'SUMA',hemi+'_tp_smoothed.1D.dset'),skiprows=5,header=None,sep=' ',names=['val'],nrows=36002)
        nodes.reset_index(inplace=True)
        nodes=nodes.drop(columns=['index'])
        nodes.to_csv('temp_write.1D',sep=' ',header=False)
        
        # Send the surface dSPM values to the volume on the pial surface
        cmd="3dSurf2Vol -surf_A {} -surf_B {} -sv {} -spec {} -sdata_1D temp_write.1D -grid_parent {} -map_func nzave -f_steps 15 -prefix {}"
        if hemi == 'rh':
            cmd=cmd.format(rh_pial,rh_smoothwm,surfvol_file,rh_specfile,surfvol_file,'temp_out.nii')
        elif hemi == 'lh':
            cmd=cmd.format(lh_pial,lh_smoothwm,surfvol_file,lh_specfile,surfvol_file,'temp_out.nii')
        subprocess.run(cmd,shell=True)

        # Move the file to the out directory
        for file in os.listdir(os.getcwd()):
            if 'temp_out.nii' in file:
                shutil.move(file,os.path.join(out_dir,file))
        
        #Expand into the volume to remove any holes from 3dSurf2Vol
        subprocess.run(("@ROI_modal_grow -input {}/temp_out.nii -outdir {}/temp_grow -niters 2 -mask {}/temp_out.nii -prefix tmp".format(out_dir,out_dir,out_dir),shell=True)
        
        # Only take the level 2 expansion; remove all other files
        for file in os.listdir(os.path.join(out_dir,'temp_grow')):
            if 'rm_rgm_02.nii.gz' in file:
                gunzip(os.path.join(out_dir,'temp_grow',file))
                os.rename(os.path.join(out_dir,'temp_grow',file[:-3]),os.path.join(out_dir,a+'_'+hemi+'.nii'))
        shutil.rmtree(os.path.join(out_dir,'temp_grow'))
        os.remove(os.path.join(fs_dir,'SUMA','temp_write.1D'))
        
        # Clean up other files
        for file in os.listdir(out_dir):
            if 'temp_out.nii' in file:
                os.remove(os.path.join(out_dir,file))
    
    # Combine masks on either side
    os.chdir(out_dir)
    cmd="3dcalc -a {} -b {} -expr '(step(a)*a)+(step(b)*b)' -prefix {}"
    lh=a+'_lh.nii'
    rh=a+'_rh.nii'
    out=a+'.nii'
    cmd=cmd.format(lh,rh,out)
    subprocess.run(cmd,shell=True)

    # Read in the smoothed files 
    lh_tp_df=pd.read_csv(os.path.join(fs_dir,'SUMA','lh_tp_smoothed.1D.dset'),header=None,delim_whitespace=True,skiprows=5,nrows=36002)
    rh_tp_df=pd.read_csv(os.path.join(fs_dir,'SUMA','rh_tp_smoothed.1D.dset'),header=None,delim_whitespace=True,skiprows=5,nrows=36002)
    
    # Combine the left and right halves to get the entire brain; take the larger value of the two 
    combined=pd.DataFrame(np.arange(0,36002))
    for x in np.arange(len(combined)):
        lh=lh_tp_df[0][x]
        rh=rh_tp_df[0][x]
        val=max(rh,lh)
        combined[0][x] = val
    percentile=np.nanpercentile([combined],[95])[0]

    # Create clusters from the volumetric dSPM mask, with the chosen percentile as the threshold
    subprocess.run(("3dclust -savemask {} -1clip {} 5 75 {}".format(os.path.join(out_dir,a+'_clust.nii'),percentile,os.path.join(out_dir,a+'.nii'))),shell=True)
