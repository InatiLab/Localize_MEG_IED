import os
import subprocess
import shutil
import pandas as pd

# Set variables
subj = '' # Set subject of interest
surfvol_file = '' # Set name of SurfVol file used in Freesurfer reconstruction

# Set directories
fs_dir = '' # Set path to all freesurfer reconstructions; should contain fsaverage
subj_fs_dir = os.path.join(fs_dir,subj) # Set path to subject's freesurfer reconstruction
mri_surf_dir = ''
scripts_dir = os.getcwd()
files_dir = os.path.join(scripts_dir,'files')
t1_dir = '' # Set path to project t1; ideally, already-axialized scan used for freesurfer
project_dir = '' # Set path to project directory

if not os.path.exists(os.path.join(t1_dir,'anat+orig.HEAD')):
    cmd="3dcopy t1.nii {}/t1+orig"
    cmd = cmd.format(project_dir)
    subprocess.run(cmd,shell=True)

#========================================================================================================================================
# Get extra freesurfer outputs and parcellations -- only run after recon-all

# Create SUMA directory
cmd = "@SUMA_Make_Spec_FS -NIFTI -sid {} -fspath {} -make_rank_dsets -extra_fs_dsets -no_ld"
cmd = cmd.format(subj,subj_fs_dir)
subprocess.run(cmd,shell=True)

# Get paths to subsequent files
surfvol_file = os.path.join(subj_fs_dir,'SUMA',subj+'_SurfVol.nii')
std_mesh = '' # Set to the standard mesh used for this project; ex. 20, 60, or 141

# Set spec files
lh_specfile = 'std.'+std_mesh+'.'+subj+'_lh.spec' # Set path to the left hemisphere spec file; make sure this is the same standard mesh as the surfaces
rh_specfile = 'std.'+std_mesh+'.'+subj+'_rh.spec' # Set path to the right hemisphere spec file; make sure this is the same standard mesh as the surfaces

# Set surface names
lh_white = 'std.'+std_mesh+'.lh.smoothwm.gii' # Set path to left hemisphere white matter surface
rh_white = 'std.'+std_mesh+'.rh.smoothwm.gii'
lh_pial = 'std.'+std_mesh+'.lh.pial.gii'
rh_pial = 'std.'+std_mesh+'.rh.pial.gii'

#========================================================================================================================================
# Get parcels

os.chdir(os.path.join(subj_fs_dir,'SUMA'))
for hemi in "lh","rh":
    # Get the label file and convert it to a 1D file
    file=os.path.join(files_dir,'std.'+std_mesh+'.'+hemi+'.Schaefer2018_400Parcels_7Networks_order.smooth3mm.lbl.niml.dset') # We are only providing the file for the 141 standard mesh. Please email if you need other files.
    convert_cmd="ConvertDset -o_1D -input {} -prefix {}/{}_schaefer_parcels.1D"
    convert_cmd=convert_cmd.format(file,os.path.join(subj_fs_dir,'SUMA'),hemi)
    subprocess.run(convert_cmd,shell=True)

    # Edit files to be AFNI-readable by removing excess rows
    filename=hemi+'_schaefer_parcels.1D.dset'
    file=pd.read_csv(os.path.join(subj_fs_dir,'SUMA',filename),delim_whitespace=True,skiprows=5,nrows=198812,header=None)
    file.to_csv(os.path.join(subj_fs_dir,'SUMA','final_'+filename),header=None,index=False,sep=' ')

    # Send the parcels into the volume using @surf_to_vol_spackle
    final_filename=os.path.join(subj_fs_dir,'SUMA','final_'+hemi+'_schaefer_parcels.1D.dset')
    surfvol_cmd="@surf_to_vol_spackle -maskset {} -spec {} -surfA {} -surfB {} -surfset {} -prefix {}"
    if hemi == 'lh':
        surfvol_cmd=surfvol_cmd.format(surfvol_file,lh_specfile,lh_white,lh_pial,final_filename,hemi+'_schaefer_parcels+orig')
    else:
        surfvol_cmd=surfvol_cmd.format(surfvol_file,rh_specfile,rh_white,rh_pial,final_filename,hemi+'_schaefer_parcels+orig')
    subprocess.run(surfvol_cmd,shell=True)

# Combine the two halves; start with right hemisphere at parcel 200 as they are both labeled 1-200 in the base file
combine_cmd="3dcalc -a temp_lh_schaefer_parcels.nii -b temp_rh_schaefer_parcels.nii -expr '(step(a)*a)+(step(b)*(b+200))' -prefix schaefer_parcels+orig"
subprocess.run(combine_cmd,shell=True)

# Clean up folder
for hemi in 'lh','rh':
    # os.remove(os.path.join(subj_fs_dir,'SUMA','final_'+hemi+'_schaefer_parcels.1D.dset'))
    # os.remove(os.path.join(subj_fs_dir,'SUMA',hemi+'_schaefer_parcels.1D.dset'))
    os.remove(os.path.join(subj_fs_dir,'SUMA','temp_'+hemi+'_schaefer_parcels.nii'))

# Write a new table that contains the information for the header
with open(os.path.join(files_dir,'rh_annot.niml.lt')) as rh_file:
    rh_lines=rh_file.readlines()
with open(os.path.join(files_dir,'lh_annot.niml.lt')) as lh_file:
    lh_lines=lh_file.readlines()

# Combine left and right hemispheres
with open('final_annot.niml.lt','w') as new_file:
    header='<VALUE_LABEL_DTABLE\nni_type="2*String"\nni_dimen="801" >'
    new_file.write(header)
    for line in rh_lines:
        if line != ' ':
            line=line.strip()
            parts=line.split(' ')
            new_line='"'+str(int(parts[0])+200)+'" "'+parts[2]+'"'
        new_file.write("\n")
        new_file.write(new_line)
    for line in lh_lines:
        if line != ' ':
            line=line.strip()
            parts=line.split(' ')
            new_line='"'+parts[0]+'" "'+parts[2]+'"'
        new_file.write("\n")
        new_file.write(new_line)
    new_file.write("\n</VALUE_LABEL_DTABLE>")

# Write the parcel key into the file header 
refit_cmd="3drefit -labeltable {} {}"
refit_cmd=refit_cmd.format('final_annot.niml.lt','schaefer_parcels+orig')
subprocess.run(refit_cmd,shell=True)

## Get parcel value for every voxel to a text file, if desired
# cmd="3dmaskdump -mask schaefer_parcels+orig schaefer_parcels+orig > schaefer_parcels.txt"
# cmd=cmd.format(surfvol_file)
# subprocess.run(cmd,shell=True)

# Clean up directory
for hemi in 'lh','rh':
    os.remove(hemi+'_annot.niml.lt')
os.remove('final_annot.niml.lt')

# Move file with parcels to project directory
for suff in '.HEAD','.BRIK':
    shutil.move(os.path.join(subj_fs_dir,'SUMA','schaefer_parcels+orig'+suff),os.path.join(project_dir,'schaefer_parcels+orig'+suff))

# Move t1 to project directory so you can align everything
shutil.copyfile(os.path.join(subj_fs_dir,'SUMA','T1.nii.gz'),os.path.join(project_dir,'schaefer_t1.nii.gz'))
