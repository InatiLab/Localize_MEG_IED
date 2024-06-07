# Input variables
if [ "$#" -eq 4 ]; then
    subj=$1 # Set patient of interest here
    fs_dir=$2 # Set directory containing subject's Freesurfer folders
    dspm_dir=$3 # Set direcotry containing subject's dSPM outputs
    trans_file=$4 # Set path to transformation matrix
else
	echo -e "\033[0;35m++ usage: $0 subj ++\033[0m"
	exit 1
fi

# Set paths
fs_dir=''  # Set directory containing subject's Freesurfer folders
dspm_dir='' # Set direcotry containing subject's dSPM outputs
trans_file='' # Set path to transformation matrix

scripts_dir=$(pwd)
utilities_dir=${scripts_dir}/utilities

#=============================================================================================
# START SCRIPT

# Get transformation matrix into 1D matrix format
cd ${dspm_dir}
if [ ! -f trans_mat.1D ]; then
    python $utilities_dir/write_trans_1D.py "$subj" "$dspm_dir" "$trans_file"
fi

# Run stc2gii_hack, written by MNE-Python. The code is publicly available at https://github.com/jbteves/stc2gii_hack
cd $dspm_dir
stc2gii_hack \
    mymodel.fif \
    full-stcs-lh.stc \
    full-stcs-rh.stc \
    myhead

# Fine-tune surfaces
sides=("lh" "rh")
for hemi in ${sides[@]}; do

    # Transform the stc surface, created above, from the tilted MEG position back to the axialized position
    ConvertSurface -ixmat_1D trans_mat.1D -i myhead-${hemi}.gii -o myhead-${hemi}_inv_shifted.gii
    
    # Gather the XYZ center coordinates of the shifted surface
    SurfaceMetrics -i myhead-${hemi}_inv_shifted.gii -coords
    centerX=`3dBrickStat -mean myhead-${hemi}_inv_shifted.gii.coord.1D.dset[1]`
    centerY=`3dBrickStat -mean myhead-${hemi}_inv_shifted.gii.coord.1D.dset[2]`
    centerZ=`3dBrickStat -mean myhead-${hemi}_inv_shifted.gii.coord.1D.dset[3]`

    # Copy smooth white matter files on the standard 60 mesh from freesurfer/SUMA to the dspm directory
    rsync -rl ${fs_dir}/SUMA/std.60.${hemi}.smoothwm.gii $dspm_dir
    
    # Get coordinates of the white matter file
    SurfaceMetrics -i  std.60.${hemi}.smoothwm.gii -coords
    
    # Get the XYZ center coordinates of the smooth white matter surface
    stdCenterX=`3dBrickStat -mean std.60.${hemi}.smoothwm.gii.coord.1D.dset[1]`
    stdCenterY=`3dBrickStat -mean std.60.${hemi}.smoothwm.gii.coord.1D.dset[2]`
    stdCenterZ=`3dBrickStat -mean std.60.${hemi}.smoothwm.gii.coord.1D.dset[3]`

    # Get the center of the smooth white matter surface 
    scaleX=`ccalc ${stdCenterX} - ${centerX}`
    scaleY=`ccalc ${stdCenterY} - ${centerY}`
    scaleZ=`ccalc ${stdCenterZ} - ${centerZ}`

    # Write each of the scales to the center_al.1D file
    echo "1 0 0 ${scaleX}" >> center_al.1D
    echo "0 1 0 ${scaleY}" >> center_al.1D
    echo "0 0 1 ${scaleZ}" >> center_al.1D

    # Shift the head file to the center of the smooth white matter
    ConvertSurface -xmat_1D center_al.1D -i myhead-${hemi}_inv_shifted.gii -o myhead-${hemi}_centered.gii
    
    # Create a new surface that's in the same space as the smooth white matter, contains the time file, and has the nodes of the standard 60 mesh
    SurfToSurf -i_gii std.60.${hemi}.smoothwm.gii -i_gii myhead-${hemi}_centered.gii -dset myhead-${hemi}.time.gii -prefix $hemi.std.60. 

    # Create a python-readable table from the results of the stc2gii_hack conversion
    cd $utilities_dir
    python mne2suma.py $subj $dspm_dir $hemi
    
    # Clean up directory
    cd $dspm_meg_dir 
    rm center_al.1D
done
