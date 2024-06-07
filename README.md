# Localize_MEG_IED

**Description**
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Source code for localization of interictal epileptiform discharges (IED) using both dSPM and ECD methods. 

**Usage**
------------------------------------------------------------------------------------------------------------------------------------------------------------------------
0. Obtain marks at IED peak from a trained epileptologist
1. do_ECD_analysis.py

   Writes ECD result; obtains earliest and peak timepoints.

   This code requires CTF software (available at https://www.ctf.com/products)
3. do_dspm_analysis.py

   Performs dSPM analysis; writes dSPM solution at first and peak timepoints.

   This code requires MNE-Python installation (available at https://mne.tools/mne-bids/stable/install.html) 
4. do_dspm_clustering.py

   Writes the dSPM solution to the surface; creates clusters using the top 5% of virtual sensors at
   the "first" and "peak" timepoints.
   
   This code requires AFNI installation (available at
   https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/main_toc.html)
   and stc2gii_hack (available at https://github.com/jbteves/stc2gii_hack).  
6. do_get_parcels.py

   Obtains the Schaefer 7Network, 400Parc parcellation from Freesurfer analysis and sends parcellation to
   the volume to relateto dSPM and ECD solutions.

   This code requires AFNI installation and previous running of the Freesurfer recon-all function.
