import pandas as pd
import os
import numpy as np
import sys

subj = sys.argv[1] # Read in subject of interest
dspm_dir = sys.argv[2] # Read in folder that contains dSPM information
hemi = sys.argv[3] # Read in hemisphere being analyzed

os.chdir(dspm_dir)

# Read in coordinates from the correctly shifted surface; drop the column that doesn't represent 
node_locs=pd.read_csv('myhead-'+hemi+'_inv_shifted.gii.coord.1D.dset',skiprows=9,header=None,delim_whitespace=True)
node_locs=node_locs.rename(columns={0:'virtual_sensor_node',1:'SUMA_X',2:'SUMA_Y',3:'SUMA_Z'})

# Write the second half of nodes if you are iterating through the second (right) hemisphere
if hemi == 'rh':
    node_locs['virtual_sensor_node']=node_locs['virtual_sensor_node']+2562

#Read in corresponding virtual sensors for each vertex
vertex_to_node=pd.read_csv(hemi+'.std.60.1D',skiprows=16,delim_whitespace=True)
vertex_to_node=vertex_to_node.rename(columns={'#0':'vertex','#1':'First node','#2':'Second node','#3':'Third node','#4':'First node weight','#5':'Second node weight','#6':'Third node weight'})

# Make a table that contains all of the information relating virtual sensors to the node
out_df=pd.DataFrame(columns=['SUMA_vertex','virtual_sensor_node','SUMA_X','SUMA_Y','SUMA_Z'])
for vertex in vertex_to_node['vertex']:
    node=float(vertex_to_node['First node'][vertex_to_node[vertex_to_node['vertex']==vertex].index.values])
    if node != -1:
        if hemi == 'rh':
            node=float(vertex_to_node['First node'][vertex_to_node[vertex_to_node['vertex']==vertex].index.values])+2562
        index=(np.where(node_locs['virtual_sensor_node']==node))[0][0]
        SUMA_X=node_locs['SUMA_X'][index]
        SUMA_Y=node_locs['SUMA_Y'][index]
        SUMA_Z=node_locs['SUMA_Z'][index]
        tmp_dict={'SUMA_vertex':vertex,'virtual_sensor_node':node,'SUMA_X':SUMA_X,'SUMA_Y':SUMA_Y,'SUMA_Z':SUMA_Z}
        tmp_df=pd.DataFrame(tmp_dict,index=[node])
        out_df=pd.concat([out_df,tmp_df])
out_df=out_df.reset_index()
out_df=out_df.drop(columns=['index'])

# Write out the table to the dspm directory
if 'SUMA_X' in locals():
    outname=subj+'_'+hemi+'_suma2mne.csv'
    out=os.path.join(dspm_dir,outname)
    out_df.to_csv(out,sep=' ',index=False)

# Remove extra files containing information about the standard mesh
os.chdir(dspm_dir)
for file in os.listdir(os.getcwd()):
    if 'std.60.' in file:
        os.remove(file)
