#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 13:07:02 2019

@author: jstout

Convert the CTF dipole fitting results from dfit to a Pandas dataframe
This routine searches for opening and closing brackets {} to define a dictionary
The dipole fit results are then stripped of a header and the dipole matrix
is converted to a dataframe

The function assemble_dipole_dframe will iterate over the different latency
marks in the data file and concatenate the results into a single dataframe

Usage:
    assemble_dipole_dframe(dipole_fit_results_filename)
    
    the dipole_fit_results_filename is an output of the dfit routine in CTF
    
Original locations [xp,yp,zp] are converted to 
    LPI format for use with afni [xLPI, yLPI, zLPI]

"""
import sys,os
import pandas as pd
import numpy as np

#Dipole column information for dataframe
dip_columns='fit_idx,//,Trial,Sample,Latency (s),xp,yp,zp,xo,yo,zo,Mom(nAm),ax,ay,az,bx,by,bz,cx,cy,cz,ox,oy,oz,conf(%),Err(%),MEG Err(%),EEG Err(%),Label'
dip_columns=dip_columns.split(',')

########################### Helper Functions #################################
def read_chunk(index, text_list):
    '''Determines if the following line is an open curly brace
    Sets the index as the current line
    Iterates until reaching and closing curly brace'''
    if len(text_list[index+1])<1:
        return index, None,None
    if text_list[index+1][0]=='{':
        key=text_list[index]
        output=''
        while text_list[index+1][0]!='}':
            output+=text_list[index]
            index+=1
        return index,key.rstrip('\n'),output
    else:
        return index,None, None
                
def increment_key(test_key, test_dict):
    '''Verify that the key will not overwrite another value
    Increment until a free key is available'''
    if test_key in test_dict.keys():
        index=1
        while test_key+str(index) in test_dict.keys():
            index+=1
        return test_key+str(index)
    return test_key
        
def dipole_results_to_dict(filename):
    '''Open dipoles results text file and return python dipole dictionary'''
    if not os.path.exists(filename):
        raise ValueError(filename + ' does not exist')
    dipole_dict={}
    fid=open(filename,'r')
    tmp=fid.readlines()
    index=0
    while index < len(tmp)-1:
        index, key,output=read_chunk(index, tmp)
        if key==None:
            index+=1
        else:
            key=increment_key(key, dipole_dict)
            dipole_dict[key]=output
    fid.close()
    return dipole_dict

def remove_fit_result_header_text(fit_text):
    '''Remove the first fifteen header lines of the string to access dipole information'''
    fit_text=fit_text.split('\n')
    return fit_text[15:]

def dipole_dict_to_dframe(dipole_dict=None, dipole_fit_key=None): 
    '''Convert the dipole results to a dataframe
    Run dipole_results_to_dict first to get the dictionary from the results 
    text file outuput from dfit (CTF software)'''
    dipole_text=remove_fit_result_header_text(dipole_dict[dipole_fit_key]) 
    tmp=pd.DataFrame([x.split('\t') for x in dipole_text], columns=dip_columns)
    return tmp.dropna(axis=0, thresh=3)  #3 was set as an arbitrary threshold - hope it doesn't cause issues

def set_column_types(dframe):    
     val_floats=['xp',
                 'yp',
                 'zp',
                 'xo',
                 'yo',
                 'zo',
                 'Mom(nAm)',
                 'ax',
                 'ay',
                 'az',
                 'bx',
                 'by',
                 'bz',
                 'cx',
                 'cy',
                 'cz',
                 'ox',
                 'oy',
                 'oz',
                 'conf(%)',
                 'Err(%)',
                 'MEG Err(%)',
                 'EEG Err(%)',
                 'Latency (s)']     
     val_ints=['Trial',
               'Sample']
     dframe[val_floats]=dframe[val_floats].astype(np.float)
     dframe[val_ints]=dframe[val_ints].astype(np.float).astype(np.int)  #Stackoverflow suggests casting as float then int.  Caused error before
     return dframe

def save_dframe(dframe, dsName):
    '''Save the output of the dataframe in pickle and csv format in a folder
    denoted dipoles in the top subject directory'''
    dipole_folder=os.path.abspath(os.path.join(dsName, 'dipoles'))
    if not os.path.exists(dipole_folder):
        os.path.mkdir(dipole_folder)
    dframe.to_csv(os.path.join(dipole_folder,'ds_name'[:-3]+'.csv'))
    dframe.to_pickle(os.path.join(dipole_folder,'ds_name'[:-3]+'.pkl'))
    
def PRI_to_LPI_scale(p,r,i):
    '''Converts the PRI CTF coordinates to LPI nifti coordinates
    Also converts cm to mm format'''
    return 10.*-r, 10.*p, 10.*i
    

############################ Main #############################################

def assemble_dipole_dframe(filename):
    '''Generates the dipole dictionary
    Loops over all fit_result keys and creates a dataframe
    Concatenates each dipole fit dataframe to the larger'''
    dipole_dict=dipole_results_to_dict(filename)
    fit_result_keys=[i for i in dipole_dict.keys() if 'Fit_Results' in i]
    dipole_dframe=dipole_dict_to_dframe(dipole_dict=dipole_dict, dipole_fit_key=fit_result_keys.pop())    
    for key in fit_result_keys:
        tmpdf=dipole_dict_to_dframe(dipole_dict=dipole_dict, dipole_fit_key=key)
        dipole_dframe=pd.concat([dipole_dframe,tmpdf])
    dipole_dframe=set_column_types(dipole_dframe)
    dipole_dframe=dipole_dframe.reset_index().drop('fit_idx', axis=1)
    dipole_dframe=dipole_dframe.rename(columns={'index':'fit_seq_idx'})
    dipole_dframe['xLPI']=10*-dipole_dframe['yp']   #PRI to LPI  >> -R,P,I
    dipole_dframe['yLPI']=10*dipole_dframe['xp']
    dipole_dframe['zLPI']=10*dipole_dframe['zp']       
    return dipole_dframe 

########### TESTING
def plot_dframe(dframe, thresh=None):
    tmp_dframe=dframe[dframe['MEG Err(%)']<thresh]
    ax=pylab.subplot(111)
    ax.xlim([-10,10])
    ax.ylim([-10,10])
    ax.scatter(tmp_dframe.xp, tmp_dframe.yp)

           
if __name__=='__main__':
    if len(sys.argv) <= 1:
        raise Error('You must input a dipole filename')
    
    file_text=dipole_results_to_dict(sys.argv[1])
    print(file_text)
